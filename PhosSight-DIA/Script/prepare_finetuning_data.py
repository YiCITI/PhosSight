#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Prepare fine-tuning data for PhosDetect (DIA version).

Outputs:
- output_train: CSV with header, 2 columns: peptide, label
  - peptide uses PhosSight encoding: lowercase s/t/y denote phosphorylation sites.
  - label: 1 if identified by DIA-NN in for-finetuning results, else 0 (sampled negatives).

- output_predict: CSV with header, 1 column: peptide

Notes:
- DIA-NN report is expected to be a Parquet file containing a 'Modified.Sequence' column.
- Only UniMod:21 (phosphorylation) is supported. Any other modifications in Modified.Sequence
  will be filtered out (to avoid feeding unsupported tokens into PhosDetect).
"""

from __future__ import annotations

import argparse
import os
import random
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Optional, Sequence, Set, Tuple

import pandas as pd


AA_RE = re.compile(r"^[A-Za-z]+$")
UNIMOD_RE = re.compile(r"\(UniMod:(\d+)\)")


def read_fasta_peptides(fasta_path: Path) -> List[Tuple[str, str]]:
    """
    Read FASTA in the same style as
    `PhosDetect/code/fasta_peptide_scoring.py::read_fasta_sequences()`.

    Returns list of (record_id, sequence). record_id is the first token after '>'.
    """
    if not fasta_path.exists():
        raise FileNotFoundError(f"FASTA not found: {fasta_path}")

    records: List[Tuple[str, str]] = []
    current_id: Optional[str] = None
    current_sequence = ""

    with fasta_path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None and current_sequence:
                    records.append((current_id, current_sequence))
                current_id = line[1:].split()[0]
                current_sequence = ""
            else:
                current_sequence += line

        if current_id is not None and current_sequence:
            records.append((current_id, current_sequence))

    return records


def diann_modified_to_phosdetect_peptide(mod_seq: str) -> Optional[str]:
    """
    Convert DIA-NN Modified.Sequence to PhosDetect peptide string.
    - Convert S(UniMod:21)/T(UniMod:21)/Y(UniMod:21) to lowercase s/t/y.
    - Reject peptides containing any UniMod other than 21.
    - Remove any remaining modification annotations like '(...)'.
    """
    if mod_seq is None:
        return None
    s = str(mod_seq)
    if s == "" or s.lower() == "nan":
        return None

    unimods = UNIMOD_RE.findall(s)
    if unimods:
        bad = [u for u in unimods if u != "21"]
        if bad:
            return None

    # Convert phospho marks to lowercase residues
    s = s.replace("S(UniMod:21)", "s")
    s = s.replace("T(UniMod:21)", "t")
    s = s.replace("Y(UniMod:21)", "y")

    # Remove any remaining (...) annotations
    s = re.sub(r"\([^)]*\)", "", s)

    # DIA-NN sometimes uses '_' as separator; remove it
    s = s.replace("_", "")

    if not AA_RE.match(s):
        return None
    return s


def find_report_parquets(result_dir: Path) -> List[Path]:
    """
    Prefer small-scale results: result_dir/for_finetuning_res/**/report.parquet
    Fallback: any report.parquet under result_dir
    """
    preferred = sorted(result_dir.glob("for_finetuning_res/**/report.parquet"))
    if preferred:
        return preferred
    return sorted(result_dir.glob("**/report.parquet"))


def extract_identified_peptides_from_reports(report_paths: Sequence[Path]) -> Set[str]:
    identified: Set[str] = set()
    for p in report_paths:
        df = pd.read_parquet(p)
        if "Modified.Sequence" not in df.columns:
            raise ValueError(
                f"Missing required column 'Modified.Sequence' in {p}. "
                f"Available columns: {list(df.columns)[:50]}"
            )
        for v in df["Modified.Sequence"].dropna().astype(str).tolist():
            pep = diann_modified_to_phosdetect_peptide(v)
            if pep:
                identified.add(pep)
    return identified


def write_train_csv(
    out_path: Path,
    positives: Sequence[str],
    candidates: Sequence[str],
    negative_ratio: float,
    seed: int,
) -> None:
    rng = random.Random(seed)
    pos = sorted(set(positives))
    cand_set = set(candidates)
    pos_set = set(pos)

    negatives_pool = sorted([p for p in cand_set if p not in pos_set])
    if not pos:
        raise ValueError("No positive peptides found. Training set cannot be created.")
    if not negatives_pool:
        raise ValueError("No negative candidates available (all candidates are positives).")

    n_neg = int(round(len(pos) * negative_ratio))
    n_neg = max(1, n_neg)
    if n_neg > len(negatives_pool):
        n_neg = len(negatives_pool)

    neg = rng.sample(negatives_pool, k=n_neg)

    df = pd.DataFrame(
        {"peptide": pos + neg, "label": [1] * len(pos) + [0] * len(neg)}
    )
    df = df.sample(frac=1.0, random_state=seed).reset_index(drop=True)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, sep="\t", index=False)


def main() -> int:
    parser = argparse.ArgumentParser(description="Prepare PhosDetect fine-tuning train/predict files from DIA-NN results.")
    parser.add_argument("--result_dir", type=Path, required=True, help="Directory containing DIA-NN results (expects report.parquet files).")
    parser.add_argument("--fasta_file", type=Path, required=True, help="FASTA of candidate peptides (sequence should already use PhosDetect encoding; s/t/y for phospho).")
    parser.add_argument("--output_train", type=Path, required=True, help="Output CSV path for training (columns: peptide, label).")
    parser.add_argument("--negative_ratio", type=float, default=1.0, help="Negatives per positive in training set (default: 1.0).")
    parser.add_argument("--seed", type=int, default=42, help="Random seed for sampling negatives/shuffling.")

    args = parser.parse_args()
    result_dir = args.result_dir
    fasta_file = args.fasta_file

    if not result_dir.exists():
        raise FileNotFoundError(f"--result_dir not found: {result_dir}")

    report_paths = find_report_parquets(result_dir)
    if not report_paths:
        raise FileNotFoundError(
            f"No report.parquet found under {result_dir}. "
            f"Expected either {result_dir}/for_finetuning_res/**/report.parquet or {result_dir}/**/report.parquet"
        )

    # Parse candidate peptides from FASTA
    fasta_records = read_fasta_peptides(fasta_file)
    candidates = []
    for rec_id, seq in fasta_records:
        seq = seq.strip()
        if not seq:
            continue
        if not AA_RE.match(seq):
            continue
        candidates.append(seq)
    candidates = sorted(set(candidates))
    if not candidates:
        raise ValueError(f"No valid peptide sequences found in FASTA: {fasta_file}")

    # Extract identified peptides from DIA-NN
    try:
        identified = extract_identified_peptides_from_reports(report_paths)
    except Exception as e:
        raise RuntimeError(
            "Failed to parse DIA-NN reports. "
            "If pandas cannot read parquet, install pyarrow (pip install pyarrow). "
            f"Underlying error: {e}"
        ) from e

    # Keep only positives that are in candidate universe (optional but usually desired)
    cand_set = set(candidates)
    positives = sorted([p for p in identified if p in cand_set])

    print(f"Found report.parquet files: {len(report_paths)}")
    print(f"Candidates from FASTA: {len(candidates)}")
    print(f"Identified (after filtering mods): {len(identified)}")
    print(f"Positives in FASTA universe: {len(positives)}")

    write_train_csv(args.output_train, positives, candidates, args.negative_ratio, args.seed)

    print(f"Wrote train TSV: {args.output_train}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

