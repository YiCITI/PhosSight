import pandas as pd
import matplotlib.pyplot as plt
import os
import argparse
from pathlib import Path
from Bio import SeqIO
from PhosSight_paper_style import color_style

# Standard style
plt.rcParams.update({
    'svg.fonttype': 'none',
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'DejaVu Sans', 'Liberation Sans'],
    'font.size': 7,
    'axes.titlesize': 7,
    'axes.labelsize': 7,
    'xtick.labelsize': 7,
    'ytick.labelsize': 7,
    'legend.fontsize': 7,
    'figure.titlesize': 7,
    'axes.linewidth': 1,
    'grid.linewidth': 0.8,
    'lines.linewidth': 1,
    'patch.linewidth': 1,
    'xtick.major.width': 1,
    'ytick.major.width': 1,
    'xtick.minor.width': 1,
    'ytick.minor.width': 1,
    'axes.unicode_minus': False
})
figure_size = (4, 3)  # inches


def to_modified_sequence(sequence: str) -> str:
    """Convert lowercase s/t/y residues to DIA-NN phospho annotation."""
    mod_map = {
        's': 'S(UniMod:21)',
        't': 'T(UniMod:21)',
        'y': 'Y(UniMod:21)'
    }
    return ''.join(mod_map.get(aa, aa) for aa in sequence)


def build_sequence_label_map(fasta_path: Path) -> dict:
    """Build mapping from modified peptide sequence to target/decoy label."""
    if not fasta_path.exists():
        raise FileNotFoundError(f"FASTA file not found: {fasta_path}")

    seq_label_map = {}
    target_count = 0
    decoy_count = 0

    for record in SeqIO.parse(str(fasta_path), "fasta"):
        record_id = str(record.id).lower()
        if 'syn' in record_id:
            label = 'target'
            target_count += 1
        elif 'dec' in record_id:
            label = 'decoy'
            decoy_count += 1
        else:
            continue

        modified_seq = to_modified_sequence(str(record.seq))
        if modified_seq in seq_label_map and seq_label_map[modified_seq] != label:
            raise ValueError(
                f"Sequence label conflict found for {modified_seq}: "
                f"{seq_label_map[modified_seq]} vs {label}"
            )
        seq_label_map[modified_seq] = label

    if not seq_label_map:
        raise ValueError(f"No target/decoy sequences were parsed from FASTA: {fasta_path}")

    print(f"Built sequence-label map from FASTA: {len(seq_label_map)} entries "
          f"(target IDs: {target_count}, decoy IDs: {decoy_count})")
    return seq_label_map


def main(
    parquet_path: Path,
    fasta_path: Path,
    output_dir: Path
):
    """Read synthetic-peptide identifications and compute site-level decoy hit ratio."""
    if parquet_path is None:
        raise ValueError("--parquet_path is required")

    try:
        df = pd.read_parquet(parquet_path)
    except FileNotFoundError:
        print(f"Error: The file {parquet_path} was not found.")
        return

    required_columns = ['Protein.Ids', 'Modified.Sequence']
    missing_columns = [col for col in required_columns if col not in df.columns]
    if missing_columns:
        raise ValueError(f"Missing required columns in parquet file: {missing_columns}")

    syn_mask = df['Protein.Ids'].astype(str).str.contains('syn_pep', regex=False)
    df_syn = df[syn_mask].copy()

    if df_syn.empty:
        raise ValueError("No rows with Protein.Ids containing 'syn_pep' were found.")

    seq_label_map = build_sequence_label_map(fasta_path)

    unimod_cound = df_syn['Modified.Sequence'].fillna('').str.count(r'\(UniMod:')
    unimod21_count = df_syn['Modified.Sequence'].fillna('').str.count(r'\(UniMod:21\)')
    no_unimod21_mask = unimod21_count == 0
    multi_unimod_mask = unimod_cound > 1
    invalid_unimod21_mask = no_unimod21_mask | multi_unimod_mask

    removed_no_unimod21 = int(no_unimod21_mask.sum())
    removed_multi_unimod21 = int(multi_unimod_mask.sum())
    if invalid_unimod21_mask.any():
        df_syn = df_syn.loc[~invalid_unimod21_mask].copy()

    if df_syn.empty:
        raise ValueError("No synthetic rows remain after removing non-single (UniMod:21) entries.")

    df_syn['target_decoy'] = df_syn['Modified.Sequence'].map(seq_label_map)

    missing_mask = df_syn['target_decoy'].isna()
    missing_count = int(missing_mask.sum())
    # if missing_mask.any():
    #     missing_sequences = df_syn.loc[missing_mask, 'Modified.Sequence'].dropna().unique().tolist()
    #     preview = ', '.join(missing_sequences[:10])
    #     print(
    #         f"Warning: {len(missing_sequences)} Modified.Sequence entries are not found in FASTA map. "
    #         f"First 10: {preview}"
    #     )

    df_syn = df_syn.loc[~missing_mask].copy()
    if df_syn.empty:
        raise ValueError("All synthetic rows are missing in FASTA map after filtering.")

    total_hits = len(df_syn)
    target_hits = (df_syn['target_decoy'] == 'target').sum()
    decoy_hits = (df_syn['target_decoy'] == 'decoy').sum()
    esti_FDP = decoy_hits / target_hits

    print(f"Removed no-(UniMod:21) rows: {removed_no_unimod21}")
    print(f"Removed multi-(UniMod:21) rows: {removed_multi_unimod21}")
    print(f"Synthetic identifications (after removing missing): {total_hits}")
    print(f"Removed missing mappings: {missing_count}")
    print(f"Target hits: {target_hits}")
    print(f"Decoy hits: {decoy_hits}")
    print(f"FDP: {esti_FDP:.6f} ({esti_FDP * 100:.3f}%)")
    print(f"Pass 1% threshold: {esti_FDP <= 0.01}")

    if output_dir is not None:
        os.makedirs(output_dir, exist_ok=True)

        fig, ax = plt.subplots(figsize=figure_size)
        ax.bar(
            ['FDP'],
            [esti_FDP],
            color=color_style.get('finetuned', '#4C72B0'),
            width=0.6
        )
        ax.axhline(
            0.01,
            color=color_style.get('original', '#DD8452'),
            linestyle='--',
            linewidth=1,
            label='1% threshold'
        )
        ax.set_ylabel('Site-level decoy ratio')
        ax.set_ylim(0, max(0.012, esti_FDP * 1.2))
        ax.grid(axis='y', linestyle='--', alpha=0.7)
        ax.legend(frameon=False)

        plt.tight_layout()
        output_file = os.path.join(output_dir, 'site_level_decoy_ratio.svg')
        plt.savefig(output_file)
        print(f"Plot saved to {output_file}")


if __name__ == "__main__":
    
    default_fasta = Path("/data1/zhiyuan/PhosSight_analysis/temp/database/syn/syn_peptides_with_decoy_sites.fasta")
    default_parquet = Path("/data1/zhiyuan/PhosSight_analysis/temp/result/syn/original/report.parquet")
    default_output_dir = Path("/data1/zhiyuan/PhosSight/PhosSight-DIA/Script/analysis/output/")
    
    parser = argparse.ArgumentParser(description="Calculate and draw site-level decoy ratio for synthetic peptides.")
    parser.add_argument(
        "--parquet_path", 
        type=Path, 
        default=default_parquet,
        help="Path to peptide-level parquet result file."
        )
    parser.add_argument(
        "--fasta_path",
        type=Path,
        default=default_fasta,
        help="Path to synthetic FASTA with decoy sites."
    )
    parser.add_argument(
        "--output_dir", 
        type=Path, 
        default=default_output_dir, 
        help="Directory to save the output plot."
        )
    args = parser.parse_args()

    main(args.parquet_path, args.fasta_path, args.output_dir)

