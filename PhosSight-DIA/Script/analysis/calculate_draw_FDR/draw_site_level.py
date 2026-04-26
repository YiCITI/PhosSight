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


def build_sequence_label_map(fasta_path: Path) -> tuple[dict, int, int]:
    """Build mapping from modified peptide sequence to target/decoy label.

    Returns
    -------
    tuple[dict, int, int]
        Sequence-label map, target record count, decoy record count.
    """
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

    if target_count == 0 or decoy_count == 0:
        raise ValueError(
            f"Invalid FASTA target/decoy counts parsed from {fasta_path}: "
            f"target={target_count}, decoy={decoy_count}"
        )

    decoy_to_target_ratio = decoy_count / target_count

    print(f"Built sequence-label map from FASTA: {len(seq_label_map)} entries "
          f"(target IDs: {target_count}, decoy IDs: {decoy_count}, "
          f"decoy/target ratio: {decoy_to_target_ratio:.6f})")
    return seq_label_map, target_count, decoy_count


def calculate_fdp(
    parquet_path: Path,
    seq_label_map: dict,
    target_count: int,
    decoy_count: int
) -> float:
    """Calculate site-level FDP for a single parquet file."""
    if parquet_path is None:
        raise ValueError("parquet_path is required")

    if target_count <= 0 or decoy_count <= 0:
        raise ValueError(
            f"Invalid target/decoy counts for correction: "
            f"target_count={target_count}, decoy_count={decoy_count}"
        )

    try:
        df = pd.read_parquet(parquet_path)
    except FileNotFoundError:
        print(f"Error: The file {parquet_path} was not found.")
        return None

    required_columns = ['Protein.Ids', 'Modified.Sequence']
    missing_columns = [col for col in required_columns if col not in df.columns]
    if missing_columns:
        raise ValueError(f"Missing required columns in parquet file {parquet_path}: {missing_columns}")

    syn_mask = df['Protein.Ids'].astype(str).str.contains('syn_pep', regex=False)
    df_syn = df[syn_mask].copy()

    if df_syn.empty:
        print(f"Warning: No rows with Protein.Ids containing 'syn_pep' were found in {parquet_path}.")
        return None

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
        print(f"Warning: No synthetic rows remain after removing non-single (UniMod:21) entries in {parquet_path}.")
        return None

    df_syn['target_decoy'] = df_syn['Modified.Sequence'].map(seq_label_map)

    missing_mask = df_syn['target_decoy'].isna()
    missing_count = int(missing_mask.sum())

    df_syn = df_syn.loc[~missing_mask].copy()
    if df_syn.empty:
        print(f"Warning: All synthetic rows are missing in FASTA map after filtering in {parquet_path}.")
        return None

    total_hits = len(df_syn)
    target_hits = (df_syn['target_decoy'] == 'target').sum()
    decoy_hits = (df_syn['target_decoy'] == 'decoy').sum()

    raw_fdp = decoy_hits / target_hits if target_hits > 0 else float('inf')
    correction_factor = target_count / decoy_count
    corrected_fdp = raw_fdp * correction_factor if target_hits > 0 else float('inf')

    print(f"\n--- Statistics for {parquet_path.name} ---")
    print(f"Removed no-(UniMod:21) rows: {removed_no_unimod21}")
    print(f"Removed multi-(UniMod:21) rows: {removed_multi_unimod21}")
    print(f"Synthetic identifications (after removing missing): {total_hits}")
    print(f"Removed missing mappings: {missing_count}")
    print(f"Target hits: {target_hits}")
    print(f"Decoy hits: {decoy_hits}")
    print(f"Library correction factor (target/decoy): {correction_factor:.6f} "
        f"({target_count}/{decoy_count})")
    print(f"Raw FDP (decoy_hits/target_hits): {raw_fdp:.6f} ({raw_fdp * 100:.3f}%)")
    print(f"Corrected FDP: {corrected_fdp:.6f} ({corrected_fdp * 100:.3f}%)")

    return corrected_fdp


def main(
    parquet_paths: list,
    labels: list,
    colors: list,
    fasta_path: Path,
    output_dir: Path
):
    """Read synthetic-peptide identifications and compute site-level FDP for multiple files."""
    
    seq_label_map, target_count, decoy_count = build_sequence_label_map(fasta_path)

    fdp_values = []
    plot_labels = []
    plot_colors = []

    for path, label, color in zip(parquet_paths, labels, colors):
        fdp = calculate_fdp(path, seq_label_map, target_count, decoy_count)
        if fdp is not None:
            fdp_values.append(fdp)
            plot_labels.append(label)
            plot_colors.append(color)

    if not fdp_values:
        print("No valid FDP values to plot.")
        return

    if output_dir is not None:
        os.makedirs(output_dir, exist_ok=True)

        fig, ax = plt.subplots(figsize=figure_size)
        
        # Plot multiple FDP bars
        ax.bar(
            plot_labels,
            fdp_values,
            color=plot_colors,
            width=0.4
        )
        ax.set_ylabel('Site-level FDP')
        ax.set_ylim(0, 0.06)
        ax.grid(axis='y', linestyle='--', alpha=0.7)

        plt.tight_layout()
        output_file = os.path.join(output_dir, 'site_level_FDP.svg')
        plt.savefig(output_file)
        print(f"\nPlot saved to {output_file}")


if __name__ == "__main__":
    
    default_fasta = Path("/data1/zhiyuan/PhosSight_analysis/temp/database/syn/syn_peptides_with_decoy_sites.fasta")
    
    # Example paths - you can modify these 3 files as needed
    default_parquets = [
        Path("/data1/zhiyuan/PhosSight_analysis/temp/result/syn/original/report.parquet"),
        Path("/data1/zhiyuan/PhosSight_analysis/temp/result/syn/pretrained_ratio_0.5/report.parquet"),
        Path("/data1/zhiyuan/PhosSight_analysis/temp/result/syn/finetuned_ratio_0.5/report.parquet")
    ]
    
    default_output_dir = Path("/data1/zhiyuan/github_repo/PhosSight/PhosSight-DIA/Script/analysis/output")
    
    parser = argparse.ArgumentParser(description="Calculate and draw site-level FDP for synthetic peptides.")
    parser.add_argument(
        "--parquet_paths", 
        type=Path, 
        nargs='+',
        default=default_parquets,
        help="Paths to multiple peptide-level parquet result files."
    )
    parser.add_argument(
        "--labels",
        type=str,
        nargs='+',
        default=["Original", "Pretrained top 50%", "Finetuned top 50%"],
        help="Labels for the parquet files in the bar plot."
    )
    parser.add_argument(
        "--colors",
        type=str,
        nargs='+',
        default=[color_style.get('original'), color_style.get('pretrained'), color_style.get('finetuned')],
        help="Colors for the bars."
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

    if len(args.parquet_paths) != len(args.labels) or len(args.parquet_paths) != len(args.colors):
        print("Warning: The number of parquet paths, labels, and colors should match.")

    main(args.parquet_paths, args.labels, args.colors, args.fasta_path, args.output_dir)

