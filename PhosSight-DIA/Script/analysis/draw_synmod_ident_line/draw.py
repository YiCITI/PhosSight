import pandas as pd
import matplotlib.pyplot as plt
import os
from PhosSight_paper_style import color_style
from pathlib import Path
import argparse

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

concentrations = ['1x', '2x', '4x', '10x', '20x']
    
def analyze_parquet_file(file_path, name):
    """
    Analyze parquet file, count the number of Modified.Sequence
    """
    print(f"Analyzing {name} file: {file_path}")
    
    # Read parquet file
    df = pd.read_parquet(file_path)
    
    # Filter Protein.Ids containing 'syn'
    if 'Protein.Ids' in df.columns:
        syn_filtered = df[df['Protein.Ids'].str.contains('syn', case=False, na=False)]
    else:
        print("Protein.Ids column not found")
        return {}
    
    # Count unique Modified.Sequence for each concentration
    results = {}
    for conc in concentrations:
        # Find records containing specific concentration in Run column
        if 'Run' in syn_filtered.columns:
            conc_filtered = syn_filtered[syn_filtered['Run'].str.contains(conc, case=False, na=False)]
            
            if 'Modified.Sequence' in conc_filtered.columns:
                unique_count = conc_filtered['Modified.Sequence'].nunique()
                results[conc] = unique_count
                print(f"Unique Modified.Sequence count at {conc} concentration: {unique_count}")
            else:
                print(f"Modified.Sequence column not found")
                results[conc] = 0
        else:
            print("Run column not found")
            results[conc] = 0
    
    return results

def plot_multi_results(results_dict, color_list, output_path):
    """
    Plot multiple comparison line charts
    results_dict: {name: {concentration: count}}
    color_list: list of colors
    output_path: filename to save
    Generate and save line chart, and return results dictionary
    """
    plt.figure(figsize=figure_size)
    for idx, (name, result) in enumerate(results_dict.items()):
        counts = [result.get(c, 0) for c in concentrations]
        plt.plot(concentrations, counts, marker='o', markersize=4, label=name, color=color_list[idx % len(color_list)])
    plt.xlabel('Concentration Ratio')
    plt.ylabel('Number of Identified Synthetic Phosphopeptides')
    plt.ylim(50, 120)
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.show()
    return results_dict


def main(
    res_dir: Path,
    output_dir: Path
):
    """Draw line chart for multiple comparison of modified peptide identification counts
    
    Args:
        res_dir: Path,
            Directory containing result parquet files
        output_dir: Path,
            Directory to save output figures
    """
    name_map = {
        "original": "Original",
        "pretrained_ratio_0.1": "Top 10%",
        "pretrained_ratio_0.2": "Top 20%",
        "pretrained_ratio_0.3": "Top 30%",
        "pretrained_ratio_0.4": "Top 40%",
        "pretrained_ratio_0.5": "Top 50%",
        "pretrained_ratio_0.6": "Top 60%",
        "pretrained_ratio_0.7": "Top 70%",
        "pretrained_ratio_0.8": "Top 80%",
        "pretrained_ratio_0.9": "Top 90%",
        "finetuned_ratio_0.1": "Top 10%",
        "finetuned_ratio_0.2": "Top 20%",
        "finetuned_ratio_0.3": "Top 30%",
        "finetuned_ratio_0.4": "Top 40%",
        "finetuned_ratio_0.5": "Top 50%",
        "finetuned_ratio_0.6": "Top 60%",
        "finetuned_ratio_0.7": "Top 70%",
        "finetuned_ratio_0.8": "Top 80%",
        "finetuned_ratio_0.9": "Top 90%"
    }
    
    original_pretrained_result_paths = [
        res_dir / "original" / "report.parquet",
        res_dir / "pretrained_ratio_0.1" / "report.parquet",
        res_dir / "pretrained_ratio_0.2" / "report.parquet",
        res_dir / "pretrained_ratio_0.3" / "report.parquet",
        res_dir / "pretrained_ratio_0.4" / "report.parquet",
        res_dir / "pretrained_ratio_0.5" / "report.parquet",
        res_dir / "pretrained_ratio_0.6" / "report.parquet",
        res_dir / "pretrained_ratio_0.7" / "report.parquet",
        res_dir / "pretrained_ratio_0.8" / "report.parquet",
        res_dir / "pretrained_ratio_0.9" / "report.parquet"
    ]
    
    original_finetuned_result_paths = [
        res_dir / "original" / "report.parquet",
        res_dir / "finetuned_ratio_0.1" / "report.parquet",
        res_dir / "finetuned_ratio_0.2" / "report.parquet",
        res_dir / "finetuned_ratio_0.3" / "report.parquet",
        res_dir / "finetuned_ratio_0.4" / "report.parquet",
        res_dir / "finetuned_ratio_0.5" / "report.parquet",
        res_dir / "finetuned_ratio_0.6" / "report.parquet",
        res_dir / "finetuned_ratio_0.7" / "report.parquet",
        res_dir / "finetuned_ratio_0.8" / "report.parquet",
        res_dir / "finetuned_ratio_0.9" / "report.parquet"
    ]

    # Define ten distinct colors
    color_list = [
        color_style['original'],
        "#ff7f0e",  # Orange
        "#1f77b4",  # Blue
        "#d62728",  # Red
        "#9467bd",  # Purple
        "#8c564b",  # Brown
        "#e377c2",  # Pink
        "#7f7f7f",  # Gray
        "#bcbd22",  # Yellow-green
        "#17becf",  # Cyan
    ]
    
    print("Start analyzing parquet files...")
    print("=" * 60)
    results_dict = {}
    for res_path in original_pretrained_result_paths:
        name = name_map[os.path.basename(os.path.dirname(str(res_path)))]
        results = analyze_parquet_file(str(res_path), name)
        results_dict[name] = results
    print("\nAnalysis results:")
    for name, result in results_dict.items():
        print(f"{name}: {result}")
    plot_multi_results(results_dict, color_list, 
                       output_path=output_dir / 'lines_comparison_pretrained.svg')
    
    print("Start analyzing parquet files...")
    print("=" * 60)
    results_dict = {}
    for res_path in original_finetuned_result_paths:
        name = name_map[os.path.basename(os.path.dirname(str(res_path)))]
        results = analyze_parquet_file(str(res_path), name)
        results_dict[name] = results
    print("\nAnalysis results:")
    for name, result in results_dict.items():
        print(f"{name}: {result}")
    plot_multi_results(results_dict, color_list, 
                       output_path=output_dir / 'lines_comparison_finetuned.svg')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Draw line chart for multiple comparison of modified peptide identification counts")
    parser.add_argument('--res_dir', type=Path, required=True, help='Directory containing result parquet files')
    parser.add_argument('--output_dir', type=Path, required=True, help='Directory to save output figures')
    args = parser.parse_args()
    
    main(
        res_dir=args.res_dir,
        output_dir=args.output_dir
    )