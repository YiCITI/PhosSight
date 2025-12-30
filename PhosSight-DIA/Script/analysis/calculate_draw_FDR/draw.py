import pandas as pd
import matplotlib.pyplot as plt
import os
import argparse
from pathlib import Path
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

def main(
    FDP_res_path: Path,
    output_dir: Path
):
    """Main function to read FDP results and plot entrapment FDP.
    
    Parameters:
        FDP_res_path (Path): Path to the CSV file containing FDP results.
    """

    # Read the CSV file into a pandas DataFrame
    try:
        df = pd.read_csv(FDP_res_path)
    except FileNotFoundError:
        print(f"Error: The file {FDP_res_path} was not found.")
        exit()

    # Define the mapping for the 'used_spectral_library' column
    library_mapping = {
        'spec_library_original': 'original',
        'spec_library_filtered_by_ratio_0.1': 'Top 10%',
        'spec_library_filtered_by_ratio_0.2': 'Top 20%',
        'spec_library_filtered_by_ratio_0.3': 'Top 30%',
        'spec_library_filtered_by_ratio_0.4': 'Top 40%',
        'spec_library_filtered_by_ratio_0.5': 'Top 50%',
        'spec_library_filtered_by_ratio_0.6': 'Top 60%',
        'spec_library_filtered_by_ratio_0.7': 'Top 70%',
        'spec_library_filtered_by_ratio_0.8': 'Top 80%',
        'spec_library_filtered_by_ratio_0.9': 'Top 90%'
    }

    # Apply the mapping to a new column
    df_use = df[df['used_spectral_library'].isin(library_mapping.keys())].copy()
    df_use['short_library_name'] = df_use['used_spectral_library'].map(library_mapping)

    # Separate the data into pretrained and finetuned
    df_use['type'] = df_use['used_spectral_library'].apply(lambda x: 'pretrained' if 'pretrained' in x else ('finetuned' if 'finetuned' in x else 'original'))
    df_use = df_use.sort_values(by=['short_library_name'])

    df_original = df_use[df_use['type'] == 'original']
    df_finetuned = pd.concat([df_use[df_use['type'] == 'finetuned'], df_original])

    fig, ax = plt.subplots(figsize=figure_size)
    colors2 = [color_style['original'] if s == 'original' else color_style['finetuned']
            for s in df_finetuned['short_library_name']]
    ax.bar(df_finetuned['short_library_name'],
        df_finetuned['entrapment_FDP'],
        color=colors2,
        width=0.6)
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", va="top", rotation_mode="anchor")

    # ax.set_xlabel('Spectral Library Filtering Threshold (Finetuned, Top N%)')
    ax.set_ylabel('Entrapment FDP')
    ax.set_ylim(0, 0.01)
    ax.grid(axis='y', linestyle='--', alpha=0.7)

    plt.tight_layout()

    # Define the output path for the plot
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, 'entrapment_FDP.svg')

    # Save the plot
    plt.savefig(output_file)

    print(f"Plot saved to {output_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Draw entrapment FDP plot.")
    parser.add_argument("--FDP_res_path", type=Path, help="Path to the CSV file containing FDP results.")
    parser.add_argument("--output_dir", type=Path, help="Directory to save the output plot.")
    args = parser.parse_args()
    
    main(args.FDP_res_path, args.output_dir)