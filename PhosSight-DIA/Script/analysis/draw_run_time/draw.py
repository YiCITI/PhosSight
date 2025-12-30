# Calculate the DIA-NN run time.
import os
import re
import matplotlib.pyplot as plt
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

def get_run_time(result_dir):
    """Get DIA-NN run time

    Parameters
    ----------
    result_dir : str
        Result directory

    Returns
    -------
    str | None
        Return run time string, or None if not found
    """
    log_file = os.path.join(result_dir, "report.log.txt")
    report_file = os.path.join(result_dir, "report.parquet")
    
    if not os.path.exists(log_file):
        print(f"Log file does not exist: {log_file}")
        return None
    
    if not os.path.exists(report_file):
        print(f"Report file does not exist: {report_file}")
        return None
    
    with open(log_file, 'r', encoding='utf-8') as f:
        log_content = f.read()
        if "Stats report saved to " in log_content:
            # Extract run time
            for line in log_content.splitlines():
                if "Stats report saved to" in line:
                    match = re.search(r'\[([^\]]+)\]', line)
                    if match:
                        time_str = match.group(1)
                        return time_str
            print("Run time not found")
            return None
        else:
            print(f"DIA-NN run failed or not completed: {result_dir}")
            return None
        
def draw_run_time(
    run_times_pretrained: dict, 
    run_times_finetuned: dict,
    output_path: str
):
    """Plot the run times of pretrained and finetuned models on the same graph, with 'original' as the intersection point"""
    # Get the key order of the two dictionaries, assuming they are the same
    keys = list(run_times_pretrained.keys())
    # 'original' is the first one
    x_labels = [str(k*100) if k != 'original' else 'original' for k in keys]
    def to_minutes(time_str):
        if time_str is None:
            return 0
        m, s = map(int, time_str.split(':'))
        return m + s / 60
    y_pre = [to_minutes(run_times_pretrained[k]) for k in keys]
    y_fine = [to_minutes(run_times_finetuned[k]) for k in keys]
    plt.figure(figsize=figure_size)
    plt.plot(x_labels, y_pre, marker='o', color=color_style['pretrained'], label='Pretrained')
    plt.plot(x_labels, y_fine, marker='o', color=color_style['finetuned'], label='Finetuned')

    # Add numerical labels for each point
    for i, k in enumerate(keys):
        # Pretrained
        if run_times_pretrained[k] is not None:
            plt.annotate(f'{y_pre[i]:.2f}', (x_labels[i], y_pre[i]), textcoords="offset points", xytext=(0,10), ha='center', color=color_style['pretrained'])
        # Finetuned
        if run_times_finetuned[k] is not None:
            plt.annotate(f'{y_fine[i]:.2f}', (x_labels[i], y_fine[i]), textcoords="offset points", xytext=(0,-15), ha='center', color=color_style['finetuned'])

    plt.xlabel('Filtering Threshold (Top N%)')
    plt.ylabel('Run Time (minutes)')
    plt.xticks(rotation=45)
    plt.legend()
    plt.tight_layout()
    plt.rcParams['svg.fonttype'] = 'none'
    plt.savefig(output_path)
    plt.close()
    print(f"Run time plot saved to: {output_path}")
    
    
def main(
    DIA_NN_result_dir: Path,
    output_path: Path
):
    pret_ratios = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    finet_ratios = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    pret_ratio_dirs = [f'pretrained_ratio_{r}' for r in pret_ratios]
    finetuned_ratio_dirs = [f'finetuned_ratio_{r}' for r in finet_ratios]
    origin_pretrained_ratio_dirs = pret_ratio_dirs + ['original']
    origin_finetuned_ratio_dirs = finetuned_ratio_dirs + ['original']
    run_times_pretrained = {}
    for i, dir in enumerate(origin_pretrained_ratio_dirs):
        result_directory = os.path.join(DIA_NN_result_dir, dir)
        time_str = get_run_time(result_directory)
        if i == len(origin_pretrained_ratio_dirs) - 1:
            run_times_pretrained['original'] = time_str
        else:
            run_times_pretrained[pret_ratios[i]] = time_str
    run_times_finetuned = {}
    for i, dir in enumerate(origin_finetuned_ratio_dirs):
        result_directory = os.path.join(DIA_NN_result_dir, dir)
        time_str = get_run_time(result_directory)
        if i == len(origin_finetuned_ratio_dirs) - 1:
            run_times_finetuned['original'] = time_str
        else:
            run_times_finetuned[finet_ratios[i]] = time_str
    # Plotting
    draw_run_time(
        run_times_pretrained,
        run_times_finetuned,
        output_path
    )

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Draw run time comparison between pretrained and finetuned models.")
    parser.add_argument("--DIA_NN_result_dir", type=Path, help="Directory containing DIA-NN results")
    parser.add_argument("--output_path", type=Path, help="Path to save the output plot")
    args = parser.parse_args()

    main(
        DIA_NN_result_dir=args.DIA_NN_result_dir,
        output_path=args.output_path
    )