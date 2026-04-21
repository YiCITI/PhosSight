#!/usr/bin/env python3
"""
Vertical stacked bar plots for Loss/Shared/Gain to highlight performance lift.

Usage:
  python plot_gain_shared_loss_with_msfragger.py [--csv1 ...] [--csv2 ...] [--output ...] [--engine-order ...]

Five search engines by default: Comet, MaxQuant, MSGF+, X!Tandem, MSFragger (same order as Metrics/plot_gain_shared_loss.py).

CSV schema examples:
  results_phosphopeptides.csv:
    Dataset,Search Engine,PhosphoRS,PhosSight,Gain,Shared,Loss,

  results_psm.csv:
    Dataset,Search Engine,PhosphoRS,PhosSight,Gain,Shared,Loss

The script ignores trailing commas and extra whitespace.
"""

import argparse
import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter


# Colors for stacked bars
GREEN = "#7F7F7F"   # Loss (gray)
BLUE = "#1F4E8C"    # Shared (higher saturation blue)
ORANGE = "#E58B00"  # Gain (higher saturation orange)
FRAME_COLOR = "#D8D8D8"
FONT_FAMILY = "sans-serif"
FONT_SIZE = 7

# Global matplotlib style configuration
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "DejaVu Sans", "Liberation Sans"],
    "font.size": 7,
    "axes.titlesize": 7,
    "axes.labelsize": 7,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 7,
    "figure.titlesize": 7,
    "axes.linewidth": 1,
    "grid.linewidth": 0.8,
    "lines.linewidth": 1,
    "patch.linewidth": 1,
    "xtick.major.width": 1,
    "ytick.major.width": 1,
    "xtick.minor.width": 1,
    "ytick.minor.width": 1,
    "axes.unicode_minus": False,
})


def _comma(x, pos):
    return f"{int(x):,}" if abs(x) >= 1 else "0"


def _ensure_csv_input(path: str, arg_name: str) -> None:
    """Reject common mistakes (e.g. passing the output .svg as --csv1)."""
    if not path or not str(path).strip():
        raise SystemExit(f"{arg_name}: path is empty.")
    path = os.path.abspath(path)
    if not os.path.isfile(path):
        raise SystemExit(f"{arg_name}: file not found:\n  {path}")
    ext = os.path.splitext(path)[1].lower()
    if ext in (".svg", ".pdf", ".png", ".jpg", ".jpeg", ".eps", ".webp"):
        raise SystemExit(
            f"{arg_name} points to a figure file ({ext}), not a table. "
            f"Use the *metrics* CSV (Gain/Shared/Loss columns), not the SVG output.\n"
            f"  You passed: {path}\n"
            f"  Example:\n"
            f"    --csv1 deeprescore2_phossight/results_peptides_label_free_deeprescore2_phossight.csv\n"
            f"    --csv2 deeprescore2_phossight/results_peptides_ucec_deeprescore2_phossight.csv\n"
            f"    --output deeprescore2_phossight/results_peptides_deeprescore2_phossight.svg"
        )
    try:
        with open(path, "rb") as f:
            head = f.read(256).lstrip()
    except OSError as e:
        raise SystemExit(f"{arg_name}: cannot read file: {path}\n{e}") from e
    if head.startswith(b"<?xml") or head.startswith(b"<svg"):
        raise SystemExit(
            f"{arg_name}: file looks like SVG/XML, not CSV. Check the path.\n  {path}"
        )


def read_csv_loose(path: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    # Drop any completely empty columns (e.g., trailing comma)
    empty_cols = [c for c in df.columns if df[c].isna().all()]
    if empty_cols:
        df = df.drop(columns=empty_cols)
    # Strip spaces from column names
    df.columns = [c.strip() for c in df.columns]
    return df


def prepare_dataframe(df: pd.DataFrame) -> pd.DataFrame:
    # Required columns
    required = {"Dataset", "Search Engine", "Gain", "Shared", "Loss"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"CSV missing columns: {missing}")
    # Keep only needed columns; preserve order of rows
    df = df[["Dataset", "Search Engine", "Gain", "Shared", "Loss"]].copy()
    # Enforce numeric
    for c in ["Gain", "Shared", "Loss"]:
        df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0).astype(int)
    return df


d = {
    'PXD023665 label free dataset': 'label free',
    'UCEC TMT dataset': 'UCEC',
}

def plot_engine_compare(engine, df_all, output_dir, ax, bar_width=0.20):
    # 挑选指定引擎数据，两行：label free 和 UCEC
    by_engine = df_all[df_all['Search Engine'] == engine].copy().reset_index(drop=True)
    if by_engine.empty:
        ax.set_title(engine, fontsize=FONT_SIZE)
        ax.text(0.5, 0.5, "No data", transform=ax.transAxes, ha="center", va="center", fontsize=FONT_SIZE)
        ax.set_xticks([])
        return
    labels = by_engine["Dataset"].map(d).tolist()
    gain = by_engine["Gain"].values
    shared = by_engine["Shared"].values
    loss = by_engine["Loss"].values
    n = len(by_engine)
    # 两根bar居中，间距gap，且不重叠
    if n == 2:
        bar_width = 0.20
        gap = 0.08
        center = 0.5
        left_x = center - (bar_width/2 + gap/2)
        right_x = center + (bar_width/2 + gap/2)
        x = np.array([left_x, right_x])
    else:
        x = np.arange(n)
    # 全部normalize为Shared为100
    total_values = shared + gain + loss
    normalized_shared = np.full(n, 100.0)
    scale_factors = 100.0 / total_values
    normalized_gain = gain * scale_factors
    normalized_loss = loss * scale_factors
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color(FRAME_COLOR)
    bar_loss = ax.bar(x, -normalized_loss, width=bar_width, color=GREEN, edgecolor="black", linewidth=0.5)
    bar_shared = ax.bar(x, normalized_shared, width=bar_width, color=BLUE, edgecolor="black", linewidth=0.5)
    bar_gain = ax.bar(x, normalized_gain, width=bar_width, bottom=normalized_shared, color=ORANGE, edgecolor="black", linewidth=0.5)
    def annotate_stack(bar_container, values, offset_factor=0.5, direction=0, color="black", fontsize=FONT_SIZE):
        for rect, val in zip(bar_container, values):
            height = rect.get_height()
            x_center = rect.get_x() + rect.get_width() / 2
            y_base = rect.get_y() + height / 2
            y_offset = direction * abs(height) * float(offset_factor)
            ax.text(
                x_center,
                y_base + y_offset,
                val,
                ha="center",
                va="center",
                fontsize=fontsize,
                fontweight="bold",
                rotation=90,
                color=color,
            )
    def annotate_outside(bar_container, values, position="top", color="black", fontsize=FONT_SIZE):
        pos_max = max((normalized_shared + normalized_gain).max(), 1)
        neg_max = max(normalized_loss.max(), 1)
        margin_top = 0.025 * max(pos_max, neg_max)
        margin_bottom = 0.035 * max(pos_max, neg_max)
        for rect, val in zip(bar_container, values):
            x_center = rect.get_x() + rect.get_width() / 2
            if position == "top":
                y = rect.get_y() + rect.get_height() + margin_top
                va = "bottom"
            else:
                y = rect.get_y() - margin_bottom * 1.2
                va = "top"
            ax.text(
                x_center,
                y,
                val,
                ha="center",
                va=va,
                fontsize=fontsize,
                fontweight="bold",
                rotation=90,
                color=color,
            )
    annotate_outside(bar_loss, loss, position="bottom", color="black")
    annotate_stack(bar_shared, shared, offset_factor=0.0, direction=0, color="white")
    annotate_outside(bar_gain, gain, position="top", color="black")
    ax.set_xticks(list(x))
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=FONT_SIZE)
    ax.yaxis.set_major_formatter(FuncFormatter(_comma))
    ax.set_ylim(-30, 170)  # 140/170
    ax.set_title(engine, fontsize=FONT_SIZE)


# Default order matches DeepRescore2 Script/Metrics/results/plot_gain_shared_loss.py (--engine-order)
DEFAULT_ENGINE_ORDER = ["Comet", "MaxQuant", "MSGF+", "X!Tandem", "MSFragger"]


def main():
    parser = argparse.ArgumentParser(
        description="Plot vertical stacked Gain/Shared/Loss bars for two datasets across search engines (includes MSFragger)."
    )
    parser.add_argument("--csv1", default='deeprescore2_phossight/results_peptides_label_free_deeprescore2_phossight.csv', help="Path to CSV file for label-free dataset")
    parser.add_argument("--csv2", default='deeprescore2_phossight/results_peptides_ucec_deeprescore2_phossight.csv', help="Path to CSV file for UCEC dataset")
    parser.add_argument("--output", default='deeprescore2_phossight/results_peptides_deeprescore2_phossight.svg', help="Output SVG path")
    parser.add_argument("--title", default="", help="(Optional) Main title for the entire figure")
    parser.add_argument(
        "--engine-order",
        default=None,
        help="Comma-separated engine order, e.g. 'Comet,MaxQuant,MSGF+,X!Tandem,MSFragger'. Default includes MSFragger.",
    )
    args = parser.parse_args()

    _ensure_csv_input(args.csv1, "--csv1")
    _ensure_csv_input(args.csv2, "--csv2")

    df1 = prepare_dataframe(read_csv_loose(args.csv1))
    df2 = prepare_dataframe(read_csv_loose(args.csv2))
    # 保证两个表有明显不同Dataset字段
    # 如果csv里名字缺乏区分，可手动加，如 df1["Dataset"] = "Label-free" ...
    df_all = pd.concat([df1, df2], ignore_index=True)

    if args.engine_order:
        engines = [s.strip() for s in args.engine_order.split(",") if s.strip()]
    else:
        engines = list(DEFAULT_ENGINE_ORDER)

    n_eng = len(engines)
    # Wider figure when adding a 5th engine (same idea as Metrics plot_gain_shared_loss.py)
    fig_w = max(4.0, 1.0 * n_eng)
    fig, axes = plt.subplots(1, n_eng, figsize=(fig_w, 4), dpi=160, sharey=True)
    if n_eng == 1:
        axes = [axes]

    for idx, (engine, ax) in enumerate(zip(engines, axes)):
        plot_engine_compare(engine, df_all, idx, ax, bar_width=0.48)

    # Add y-axis label to indicate percentage
    axes[0].set_ylabel("Percentage (%)")

    # Create shared legend for all subplots, placed above the figure
    legend_handles = [
        plt.Rectangle((0, 0), 1, 1, color=GREEN, ec="black", lw=0.5, label="Loss"),
        plt.Rectangle((0, 0), 1, 1, color=BLUE, ec="black", lw=0.5, label="Shared"),
        plt.Rectangle((0, 0), 1, 1, color=ORANGE, ec="black", lw=0.5, label="Gain"),
    ]
    
    if args.title:
        fig.suptitle(args.title, fontsize=FONT_SIZE, fontweight="bold", y=0.98)
    
    # Adjust layout first, then add legend
    plt.tight_layout(rect=(0, 0, 1, 0.88))  # Leave space at top for legend
    fig.subplots_adjust(wspace=0.05, top=0.82)  # 横向间距调小，为图例留出空间
    
    # Add shared legend above all subplots
    fig.legend(handles=legend_handles, loc="upper center", frameon=True, ncol=3, 
               fontsize=FONT_SIZE, bbox_to_anchor=(0.5, 0.95))
    out = args.output or "comparison_5engines.svg"
    fig.savefig(out)
    print(f"Saved combined comparison figure to: {out}")

if __name__ == "__main__":
    main()


