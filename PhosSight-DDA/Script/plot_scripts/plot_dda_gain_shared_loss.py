#!/usr/bin/env python3
"""
Vertical stacked bar plots for Loss/Shared/Gain to highlight performance lift.

Usage:
  python plot_gain_shared_loss.py --csv <path.csv> --title "Figure title" --output <out.pdf>

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


GREEN = "#3BA44D"   # Loss
BLUE = "#2C5AA0"    # Shared
ORANGE = "#F2A124"  # Gain
FRAME_COLOR = "#D8D8D8"


def _comma(x, pos):
    return f"{int(x):,}" if abs(x) >= 1 else "0"


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


def plot_vertical_bars(df: pd.DataFrame, title: str, output: str):
    # Build x labels: "Dataset - Engine" shortened as needed
    labels = []
    for _, row in df.iterrows():
        ds = row["Dataset"]
        eng = row["Search Engine"]
        labels.append(f"{eng}")

    gain = df["Gain"].values
    shared = df["Shared"].values
    loss = df["Loss"].values

    n = len(df)
    x = range(n)

    plt.figure(figsize=(4, 8), dpi=150)
    ax = plt.gca()

    # Draw a light frame background similar to reference
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color(FRAME_COLOR)

    # Normalize shared values to 0-100 range for each search engine
    # Calculate total for each engine (shared + gain + loss)
    total_values = shared + gain + loss
    
    # Normalize shared to 0-100, then scale gain and loss proportionally
    normalized_shared = np.full(n, 100.0)  # All shared bars are 100
    scale_factors = 100.0 / total_values  # Scale factor for each engine
    
    normalized_gain = gain * scale_factors
    normalized_loss = loss * scale_factors

    # Plot stacked bars: loss (down), shared (middle), gain (top)
    # Implemented by drawing loss as negative, shared at 0, gain on top of shared
    # Narrower bars to match reference style
    bar_width = 0.6
    bar_loss = ax.bar(x, -normalized_loss, width=bar_width, color=GREEN, edgecolor="black", linewidth=0.5)
    bar_shared = ax.bar(x, normalized_shared, width=bar_width, color=BLUE, edgecolor="black", linewidth=0.5)
    bar_gain = ax.bar(x, normalized_gain, width=bar_width, bottom=normalized_shared, color=ORANGE, edgecolor="black", linewidth=0.5)

    # Annotate numbers centered within a section (used for Shared)
    def annotate_stack(bar_container, values, offset_factor=0.5, direction=0, color="black", fontsize=12):
        for rect, val in zip(bar_container, values):
            height = rect.get_height()
            x_center = rect.get_x() + rect.get_width() / 2
            y_base = rect.get_y() + height / 2
            # Offset relative to this bar's absolute height
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

    # Annotate outside the bar: position='top' or 'bottom'
    def annotate_outside(bar_container, values, position="top", color="black", fontsize=12):
        # margin based on overall scale
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
                y = rect.get_y() - margin_bottom * 1.1
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

    # Labels: Loss below entire bar, Gain above entire bar, Shared centered
    annotate_outside(bar_loss, loss, position="bottom", color="black")
    annotate_stack(bar_shared, shared, offset_factor=0.0, direction=0, color="white")
    annotate_outside(bar_gain, gain, position="top", color="black")

    # Axis formatting
    ax.set_xticks(list(x))
    ax.set_xticklabels(labels, rotation=45, ha="right")
    ax.yaxis.set_major_formatter(FuncFormatter(_comma))

    # Fixed Y-axis range: -50 to 200
    ax.set_ylim(-30, 150)  # 150 / 180

    # Title and legend
    ax.set_title(title, fontsize=12, fontweight="bold")
    legend_handles = [
        plt.Rectangle((0, 0), 1, 1, color=GREEN, ec="black", lw=0.5, label="Loss"),
        plt.Rectangle((0, 0), 1, 1, color=BLUE, ec="black", lw=0.5, label="Shared"),
        plt.Rectangle((0, 0), 1, 1, color=ORANGE, ec="black", lw=0.5, label="Gain"),
    ]
    ax.legend(handles=legend_handles, loc="upper left", frameon=True, ncol=3)

    plt.tight_layout()
    out = output or os.path.splitext(os.path.abspath(args.csv))[0] + "_bars.svg"
    plt.savefig(out)
    print(f"Saved figure to: {out}")


def main():
    parser = argparse.ArgumentParser(description="Plot vertical stacked Gain/Shared/Loss bars from CSV")
    parser.add_argument("--csv", required=True, help="Path to CSV file")
    parser.add_argument("--title", required=True, help="Figure title")
    parser.add_argument("--output", default=None, help="Output path (pdf/png)")
    parser.add_argument("--filter-dataset", default=None, help="Optional dataset name filter")
    parser.add_argument("--engine-order", default=None, help="Comma-separated engine order, e.g., 'Comet,MaxQuant,MSGF+,X!Tandem'")
    global args
    args = parser.parse_args()

    df = read_csv_loose(args.csv)
    df = prepare_dataframe(df)

    if args.filter_dataset:
        df = df[df["Dataset"].astype(str) == args.filter_dataset]
        if df.empty:
            raise SystemExit(f"No rows match Dataset == '{args.filter_dataset}'")

    if args.engine_order:
        order = [s.strip() for s in args.engine_order.split(",") if s.strip()]
        cat = pd.Categorical(df["Search Engine"], categories=order, ordered=True)
        df = df.assign(_ord=cat).sort_values(["Dataset", "_ord"]).drop(columns=["_ord"]).reset_index(drop=True)

    plot_vertical_bars(df, args.title, args.output)


if __name__ == "__main__":
    main()


