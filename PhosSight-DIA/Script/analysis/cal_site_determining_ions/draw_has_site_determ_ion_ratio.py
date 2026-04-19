import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


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


FIGSIZE = (5.0, 3)
PIE_RADIUS = 0.9
PCT_DISTANCE = 0.68


def _format_autopct(pct: float, total: float) -> str:
	if total <= 0 or pd.isna(pct):
		return "0.0%\n(n=0)"
	count = int(round(pct * total / 100.0))
	return f"{pct:.1f}%\n(n={count})"


def _make_autopct(values):
	total = float(sum(values))
	return lambda pct: _format_autopct(pct, total)


def _validate_columns(df: pd.DataFrame) -> None:
	required = {'site_determining_count_a', 'site_determining_count_b'}
	missing = required - set(df.columns)
	if missing:
		missing_cols = ', '.join(sorted(missing))
		raise ValueError(f"CSV is missing required columns: {missing_cols}")


def _normalize_binary(series: pd.Series) -> pd.Series:
	if series.isna().any():
		raise ValueError(f"Column {series.name} contains missing values. Please clean the data before running.")
	# Convert values to a 0/1 presence indicator.
	numeric = pd.to_numeric(series, errors='raise')
	return (numeric > 0).astype(int)


def draw_b_presence_pie(b_binary: pd.Series, output_path: Path) -> None:
	present_count = int((b_binary == 1).sum())
	absent_count = int((b_binary == 0).sum())

	labels = ['B absent', 'B present']
	sizes = [absent_count, present_count]

	# Keep the same color families as the 4-combination pie: blue for B absent, orange for B present.
	colors = ['#5B8FD1', '#E89A5B']

	fig, ax = plt.subplots(figsize=FIGSIZE)
	wedges, _, _ = ax.pie(
		sizes,
		labels=None,
		colors=colors,
		autopct=_make_autopct(sizes),
		pctdistance=PCT_DISTANCE,
		radius=PIE_RADIUS,
		startangle=90,
		counterclock=False,
		wedgeprops={'edgecolor': 'white', 'linewidth': 1},
	)
	ax.axis('equal')
	ax.legend(
		wedges,
		labels,
		loc='center left',
		bbox_to_anchor=(1.0, 0.5),
		frameon=False,
	)
	fig.subplots_adjust(left=0.06, right=0.72, top=0.95, bottom=0.08)
	fig.savefig(output_path)
	plt.close(fig)


def draw_ab_combination_pie(a_binary: pd.Series, b_binary: pd.Series, output_path: Path) -> None:
	combo_counts = {
		'A absent / B absent': int(((a_binary == 0) & (b_binary == 0)).sum()),
		'A present / B absent': int(((a_binary == 1) & (b_binary == 0)).sum()),
		'A absent / B present': int(((a_binary == 0) & (b_binary == 1)).sum()),
		'A present / B present': int(((a_binary == 1) & (b_binary == 1)).sum()),
	}

	labels = list(combo_counts.keys())
	sizes = list(combo_counts.values())

	# Use blue shades for B=0 and orange shades for B=1; same families as the first pie but not identical colors.
	colors = ['#9BBCE5', '#3F73B5', '#F2C49B', '#D6762C']

	fig, ax = plt.subplots(figsize=FIGSIZE)
	wedges, _, _ = ax.pie(
		sizes,
		labels=None,
		colors=colors,
		autopct=_make_autopct(sizes),
		pctdistance=PCT_DISTANCE,
		radius=PIE_RADIUS,
		startangle=90,
		counterclock=False,
		wedgeprops={'edgecolor': 'white', 'linewidth': 1},
	)
	ax.axis('equal')
	ax.legend(
		wedges,
		labels,
		loc='center left',
		bbox_to_anchor=(1.0, 0.5),
		frameon=False,
	)
	fig.subplots_adjust(left=0.06, right=0.72, top=0.95, bottom=0.08)
	fig.savefig(output_path)
	plt.close(fig)


def main(input_csv: Path) -> None:
	if not input_csv.exists() or not input_csv.is_file():
		raise FileNotFoundError(f"Input CSV does not exist or is not a file: {input_csv}")

	output_dir = input_csv.parent

	df = pd.read_csv(input_csv)
	_validate_columns(df)

	a_binary = _normalize_binary(df['site_determining_count_a'])
	b_binary = _normalize_binary(df['site_determining_count_b'])

	stem = input_csv.stem
	b_presence_output = output_dir / f"{stem}_b_presence_ratio_pie.svg"
	ab_combo_output = output_dir / f"{stem}_ab_01_combination_ratio_pie.svg"

	draw_b_presence_pie(b_binary, b_presence_output)
	draw_ab_combination_pie(a_binary, b_binary, ab_combo_output)

	print(f"Saved: {b_presence_output}")
	print(f"Saved: {ab_combo_output}")


if __name__ == '__main__':
	parser = argparse.ArgumentParser(
		description='Draw pie charts for site_determining_count_a/b statistics from CSV.'
	)
	parser.add_argument('input_csv', type=Path, help='Input CSV file path')
	args = parser.parse_args()

	main(input_csv=args.input_csv)
