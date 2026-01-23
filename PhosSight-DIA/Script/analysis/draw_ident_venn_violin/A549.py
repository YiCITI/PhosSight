import pandas as pd
from pathlib import Path
import pickle
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from itertools import combinations
from scipy.stats import ranksums
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
from PhosSight_paper_style import color_style
from tqdm import tqdm
import argparse
from matplotlib_venn import venn2, _venn2
import seaborn as sns

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
venn_alpha = 0.7

class PSMUpSetPlotAnalyzer_all_and_phos:
    def __init__(self):
        self.all_PSMs = {}
        self.all_peptides = {}
        self.phos_PSMs = {}
        self.phos_peptides = {}


    def load_and_filter(self, thres_report_path_dict: dict, cache_path=None, use_cache=True):
        """
        Load spectra, return all PSMs and phosphopeptide PSMs
        thres_report_path_dict: {threshold_name: parquet_path}
        """
        import pickle
        if cache_path is not None and use_cache and Path(cache_path).exists():
            print(f"Cache file detected, reading directly: {cache_path}")
            with open(cache_path, 'rb') as f:
                cache_data = pickle.load(f)
                self.all_PSMs = cache_data.get('all_PSMs', {})
                self.all_peptides = cache_data.get('all_peptides', {})
                self.phos_PSMs = cache_data.get('phos_PSMs', {})
                self.phos_peptides = cache_data.get('phos_peptides', {})
            for name in self.all_PSMs:
                print(f"{name} PSM count after filtering: {len(self.all_PSMs[name])}")
            for name in self.all_peptides:
                print(f"{name} Peptide count after filtering: {len(self.all_peptides[name])}")
            for name in self.phos_PSMs:
                print(f"{name} Phosphopeptide PSM count after filtering: {len(self.phos_PSMs[name])}")
            for name in self.phos_peptides:
                print(f"{name} Phosphopeptide count after filtering: {len(self.phos_peptides[name])}")
            return
            
        for name, path in thres_report_path_dict.items():
            print(f"Loading and filtering: {name} -> {path}")
            df = pd.read_parquet(path)
            if 'Run' not in df.columns or 'Precursor.Id' not in df.columns:
                raise ValueError(f"{path} missing 'Run' or 'Precursor.Id' column")
            # Combine Run and Precursor.Id as unique PSM id
            df['PSM_ID'] = df['Run'].astype(str) + "|" + df['Precursor.Id'].astype(str)
            self.all_PSMs[name] = set(df['PSM_ID'])
            self.all_peptides[name] = set(df['Modified.Sequence'])
            # Check if required columns exist
            if 'Modified.Sequence' not in df.columns:
                raise ValueError(f"{path} missing 'Modified.Sequence' column")
            self.phos_PSMs[name] = set(df[df['Modified.Sequence'].str.contains("UniMod:21", na=False)]['PSM_ID'])
            self.phos_peptides[name] = set(df[df['Modified.Sequence'].str.contains("UniMod:21", na=False)]['Modified.Sequence'])
            print(f"{name} PSM count after filtering: {len(self.all_PSMs[name])}")
            print(f"{name} Phosphopeptide PSM count after filtering: {len(self.phos_PSMs[name])}")

        if cache_path is not None:
            with open(cache_path, 'wb') as f:
                pickle.dump({'all_PSMs': self.all_PSMs, 'all_peptides': self.all_peptides, 'phos_PSMs': self.phos_PSMs, 'phos_peptides': self.phos_peptides}, f)
            print(f"Filtered results saved to cache: {cache_path}")
        
        
    def plot_venn2(
        self,
        key2, 
        colors, 
        output_file_PSMs,
        output_file_peptides,
        is_phos=True
    ):
        """
        Plot venn2 and save
        
        key2: List of two keys to plot venn2
        Only supports 2 sets
        Example: ["original", "finetuned_top_50%"]
        Ensure self.phos_peptides and self.phos_PSMs contain these two keys
        Example: self.all_PSMs = {"original": set(...), "finetuned_top_50%": set(...)}
        Example: self.phos_PSMs = {"original": set(...), "finetuned_top_50%": set(...)}
        """
        if len(key2) != 2:
            raise ValueError("venn2 only supports 2 sets")
        if is_phos:
            sets = [self.phos_PSMs[label] for label in key2]
        else:
            sets = [self.all_PSMs[label] for label in key2]
        plt.figure(figsize=figure_size)
        venn2(
        # venn2_unweighted(
            sets,
            set_labels=key2,
            set_colors=colors,
            alpha=venn_alpha
        )
        plt.rcParams['svg.fonttype'] = 'none'
        plt.savefig(output_file_PSMs)
        plt.close()
        print(f"PSM level Venn diagram saved to: {output_file_PSMs}")
        
        if is_phos:
            sets = [self.phos_peptides[label] for label in key2]
        else:
            sets = [self.all_peptides[label] for label in key2]
        plt.figure(figsize=figure_size)
        venn2(
        # venn2_unweighted(
            sets,
            set_labels=key2,
            set_colors=colors,
            alpha=venn_alpha
        )
        plt.rcParams['svg.fonttype'] = 'none'
        plt.savefig(output_file_peptides)
        plt.close()
        print(f"Peptide level Venn diagram saved to: {output_file_peptides}")
        
    
    def _read_score_file(
        self,
        score_file,
        id_column,
        score_column
    ):
        """Read score file, return a dictionary {PSM_ID: score_value}"""
        score_path = Path(score_file)
        cache_path = score_path.with_suffix('.pkl')

        if cache_path.exists():
            try:
                with open(cache_path, 'rb') as f:
                    pep_score_dict = pickle.load(f)
                print(f"Detected {score_column} cache, reading directly: {cache_path}")
                return pep_score_dict
            except Exception as exc:
                print(f"Failed to read {score_column} cache, regenerating: {exc}")

        from collections import defaultdict
        pep_score_dict = defaultdict(list)
        df = pd.read_csv(score_path, sep="\t")
        if id_column not in df.columns or score_column not in df.columns:
            raise ValueError(f"{score_path} missing '{id_column}' or '{score_column}' column")
        for _, row in tqdm(df.iterrows(), total=len(df), desc="Reading score file"):
            pep_id = row[id_column]
            score = row[score_column]
            pep_score_dict[pep_id].append(score)
        pep_score_dict = dict(pep_score_dict)
        try:
            with open(cache_path, 'wb') as f:
                pickle.dump(pep_score_dict, f)
            print(f"Score cache saved: {cache_path}")
        except Exception as exc:
            print(f"Failed to save score cache: {exc}")
        return pep_score_dict
        
        
    def get_venn2_groups(
        self,
        key2,
        is_phos,
        peptide_str_converters=None
    ):
        """
        Get peptide ID sets for the three parts of venn2

        key2: List of two keys to plot venn2
        Example: ["original", "finetuned_top_50%"]
        """
        if len(key2) != 2:
            raise ValueError("get_venn2_groups currently only supports two sets")
        
        if is_phos:
            target_dict = self.phos_peptides
        else:
            target_dict = self.all_peptides

        missing_keys = [label for label in key2 if label not in target_dict]
        if missing_keys:
            raise KeyError(f"Keys not found in { 'phos_peptides' if is_phos else 'all_peptides' }: {', '.join(missing_keys)}")
        set_a = target_dict[key2[0]]
        set_b = target_dict[key2[1]]
        
        peptide_str_converters = peptide_str_converters or []
        def apply_converters(pep):
            for converter in peptide_str_converters:
                pep = converter(pep)
            return pep
        
        # _convert_modified_pep_to_sty_pep = lambda mod_pep: mod_pep.replace("S(UniMod:21)", "s").replace("T(UniMod:21)", "t").replace("Y(UniMod:21)", "y").replace("C(UniMod:4)", "C")
        
        set_a = {apply_converters(pep) for pep in set_a}
        set_b = {apply_converters(pep) for pep in set_b}

        group_labels = [f"{key2[0]} only", "Shared", f"{key2[1]} only"]
        group_pep_ids = [
            set_a - set_b,
            set_a & set_b,
            set_b - set_a
        ]
        
        return group_labels, group_pep_ids
    
        
    def plot_violin(
        self,
        key2,
        is_phos,
        colors,
        output_file,
        score_file,
        id_column='sequence',
        score_column='detectability_score',
        y_label=None,
        peptide_str_converters=None
    ):
        """
        Plot violin plot of score distribution for the three parts of venn diagram and save

        key2: List of two keys to plot violin plot
        Example: ["original", "finetuned_top_50%"]
        score_file: Path to the file containing score values
        """
        if score_file is None:
            raise ValueError("plot_violin requires score_file path")

        group_labels, group_pep_ids = self.get_venn2_groups(key2, is_phos, peptide_str_converters)

        pep_score_dict = self._read_score_file(score_file, id_column, score_column)

        # Checking missing peptide IDs in score dictionary
        missing_entries = {
            label: [pep_id for pep_id in ids if pep_id not in pep_score_dict]
            for label, ids in zip(group_labels, group_pep_ids)
        }
        for label, missing in missing_entries.items():
            if missing:
                print(
                    f"Warning: {label} has {len(missing)} peptides missing {score_column}, removed from visualization"
                )

        group_pep_ids = [
            {pep_id for pep_id in ids if pep_id in pep_score_dict}
            for ids in group_pep_ids
        ]
        # TODO: not skip these missing entries
        # missing_details = {
        #     label: missing
        #     for label, missing in missing_entries.items()
        #     if missing
        # }
        # if missing_details:
        #     detail_msg = "; ".join(
        #         f"{label} missing {len(missing)} items" for label, missing in missing_details.items()
        #     )
        #     raise KeyError(
        #         f"Peptide IDs not found in {score_column} file: {detail_msg}"
        #     )

        data = [
            [score for pep_id in ids for score in pep_score_dict[pep_id]]
            for ids in group_pep_ids
        ]

        if any(len(values) == 0 for values in data):
            empty_groups = ", ".join(
                label for label, values in zip(group_labels, data) if len(values) == 0
            )
            raise ValueError(
                f"Detected missing {score_column} data for the following groups: {empty_groups}"
            )

        colors_to_use = _venn2._compute_colors(color_a = colors[0], color_b = colors[1])

        fig, ax = plt.subplots(figsize=figure_size)

        records = [
            {"Score": value, "Type": label}
            for label, values in zip(group_labels, data)
            for value in values
        ]
        df_plot = pd.DataFrame(records)
        
        # Rearranging to have "Shared" in the middle
        colors_to_use = [colors_to_use[0], colors_to_use[2], colors_to_use[1]]

        ax = sns.violinplot(
            data=df_plot,
            x='Type',
            y='Score',
            ax=ax,
            palette=colors_to_use,
            inner=None
        )

        poly_collections = [
            coll for coll in ax.collections if isinstance(coll, PolyCollection)
        ]
        for idx, coll in enumerate(poly_collections[:len(group_labels)]):
            coll.set_alpha(venn_alpha)
            coll.set_edgecolor(colors_to_use[idx])

        violin_positions = list(range(len(group_labels)))
        bp = ax.boxplot(
            data,
            positions=violin_positions,
            patch_artist=True,
            widths=0.25,
            boxprops=dict(facecolor='white'),
            medianprops=dict(color='black'),
            whiskerprops=dict(color='black'),
            capprops=dict(color='black'),
            flierprops=dict(marker='', markersize=0)
        )

        for idx, patch in enumerate(bp['boxes']):
            patch.set_facecolor(colors_to_use[idx])
            patch.set_alpha(venn_alpha)
            patch.set_edgecolor(colors_to_use[idx])

        for idx, group_data in enumerate(data):
            if len(group_data) > 0:
                median_val = pd.Series(group_data).median()
                ax.text(violin_positions[idx] + 0.15, median_val, f'{median_val:.2f}', 
                        ha='left', va='center', fontsize=6, color='black')

        y_values = [value for group in data for value in group]
        y_min, y_max = min(y_values), max(y_values)
        y_range = (y_max - y_min) or 1
        y_start = y_max + 0.05 * y_range
        y_step = 0.08 * y_range

        if len(group_labels) == 3:
            comparison_layout = [
                ((0, 1), 0),
                ((1, 2), 0),
                ((0, 2), 1),
            ]
        else:
            comparison_layout = [
                (pair, idx)
                for idx, pair in enumerate(combinations(range(len(group_labels)), 2))
            ]

        max_level = 0
        for (idx_a, idx_b), level in comparison_layout:
            if idx_a >= len(data) or idx_b >= len(data):
                continue
            group_a = data[idx_a]
            group_b = data[idx_b]
            _, p_val = ranksums(group_a, group_b)
            line_y = y_start + level * y_step
            cap_height = line_y + 0.02 * y_range
            ax.plot([idx_a, idx_a, idx_b, idx_b], [line_y, cap_height, cap_height, line_y], color='black', linewidth=1)
            ax.text((idx_a + idx_b) / 2, cap_height + 0.01 * y_range, f"p = {p_val:.2e}", ha='center', va='bottom', fontsize=6)
            max_level = max(max_level, level)

        ax.set_ylim(top=y_start + (max_level + 1) * y_step)

        ax.set_ylabel(y_label)
        ax.set_xlabel("")
        ax.set_xticklabels(group_labels)
        plt.rcParams['svg.fonttype'] = 'none'
        plt.savefig(output_file)
        plt.close()
        print(f"Violin plot saved to: {output_file}")
        
        
    def calculate_delta_rt_diann_scores(
        self,
        thres_report_path_dict: dict,
        output_score_file: str
    ):
        """
        Calculate Delta RT and save as score file format, using DIA-NN predicted RT
        """
        print(f"Calculating Delta RT and saving to: {output_score_file}")
        
        all_deltas = []
        for name, path in thres_report_path_dict.items():
            print(f"Processing: {name} -> {path}")
            df = pd.read_parquet(path)
            if "Modified.Sequence" not in df.columns or "RT" not in df.columns or "Predicted.RT" not in df.columns:
                print(f"Warning: {path} missing 'Modified.Sequence', 'RT' or 'Predicted.RT' column, skipping")
                continue
            
            for mod_seq, rt, pred_rt in tqdm(zip(df["Modified.Sequence"], df["RT"], df["Predicted.RT"]), total=len(df), desc=f"Calculating Delta RT for {name}"):
                if pd.isna(mod_seq) or pd.isna(rt) or pd.isna(pred_rt):
                    continue
                delta = abs(float(rt) - float(pred_rt))
                all_deltas.append({"peptide": str(mod_seq), "delta_rt": delta})
                
        df_out = pd.DataFrame(all_deltas)
        # Ensure directory exists
        Path(output_score_file).parent.mkdir(parents=True, exist_ok=True)
        df_out.to_csv(output_score_file, sep="\t", index=False)
        print(f"Delta RT calculation completed, {len(df_out)} records, saved to {output_score_file}")
        
        
def main(
    res_dir: Path,
    detectability_score_path: Path,
    output_dir: Path
):    
    analyzer = PSMUpSetPlotAnalyzer_all_and_phos()
    thres_report_path_dict = {
        "Original": res_dir / "original" / "report.parquet",
        "Top 50%": res_dir / "finetuned_ratio_0.5" / "report.parquet"
    }
    analyzer.load_and_filter(thres_report_path_dict, cache_path=output_dir / "venn_cache.pkl")
    
    analyzer.plot_venn2(
        ["Original", "Top 50%"],
        [color_style["original"], color_style["finetuned"]],
        output_dir / "psm_venn2_phos_202503_A549_0h_24h.svg",
        output_dir / "peptide_venn2_phos_202503_A549_0h_24h.svg",
        is_phos=True
    )
    
    analyzer.plot_venn2(
        ["Original", "Top 50%"],
        [color_style["original"], color_style["finetuned"]],
        output_dir / "psm_venn2_all_202503_A549_0h_24h.svg",
        output_dir / "peptide_venn2_all_202503_A549_0h_24h.svg",
        is_phos=False
    )
    
    analyzer.plot_violin(
        ["Original", "Top 50%"],
        True,
        [color_style["original"], color_style["finetuned"]],
        output_dir / "peptide_violin_detectability_phos_202503_A549_0h_24h.svg",
        score_file=detectability_score_path,
        id_column='Sequence',
        score_column='Score',
        y_label="Detectability",
        peptide_str_converters = [lambda s: s.replace("S(UniMod:21)", "s").replace("T(UniMod:21)", "t").replace("Y(UniMod:21)", "y").replace("C(UniMod:4)", "C")]
    )
    
    analyzer.plot_violin(
        ["Original", "Top 50%"],
        False,
        [color_style["original"], color_style["finetuned"]],
        output_dir / "peptide_violin_detectability_all_202503_A549_0h_24h.svg",
        score_file=detectability_score_path,
        id_column='Sequence',
        score_column='Score',
        y_label="Detectability",
        peptide_str_converters = [lambda s: s.replace("S(UniMod:21)", "s").replace("T(UniMod:21)", "t").replace("Y(UniMod:21)", "y").replace("C(UniMod:4)", "C")]
    )
    
    delta_rt_diann_score_file = output_dir / "delta_rt_diann_scores.txt"
    analyzer.calculate_delta_rt_diann_scores(
        thres_report_path_dict,
        output_score_file = delta_rt_diann_score_file
    )
    
    analyzer.plot_violin(
        ["Original", "Top 50%"],
        True,
        [color_style["original"], color_style["finetuned"]],
        output_dir / "peptide_violin_delta_rt_diann_phos_202503_A549_0h_24h.svg",
        score_file=delta_rt_diann_score_file,
        id_column='peptide',
        score_column='delta_rt',
        y_label="Delta RT"
    )
    
    analyzer.plot_violin(
        ["Original", "Top 50%"],
        False,
        [color_style["original"], color_style["finetuned"]],
        output_dir / "peptide_violin_delta_rt_diann_all_202503_A549_0h_24h.svg",
        score_file=delta_rt_diann_score_file,
        id_column='peptide',
        score_column='delta_rt',
        y_label="Delta RT"
    )

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="A549 DIA-NN result venn and violin plot analysis")
    parser.add_argument("--res_dir", type=Path, required=True, help="Directory containing DIA-NN results")
    parser.add_argument("--detect_path", type=Path, required=True, help="Path to detectability score file")
    parser.add_argument("--output_dir", type=Path, required=True, help="Directory to save output plots")
    args = parser.parse_args()
    
    args.output_dir.mkdir(parents=True, exist_ok=True)
    
    main(
        res_dir=args.res_dir,
        detectability_score_path=args.detect_path,
        output_dir=args.output_dir
    )