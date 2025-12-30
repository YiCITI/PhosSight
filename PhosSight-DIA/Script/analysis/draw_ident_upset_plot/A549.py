import pandas as pd
import argparse
from pathlib import Path
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

# Do not draw in this script anymore, just prepare the data for R script from scy

class PSMUpSetPlotAnalyzer_without_filter:
    def __init__(self):
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
                # self.all_PSMs = cache_data.get('all_PSMs', {})
                self.phos_PSMs = cache_data.get('phos_PSMs', {})
                self.phos_peptides = cache_data.get('phos_peptides', {})
            # for name in self.all_PSMs:
            #     print(f"{name} PSM count after filtering: {len(self.all_PSMs[name])}")
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
            # self.all_PSMs[name] = set(df['PSM_ID'])
            # Check if required columns exist
            if 'Modified.Sequence' not in df.columns:
                raise ValueError(f"{path} missing 'Modified.Sequence' column")
            self.phos_PSMs[name] = set(df[df['Modified.Sequence'].str.contains("UniMod:21", na=False)]['PSM_ID'])
            self.phos_peptides[name] = set(df[df['Modified.Sequence'].str.contains("UniMod:21", na=False)]['Modified.Sequence'])
            # print(f"{name} PSM count after filtering: {len(self.all_PSMs[name])}")
            print(f"{name} Phosphopeptide PSM count after filtering: {len(self.phos_PSMs[name])}")

        if cache_path is not None:
            with open(cache_path, 'wb') as f:
                pickle.dump({
                    # 'all_PSMs': self.all_PSMs, 
                    'phos_PSMs': self.phos_PSMs,
                    'phos_peptides': self.phos_peptides
                }, f)
            print(f"Filtered results saved to cache: {cache_path}")
            # Serialize to matrix and save as csv
            def save_matrix(ident_dict, out_path, level):
                all_identifiers = set()
                for s in ident_dict.values():
                    all_identifiers.update(s)
                all_identifiers = sorted(all_identifiers)
                matrix = pd.DataFrame(0, index=all_identifiers, columns=list(ident_dict.keys()))
                for name, s in ident_dict.items():
                    matrix.loc[list(s), name] = 1
                matrix.to_csv(out_path)
                print(f"Saved {level} matrix to: {out_path}")

            if cache_path is not None:
                base = str(cache_path).rsplit('.', 1)[0]
                # save_psm_matrix(self.all_PSMs, base + '_all.csv')
                save_matrix(self.phos_PSMs, base + '_phos_PSMs.csv', "PSM")
                save_matrix(self.phos_peptides, base + '_phos_peptides.csv', "peptide")


    # def plot_upset(
    #     self, 
    #     # output_file_all, 
    #     output_file_phos
    # ):
    #     """
    #     Draw upset plot and save
    #     """
    #     # # Construct upsetplot input (using from_contents)
    #     # fig1 = plt.figure(figsize=figure_size)
    #     # df_upset = from_contents({name: self.all_PSMs[name] for name in self.all_PSMs})
    #     # top_combos = df_upset.index.value_counts().nlargest(15).index
    #     # df_upset_top = df_upset[df_upset.index.isin(top_combos)]
    #     # upset = UpSet(
    #     #     df_upset_top,
    #     #     sort_by="cardinality",
    #     #     show_counts=True,
    #     #     element_size=40
    #     # )
    #     # upset.plot(fig=fig1)
    #     # for ax in fig1.get_axes():
    #     #     ax.tick_params(labelsize=font_size)
    #     #     ax.xaxis.label.set_size(font_size)
    #     #     ax.yaxis.label.set_size(font_size)
    #     #     for label in ax.get_xticklabels() + ax.get_yticklabels():
    #     #         label.set_fontsize(font_size)
    #     # # plt.title("Overlap of Identified PSMs Across Filtering Methods")
    #     # # plt.tight_layout()
    #     # plt.savefig(output_file_all, bbox_inches="tight")
    #     # plt.close()
    #     # print(f"UpSet plot saved to: {output_file_all}")
        
    #     fig2 = plt.figure(figsize=figure_size)
    #     df_upset = from_contents({name: self.phos_PSMs[name] for name in self.phos_PSMs})
    #     top_combos = df_upset.index.value_counts().nlargest(15).index
    #     df_upset_top = df_upset[df_upset.index.isin(top_combos)]
    #     upset = UpSet(
    #         df_upset_top,
    #         sort_by="cardinality",
    #         show_counts=True,
    #         element_size=40
    #     )
    #     upset.plot(fig=fig2)
    #     for ax in fig2.get_axes():
    #         ax.tick_params(labelsize=font_size)
    #         ax.xaxis.label.set_size(font_size)
    #         ax.yaxis.label.set_size(font_size)
    #         for label in ax.get_xticklabels() + ax.get_yticklabels():
    #             label.set_fontsize(font_size)
    #     # plt.title("Overlap of Identified PSMs Across Filtering Methods")
    #     # plt.tight_layout()
    #     plt.savefig(output_file_phos, bbox_inches="tight")
    #     plt.close()
    #     print(f"UpSet plot saved to: {output_file_phos}")


def main(
    res_dir: Path,
    output_dir: Path
):
    """Main function to run the PSM UpSet plot analysis.
    
    Args:
        res_dir (Path): Directory containing result parquet files.
        output_dir (Path): Directory to save output files.
    """
    thres_report_path_dict = {
        "Original": res_dir / "original" / "report.parquet",
        "Top 10%": res_dir / "finetuned_ratio_0.1" / "report.parquet",
        "Top 20%": res_dir / "finetuned_ratio_0.2" / "report.parquet",
        "Top 30%": res_dir / "finetuned_ratio_0.3" / "report.parquet",
        "Top 40%": res_dir / "finetuned_ratio_0.4" / "report.parquet",
        "Top 50%": res_dir / "finetuned_ratio_0.5" / "report.parquet",
        "Top 60%": res_dir / "finetuned_ratio_0.6" / "report.parquet",
        "Top 70%": res_dir / "finetuned_ratio_0.7" / "report.parquet",
        "Top 80%": res_dir / "finetuned_ratio_0.8" / "report.parquet",
        "Top 90%": res_dir / "finetuned_ratio_0.9" / "report.parquet"
    }
    analyzer = PSMUpSetPlotAnalyzer_without_filter()
    analyzer.load_and_filter(thres_report_path_dict, cache_path=output_dir / "upset_plot_202503_A549_finetuned.pkl")
    # Get two .csv files for PSM and peptide matrix
    
    # Do not use this plot for now, use R code from scy 
    # analyzer.plot_upset(
    #     output_dir / "psm_upset_plot_phos_202503_A549_0h_24h_finetuned.svg"
    # )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='PSM UpSet Plot Analysis for A549 Dataset')
    parser.add_argument('--res_dir', type=str, required=True, help='Result directory path containing parquet result files')
    parser.add_argument('--output_dir', type=str, required=True, help='Output directory path')
    args = parser.parse_args()
    
    main(
        Path(args.res_dir), 
        Path(args.output_dir)
    )