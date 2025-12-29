import pandas as pd
# import matplotlib.pyplot as plt
from protein_species_mapper import ProteinSpeciesMapper
# from upsetplot import UpSet, from_contents
from pathlib import Path
import argparse
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

# Standard style
font_size = 7
figure_size = (40, 30)  # inches

class PSMUpSetPlotAnalyzer_with_SpeciesFilter:
    def __init__(self, species_fasta_dict, target_species=None):
        if target_species is None:
            target_species = ["yeast", "synthetic"]
        self.mapper = ProteinSpeciesMapper(species_fasta_dict)
        self.target_species = target_species
        self.PSM_results = {}
        self.pep_results = {}


    def load_and_filter(self, thres_report_path_dict: dict, cache_path=None, use_cache=True):
        """
        thres_report_path_dict: {threshold_name: parquet_path}
        Load and filter all report files, keeping only PSMs of target species
        """
        import pickle
        if cache_path is not None and use_cache and Path(cache_path).exists():
            print(f"Cache file detected, reading directly: {cache_path}")
            with open(cache_path, 'rb') as f:
                cached = pickle.load(f)

            self.PSM_results = cached['PSM_results']
            self.pep_results = cached['pep_results']

            if self.PSM_results:
                for name in self.PSM_results:
                    print(f"{name} PSM count after filtering: {len(self.PSM_results[name])}")
                return
            if self.pep_results:
                for name in self.pep_results:
                    print(f"{name} Peptide count after filtering: {len(self.pep_results[name])}")
                return
            
        for name, path in thres_report_path_dict.items():
            print(f"Loading and filtering: {name} -> {path}")
            df = pd.read_parquet(path)
            # Species filtering
            df['Species'] = df['Protein.Ids'].apply(self.mapper.map_ids_to_species)
            df_target = df[df['Species'].isin(self.target_species)].copy()
            # Combine Run and Precursor.Id as unique PSM id
            df_target['PSM_ID'] = df_target['Run'].astype(str) + "|" + df_target['Precursor.Id'].astype(str)
            self.PSM_results[name] = set(df_target['PSM_ID'])
            self.pep_results[name] = set(df_target['Modified.Sequence'])
            print(f"{name} PSM count after filtering: {len(self.PSM_results[name])}")
            print(f"{name} Peptide count after filtering: {len(self.pep_results[name])}")
            
        if cache_path is not None:
            with open(cache_path, 'wb') as f:
                pickle.dump({
                    'PSM_results': self.PSM_results,
                    'pep_results': self.pep_results
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
                if cache_path is not None:
                    matrix.to_csv(out_path)
                print(f"Saved {level} matrix to: {out_path}")
                
            if cache_path is not None:
                base = str(cache_path).rsplit('.', 1)[0]
                save_matrix(self.PSM_results, base + '_PSMs.csv', "PSM")
                save_matrix(self.pep_results, base + '_peptides.csv', "peptide")


    # def plot_upset(self, output_file="psm_upset_plot.svg"):
    #     """
    #     Draw upset plot and save
    #     """
    #     # Construct upsetplot input (using from_contents)
    #     # fig = plt.figure(figsize=figure_size)
    #     fig=plt.figure(figsize=figure_size)
    #     df_upset = from_contents({name: self.PSM_results[name] for name in self.PSM_results})
    #     top_combos = df_upset.index.value_counts().nlargest(15).index
    #     df_upset_top = df_upset[df_upset.index.isin(top_combos)]
    #     upset = UpSet(
    #         df_upset_top,
    #         sort_by="cardinality",
    #         show_counts=True,
    #         element_size=40
    #     )
    #     upset.plot(fig=fig)
    #     plt.savefig(output_file, dpi=300, bbox_inches="tight")
    #     plt.close()
    #     print(f"UpSet plot saved to: {output_file}")
        

def main(
    res_dir: Path,
    fasta_dir: Path,
    output_dir: Path
):
    """Main function to run the PSM UpSet plot analysis.
    
    Args:
        res_dir (Path): Directory containing result parquet files.
        fasta_dir (Path): Directory containing FASTA files.
        output_dir (Path): Directory to save output files.
    """
    fasta_file_synthetic = fasta_dir / "syn_pep_strip_with_protein.fasta"
    fasta_file_yeast = fasta_dir / "uniprotkb_proteome_UP000002311_2025_07_07.fasta"
    fasta_file_ecoli = fasta_dir / "uniprotkb_proteome_UP000000625_2025_09_10.fasta"
    fasta_file_castor = fasta_dir / "uniprotkb_proteome_UP000008311_2025_09_18.fasta"
    # Create batch processor
    species_fasta_dict = {
        "synthetic": str(fasta_file_synthetic),
        "yeast": str(fasta_file_yeast),
        "ecoli": str(fasta_file_ecoli),
        "castor": str(fasta_file_castor)
    }
    
    target_species = ["yeast", "synthetic"]
    
    
    analyzer = PSMUpSetPlotAnalyzer_with_SpeciesFilter(species_fasta_dict, target_species)
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
    analyzer.load_and_filter(thres_report_path_dict, cache_path=output_dir / "upset_plot_JPST000859_finetuned.pkl")
    # Get two .csv files for PSM and peptide matrix
    
    # Do not use this plot for now, use R code from scy 
    # analyzer.plot_upset(output_dir / "psm_upset_plot_JPST000859_finetuned.svg")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='PSM UpSet Plot Analysis for syn yeast Dataset')
    parser.add_argument('--res_dir', type=str, required=True, help='Result directory path containing parquet result files')
    parser.add_argument('--fasta_dir', type=str, required=True, help='FASTA directory path')
    parser.add_argument('--output_dir', type=str, required=True, help='Output directory path')
    args = parser.parse_args()
    
    main(
        Path(args.res_dir), 
        Path(args.fasta_dir), 
        Path(args.output_dir)
    )