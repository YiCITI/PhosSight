# Entrapment FDP calculator for batch processing of search results

import os
import argparse
import pandas as pd
from pathlib import Path
from typing import Dict, Tuple, List, Optional
import logging
from tqdm import tqdm
from protein_species_mapper import ProteinSpeciesMapper

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


class EntrapmentFDPCalculator:
    """Calculator class for calculating entrapment FDP of search results"""

    def __init__(self, 
                 species_fasta_dict: Dict[str, str], 
                 target_species: List[str] = ["yeast", "synthetic"]
        ):
        """
        Initialize the calculator
        
        species_fasta_dict: Dictionary mapping species name (str) -> fasta file path (str)
        target_species: List of target species for ratio calculation, default is ["yeast", "synthetic"]
        """
        self.mapper = ProteinSpeciesMapper(species_fasta_dict)
        self.target_species = target_species

    def _read_result_file(self, result_file: str) -> Optional[pd.DataFrame]:
        """Read DIA-NN result file and map protein IDs to species"""
        if not os.path.exists(result_file):
            logger.error(f"Result file does not exist: {result_file}")
            return None

        logger.info(f"Reading result file: {result_file}")
        df = pd.read_parquet(result_file)
        try:
            df['Species'] = list(tqdm(
                (self.mapper.map_ids_to_species(pid) for pid in df['Protein.Ids']),
                total=len(df),
                desc='Mapping protein IDs to species'
            ))
        except ValueError as e:
            logger.error(f"Error occurred while mapping protein IDs to species: {e}")
            raise
        
        return df
    
    
    def _calculate_target_entrapment_ratio(self, spectral_library_file: str) -> Optional[float]:
        """Calculate the ratio of total peptides for yeast and ecoli in the spectral library"""
        if not os.path.exists(spectral_library_file):
            logger.error(f"Spectral library file does not exist: {spectral_library_file}")
            return None
        
        logger.info(f"Reading spectral library file: {spectral_library_file}")
        df_spec_lib = pd.read_parquet(spectral_library_file)
        
        df_spec_lib['Species'] = list(tqdm(
            (self.mapper.map_ids_to_species(pid) for pid in df_spec_lib['Protein.Ids']),
            total=len(df_spec_lib),
            desc='Mapping spectral library protein IDs to species'
        ))
        
        target_count = sum(df_spec_lib['Species'].isin(self.target_species))
        multi_species_count = df_spec_lib[df_spec_lib['Species'].str.contains(';')].shape[0]
        entrapment_count = df_spec_lib.shape[0] - target_count - multi_species_count

        if entrapment_count == 0:
            logger.error("No entrapment peptides in spectral library, cannot calculate ratio")
            return None
        
        ratio = target_count / entrapment_count
        logger.info(f"target/entrapment ratio: {ratio:.4f} ({target_count}/{entrapment_count})")
        return ratio


    def _calculate_entrapment_FDP(self, result_df: pd.DataFrame, target_entrapment_ratio: float) -> Tuple[float, float, Dict[str, int]]:
        """Calculate entrapment FDP of search results
        
        Parameters
        ----------
        result_df : pd.DataFrame
            DataFrame containing search results, must include 'Species' column
        target_entrapment_ratio : float
            Ratio of target to entrapment peptides in the spectral library
            
        Returns
        -------
        Tuple[float, float, Dict[str, int]]
            entrapment_FDP : float
                Entrapment FDP calculated
            dict_nums : Dict[str, int]
                Dictionary containing identification counts for each category, keys include:
                - 'num_yeast': Number of identified yeast
                - 'num_synthetic': Number of identified synthetic
                - 'num_target': Number of identified target species
                - 'num_multi_species': Number of identified multi-species peptides
                - 'num_entrapment': Number of identified entrapment peptides
        """
        ident_yeast_num = result_df[result_df['Species'] == 'yeast'].shape[0]
        ident_synthetic_num = result_df[result_df['Species'] == 'synthetic'].shape[0]
        ident_target_num = sum(result_df['Species'].isin(self.target_species))
        ident_multi_species_num = result_df[result_df['Species'].str.contains(';')].shape[0]
        ident_entrapment_num = result_df.shape[0] - ident_target_num - ident_multi_species_num
        
        logger.info(f"Identification results statistics - yeast: {ident_yeast_num}, entrapment: {ident_entrapment_num}, synthetic: {ident_synthetic_num}")
        
        # Method:
        # (Number of identified entrapment + Total target candidates / Total entrapment candidates * Number of identified entrapment)
        #                                   /
        # (Number of identified target and entrapment peptides)
        if (ident_yeast_num + ident_synthetic_num + ident_entrapment_num) == 0:
            raise ValueError("Sum of identified yeast, synthetic peptides, and entrapment peptides is 0, cannot calculate entrapment FDP")
        else:
            entrapment_FDP = (ident_entrapment_num + target_entrapment_ratio * ident_entrapment_num) / (ident_target_num + ident_entrapment_num)

        return entrapment_FDP, {
            "num_yeast": ident_yeast_num,
            "num_synthetic": ident_synthetic_num,
            "num_target": ident_target_num,
            "num_multi_species": ident_multi_species_num,
            "num_entrapment": ident_entrapment_num
        }


    def calculate_for_single_result(self, result_file: str, spectral_library_file: str) -> Tuple[float, float, Dict[str, int]]:
        """Calculate entrapment FDP for a single result file"""
        logger.info(f"Start processing result file: {result_file}")
        
        # Calculate ratio of target and entrapment peptides
        target_entrapment_ratio = self._calculate_target_entrapment_ratio(spectral_library_file)
        if target_entrapment_ratio is None:
            raise ValueError("Failed to calculate target and entrapment peptide ratio")
        
        # Read result file
        df_result = self._read_result_file(result_file)
        if df_result is None:
            raise ValueError("Failed to read result file")
        # Save df_result to file for checking Species column content
        output_path = Path(result_file).with_suffix('.with_species.csv')
        df_result.to_csv(output_path, index=False, encoding='utf-8-sig')
        logger.info(f"Saved result with Species column to: {output_path}")

        # Calculate entrapment FDP
        entrapment_FDP, dict_nums = self._calculate_entrapment_FDP(df_result, target_entrapment_ratio)
        
        logger.info(f"Entrapment FDP: {entrapment_FDP:.6f}")
        
        return entrapment_FDP, dict_nums
    

class BatchEntrapmentFDPProcessor:
    """Processor class for batch processing entrapment FDP calculation"""
    
    def __init__(self, base_dir: str, 
                 species_fasta_dict: Dict[str, str]):
        """
        Initialize batch processor
        
        Parameters
        ----------
        base_dir : str
            Base directory containing all result folders
        species_fasta_dict : Dict[str, str]
            Dictionary mapping species name (str) -> fasta file path (str)
            Example: {
                "synthetic": "/path/to/synthetic.fasta",
                "yeast": "/path/to/yeast.fasta",
                "ecoli": "/path/to/ecoli.fasta",
                "castor": "/path/to/castor.fasta"
            }
        """
        self.base_dir = Path(base_dir)
        self.calculator = EntrapmentFDPCalculator(
            species_fasta_dict=species_fasta_dict
        )

    def process_all_results_spec_libs(self, res_spec_lib_path_dict: Dict[str, str]) -> pd.DataFrame:
        """
        Batch process all result folders
        
        Parameters
        ----------
        res_spec_lib_path_dict : Dict[str, str]
            Mapping of result file paths to spectral library file paths
            Example: {
                "/path/to/original/report.parquet": "/path/to/original/spec_library_original.parquet",
                "/path/to/pretrained_ratio_0.5/report.parquet": "/path/to/pretrained_ratio_0.5/spec_library_filtered_by_ratio_0.5.parquet",
                ...
            }
        Returns
        -------
        pd.DataFrame
            DataFrame containing all results
        """
        
        results = []

        for result_file, spectral_library_file in res_spec_lib_path_dict.items():
            try:
                entrapment_FDP, dict_nums = self.calculator.calculate_for_single_result(
                    str(result_file), 
                    str(spectral_library_file)
                )

                results.append({
                    'used_spectral_library': Path(spectral_library_file).stem,
                    'entrapment_FDP': entrapment_FDP,
                    **dict_nums
                })
                
            except Exception as e:
                logger.error(f"Error processing file {result_file}: {e}")
                results.append({
                    'used_spectral_library': Path(spectral_library_file).stem,
                    'entrapment_FDP': None,
                    'error': str(e)
                })
        
        df_results = pd.DataFrame(results)
        logger.info(f"Batch processing completed, processed {len(results)} files")

        return df_results
    
    def save_results_to_csv(self, results_df: pd.DataFrame, output_file: str = None):
        """Save results to CSV file"""
        if output_file is None:
            output_file = self.base_dir / "entrapment_FDP_results.csv"
        
        results_df.to_csv(output_file, index=False, encoding='utf-8-sig')
        logger.info(f"Results saved to: {output_file}")


def main(
    fasta_dir: Path,
    DIA_NN_result_dir: Path,
    spec_lib_dir: Path,
    output_dir: Path
):
    # Ensure output directory exists
    output_dir.mkdir(parents=True, exist_ok=True)

    # Configure file logging
    log_file = output_dir / "entrapment_FDP_calculation.log"
    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    file_handler.setFormatter(formatter)
    logging.getLogger().addHandler(file_handler)
    logger.info(f"Logging initialized. Log file: {log_file}")
    
    # FASTA file paths
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
    
    processor = BatchEntrapmentFDPProcessor(
        base_dir=str(output_dir),
        species_fasta_dict=species_fasta_dict
    )
    res_spec_lib_path_dict = {
        DIA_NN_result_dir / "original" / "report.parquet": spec_lib_dir / "spec_library_original.parquet",
        DIA_NN_result_dir / "pretrained_ratio_0.1" / "report.parquet": spec_lib_dir / "pretrained_filtered" / "spec_library_filtered_by_ratio_0.1.parquet",
        DIA_NN_result_dir / "pretrained_ratio_0.2" / "report.parquet": spec_lib_dir / "pretrained_filtered" / "spec_library_filtered_by_ratio_0.2.parquet",
        DIA_NN_result_dir / "pretrained_ratio_0.3" / "report.parquet": spec_lib_dir / "pretrained_filtered" / "spec_library_filtered_by_ratio_0.3.parquet",
        DIA_NN_result_dir / "pretrained_ratio_0.4" / "report.parquet": spec_lib_dir / "pretrained_filtered" / "spec_library_filtered_by_ratio_0.4.parquet",
        DIA_NN_result_dir / "pretrained_ratio_0.5" / "report.parquet": spec_lib_dir / "pretrained_filtered" / "spec_library_filtered_by_ratio_0.5.parquet",
        DIA_NN_result_dir / "pretrained_ratio_0.6" / "report.parquet": spec_lib_dir / "pretrained_filtered" / "spec_library_filtered_by_ratio_0.6.parquet",
        DIA_NN_result_dir / "pretrained_ratio_0.7" / "report.parquet": spec_lib_dir / "pretrained_filtered" / "spec_library_filtered_by_ratio_0.7.parquet",
        DIA_NN_result_dir / "pretrained_ratio_0.8" / "report.parquet": spec_lib_dir / "pretrained_filtered" / "spec_library_filtered_by_ratio_0.8.parquet",
        DIA_NN_result_dir / "pretrained_ratio_0.9" / "report.parquet": spec_lib_dir / "pretrained_filtered" / "spec_library_filtered_by_ratio_0.9.parquet",
        DIA_NN_result_dir / "finetuned_ratio_0.1" / "report.parquet": spec_lib_dir / "finetuned_filtered" / "spec_library_filtered_by_ratio_0.1.parquet",
        DIA_NN_result_dir / "finetuned_ratio_0.2" / "report.parquet": spec_lib_dir / "finetuned_filtered" / "spec_library_filtered_by_ratio_0.2.parquet",
        DIA_NN_result_dir / "finetuned_ratio_0.3" / "report.parquet": spec_lib_dir / "finetuned_filtered" / "spec_library_filtered_by_ratio_0.3.parquet",
        DIA_NN_result_dir / "finetuned_ratio_0.4" / "report.parquet": spec_lib_dir / "finetuned_filtered" / "spec_library_filtered_by_ratio_0.4.parquet",
        DIA_NN_result_dir / "finetuned_ratio_0.5" / "report.parquet": spec_lib_dir / "finetuned_filtered" / "spec_library_filtered_by_ratio_0.5.parquet",
        DIA_NN_result_dir / "finetuned_ratio_0.6" / "report.parquet": spec_lib_dir / "finetuned_filtered" / "spec_library_filtered_by_ratio_0.6.parquet",
        DIA_NN_result_dir / "finetuned_ratio_0.7" / "report.parquet": spec_lib_dir / "finetuned_filtered" / "spec_library_filtered_by_ratio_0.7.parquet",
        DIA_NN_result_dir / "finetuned_ratio_0.8" / "report.parquet": spec_lib_dir / "finetuned_filtered" / "spec_library_filtered_by_ratio_0.8.parquet",
        DIA_NN_result_dir / "finetuned_ratio_0.9" / "report.parquet": spec_lib_dir / "finetuned_filtered" / "spec_library_filtered_by_ratio_0.9.parquet"
    }
    
    # Process all folders containing result files
    logger.info("Starting batch processing of all result folders...")
    results_df = processor.process_all_results_spec_libs(res_spec_lib_path_dict)

    # Display results
    print("\n=== Entrapment FDP Calculation Results ===")
    print(results_df.to_string(index=False))
    
    # Save results to CSV
    processor.save_results_to_csv(results_df)
    
    # Display statistics
    if 'entrapment_FDP' in results_df.columns:
        successful_results = results_df.dropna(subset=['entrapment_FDP'])
        if not successful_results.empty:
            print(f"\n=== Statistics ===")
            print(f"Number of successfully processed folders: {len(successful_results)}")
            print(f"FDP range: {successful_results['entrapment_FDP'].min():.6f} - {successful_results['entrapment_FDP'].max():.6f}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate Entrapment FDP")
    parser.add_argument("--fasta_dir", type=Path, required=True, help="Path to FASTA directory")
    parser.add_argument("--DIA_NN_result_dir", type=Path, required=True, help="Path to DIA-NN result directory")
    parser.add_argument("--spec_lib_dir", type=Path, required=True, help="Path to spectral library directory")
    parser.add_argument("--output_dir", type=Path, required=True, help="Path to output directory")
    
    args = parser.parse_args()

    main(
        fasta_dir=args.fasta_dir,
        DIA_NN_result_dir=args.DIA_NN_result_dir,
        spec_lib_dir=args.spec_lib_dir,
        output_dir=args.output_dir
    )
