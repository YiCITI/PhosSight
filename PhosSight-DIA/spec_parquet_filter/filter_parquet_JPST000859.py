
from pathlib import Path
import os

# TODO: delete these comments
# spec_lib_dir = Path("~/PhosSight_analysis/spec_lib/JPST000859").expanduser()
# fasta_dir = Path("~/PhosSight_analysis/data/derivedfasta/JPST000859").expanduser()

def step1_generate_original_parquet(spec_lib_dir: Path):
    """Generate the original parquet file from the input data.

    Args:
        spec_lib_dir: Path to the spectral library directory
    """
    from filter_parquet.filter_parquet_using_pep_list import filter_syn_pep_in_parquet_by_peptide_list

    txt_file_syn_pep = Path(spec_lib_dir) / "syn_pep_STY.txt"
    parquet_in = Path(spec_lib_dir) / "spectral-library-all.parquet"
    parquet_out_original = Path(spec_lib_dir) / "spec_library_original.parquet"
    filter_syn_pep_in_parquet_by_peptide_list(
        txt_syn_pep_file_path=txt_file_syn_pep,
        input_parquet_path=parquet_in,
        output_parquet_path=parquet_out_original
    )
    

def step2_generate_filtered_parquet_using_pretrained_model(spec_lib_dir: Path, fasta_dir: Path):
    """Generate filtered parquet files using pretrained model.

    Args:
        spec_lib_dir: Path to the spectral library directory
        fasta_dir: Path to the fasta/peptide lists directory
    """

    from filter_parquet.filter_parquet_using_pep_list import filter_parquet_by_peptide_list

    ratio_suffix = [f'ratio_{i}' for i in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]]
    score_suffix = [f'score_{i}' for i in [0.55, 0.6, 0.65]]
    output_dir = Path(spec_lib_dir) / "pretrained_filtered"

    os.makedirs(output_dir, exist_ok=True)
    suffix = ratio_suffix + score_suffix

    for s in suffix:
        txt_file = Path(fasta_dir) / "for_filter_spec_lib_parquet_pretrained" / f"filtered_peptides_pretrained_{s}.txt"
        parquet_in = Path(spec_lib_dir) / "spectral-library-all.parquet"
        parquet_out = output_dir / f"spec_library_filtered_by_{s}.parquet"
        filter_parquet_by_peptide_list(
            txt_file,
            parquet_in,
            parquet_out
        )


def step3_generate_filtered_parquet_using_finetuned_model(spec_lib_dir: Path, fasta_dir: Path):
    """Generate filtered parquet files using fine-tuned model.

    Args:
        spec_lib_dir: Path to the spectral library directory
        fasta_dir: Path to the fasta/peptide lists directory
    """
    from filter_parquet.filter_parquet_using_pep_list import filter_parquet_by_peptide_list

    ratio_suffix = [f'ratio_{i}' for i in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]]
    score_suffix = [f'score_{i}' for i in [0.1, 0.2, 0.5, 0.7, 0.9]]
    output_dir = Path(spec_lib_dir) / "finetuned_filtered"

    os.makedirs(output_dir, exist_ok=True)
    suffix = ratio_suffix + score_suffix

    for s in suffix:
        txt_file = Path(fasta_dir) / "for_filter_spec_lib_parquet_finetuned" / f"filtered_peptides_finetuned_{s}.txt"
        parquet_in = Path(spec_lib_dir) / "spectral-library-all.parquet"
        parquet_out = output_dir / f"spec_library_filtered_by_{s}.parquet"
        filter_parquet_by_peptide_list(
            txt_file,
            parquet_in,
            parquet_out
        )

if __name__ == "__main__":
    import sys
    import argparse

    parser = argparse.ArgumentParser(description='Process parquet files with different steps')
    parser.add_argument('--step', type=int, choices=[1, 2, 3], required=True,
                       help='Step to run: 1=generate original parquet, 2=filter using pretrained model, 3=filter using finetuned model')
    parser.add_argument('--spec-lib-dir', type=str, required=True,
                        help='Path to spec lib dir (required)')
    parser.add_argument('--fasta-dir', type=str, required=False,
                        help='Path to fasta dir (required for steps 2 and 3)')

    args = parser.parse_args()

    # Convert to Path and expand ~. fasta_dir may be optional for step 1.
    spec_lib_dir = Path(args.spec_lib_dir).expanduser()
    fasta_dir = Path(args.fasta_dir).expanduser() if args.fasta_dir else None
    print(f"Using spec_lib_dir={spec_lib_dir}")
    if fasta_dir:
        print(f"Using fasta_dir={fasta_dir}")
    else:
        print("No fasta_dir provided (this is OK for step 1)")

    # Validate required inputs for each step
    if args.step in (2, 3) and not fasta_dir:
        print("Error: --fasta-dir is required for steps 2 and 3.")
        sys.exit(2)

    if args.step == 1:
        print("Running Step 1: Generate original parquet file")
        step1_generate_original_parquet(spec_lib_dir)
        print("Step 1 completed.")
    elif args.step == 2:
        print("Running Step 2: Generate filtered parquet using pretrained model")
        step2_generate_filtered_parquet_using_pretrained_model(spec_lib_dir, fasta_dir)
        print("Step 2 completed.")
    elif args.step == 3:
        print("Running Step 3: Generate filtered parquet using finetuned model")
        step3_generate_filtered_parquet_using_finetuned_model(spec_lib_dir, fasta_dir)
        print("Step 3 completed.")
