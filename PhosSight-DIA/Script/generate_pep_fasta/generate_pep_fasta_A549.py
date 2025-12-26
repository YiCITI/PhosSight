# Generate fasta databases for Zhu lab's human phosphopeptides

from process_fasta.generate_samples import generate_modified_pep_fasta_file
from process_fasta.txt2fasta import txt_to_fasta
from process_fasta.combine_fasta import combine_fasta
from process_fasta.filter_pep_fasta import filter_peptides_by_ratio, filter_peptides_above_thres
import argparse
from os import path, makedirs

# TODO: Update work_dir to the correct path before running
# work_dir = path.join(path.dirname(path.abspath(__file__)), '202503_A549_0h_24h')


def step1_prepare_all_peptides(work_dir):
    """
    Step 1: Prepare all possible peptides
    Generate peptide files from human fasta
    """
    print("=== Step 1: Prepare all possible peptides ===")
    
    print("Generating modified peptide FASTA file...")
    generate_modified_pep_fasta_file(
        file_path=path.join(work_dir, 'uniprotkb_proteome_UP000005640_2025_09_25.fasta'),
        output_file=path.join(work_dir, 'human_all_peptides_2_7_46.fasta'),
        miss_cleavage=2,
        min_length=7,
        max_length=46,
        is_replace_STY_to_BJX=False
    )
    
    print("Step 1 completed! Please check the generated files, then run step 2.")


def step2_filter_by_pretrained_model(work_dir):
    """
    Step 2: Pre-trained model peptide filtering
    Filter peptides based on pre-trained model scores
    """
    print("=== Step 2: Pre-trained model peptide filtering ===")
    
    ratio_set = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    if not path.exists(path.join(work_dir, 'for_filter_spec_lib_parquet_pretrained')):
        makedirs(path.join(work_dir, 'for_filter_spec_lib_parquet_pretrained'))
    for ratio in ratio_set:
        print(f"Filtering database peptides based on the provided list with ratio {ratio}...")
        filter_peptides_by_ratio(
            peptides_list_file_for_filtering=path.join(work_dir, 'peptide_scores_pretrained.txt'),
            all_peptides_fasta_file=path.join(work_dir, 'human_all_peptides_2_7_46.fasta'),
            output_txt_file=path.join(work_dir, 'for_filter_spec_lib_parquet_pretrained', f'filtered_peptides_pretrained_ratio_{ratio}.txt'),
            ratio_to_keep=ratio
        )
    
    print("Step 2 completed! Please check the generated pre-trained model filtering results, then run step 3.")


def step3_filter_by_finetuned_model(work_dir):
    """
    Step 3: Fine-tuned model peptide filtering
    Filter peptides based on fine-tuned model scores
    """
    print("=== Step 3: Fine-tuned model peptide filtering ===")
    
    ratio_set = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    if not path.exists(path.join(work_dir, 'for_filter_spec_lib_parquet_finetuned')):
        makedirs(path.join(work_dir, 'for_filter_spec_lib_parquet_finetuned'))
    for ratio in ratio_set:
        print(f"Filtering database peptides based on the provided list with ratio {ratio}...")
        filter_peptides_by_ratio(
            peptides_list_file_for_filtering=path.join(work_dir, 'peptide_scores_finetuned.txt'),
            all_peptides_fasta_file=path.join(work_dir, 'human_all_peptides_2_7_46.fasta'),
            output_txt_file=path.join(work_dir, 'for_filter_spec_lib_parquet_finetuned', f'filtered_peptides_finetuned_ratio_{ratio}.txt'),
            ratio_to_keep=ratio
        )
    
    print("Step 3 completed! All steps have been completed.")


def main():
    # """
    # 主函数 - 选择要运行的步骤
    # """
    # print("JPST000859 Database Generation Tool")
    # print("Available steps:")
    # print("1. Prepare all possible peptides (generate and combine from proteomes)")
    # print("2. Pre-trained model peptide filtering")
    # print("3. Fine-tuned model peptide filtering")
    # print("\nPlease call the corresponding functions as needed:")
    # print("step1_prepare_all_peptides()")
    # print("step2_filter_pretrained_model()")
    # print("step3_filter_finetuned_model()")
    # print("\nOr uncomment the code below to run specific steps.")
    
    parser = argparse.ArgumentParser(description="JPST000859 Database Generation Tool")
    parser.add_argument('--step', type=int, choices=[1, 2, 3], required=True,
                        help="Choose the step to run: 1=prepare peptides, 2=pretrained model filtering, 3=fine-tuned model filtering")
    parser.add_argument('--work_dir', type=str, required=True,
                        help="Path to the working directory")
    args = parser.parse_args()

    if args.step == 1:
        step1_prepare_all_peptides(args.work_dir)
    elif args.step == 2:
        step2_filter_by_pretrained_model(args.work_dir)
    elif args.step == 3:
        step3_filter_by_finetuned_model(args.work_dir)

if __name__ == "__main__":
    main()