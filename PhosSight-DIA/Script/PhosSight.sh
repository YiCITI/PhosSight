#!/bin/bash
set -e

export PYTHONUTF8=1
export PYTHONIOENCODING=utf-8
export LANG=en_US.UTF-8
export LC_ALL=en_US.UTF-8
CONDA_NO_PLUGINS=true

if ! command -v singularity &> /dev/null; then
    echo "Singularity not found. Please install or load Singularity."
    exit 1
fi

# Set pathes
anacondaPath=~/miniconda3
diann_singularity_img_path=~/diann-2.2.0/diann-2.2.0.img
diann_executable_path=~/diann-2.2.0/diann-linux
syn_raw_dir=~/PhosSight_analysis/temp/dataset/syn
syn_spec_lib_dir=~/PhosSight_analysis/temp/spec_lib/syn
syn_fasta_dir=~/PhosSight_analysis/temp/database/syn
syn_result_dir=~/PhosSight_analysis/temp/result/syn
A549_raw_dir=~/PhosSight_analysis/temp/dataset/A549
A549_spec_lib_dir=~/PhosSight_analysis/temp/spec_lib/A549
A549_fasta_dir=~/PhosSight_analysis/temp/database/A549
A549_result_dir=~/PhosSight_analysis/temp/result/A549
PhosSight_DIA_dir=~/github_repo/PhosSight/PhosSight-DIA/
analysis_dir=~/github_repo/PhosSight/PhosSight-DIA/Script/analysis/

# Create environment and install dependencies
source $anacondaPath/etc/profile.d/conda.sh
# conda create -n PhosSight_DIA python=3.13.2 -y
conda activate PhosSight_DIA
# pip install -r $PhosSight_DIA_dir/Install/requirements.txt

# ============================================================================================
# Predict theoretical spectra using DIA-NN's built-in theoretical spectrum prediction function
# ============================================================================================
# syn
# singularity exec $diann_singularity_img_path $diann_executable_path --lib  --threads 32 --verbose 1 --out $syn_spec_lib_dir/report.parquet --qvalue 0.01 --matrices --temp $syn_spec_lib_dir --out-lib $syn_spec_lib_dir/report-lib.parquet --gen-spec-lib --predictor --fasta $syn_fasta_dir/syn_pep_strip_with_protein.fasta --fasta $syn_fasta_dir/uniprotkb_proteome_UP000002311_2025_07_07.fasta --fasta $syn_fasta_dir/uniprotkb_proteome_UP000000625_2025_09_10.fasta --fasta $syn_fasta_dir/uniprotkb_proteome_UP000008311_2025_09_18.fasta --fasta-search --min-fr-mz 110 --max-fr-mz 1600 --met-excision --min-pep-len 7 --max-pep-len 46 --min-pr-mz 400 --max-pr-mz 1250 --min-pr-charge 1 --max-pr-charge 4 --cut K*,R* --missed-cleavages 1 --unimod4 --var-mods 1 --var-mod UniMod:21,79.966331,STY --window 7 --mass-acc 8 --mass-acc-ms1 10 --peptidoforms --rt-profiling
# singularity exec $diann_singularity_img_path $diann_executable_path --lib $syn_spec_lib_dir/report-lib.predicted.speclib --threads 32 --verbose 1 --out $syn_spec_lib_dir/report.parquet --qvalue 0.01 --matrices --temp $syn_spec_lib_dir --out-lib $syn_spec_lib_dir/spectral-library-all.parquet --gen-spec-lib --unimod4 --var-mods 1 --var-mod UniMod:21,79.966331,STY --window 7 --mass-acc 8 --mass-acc-ms1 10 --peptidoforms --rt-profiling
# # A549
# singularity exec $diann_singularity_img_path $diann_executable_path --lib  --threads 32 --verbose 1 --out $A549_spec_lib_dir/report.parquet --qvalue 0.01 --matrices --temp $A549_spec_lib_dir --out-lib $A549_spec_lib_dir/report-lib.parquet --gen-spec-lib --predictor --fasta $A549_fasta_dir/uniprotkb_proteome_UP000005640_2025_09_25.fasta --fasta-search --min-fr-mz 200 --max-fr-mz 1800 --met-excision --min-pep-len 7 --max-pep-len 46 --min-pr-mz 350 --max-pr-mz 1650 --min-pr-charge 1 --max-pr-charge 4 --cut K*,R* --missed-cleavages 1 --unimod4 --var-mods 1 --var-mod UniMod:21,79.966331,STY --window 10 --mass-acc 10 --mass-acc-ms1 4 --peptidoforms --rt-profiling
# singularity exec $diann_singularity_img_path $diann_executable_path --lib $A549_spec_lib_dir/report-lib.predicted.speclib --threads 32 --verbose 1 --out $A549_spec_lib_dir/report.parquet --qvalue 0.01 --matrices --temp $A549_spec_lib_dir --out-lib $A549_spec_lib_dir/spectral-library-all.parquet --gen-spec-lib --unimod4 --var-mods 1 --var-mod UniMod:21,79.966331,STY --window 10 --mass-acc 10 --mass-acc-ms1 4 --peptidoforms --rt-profiling

# # ============================================================================================
# # Generate peptide list from fasta for PhosDetect to (filter original spectral library and) predict detectability
# # ============================================================================================
# python $PhosSight_DIA_dir/Script/generate_pep_fasta/generate_pep_fasta_syn_4species.py --step 1 --work_dir $syn_fasta_dir
# python $PhosSight_DIA_dir/Script/generate_pep_fasta/generate_pep_fasta_A549.py --step 1 --work_dir $A549_fasta_dir

# # ============================================================================================
# # Filter (or check) the original spectral library
# # ============================================================================================
# python $PhosSight_DIA_dir/Script/spec_parquet_filter/filter_parquet_syn.py --step 1 --spec-lib-dir $syn_spec_lib_dir --fasta-dir $syn_fasta_dir
# python $PhosSight_DIA_dir/Script/spec_parquet_filter/filter_parquet_A549.py --step 1 --spec-lib-dir $A549_spec_lib_dir --fasta-dir $A549_fasta_dir
# # Or just copy the original spectral library

# # ============================================================================================
# # Run DIA-NN spectral library search using original spectral libraries to generate baseline results
# # ============================================================================================
# python $PhosSight_DIA_dir/Script/run_diann_syn.py --step original --diann-cmd-prefix "singularity exec $diann_singularity_img_path $diann_executable_path" --output_dir $syn_result_dir --raw_dir $syn_raw_dir --spec_lib_dir $syn_spec_lib_dir
# python $PhosSight_DIA_dir/Script/run_diann_A549.py --step original --diann_cmd_prefix "singularity exec $diann_singularity_img_path $diann_executable_path" --output_dir $A549_result_dir --raw_dir $A549_raw_dir --spec_lib_dir $A549_spec_lib_dir

# # ==============================================================================================
# # Run small-scale library search experiments to generate training data for fine-tuning
# # ==============================================================================================
# python $PhosSight_DIA_dir/Script/run_diann_syn.py --step three_rep --diann-cmd-prefix "singularity exec $diann_singularity_img_path $diann_executable_path" --output_dir $syn_result_dir --raw_dir $syn_raw_dir --spec_lib_dir $syn_spec_lib_dir
# python $PhosSight_DIA_dir/Script/run_diann_A549.py --step two_rep --diann_cmd_prefix "singularity exec $diann_singularity_img_path $diann_executable_path" --output_dir $A549_result_dir --raw_dir $A549_raw_dir --spec_lib_dir $A549_spec_lib_dir

# # =============================================================================================
# # Score detectability of theoretical spectra using PhosDetect
# # =============================================================================================


# # =============================================================================================
# # Filter the spectral library to generate spectral library for DIA analysis
# # =============================================================================================
# # For pretrained model
# python $PhosSight_DIA_dir/Script/generate_pep_fasta/generate_pep_fasta_syn_4species.py --step 2 --work_dir $syn_fasta_dir
# python $PhosSight_DIA_dir/Script/generate_pep_fasta/generate_pep_fasta_A549.py --step 2 --work_dir $A549_fasta_dir
# python $PhosSight_DIA_dir/Script/spec_parquet_filter/filter_parquet_syn.py --step 2 --spec-lib-dir $syn_spec_lib_dir --fasta-dir $syn_fasta_dir
# python $PhosSight_DIA_dir/Script/spec_parquet_filter/filter_parquet_A549.py --step 2 --spec-lib-dir $A549_spec_lib_dir --fasta-dir $A549_fasta_dir
# # For fine-tuned model
# python $PhosSight_DIA_dir/Script/generate_pep_fasta/generate_pep_fasta_syn_4species.py --step 3 --work_dir $syn_fasta_dir
# python $PhosSight_DIA_dir/Script/generate_pep_fasta/generate_pep_fasta_A549.py --step 3 --work_dir $A549_fasta_dir
# python $PhosSight_DIA_dir/Script/spec_parquet_filter/filter_parquet_syn.py --step 3 --spec-lib-dir $syn_spec_lib_dir --fasta-dir $syn_fasta_dir
# python $PhosSight_DIA_dir/Script/spec_parquet_filter/filter_parquet_A549.py --step 3 --spec-lib-dir $A549_spec_lib_dir --fasta-dir $A549_fasta_dir

# # =============================================================================================
# # Run DIA-NN spectral library search (spectral libraries filtered by pretrained model and fine-tuned model)
# # =============================================================================================
# python $PhosSight_DIA_dir/Script/run_diann_syn.py --step pretrained --diann-cmd-prefix "singularity exec $diann_singularity_img_path $diann_executable_path" --output_dir $syn_result_dir --raw_dir $syn_raw_dir --spec_lib_dir $syn_spec_lib_dir
# python $PhosSight_DIA_dir/Script/run_diann_A549.py --step pretrained --diann_cmd_prefix "singularity exec $diann_singularity_img_path $diann_executable_path" --output_dir $A549_result_dir --raw_dir $A549_raw_dir --spec_lib_dir $A549_spec_lib_dir
# python $PhosSight_DIA_dir/Script/run_diann_syn.py --step finetuned --diann-cmd-prefix "singularity exec $diann_singularity_img_path $diann_executable_path" --output_dir $syn_result_dir --raw_dir $syn_raw_dir --spec_lib_dir $syn_spec_lib_dir
# python $PhosSight_DIA_dir/Script/run_diann_A549.py --step finetuned --diann_cmd_prefix "singularity exec $diann_singularity_img_path $diann_executable_path" --output_dir $A549_result_dir --raw_dir $A549_raw_dir --spec_lib_dir $A549_spec_lib_dir

# # =============================================================================================
# # Analyze DIA-NN results to evaluate the effectiveness of PhosDetect
# # =============================================================================================

# pip install -e $analysis_dir
# pip install matplotlib==3.10.6
# pip install scipy==1.16.2
# pip install matplotlib-venn==1.1.2
# pip install seaborn==0.13.2

# # Entrapment FDR calculation and plotting for syn dataset
# python $analysis_dir/calculate_draw_FDR/entrapment_fdr_calculator.py --fasta_dir $syn_fasta_dir --DIA_NN_result_dir $syn_result_dir --spec_lib_dir $syn_spec_lib_dir --output_dir $analysis_dir/output/
# python $analysis_dir/calculate_draw_FDR/draw.py --FDP_res_path $analysis_dir/output/entrapment_FDP_results.csv --output_dir $analysis_dir/output/

# Upset plot for identification overlap in syn dataset
# python $analysis_dir/draw_ident_upset_plot/syn_yeast.py --res_dir $syn_result_dir --fasta_dir $syn_fasta_dir --output_dir $analysis_dir/output/
# # TODO: Add R script from scy to generate upset plot

# Upset plot for identification overlap in A549 dataset
# python $analysis_dir/draw_ident_upset_plot/A549.py --res_dir $A549_result_dir --output_dir $analysis_dir/output/
# # TODO: Add R script from scy to generate upset plot

# Venn diagram and violin plot for identification overlap and delta RT in A549 dataset
# python $analysis_dir/draw_ident_venn_violin/A549.py --res_dir $A549_result_dir --detect_path $A549_fasta_dir/peptide_scores_finetuned.txt --output_dir $analysis_dir/output/

# Draw run time comparison for syn dataset
# python $analysis_dir/draw_run_time/draw.py --DIA_NN_result_dir $syn_result_dir --output_path $analysis_dir/output/syn_run_time_comparison.svg
# # Draw run time comparison for A549 dataset
# python $analysis_dir/draw_run_time/draw.py --DIA_NN_result_dir $A549_result_dir --output_path $analysis_dir/output/A549_run_time_comparison.svg

# Draw lines for modified synthetic peptide identifications
python $analysis_dir/draw_site_ident_line/draw_count_modpep_ident.py --res_dir $syn_result_dir --output_dir $analysis_dir/output/

conda deactivate