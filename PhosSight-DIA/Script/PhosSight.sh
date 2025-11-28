#!/bin/bash

export PYTHONUTF8=1
export PYTHONIOENCODING=utf-8
export LANG=en_US.UTF-8
export LC_ALL=en_US.UTF-8
CONDA_NO_PLUGINS=true

# Set pathes

anacondaPath=~/miniconda3
syn_spec_lib_dir=~/PhosSight_analysis/temp/spec_lib/syn
syn_fasta_dir=~/PhosSight_analysis/temp/data/derivedfasta/syn
A549_spec_lib_dir=~/PhosSight_analysis/temp/spec_lib/A549
A549_fasta_dir=~/PhosSight_analysis/temp/data/derivedfasta/A549

# Create environment and install dependencies
source $anacondaPath/etc/profile.d/conda.sh
conda create -n PhosSight_DIA python=3.13.2 -y
conda activate PhosSight_DIA
pip install -r ../Install/requirements.txt

# Predict theoretical spectra using DIA-NN's built-in theoretical spectrum prediction function

# Generate peptide list from fasta for PhosDetect to predict detectability

# Filter (or check) the original spectral library
python ../spec_parquet_filter/filter_parquet_syn.py --step 1 --spec-lib-dir $syn_spec_lib_dir --fasta-dir $syn_fasta_dir
python ../spec_parquet_filter/filter_parquet_A549.py --step 1 --spec-lib-dir $A549_spec_lib_dir --fasta-dir $A549_fasta_dir

# Run small-scale library search experiments to generate training data for fine-tuning

# Score detectability of theoretical spectra using PhosDetect

# Filter the spectral library to generate spectral library for DIA analysis

# Run DIA-NN spectral library search (spectral libraries filtered by pretrained model and fine-tuned model)

# Analyze DIA-NN results to evaluate the effectiveness of PhosDetect

conda deactivate