#!/bin/bash

export PYTHONUTF8=1
export PYTHONIOENCODING=utf-8
export LANG=en_US.UTF-8
export LC_ALL=en_US.UTF-8
CONDA_NO_PLUGINS=true

# Set pathes

anacondaPath=~/miniconda3
syn_spec_lib_dir=~/PhosSight_analysis/spec_lib/JPST000859
syn_fasta_dir=~/PhosSight_analysis/data/derivedfasta/JPST000859
A549_spec_lib_dir=~/PhosSight_analysis/spec_lib/202503_A549_0h_24h
A549_fasta_dir=~/PhosSight_analysis/data/derivedfasta/202503_A549_0h_24h

# Create environment and install dependencies
source $anacondaPath/etc/profile.d/conda.sh
conda create -n PhosSight_DIA python=3.13.2 -y
conda activate PhosSight_DIA
pip install -r ../Install/requirements.txt

# Predict theoretical spectra using DIA-NN's built-in theoretical spectrum prediction function

# Generate peptide list from fasta for PhosDetect to predict detectability

# Filter (or check) the original spectral library

# Run small-scale library search experiments to generate training data for fine-tuning

# Score detectability of theoretical spectra using PhosDetect

# Filter the spectral library to generate spectral library for DIA analysis

# Run DIA-NN spectral library search (spectral libraries filtered by pretrained model and fine-tuned model)

# Analyze DIA-NN results to evaluate the effectiveness of PhosDetect

conda deactivate