#!/bin/bash

echo "Installing..."

##Parameter path
#DeepRescore2Path="E:/Project/DeepRescore2/Github/DeepRescore2"
#anacondaPath="/C/ProgramData/anaconda3"
#scriptPath="$DeepRescore2Path/Script"


if [ -n "$1" ]; then
    DeepRescore2Path="$1"
    echo "DeepRescore2Path: $DeepRescore2Path"
else
    echo "Please give DeepRescore2Path"
    exit 1
fi

if [ -n "$2" ]; then
    anacondaPath="$2"
    echo "anacondaPath: $anacondaPath"
else
    echo "Please give anacondaPath"
    exit 2
fi

scriptPath="$DeepRescore2Path/Script"
installPath="$DeepRescore2Path/Install"

require_yaml() {
    yaml_path="$1"
    env_name="$2"

    if [ ! -f "$yaml_path" ]; then
        echo "YAML file for $env_name not found at $yaml_path"
        exit 10
    fi
}

env_exists() {
    conda env list | awk '{print $1}' | grep -Fxq "$1"
}

# Prepare conda for non-interactive usage and prefer conda-forge to avoid Anaconda TOS prompts
if [ -f "$anacondaPath/etc/profile.d/conda.sh" ]; then
    # shellcheck source=/dev/null
    source "$anacondaPath/etc/profile.d/conda.sh"
    conda config --system --set always_yes yes || true
    conda config --system --set channel_priority strict || true
fi

#====================================Download==========================================#
##1. AutoRT

if [ ! -d "$DeepRescore2Path/Script/AutoRT" ]; then
    git clone https://github.com/bzhanglab/AutoRT.git
    mv AutoRT $DeepRescore2Path/Script
fi

##2. pDeep3
if [ ! -d "$DeepRescore2Path/Script/pDeep3" ]; then
    git clone https://github.com/pFindStudio/pDeep3.git
    mv pDeep3 $DeepRescore2Path/Script/pDeep3
fi

##3. PhosphoRS
if [ ! -d "$DeepRescore2Path/Script/PhosphoRS" ]; then
    curl -o phosphoRS-cli.zip -LJ https://github.com/lmsac/phosphoRS-cli/releases/download/v1.0.0/phosphoRS-cli.zip
    unzip phosphoRS-cli.zip -d phosphoRS-cli
    mv phosphoRS-cli $DeepRescore2Path/Script/PhosphoRS
    rm -f phosphoRS-cli.zip
fi

##4. SpectralEntropy
if [ ! -d "$DeepRescore2Path/Script/SpectralEntropy" ]; then
    git clone https://github.com/YuanyueLi/SpectralEntropy.git
    mv SpectralEntropy $DeepRescore2Path/Script/SpectralEntropy
fi

#=====================================Install==========================================#
##1.AutoRT

source $anacondaPath/etc/profile.d/conda.sh
autort_env="AutoRT"
autort_yaml="$installPath/AutoRT.yml"
require_yaml "$autort_yaml" "$autort_env"

if env_exists "$autort_env"; then
    echo "AutoRT environment already exists, updating from $autort_yaml."
    conda env update -n "$autort_env" -f "$autort_yaml" --prune
else
    echo "Creating AutoRT environment from $autort_yaml..."
    conda env create -f "$autort_yaml"
    echo "AutoRT environment created."
fi

conda activate AutoRT
python --version
conda deactivate

##2.pDeep3

source $anacondaPath/etc/profile.d/conda.sh
pdeep3_env="pDeep3"
pdeep3_yaml="$installPath/pDeep3.yml"
require_yaml "$pdeep3_yaml" "$pdeep3_env"

if env_exists "$pdeep3_env"; then
    echo "pDeep3 environment already exists, updating from $pdeep3_yaml."
    conda env update -n "$pdeep3_env" -f "$pdeep3_yaml" --prune
else
    echo "Creating pDeep3 environment from $pdeep3_yaml..."
    conda env create -f "$pdeep3_yaml"
    echo "pDeep3 environment created."
fi

conda activate pDeep3
python --version
pip install -e $DeepRescore2Path/Script/pDeep3/pDeep3
conda deactivate

##3.Phossight

source $anacondaPath/etc/profile.d/conda.sh
phossight_env="phossight"
phossight_yaml="$installPath/phossight.yml"
require_yaml "$phossight_yaml" "$phossight_env"

if env_exists "$phossight_env"; then
    echo "phossight environment already exists, updating from $phossight_yaml."
    conda env update -n "$phossight_env" -f "$phossight_yaml" --prune
else
    echo "Creating phossight environment from $phossight_yaml..."
    conda env create -f "$phossight_yaml"
    echo "phossight environment created."
fi

conda activate phossight
python --version
conda deactivate

sourceDir="$scriptPath/pDeep3/SpectralEntropyScripts/"
destinationDir="$scriptPath/pDeep3/SpectralEntropy/"

if [ -d "$sourceDir" ] && [ -d "$destinationDir" ]; then
    for file in "$sourceDir"*
    do
        if [ -f "$file" ]; then
            baseName="$(basename $file)"
            destinationPath="$destinationDir$baseName"
            
            if [ ! -f "$destinationPath" ]; then
                cp "$file" "$destinationPath"
                echo "Copied $baseName to $destinationDir"
            else
                echo "File $baseName already exists in $destinationDir"
            fi
        fi
    done
else
    echo "Source and/or destination directory does not exist."
fi

if [ ! -f "$scriptPath/pDeep3/pDeep3/pDeep/cmd/tune_and_predict.py" ]; then
    sourceDir="$scriptPath/pDeep3/pDeep3Scripts/tune_and_predict.py"
    destinationDir="$scriptPath/pDeep3/pDeep3/pDeep/cmd/"
    cp "$sourceDir" "$destinationDir"
    echo "Copied tune_and_predict.py to $destinationDir"
else
    echo "tune_and_predict.py already exists in $destinationDir"
fi

if [ ! -f "$scriptPath/pDeep3/pDeep3/pDeep/load_data.py" ]; then
    sourceDir="$scriptPath/pDeep3/pDeep3Scripts/load_data.py"
    destinationDir="$scriptPath/pDeep3/pDeep3/pDeep/"
    cp "$sourceDir" "$destinationDir"
    echo "Copied load_data.py to $destinationDir"
else
    echo "load_data.py already exists in $destinationDir"
fi

##4.R environment

source $anacondaPath/etc/profile.d/conda.sh
r_env="R_env"
r_yaml="$installPath/R_env.yml"
require_yaml "$r_yaml" "$r_env"

if env_exists "$r_env"; then
    echo "R_env environment already exists, updating from $r_yaml."
    conda env update -n "$r_env" -f "$r_yaml" --prune
else
    echo "Creating R_env environment from $r_yaml..."
    conda env create -f "$r_yaml"
    echo "R_env environment created."
fi

