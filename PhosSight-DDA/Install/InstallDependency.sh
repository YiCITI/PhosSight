#!/bin/bash

echo "Installing..."

##Parameter path
#PhosSightPath="E:/Project/PhosSight/Github/PhosSight"
#anacondaPath="/C/ProgramData/anaconda3"
#scriptPath="$PhosSightPath/Script"


if [ -n "$1" ]; then
    PhosSightPath="$1"
    echo "PhosSightPath: $PhosSightPath"
else
    echo "Please give PhosSightPath"
    exit 1
fi

if [ -n "$2" ]; then
    anacondaPath="$2"
    echo "anacondaPath: $anacondaPath"
else
    echo "Please give anacondaPath"
    exit 2
fi

scriptPath="$PhosSightPath/Script"

# Prepare conda for non-interactive usage and prefer conda-forge to avoid Anaconda TOS prompts
if [ -f "$anacondaPath/etc/profile.d/conda.sh" ]; then
    # shellcheck source=/dev/null
    source "$anacondaPath/etc/profile.d/conda.sh"
    conda config --system --set always_yes yes || true
    conda config --system --set channel_priority strict || true
fi

#====================================Download==========================================#
##1. AutoRT

if [ ! -d "$PhosSightPath/Script/AutoRT" ]; then
    git clone https://github.com/bzhanglab/AutoRT.git
    mv AutoRT $PhosSightPath/Script
fi

##2. pDeep3
if [ ! -d "$PhosSightPath/Script/pDeep3" ]; then
    git clone https://github.com/pFindStudio/pDeep3.git
    mv pDeep3 $PhosSightPath/Script/pDeep3
fi

##3. PhosphoRS
if [ ! -d "$PhosSightPath/Script/PhosphoRS" ]; then
    curl -o phosphoRS-cli.zip -LJ https://github.com/lmsac/phosphoRS-cli/releases/download/v1.0.0/phosphoRS-cli.zip
    unzip phosphoRS-cli.zip -d phosphoRS-cli
    mv phosphoRS-cli $PhosSightPath/Script/PhosphoRS
    rm -f phosphoRS-cli.zip
fi

##4. SpectralEntropy
if [ ! -d "$PhosSightPath/Script/SpectralEntropy" ]; then
    git clone https://github.com/YuanyueLi/SpectralEntropy.git
    mv SpectralEntropy $PhosSightPath/Script/SpectralEntropy
fi

#=====================================Install==========================================#
##1.AutoRT

source $anacondaPath/etc/profile.d/conda.sh
autort_env="AutoRT"

if conda env list | grep -q $autort_env; then
    echo "AutoRT environment already exists."
else
    echo "Creating AutoRT environment..."
    conda create -n AutoRT --override-channels -c conda-forge python=3.8 -y
    echo "AutoRT environment created and activated."
fi

conda activate AutoRT
python --version

# conda install --override-channels -c conda-forge tensorflow -y
#git clone https://github.com/bzhanglab/AutoRT
pip install -r $scriptPath/AutoRT/requirements.txt

conda deactivate

##2.pDeep3

source $anacondaPath/etc/profile.d/conda.sh
pdeep3_env="pDeep3"

if conda env list | grep -q $pdeep3_env; then
    echo "pDeep3 environment already exists."
else
    echo "Creating pDeep3 environment..."
    conda create -n pDeep3 --override-channels -c conda-forge python=3.6 -y

    echo "pDeep3 environment created and activated."
fi

conda activate pDeep3
python --version

# Use pip for legacy TF 1.13.1 to avoid Anaconda channels
pip install tensorflow==1.13.1
pip install $scriptPath/pDeep3/pDeep3/.
pip install -e $scriptPath/pDeep3/pDeep3/.
pip install Cython

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

##3.R environment

source $anacondaPath/etc/profile.d/conda.sh
r_env="R_env"

if conda env list | grep -q $r_env; then
    echo "R_env environment already exists."
else
    echo "Creating R_env environment..."
    conda env create -f $PhosSightPath/Install/environment_R.yml --prefix $anacondaPath/envs/R_env

    echo "R_env environment created and activated."
fi

