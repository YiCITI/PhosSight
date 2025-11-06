#!/bin/bash

export PYTHONUTF8=1
export PYTHONIOENCODING=utf-8
export LANG=en_US.UTF-8
export LC_ALL=en_US.UTF-8
CONDA_NO_PLUGINS=true


if [ -n "$1" ]; then
    param_path="$1"
else
    echo "Please give parameter path"
    exit 1
fi

#=======================================Step0: Preparation===========================#

echo "Step0: Preparation"

##Read parameter file
###Parameters used for identification

while IFS=$'\t' read -r Name Value || [ -n "$Name" ]; do
  
    # echo "Name:$Name"
    # echo "Value:$Value"
    eval "$Name=\"$Value\""

done < "$param_path"

echo "PhosSightPath: $PhosSightPath"
echo "anacondaPath: $anacondaPath"
echo "phosSightNegativePath: $phosSightNegativePath"
echo "phosSightPretrainedModelPath: $phosSightPretrainedModelPath"

# Judge input data exist or not
missing_dirs=""

if [ ! -d "$PhosSightPath/Script" ]; then
    missing_dirs+="Script, "
fi

if [ ! -d "$PhosSightPath/Script/AutoRT" ]; then
    missing_dirs+="AutoRT, "
fi

if [ ! -d "$PhosSightPath/Script/DeepRelocalization" ]; then
    missing_dirs+="DeepRelocalization, "
fi

if [ ! -d "$PhosSightPath/Script/Features" ]; then
    missing_dirs+="Features, "
fi

if [ ! -d "$PhosSightPath/Script/generate_train_prediction" ]; then
    missing_dirs+="generate_train_prediction, "
fi

if [ ! -d "$PhosSightPath/Script/GenerateFeatureMatrix" ]; then
    missing_dirs+="GenerateFeatureMatrix, "
fi

if [ ! -d "$PhosSightPath/Script/KinaseActivityScoreInference" ]; then
    missing_dirs+="KinaseActivityScoreInference, "
fi

if [ ! -d "$PhosSightPath/Script/pDeep3/pDeep3" ]; then
    missing_dirs+="pDeep3_pDeep3, "
fi

if [ ! -d "$PhosSightPath/Script/pDeep3/SpectralEntropy" ]; then
    missing_dirs+="SpectralEntropy, "
fi

if [ ! -d "$PhosSightPath/Script/pDeep3" ]; then
    missing_dirs+="pDeep3, "
fi

if [ ! -d "$PhosSightPath/Script/Percolator" ]; then
    missing_dirs+="Percolator, "
fi

if [ ! -d "$PhosSightPath/Script/PGA" ]; then
    missing_dirs+="PGA, "
fi

if [ ! -d "$PhosSightPath/Script/PhosphoRS" ]; then
    missing_dirs+="PhosphoRS, "
fi

if [ ! -d "$PhosSightPath/Script/PhosphoRS/phosphoRS-cli" ]; then
    missing_dirs+="phosphoRS-cli, "
fi

if [ ! -d "$PhosSightPath/Script/TMTQuantification" ]; then
    missing_dirs+="TMTQuantification, "
fi

if [ -n "$missing_dirs" ]; then
    missing_dirs="${missing_dirs%, }" 
    echo "PhosSightPath is missing the following subdirectories: $missing_dirs"
    exit 2
else
    echo "PhosSightPath contains Script and all subdirectories."
fi

if ls "$rawSpectraPath"/*.raw 1> /dev/null 2>&1; then
  echo "RAW files provided"
else
  echo "no RAW files, please provide!"
  exit 3 
fi

if ls "$spectraPath"/*.mgf 1> /dev/null 2>&1; then
  echo "MGF files provided"
else
  echo "no MGF files, please provide!"
  exit 4 
fi

if [ -f "$inputFeaturePath" ]; then
  echo "feature file provided"
else
  echo "no feature file, please provide"
  exit 5
fi


scriptPath="$PhosSightPath/Script"
phosphoRSPath="$scriptPath/PhosphoRS/phosphoRS-cli/phosphoRS.exe"
pDeep3_modelPath="$scriptPath/pDeep3/PreTrainedPhosphoModel/transfer-phos-wb-QE.ckpt"
##New generated data
###Path to the PhosphoRS results file
phosphoRSResultsPath="$outputPath/PhosphoRS" 
TXTPath="$phosphoRSResultsPath/TXT"
xmlPath="$phosphoRSResultsPath/xml"
ResultsPath="$phosphoRSResultsPath/Results"
ResultsAddIsoformSequencePath="$phosphoRSResultsPath/Results_AddIsoformSequence"
## Path to the generated feature file
featurePath="$outputPath/Features"
## Path to the PGA filtering results file
PGAPath="$outputPath/PGA"
PGA_peptide_level_Path="$PGAPath/peptide_level"
PGA_psm_level_Path="$PGAPath/psm_level"
## Path for AutoRT and pDeep3 training and prediction data
dataPath="$outputPath/generate_train_prediction"
autoRT_trainPath="$dataPath/autoRT_train"
autoRT_predictionPath="$dataPath/autoRT_prediction"
pDeep3_trainPath="$outputPath/generate_train_prediction/pDeep3_train"
pDeep3_predictionPath="$outputPath/generate_train_prediction/pDeep3_prediction"
phosSight_trainPath="$outputPath/generate_train_prediction/phosSight_train"
phosSight_predictionPath="$outputPath/generate_train_prediction/phosSight_prediction"
autoRT_resultsPath="$outputPath/autoRT_Results"
tf_modelPath="$autoRT_resultsPath/tf_model"
tf_predictionPath="$autoRT_resultsPath/tf_prediction"
pDeep3_resultsPath="$outputPath/pDeep3_Results"
pLabelPath="$pDeep3_resultsPath/pLabel"
PercolatorPath="$outputPath/Percolator"
Method1ResultsPath=$outputPath
phosSight_resultsPath="$outputPath/phosSight_Results"
pth_modelPath="$phosSight_resultsPath/pth_model"
pth_predictionPath="$phosSight_resultsPath/pth_prediction"
PhosSightResultsPath="$PercolatorPath/PhosSight"
###Modifications
if [ "$FixedMods" != "null" ]; then
  Mods="$VariableMods;$FixedMods"
else
  Mods="$VariableMods"
fi
###Judge folder exit or not if not build it
folders=($phosphoRSResultsPath $TXTPath $xmlPath $ResultsPath $ResultsAddIsoformSequencePath
         $featurePath $PGAPath $PGA_peptide_level_Path $PGA_psm_level_Path
         $dataPath $autoRT_trainPath $autoRT_predictionPath $pDeep3_trainPath $pDeep3_predictionPath
         $phosSight_trainPath $phosSight_predictionPath $autoRT_resultsPath $tf_modelPath $tf_predictionPath $pDeep3_resultsPath $pLabelPath $PercolatorPath
         $Method1ResultsPath $PhosSightResultsPath $phosSight_resultsPath $pth_modelPath $pth_predictionPath)

for folder in "${folders[@]}"; do
  if [ ! -d "$folder" ]; then
    mkdir -p "$folder"
    #echo "Folder $folder created."
  fi
done

cat "$spectraPath"/*.mgf > "$outputPath/Combined.mgf"

source $anacondaPath/etc/profile.d/conda.sh
conda activate R_env

####MaxQuant
if [ "$searchEngine" = "maxquant" ]; then
  AddModificationAdjustChargeCommand="python $scriptPath/Features/AddModifedSequenceAdjustCharge.py \
    \"$inputFeaturePath\" \
    \"$spectraPath\" \
    \"$featurePath/features.txt\" \
    \"$Mods\" \
    \"$ModsReplace\""
  eval "$AddModificationAdjustChargeCommand"
fi
####Comet, MSGF, X!Tandem
if [ "$searchEngine" = "comet" ] || [ "$searchEngine" = "msgf" ] || [ "$searchEngine" = "xtandem" ]; then
  AddModifedSequenceCommand="python $scriptPath/Features/AddModifedSequence.py \
    \"$inputFeaturePath\" \
    \"$featurePath/features.txt\" \
    \"$Mods\" \
    \"$ModsReplace\""
  eval "$AddModifedSequenceCommand"
fi

#====================================================================================#

#=======================================Step1: Phosphosite localization using PhosphoRS===========================#

echo "Step1: Phosphosite localization using PhosphoRS"
echo "Step1.1: Generate PhosphoRS input"

GeneratePhosphoRSInput_Command1="python $scriptPath/PhosphoRS/GeneratePhosphoRSCSVFile.py $featurePath/features.txt $TXTPath"

$GeneratePhosphoRSInput_Command1

GeneratePhosphoRSInput_Command2="Rscript $scriptPath/PhosphoRS/generate_phosphoRS_input_xml_folder.R \
                                \"$TXTPath\" \
                                \"$spectraPath\" \
                                \"$xmlPath\" \
                                \"$Mods\""

# eval "$GeneratePhosphoRSInput_Command2"
echo ">>> Skip Step1.1..."

echo "Step1.2: Run PhosphoRS"

files=$(ls $xmlPath)

for file in $files; do
  file_path="$xmlPath/$file"
  name="${file/.xml/.csv}"
  output_path="$ResultsPath/$name"
  if [ ! -f "$output_path" ]; then
    RunPhosphoRS_Command="$phosphoRSPath -i $file_path -o $output_path"
    $RunPhosphoRS_Command
    echo "Processing File: $name"
  else
    echo "File $name already exists in $ResultsPath, skip!"
  fi

done

echo "Step1.3: Add Isoform Sequence"

AddIsoformSequence_Command="python $scriptPath/PhosphoRS/AddIsoformSequenceForPhosphoRSResults.py \
\"$TXTPath\" \
\"$ResultsPath\" \
\"$featurePath/features.txt\" \
\"$ResultsAddIsoformSequencePath\" \
\"$featurePath/features2.txt\" \
\"$Mods\""

eval "$AddIsoformSequence_Command"
# echo ">>> Skip Step1.3..."

echo "Step1.4: Combine PhosphoRS Results"
CombinePhosphoRSResults_Command="python $scriptPath/PhosphoRS/CombinePhosphoRSResults.py \
$ResultsAddIsoformSequencePath $phosphoRSResultsPath/PhosphoRS.txt"

eval "$CombinePhosphoRSResults_Command"
# echo ">>> Skip Step1.4..."

echo "Step1.5: Add PhosphoRS Results To Features"
AddPhosphoRSToFeatures_Command="Rscript $scriptPath/PhosphoRS/combine_features_withlocalization.R \
$featurePath/features2.txt \
$phosphoRSResultsPath/PhosphoRS.txt \
$featurePath/features.PhosphoRS.txt"

eval "$AddPhosphoRSToFeatures_Command"
# echo ">>> Skip Step1.5..."

echo "Step1.6: Add PhosphoRS Probability"
AddPhosphoRSProbability_Command="python $scriptPath/PhosphoRS/GetPhosphoRSSiteProbability.py \
$featurePath/features.PhosphoRS.txt \
$featurePath/features.PhosphoRS.txt"

eval "$AddPhosphoRSProbability_Command"
# echo ">>> Skip Step1.6..."

#====================================================================================#

#=======================================Step2: Sequence quality control using PGA===========================#

echo "Step 2: Sequence quality control using PGA"
echo "Step2.1: Generate PGA input"

GeneratePGAInput_Command="Rscript $scriptPath/PGA/got_pga_input.R \
$featurePath/features.PhosphoRS.txt \
$searchEngine \
$PGAPath/pga-rawPSMs.txt"

eval "$GeneratePGAInput_Command"
# echo ">>> Skip Step2.1..."

conda deactivate

echo "Step2.2: Docker run PGA"
docker_Command="docker run --rm -v $outputPath:/opt/ -t proteomics/pga:latest Rscript"
cp "$scriptPath/PGA/calculate_fdr.R" "$PGAPath/"
cd "$outputPath"
Calculate_FDR_Command="$docker_Command ./PGA/calculate_fdr.R \
\"./PGA/\" \
\"pga\" \
\"$decoyPrefix\" \
\"FALSE\""

eval "$Calculate_FDR_Command"
# echo ">>> Skip Step2.2..."

#====================================================================================#

#=======================================Step 3: Generate train and prediction datasets===========================#

echo "Step 3: Generate train and prediction datasets"

source $anacondaPath/etc/profile.d/conda.sh
conda activate R_env

generate_train_prediction_Command="Rscript $scriptPath/generate_train_prediction/got_train_prediction.R \
\"$PGAPath/peptide_level/pga-peptideSummary.txt\" \
\"$PGAPath/psm_level/pga-peptideSummary.txt\" \
\"$featurePath/features.PhosphoRS.txt\" \
\"$autoRT_trainPath/\" \
\"$autoRT_predictionPath/\" \
\"$pDeep3_trainPath/\" \
\"$pDeep3_predictionPath/\" \
\"$PGAPath/pga-rawPSMs.txt\" \
\"$phosphoRSResultsPath/PhosphoRS.txt\" \
\"$VariableMods\" \
\"$FixedMods\" \
\"$phosSightNegativePath\" \
\"$phosSight_trainPath/\" \
\"$phosSight_predictionPath/\""

eval "$generate_train_prediction_Command"
# echo ">>> Skip Step3..."

conda deactivate

#====================================================================================#

#=======================================Step 4: RT prediction using AutoRT===========================#

echo "Step 4: RT prediction using AutoRT"

source $anacondaPath/etc/profile.d/conda.sh
conda activate AutoRT
# conda activate pDeep3

echo "Step 4.1: AutoRT Train"
autoRT_train_Command="python $scriptPath/AutoRT/autort.py train -i $autoRT_trainPath/auto_rt_train.txt -o $tf_modelPath/ -e 40 -b 64 -u m -m $scriptPath/AutoRT/models/ptm_base_model/phosphorylation_sty/model.json -rlr -n 10"

eval "$autoRT_train_Command"
# echo ">>> Skip Step4.1..."

echo "Step 4.2: AutoRT Predict Phospho"
autoRT_predict_phospho_Command="python $scriptPath/AutoRT/autort.py predict -t $autoRT_predictionPath/auto_rt_prediction.Phospho.txt -s $tf_modelPath/model.json -o $tf_predictionPath/ -p phospho.prediction"

eval "$autoRT_predict_phospho_Command"
# echo ">>> Skip Step4.2..."

echo "Step 4.3: AutoRT Predict nonPhospho"
autoRT_predict_nonPhospho_Command="python $scriptPath/AutoRT/autort.py predict -t $autoRT_predictionPath/auto_rt_prediction.nonPhospho.txt -s $tf_modelPath/model.json -o $tf_predictionPath/ -p nonPhospho.prediction"

eval "$autoRT_predict_nonPhospho_Command"
# echo ">>> Skip Step4.3..."

conda deactivate

#====================================================================================#

#=======================================Step 5: Spectrum prediction using pDeep3===========================#

echo "Step 5: Spectrum prediction using pDeep3"

source $anacondaPath/etc/profile.d/conda.sh
conda activate pDeep3

echo "Step 5.1: Generate pDeep3 parameters"
generate_pLabel_parameters_Command="python $scriptPath/pDeep3/generate_pLabel_parameters.py $pDeep3_trainPath/pdeep3_train.txt $rawSpectraPath $pLabelPath $pDeep3_resultsPath/pLabelParams.cfg $pDeep3_resultsPath/Train.Phospho.cfg $pDeep3_predictionPath/pdeep3_prediction.Phospho.txt $pDeep3_resultsPath/TrainingData.psmlabel $pDeep3_resultsPath/TrainingData.psmlabel $pDeep3_modelPath $pDeep3_resultsPath/Train.nonPhospho.cfg $pDeep3_predictionPath/pdeep3_prediction.nonPhospho.txt"

eval "$generate_pLabel_parameters_Command"
# echo ">>> Skip Step5.1..."

echo "Step 5.2: Run psmLabel"
cd "$scriptPath/pDeep3/pDeep3/pDeep/psmLabel/"
psmLabel_Command="./psmLabel.exe $pDeep3_resultsPath/pLabelParams.cfg"

eval "$psmLabel_Command"
# echo ">>> Skip Step5.2..."

echo "Step 5.3: Combine pLabel"
combine_pLabel_Command="python $scriptPath/pDeep3/CombinepLabelFiles.py $pLabelPath $pDeep3_resultsPath/TrainingData.psmlabel"

eval "$combine_pLabel_Command"
# echo ">>> Skip Step5.3..."

echo "Step 5.4: pDeep3 for Phospho"

run_Phospho_Command="python $scriptPath/pDeep3/Run/run.py $pDeep3_resultsPath/Train.Phospho.cfg $pDeep3_resultsPath/pDeep3_Predict.Phospho.txt"

eval "$run_Phospho_Command"
# echo ">>> Skip Step5.4..."

echo "Step 5.5: pDeep3 for nonPhospho"
run_nonPhospho_Command="python $scriptPath/pDeep3/Run/run.py $pDeep3_resultsPath/Train.nonPhospho.cfg $pDeep3_resultsPath/pDeep3_Predict.nonPhospho.txt"

eval "$run_nonPhospho_Command"
# echo ">>> Skip Step5.5..."

echo "Step 5.6: Build SpectralEntropy"
cd "$scriptPath/pDeep3/SpectralEntropy/"
build_SpectralEntropy_Command="python setup.py build_ext --inplace"

eval "$build_SpectralEntropy_Command"
# echo ">>> Skip Step5.6..."

echo "Step 5.7: Run SpectralEntropy/program.py - Phospho"

run_Phospho_SpectralEntropy_Command="python $scriptPath/pDeep3/SpectralEntropy/program.py $pDeep3_predictionPath/pdeep3_prediction.Phospho.txt $outputPath/Combined.mgf $pDeep3_resultsPath/pDeep3_Predict.Phospho.txt $pDeep3_resultsPath/pDeep3PredictionResults.Phospho.txt"

eval "$run_Phospho_SpectralEntropy_Command"
# echo ">>> Skip Step5.7..."

echo "Step 5.8: Run SpectralEntropy/program.py - nonPhospho"

run_nonPhospho_SpectralEntropy_Command="python $scriptPath/pDeep3/SpectralEntropy/program.py $pDeep3_predictionPath/pdeep3_prediction.nonPhospho.txt $outputPath/Combined.mgf $pDeep3_resultsPath/pDeep3_Predict.nonPhospho.txt $pDeep3_resultsPath/pDeep3PredictionResults.nonPhospho.txt"

eval "$run_nonPhospho_SpectralEntropy_Command"
# echo ">>> Skip Step5.8..."

conda deactivate

#=======================================Step 6: PhosSight: Peptide Detectability===========================#

source $anacondaPath/etc/profile.d/conda.sh
conda activate cv

echo "Step 6: Run PhosSight"

max_len=53
echo "Step 6.1: Train PhosSight model"
phosSight_train_Command="python $scriptPath/PhosDetect/code/program_train.py -w $phosSightPretrainedModelPath -m $pth_modelPath/phossight_finetune.pth -p $phosSight_trainPath/phossight_train.txt -e 5 --patience 5 -lr 1e-4 --device cuda --max_len $max_len"

eval "$phosSight_train_Command"
# echo ">>> Skip Step6.1..."

echo "Step 6.2: PhosSight for nonPhospho"

phosSight_predict_Command="python $scriptPath/PhosDetect/code/program_predict.py -m $pth_modelPath/phossight_finetune.pth -p $phosSight_predictionPath/phossight_prediction.nonPhospho.txt -o $pth_predictionPath/PhosSight.Predict.nonPhospho.txt --max_len $max_len"

eval "$phosSight_predict_Command"
# echo ">>> Skip Step6.2..."

echo "Step 6.3: PhosSight for Phospho"
phosSight_predict_Command="python $scriptPath/PhosDetect/code/program_predict.py -m $pth_modelPath/phossight_finetune.pth -p $phosSight_predictionPath/phossight_prediction.Phospho.txt -o $pth_predictionPath/PhosSight.Predict.Phospho.txt --max_len $max_len" 

eval "$phosSight_predict_Command"
# echo ">>> Skip Step6.3..."

#====================================================================================#

#=======================================Step 7: Deep-relocalization===========================#

echo "Step 7: Deep-relocalization"

source $anacondaPath/etc/profile.d/conda.sh
conda activate R_env

run_DeepLocalization_Command="Rscript $scriptPath/DeepRelocalization/calculate_localization_probability_entropy.R \
\"$phosphoRSResultsPath/PhosphoRS.txt\" \
\"$autoRT_resultsPath/tf_prediction/phospho.prediction.tsv\" \
\"$pDeep3_resultsPath/pDeep3PredictionResults.Phospho.txt\" \
\"$featurePath/features.PhosphoRS.txt\" \
\"$featurePath/Features.Localization.entropy.txt\" \
\"$VariableMods\" \
\"$FixedMods\" \
\"$pth_predictionPath/PhosSight.Predict.Phospho.txt\""


eval "$run_DeepLocalization_Command"
# echo ">>> Skip Step7..."

#====================================================================================#

#=======================================Step 8: Rescoring using Percolator===========================#

echo "Step 8: Rescoring using Percolator"
echo "Step 8.1: Method1 Results"

Method1Results_Command="Rscript $scriptPath/Percolator/PhosphoRSResults.R $PGAPath/peptide_level/pga-peptideSummary.txt $PGAPath/psm_level/pga-peptideSummary.txt $featurePath/Features.Localization.entropy.txt $Method1ResultsPath/Method1Results.txt"

eval "$Method1Results_Command"
# echo ">>> Skip 8.1..."

echo "Step 8.2: Generate Percolator Input"
GeneratePercolatorInput_Command="Rscript $scriptPath/Percolator/format_percolator_input_PhosSight.R \
\"$featurePath/Features.Localization.entropy.txt\" \
\"$PGAPath/pga-rawPSMs.txt\" \
\"$PercolatorPath/PhosSight.pin\" \
\"$searchEngine\" \
\"$phosphoRSResultsPath/PhosphoRS.txt\" \
\"$pDeep3_resultsPath/pDeep3PredictionResults.Phospho.txt\" \
\"$pDeep3_resultsPath/pDeep3PredictionResults.nonPhospho.txt\" \
\"$autoRT_resultsPath/tf_prediction/phospho.prediction.tsv\" \
\"$autoRT_resultsPath/tf_prediction/nonPhospho.prediction.tsv\" \
\"$VariableMods\" \
\"$FixedMods\" \
\"$pth_predictionPath/PhosSight.Predict.Phospho.txt\" \
\"$pth_predictionPath/PhosSight.Predict.nonPhospho.txt\""

eval "$GeneratePercolatorInput_Command"
# echo ">>> Skip Step 8.2..."

conda deactivate

echo "Step 8.3: Run Percolator"
docker_Command="docker run --rm -u 0:0 -v $outputPath/:/data/ -t bzhanglab/percolator:3.4 percolator"
Percolator_Command="$docker_Command ./Percolator/PhosSight.pin -r ./Percolator/PhosSight/PhosSight.pep.txt -m ./Percolator/PhosSight/PhosSight.psms.txt -w ./Percolator/PhosSight/PhosSight.weights.txt -M ./Percolator/PhosSight/PhosSight.decoy.psms.txt"

eval "$Percolator_Command"
# echo ">>> Skip Step 8.3..."

#====================================================================================#

#=======================================Step 9: Get PhosSight Results===========================#

echo "Step 9: PhosSight Results"

source $anacondaPath/etc/profile.d/conda.sh
conda activate R_env

PhosSightResults_Command="Rscript $scriptPath/Percolator/GetPhosSightResults_v2.R $featurePath $PhosSightResultsPath $outputPath"

eval "$PhosSightResults_Command"
# echo ">>> Skip Step 9..."

conda deactivate

#====================================================================================#