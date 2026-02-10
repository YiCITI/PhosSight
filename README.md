
## Contents
  - [Contents](#contents)
  - [Phosdetect](#phosdetect)
    - [Files Structure](#files-structure)
    - [Model Description](#model-description)
    - [How to Train](#how-to-train)
    - [How to Test](#how-to-test)
    - [How to Infer](#how-to-infer)
  - [DDA](#dda)
    - [Files Structure](#files-structure-1)
    - [How to Use](#how-to-use)
    - [Data type](#data-type)
    - [Download example data](#download-example-data)
    - [Directory structure of input](#directory-structure-of-input)
    - [Parameters of PhosSight](#parameters-of-phossight)
    - [Run PhosSight](#run-phossight)
    - [Output](#output)
    - [Other functions](#other-functions)
      - [Quantification for TMT dataset](#quantification-for-tmt-dataset)
      - [Kinase activity score inference](#kinase-activity-score-inference)
  - [DIA](#dia)
    - [Files Structure](#files-structure-2)
  - [TMT Quantification Analysis](#tmt-quantification-analysis)
    - [Files Structure](#files-structure-3)
    - [How to Use](#how-to-use-1)
  - [Contact](#contact)
  - [Acknowledgements](#acknowledgements)
  - [References](#references)

## Phosdetect

Phosdetect is a deep learning model for phosphopeptide identification based on BiGRU architecture with physicochemical property features. Both DDA and DIA workflows are based on Phosdetect for phosphopeptide identification and analysis.

### Files Structure

```
Phosdetect
|---model.py              # Model definition (includes V2 improvements)
|---train.py              # Training script
|---test.py               # Testing script
|---infer.py              # Inference script for custom sequences
|---run_training.py       # Convenient training wrapper script
|---run_testing.py        # Convenient testing wrapper script
|---README.md             # Phosdetect documentation
|---data/                 # Data directory
|---|---train/            # Training dataset
|---|---|---balanced_dataset_1.csv
|---|---test/             # Test datasets
|---|---|---DeepDetect_ecoli.csv
|---|---|---DeepDetect_human.csv
|---|---|---DeepDetect_mouse.csv
|---|---|---DeepDetect_yeast.csv
|---|---|---DeepRescore2_HCC.csv
|---|---|---DeepRescore2_label_free.csv
|---|---|---DeepRescore2_UCEC.csv
|---logs/                 # Training and testing logs directory
|---|---20250810_164602/
|---|---20250810_165125/
|---models/               # Model weights directory
|---|---20250810_164602/
|---|---20250810_165125/
|---|---|---best_model.pth
```

- **model.py** implements the BiGRU-based models for phosphopeptide detection, including the improved V2 version with physicochemical property features.
- **train.py** and **run_training.py** provide training functionality for the Phosdetect model.
- **test.py** and **run_testing.py** provide testing functionality for evaluating model performance.
- **infer.py** allows inference on custom peptide sequences.
- **data/** contains the training and test datasets:
  - **data/train/** stores the training dataset (`balanced_dataset_1.csv`)
  - **data/test/** stores seven test datasets for model evaluation:
    - `DeepDetect_ecoli.csv` - E. coli dataset from DeepDetect
    - `DeepDetect_human.csv` - Human dataset from DeepDetect
    - `DeepDetect_mouse.csv` - Mouse dataset from DeepDetect
    - `DeepDetect_yeast.csv` - Yeast dataset from DeepDetect
    - `DeepRescore2_HCC.csv` - HCC dataset from DeepRescore2
    - `DeepRescore2_label_free.csv` - Label-free dataset from DeepRescore2
    - `DeepRescore2_UCEC.csv` - UCEC dataset from DeepRescore2
- **logs/** stores training and testing logs with timestamps.
- **models/** stores trained model weights with timestamps.

### Model Description

Phosdetect uses a bidirectional GRU (BiGRU) architecture to identify phosphopeptides from peptide sequences. The V2 improved version incorporates physicochemical property features including hydrophobicity, charge, and polarity for each amino acid and phosphorylated residue.

The model supports the following amino acids and modifications:
- Standard 20 amino acids (A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y)
- Phosphorylated residues (s for pSer, t for pThr, y for pTyr)
- Padding token (Z)

### How to Train

To train the Phosdetect model, you can use either the convenient wrapper script or the training script directly.

#### Using the convenient wrapper script:

```bash
cd Phosdetect
python run_training.py --data_path ./data/train/balanced_dataset_1.csv
```

#### Using the training script directly:

```bash
cd Phosdetect
python train.py -p ./data/train/balanced_dataset_1.csv \
                -m ./models/best_model.pth \
                --model_type bigru_improved_v2 \
                -b 128 -e 100 -lr 0.0005 -d 0.3
```

#### Training Parameters:

- `--data_path` or `-p`: Path to training data (CSV format: sequence, label)
- `--batch_size` or `-b`: Batch size (default: 128)
- `--epochs` or `-e`: Number of training epochs (default: 100)
- `--learning_rate` or `-lr`: Learning rate (default: 0.0005)
- `--dropout` or `-d`: Dropout rate (default: 0.3)
- `--model_type`: Model type (default: bigru_improved_v2)
- `--patience`: Early stopping patience (default: 15)
- `--optimizer`: Optimizer type - adam, adamw, or sgd (default: adam)
- `--scheduler`: Learning rate scheduler - plateau, cosine, or step (default: plateau)
- `--in_features`: Input feature dimension (default: 10)
- `--out_features`: Output feature dimension (default: 20)
- `--num_layers`: Number of GRU layers (default: 2)
- `--log_file`: Path to log file (optional)
- `--verbose`: Enable verbose output

The training process will automatically:
- Create timestamped directories for logs and model weights
- Implement early stopping to prevent overfitting
- Save the best model based on validation performance
- Log training metrics to both console and log file

### How to Test

To test the trained Phosdetect model, you can use either the convenient wrapper script or the testing script directly.

#### Using the convenient wrapper script:

```bash
cd Phosdetect
python run_testing.py --model_path ./models/20250810_165125/best_model.pth
```

The wrapper script will automatically test on multiple datasets:
- DeepDetect_ecoli
- DeepDetect_human
- DeepDetect_mouse
- DeepDetect_yeast
- DeepRescore2_HCC
- DeepRescore2_label_free
- DeepRescore2_UCEC

#### Using the testing script directly:

```bash
cd Phosdetect
python test.py -p ./data/test/DeepDetect_human.csv \
               -m ./models/best_model.pth \
               --model_type bigru_improved_v2
```

#### Testing Parameters:

- `-p` or `--data_path`: Path to test data (CSV format: sequence, label)
- `-m` or `--model_path`: Path to trained model weights
- `-b` or `--batch_size`: Batch size (default: 128)
- `--model_type`: Model type (default: bigru_improved_v2)
- `--log_file`: Path to log file (optional)
- `--output_file`: Path to save results CSV (optional)
- `--in_features`: Input feature dimension (must match training, default: 10)
- `--out_features`: Output feature dimension (must match training, default: 20)
- `--num_layers`: Number of GRU layers (must match training, default: 2)
- `--dropout`: Dropout rate (must match training, default: 0.3)
- `--verbose`: Enable verbose output

The testing process will:
- Calculate performance metrics (AUC, Accuracy, F1-score, Precision, Recall)
- Save results to a CSV file if specified
- Log testing metrics to both console and log file

### How to Infer

To perform inference on custom peptide sequences:

```bash
cd Phosdetect
python infer.py --model_path ./models/20250810_165125/best_model.pth \
                --seq_file sequences.txt \
                --output_csv results.csv
```

Or provide sequences directly via command line:

```bash
cd Phosdetect
python infer.py --model_path ./models/20250810_165125/best_model.pth \
                --seq "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWSTPSELGHAGLNGDILVWNPVLEDAFELSSMGIRVDADTLKHQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQAGVWPAAVRESVPSLL" \
                --output_csv results.csv
```

#### Inference Parameters:

- `--model_path`: Path to trained model weights (required)
- `--seq`: Peptide sequences provided directly on command line (space-separated)
- `--seq_file`: Path to text file with sequences (one per line)
- `--batch_size`: Batch size for inference (default: 128)
- `--in_features`: Input feature dimension (must match training, default: 10)
- `--out_features`: Output feature dimension (must match training, default: 20)
- `--num_layers`: Number of GRU layers (must match training, default: 2)
- `--dropout`: Dropout rate (must match training, default: 0.3)
- `--output_csv`: Path to save results CSV (optional)

The inference will output a CSV file with columns: `sequence` and `probability`, where probability indicates the model's confidence that the sequence is a phosphopeptide.

## DDA

### Files Structure

```
PhosSight-DDA
|---Script
|---|---DeepRelocalization
|---|---Features
|---|---generate_train_prediction
|---|---GenerateFeatureMatrix
|---|---pDeep3
|---|---Percolator
|---|---PGA
|---|---PhosphoRS
|---|---PhosSight.sh
|---Parameters
|---|---PXD000138_maxquant.param
|---|---PXD023665_comet.param
|---|---PXD023665_maxquant.param
|---|---PXD023665_msgf.param
|---|---PXD023665_xtandem.param
|---|---UCEC_comet.param
|---|---UCEC_maxquant.param
|---|---UCEC_msgf.param
|---|---UCEC_xtandem.param
|---Install
|---|---InstallDependency.sh
|---|---environment_R.yml
```

- **Script** implements PhosSight to improve phosphopeptide identification and phosphosite localization.
- **Parameters** include 9 parameter files for the three test datasets of four search engines used in our manuscript, including label free dataset (PRIDE ID: PXD000138 and PXD023665) and UCEC TMT dataset, respectively.
- **Install** includes scripts for installing PhosSight, e.g., configuring the deep learning conda environment, and installing R packages.

### How to Use

1. Install [Git](https://git-scm.com/downloads), [Docker](https://docs.docker.com/install/) and [Anaconda](https://www.anaconda.com/download) on Windows Platform.
2. Clone the repository:
   ```shell
   git clone https://github.com/YiCITI/PhosSight.git
   cd PhosSight
   ```
3. Install the required packages and dependencies:
   ```shell
   cd PhosSight-DDA/Install
   # use Git Bash to run the script below!!!
   bash InstallDependency.sh  # run the script to install all dependencies or run step by step manually
   ```
4. Download the datasets from [Zenodo](https://zenodo.org/records/10049730)
5. Modify the parameter files in `Parameters/` and run the `Script/PhosSight.sh` script:
   ```shell
    cd PhosSight-DDA/Script
    # use Git Bash to run the script below!!!
    # For example:
    bash PhosSight.sh ../Parameters/PXD000138_maxquant.param  # run the script or execute step by step manually
   ```

### Data type

- Instrument Type: The proposed approach, PhosSight, is designed to be applicable to mass spectrometry-based proteomics data obtained from various types of instruments, including but not limited to Orbitrap, Q-TOF, and ion trap instruments.
- Peptide Type: PhosSight is applicable to different types of peptides, including both labeled (e.g., TMT-labeled) and unlabeled peptides. However, it is important to note that PhosSight is specifically designed to handle peptides with phosphorylation modification.
- The current version supports four search engines, [MS-GF+ (v2019.02.28)](https://github.com/MSGFPlus/msgfplus), [Comet (2018.01 rev.4)](http://comet-ms.sourceforge.net/), [X!Tandem (v2017.2.1.2)](https://www.thegpm.org/TANDEM/), and [MaxQuant (v1.6.5.0)](https://maxquant.org/).
- Computational Requirements: Currently, PhosSight only supports running on Windows systems. The computational requirements for running PhosSight depend on the size of the dataset and the specific hardware configuration. PhosSight utilizes deep learning models, and the computational demands may increase with larger datasets. We recommend running PhosSight on a machine with sufficient computational resources, such as a multi-core CPU and a GPU, to ensure efficient processing.

### Download example data

- Please go to https://zenodo.org/records/10049730 and download **ExampleData1.zip (Synthetic dataset, PXD000138)**, **ExampleData2.zip (Label free dataset, PXD023665)**, **ExampleData3.zip (TMT dataset, UCEC)** used in our manuscript. Unzip these files as the input for PhosSight.

### Directory structure of input

In order to perform PhosSight, the input dataset for PhosSight must be prepared as follows.

```
|---Raw_input_directory
|---|---MGF
|---|---|---Spectra1.mgf
|---|---|---Spectra2.mgf
             ...
|---|---|---SpectraN.mgf
|---|---RAW
|---|---|---Spectra1.raw
|---|---|---Spectra2.raw
             ...
|---|---|---SpectraN.raw
|---|---features_matrix.txt
```

- **MGF** includes the MS/MS spectra (MGF format).
- **RAW** includes the MS/MS spectra (RAW format).
- **features_matrix.txt** is the path to the feature matrix which contains all the necessary features as follows:

<table>
  <tr>
    <th rowspan="1">Feature groups</th>
    <th>Feature name</th>
    <th>Feature description</th>
  </tr>
  <tr>
    <td rowspan="2">Features based on deep learning</td>
    <td>RT Ratio</td>
    <td>RT ratio  between observed RT and predicted RT</td>
  </tr>
  <tr>
    <td>Spectrum similarity</td>
    <td>The spectral similarity characterized by entropy distance between predicted MS/MS spectrum and experimental MS/MS spectrum of a peptide</td>
  </tr>
  <tr>
    <td>Peptide detectability</td>
    <td>The detectability of a peptide predicted by PhosDetect</td>
  </tr>
  <tr>
    <td rowspan="7">Search engine independent features</td>
    <td>Mass_Error</td>
    <td>Difference between theoretical and experimental mass</td>
  </tr>
  <tr>
    <td>Charge</td>
    <td>Peptide charge</td>
  </tr>
  <tr>
    <td>Abs_Mass_Error</td>
    <td>Absolute value of the difference between theoretical and experimental mass</td>
  </tr>
  <tr>
    <td>Ln_Total_Intensity</td>
    <td>Total intensity, natural logarithm transformed</td>
  </tr>
  <tr>
    <td>Match_Ions_Intensity</td>
    <td>Total intensity of matched ions, natural logarithm transformed</td>
  </tr>
  <tr>
    <td>Max_Match_Ion_Intensity</td>
    <td>Max intensity of matched fragment ions</td>
  </tr>
  <tr>
    <td>Rel_Match_Ions_Intensity</td>
    <td>The total intensity of all matched ions divided by the total intensity of the spectrum</td>
  </tr>
  <tr>
    <td rowspan="5">Search engine specific features (Comet (2018.01 rev.4))</td>
    <td>xcorr</td>
    <td>Cross-correlation of the experimental and theoretical spectra</td>
  </tr>
  <tr>
    <td>deltacn</td>
    <td>The normalized difference of XCorr values between the best sequence and the next best sequence</td>
  </tr>
  <tr>
    <td>spscore</td>
    <td>The spscore of Comet</td>
  </tr>
  <tr>
    <td>sprank</td>
    <td>The sprank score of Comet</td>
  </tr>
  <tr>
    <td>Ln_expect</td>
    <td>Comet  Evalue, natural logarithm transformed</td>
  </tr>
  <tr>
    <td rowspan="3">Search engine specific features (MaxQuant (v1.6.5.0))</td>
    <td>Score</td>
    <td>Andromeda score</td>
  </tr>
  <tr>
    <td>Ln-PEP</td>
    <td>Posterior Error Probability of the identification, natural logarithm transformed</td>
  </tr>
  <tr>
    <td>Delta_Score</td>
    <td>Score difference to the second best identified peptide</td>
  </tr>
  <tr>
    <td rowspan="4">Search engine specific features (MS-GF+ (v2019.02.28))</td>
    <td>MS-GF:RawScore</td>
    <td>Raw match score of MS-GF+</td>
  </tr>
  <tr>
    <td>MS-GF:DeNovoScore</td>
    <td>Maximum possible raw match score to this spectrum</td>
  </tr>
  <tr>
    <td>MS-GF:SpecEValue</td>
    <td>Negative MS-GF+ Spectral E Value, logged</td>
  </tr>
  <tr>
    <td>Ln-MS-GF:EValue</td>
    <td>Negative MS-GF+ E value, logged</td>
  </tr>
  <tr>
    <td rowspan="2">Search engine specific features (X!Tandem (v2017.2.1.2))</td>
    <td>Ln-X!Tandem:expect</td>
    <td>X!Tandem  Evalue, natural logarithm transformed</td>
  </tr>
  <tr>
    <td>X!Tandem:hyperscore</td>
    <td>X!Tandem hyperscore</td>
  </tr>
</table>

We used PDV (PDV-1.6.1.beta.features-jar-with-dependencies.jar) attached under the 'Script/GenerateFeatureMatrix' folder to generate feature matrix. The script to run this jar file based on the Comet (2018.01 rev.4) identifications is as follows:

```sh
java -Xmx100g -jar ./Script/GenerateFeatureMatrix/PDV-1.6.1.beta.features-jar-with-dependencies.jar \
  -r ./ExampleData/PXD023665/Comet.pep.xml \
  -rt 2 \
  -s ./ExampleData/Combined.mgf \
  -st 1 \
  -i * \
  -k s \
  -o . \
  -a 0.02 \
  -c 0 \
  -decoy REV_ \
  -ft pdf \
  --features

```

### Parameters of PhosSight

Each column of the parameter file is described as follows (Please change the 'Value' column based on your data):

| Name             | Value         | Description                                                                                                                                         |
| ---------------- | ------------- | --------------------------------------------------------------------------------------------------------------------------------------------------- |
| PhosSightPath    | PhosSight_DIR | PhosSight directory                                                                                                                                 |
| anacondaPath     | ANACONDA_DIR  | Anaconda directory. Default is /C/ProgramData/anaconda3                                                                                             |
| decoyPrefix      | DECOY_PREFIX  | Decoy prefix used for searching. Default is XXX_                                                                                                    |
| searchEngine     | SEARCH_ENGINE | Four search engines, msgf, comet, xtandem, maxquant, are supported                                                                                  |
| rawSpectraPath   | RAW_DIR       | Path to the MS/MS spectra (RAW) directory                                                                                                           |
| spectraPath      | MGF_DIR       | Path to the MS/MS spectra (MGF) directory                                                                                                           |
| inputFeaturePath | FEATURE_DIR   | Path to the feature matrix                                                                                                                          |
| outputPath       | OUT_DIR       | Output directory                                                                                                                                    |
| VariableMods     | VAR_MOD       | Variable modifications used for searching, e.g. '1,Oxidation,M,15.994919,1;2,Phospho,S,79.966331,2;3,Phospho,T,79.966331,2;4,Phospho,Y,79.966331,2' |
| FixedMods        | Fix_MOD       | Fixed modifications used for searching, e.g. '5,Carbamidomethyl,C,57.021464,3'. If null, use 'null'                                                 |
| ModsReplace      | RENAME_MOD    | Some modifications need to rename, e.g. '[79.966331],Phospho'. If null, use 'null'                                                                  |

As a reference, we prepared 9 parameter files for the three test datasets of four search engines used in our manuscript, including label free dataset (PRIDE ID: PXD000138 and PXD023665) and UCEC TMT dataset, respectively. Please check the 'PhosSight/Parameters' folder.

### Run PhosSight

- Open the Docker Desktop.
- Open the PhosSight Parameters folder and edit the parameters.
- Open the PhosSight Script folder and run PhosSight.

```
$ cd PhosSight/Script
$ bash PhosSight.sh $param_path
```

### Output

PhosSight will output results of each step, including

* Features
* PhosphoRS
* PGA
* generate_train_prediction
* autoRT_Results
* pDeep3_Results
* Percolator

PhosSight also output two tables as the final results:

* File named 'Method1Results.txt' which is filtered using both PGA FDR < 1% and PhosphoRS localization probability > 0.75.
* File named 'PhosSightResults.txt' which is filtered using both q-value < 1% and DeepLocalization probability > 0.75.




## DIA

### Files Structure

```
PhosSight-DIA
|--- spec_parquet_filter    # Filter the spectral library in parquet format
|--- Script                 # Scripts for PhosSight-DIAs
|     |--- PhosSight.sh     # Main script for PhosSight-DIA pipeline
```

### Prerequisites

Before using PhosSight-DIA, please install Singularity and set up a Singularity container as recommended by DIA-NN. For detailed instructions, you can refer to the following official resources:

DIA-NN's official documentation: [Official Setup Guide](https://github.com/vdemichev/DiaNN/?tab=readme-ov-file#installation)

Recommended community guide: [Community Guide](https://github.com/vdemichev/DiaNN/issues/1202#issuecomment-2417108281)

We recommend reviewing these links to ensure a correct installation.



Coming soon...

## TMT Quantification Analysis

TMT Quantification Analysis provides a complete workflow for TMT-labeled phosphoproteomic data quantification and analysis, including MASCI-based ion intensity extraction, site-level quantification, and comparative analysis between PhosSight and PhosphoRS methods.In our manuscript, we used [MASCI](https://github.com/PNNL-Comp-Mass-Spec/MASIC) to perform the TMT quantification for TMT10 (UCEC)  datasets. We prepared the original scripts we used for the quantification. You can change the input data path and parameters used for MASCI following our scripts to do the TMT quantification.


### Files Structure

```
TMTquantificationanalysis
|---code/                          # Main analysis scripts
|---|---generate_MASCIParameters.py # Generate MASCI parameter files
|---|---MASCI_param/               # MASCI parameter files
|---|---|---TMT10_LTQ-FT_10ppm_ReporterTol0.003Da_2014-08-06.xml
|---|---|---TMT11_LTQ-FT_10ppm_ReporterTol0.003Da_2017-03-17.xml
|---|---run_MASCI_simple.bat        # Run MASCI quantification (Windows)
|---|---run_MASCI_simple.sh          # Run MASCI quantification (Linux)
|---|---run_complete_workflow_full.bat # Complete workflow automation
|---|---step1_GetIonIntensity_PhosSight_17plex.R    # Extract ion intensities for PhosSight
|---|---step1_GetIonIntensity_PhosphoRS_17plex.R     # Extract ion intensities for PhosphoRS
|---|---step2_sitequant_full_17plex.py                # Site-level quantification
|---|---step3_AdjustSiteTable_17plex.py               # Adjust site table
|---|---step4_MergeSiteTables_17plex.py               # Merge site tables
|---|---step5_GetUniprotIDGeneName_17plex.py           # Add Uniprot ID and gene names
|---|---README.md                    # Detailed workflow documentation
|---drawfig7/                       # Scripts for Figure 7 generation
|---|---Create_QuantifiableSites_100Samples.R
|---|---Create_PhosSight_vs_PhosphoRS_GainSharedLoss.R
|---|---Create_MissingValue_SiteCount_Comparison.R
|---|---README.md
|---drawsuppfigs8/                  # Scripts for Supplementary Figure 8 generation
|---|---Create_PhosSight_Boxplot_AllSamples.R
|---|---Create_PhosSight_Density_Plot.R
|---|---Create_PhosSight_PCA_Plot.R
|---|---README.md
```

- **code/** contains the main analysis scripts for TMT quantification workflow:
  - **generate_MASCIParameters.py** generates MASCI parameter files for batch processing
  - **MASCI_param/** stores MASCI parameter files for TMT10 and TMT11 labeling
  - **run_MASCI_simple.bat/sh** runs MASCI quantification on RAW files
  - **step1_GetIonIntensity_*.R** extracts reporter ion intensities from MASCI results and matches with identification results
  - **step2_sitequant_full_17plex.py** performs site-level quantification using full sitequant workflow
  - **step3-5** scripts perform post-processing: table adjustment, merging, and annotation
- **drawfig7/** contains R scripts for generating comparison figures between PhosSight and PhosphoRS methods
- **drawsuppfigs8/** contains R scripts for generating quality control and analysis plots

### Data Requirements

The TMT Quantification Analysis workflow requires the following data files:

#### TMT01-TMT17 Dataset

- **TMT批次**: TMT01-TMT17 (17 batches total)
- **TMT通道**: 170 channels (17 × 10 channels per batch)
- **样本通道**: 153 channels (excluding 17 reference channels 126)

#### Required Data Files

1. **MASCI Quantification Results**
   - Location: MASCI output directories for each TMT batch (e.g., `QuantificationResults01` to `QuantificationResults17`)
   - Each folder should contain MASCI output files for the corresponding TMT batch
   - Required files: ReporterIons and SICstats files
   - **Note**: Update the path in the R scripts (`step1_GetIonIntensity_*.R`) to match your MASCI output directory

2. **Identification Results**
   - PhosSight identification results for each TMT batch (TMT01-TMT17)
   - PhosphoRS identification results for each TMT batch (TMT01-TMT17)
   - These will be merged into: `PhosSightResults_TMT01_17.txt` and `PhosphoRSResults_TMT01_17.txt`

3. **FASTA File**
   - Protein sequence database file (e.g., `swiss_prot_human_20190214_target_conts_96Libraries.fasta`)
   - Required for site-level quantification and precise site localization

4. **RAW Files** (for MASCI quantification)
   - Original RAW mass spectrometry files for each TMT batch
   - Used as input for MASCI to extract reporter ion intensities

**Note**: The TMT01-TMT17 dataset used in our manuscript can be downloaded via （https://pdc.cancer.gov/pdc/study/PDC000125） or obtained from the corresponding authors. Please ensure all data files are properly organized before running the workflow.

### How to Use

#### Prerequisites

1. Install [MASCI](https://github.com/PNNL-Comp-Mass-Spec/MASIC) for TMT reporter ion quantification
2. Ensure R environment is set up with required packages (tidyverse, data.table, ggplot2)
3. Ensure Python environment is set up with required packages (pandas, numpy)
4. Prepare the TMT01-TMT17 dataset and organize files according to the data requirements above

#### Step 1: MASCI Quantification

Run MASCI to extract reporter ion intensities from RAW files:

```bash
cd TMTquantificationanalysis/code
# Windows
run_MASCI_simple.bat

# Linux
bash run_MASCI_simple.sh
```

Or use the parameter generator:

```bash
python generate_MASCIParameters.py <input_dir> <output_script> <output_dir> <TMTType>
```

#### Step 2: Extract Ion Intensities

Match MASCI quantification results with identification results:

```bash
# For PhosSight results
Rscript step1_GetIonIntensity_PhosSight_17plex.R

# For PhosphoRS results
Rscript step1_GetIonIntensity_PhosphoRS_17plex.R
```

#### Step 3: Site-Level Quantification

Perform site-level quantification using the full sitequant workflow:

```bash
python step2_sitequant_full_17plex.py -i <intensity_file> -f <fasta_file> -o <output_file>
```

#### Step 4-5: Post-Processing

Adjust, merge, and annotate site tables:

```bash
python step3_AdjustSiteTable_17plex.py -i <input> -o <output>
python step4_MergeSiteTables_17plex.py -i <input> -o <output>
python step5_GetUniprotIDGeneName_17plex.py -i <input> -o <output>
```

#### Complete Workflow

Run the complete workflow automatically:

```bash
cd TMTquantificationanalysis/code
run_complete_workflow_full.bat
```

#### Generate Figures

Generate comparison and analysis figures:

```bash
cd TMTquantificationanalysis/drawfig7
Rscript Create_QuantifiableSites_100Samples.R
Rscript Create_PhosSight_vs_PhosphoRS_GainSharedLoss.R
Rscript Create_MissingValue_SiteCount_Comparison.R

cd ../drawsuppfigs8
Rscript Create_PhosSight_Boxplot_AllSamples.R
Rscript Create_PhosSight_Density_Plot.R
Rscript Create_PhosSight_PCA_Plot.R
```

**Note**: Large data files (>100MB) required for these scripts are not included in the repository. Please refer to the README files in each folder for data file requirements.

## Contact

Xinpei Yi - [@yixinpei](https://scholar.google.com/citations?user=Z4lICl8AAAAJ&hl=en) - yixinpei13@gmail.com
<br/>Lab Website: [YiCITI](https://github.com/YiCITI/)


## Acknowledgements

We acknowledge the following open-source software and tools used in this project:

### Deep Learning and Rescoring
- [DeepRescore2](https://github.com/bzhanglab/DeepRescore2) - Deep learning-based peptide rescoring
- [AutoRT](https://github.com/bzhanglab/AutoRT) - Retention time prediction
- [pDeep3](https://github.com/pFindStudio/pDeep3) - MS/MS spectrum prediction

### Phosphosite Localization
- [PhosphoRS](https://github.com/lmsac/phosphoRS-cli) - Phosphosite localization probability calculation

### Spectral Analysis
- [SpectralEntropy](https://github.com/YuanyueLi/SpectralEntropy) - Spectral similarity calculation using entropy distance

### TMT Quantification
- [MASCI](https://github.com/PNNL-Comp-Mass-Spec/MASIC) - Mass spectrometry analysis and quantification for TMT-labeled datasets

### Search Engines
- [MS-GF+](https://github.com/MSGFPlus/msgfplus) - Database search engine
- [Comet](http://comet-ms.sourceforge.net/) - Database search engine
- [X!Tandem](https://www.thegpm.org/TANDEM/) - Database search engine
- [MaxQuant](https://maxquant.org/) - Database search and quantification platform

### DIA Analysis
- [DIA-NN](https://github.com/vdemichev/DiaNN) - Data-independent acquisition data analysis

## References
Ben Wang, Zhiyuan Cheng, Chengying She, Jiahui Zhang, Lin Lv, Hongwen Zhu, Lizhuang Liu, Yan Fu, Xinpei Yi, **PhosSight: A Unified Deep Learning Framework Boosting and Accelerating Phosphoproteomic Identification to Enable Biological Discoveries**. *Under review* (2026).

