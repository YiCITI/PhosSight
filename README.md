# Title of Paper

## Contents

- [Title of Paper](#title-of-paper)
  - [Contents](#contents)
  - [DDA](#dda)
    - [Files Structure](#files-structure)
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
  - [Contact](#contact)
  - [Acknowledgements](#acknowledgements)
  - [References](#references)

## DDA

### Files Structure

```
PhosSight-DDA
|---Script
|---|---DeepRelocalization
|---|---Features
|---|---generate_train_prediction
|---|---GenerateFeatureMatrix
|---|---KinaseActivityScoreInference
|---|---pDeep3
|---|---Percolator
|---|---PGA
|---|---PhosphoRS
|---|---TMTQuantification
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

| Name             | Value            | Description                                                                                                                                         |
| ---------------- | ---------------- | --------------------------------------------------------------------------------------------------------------------------------------------------- |
| PhosSightPath | PhosSight_DIR | PhosSight directory                                                                                                                              |
| anacondaPath     | ANACONDA_DIR     | Anaconda directory. Default is /C/ProgramData/anaconda3                                                                                             |
| decoyPrefix      | DECOY_PREFIX     | Decoy prefix used for searching. Default is XXX_                                                                                                    |
| searchEngine     | SEARCH_ENGINE    | Four search engines, msgf, comet, xtandem, maxquant, are supported                                                                                  |
| rawSpectraPath   | RAW_DIR          | Path to the MS/MS spectra (RAW) directory                                                                                                           |
| spectraPath      | MGF_DIR          | Path to the MS/MS spectra (MGF) directory                                                                                                           |
| inputFeaturePath | FEATURE_DIR      | Path to the feature matrix                                                                                                                          |
| outputPath       | OUT_DIR          | Output directory                                                                                                                                    |
| VariableMods     | VAR_MOD          | Variable modifications used for searching, e.g. '1,Oxidation,M,15.994919,1;2,Phospho,S,79.966331,2;3,Phospho,T,79.966331,2;4,Phospho,Y,79.966331,2' |
| FixedMods        | Fix_MOD          | Fixed modifications used for searching, e.g. '5,Carbamidomethyl,C,57.021464,3'. If null, use 'null'                                                 |
| ModsReplace      | RENAME_MOD       | Some modifications need to rename, e.g. '[79.966331],Phospho'. If null, use 'null'                                                                  |

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

### Other functions

#### Quantification for TMT dataset

In our manuscript, we used [MASCI](https://github.com/PNNL-Comp-Mass-Spec/MASIC) to perform the TMT quantification for both TMT10 (UCEC) and TMT11 (HCC) datasets. We prepared the original scripts we used for the quantification under the 'Script/TMTQuantification' folder. You can change the input data path and parameters used for MASCI following our scripts to do the TMT quantification.

#### Kinase activity score inference

In our manuscript, we performed kinase activity score inference for the HCC datasets. We prepared the original scripts we used under the 'Script/KinaseActivityScoreInference' folder. The excel file ('mmc4.xlsx') contains the list of known targets that we used for the inference. You can change the input data path and parameters to do the kinase activity score inference.

## DIA

Coming soon...

## Contact

Xinpei Yi - [@yixinpei](https://twitter.com/yixinpei) - yixinpei13@gmail.com

## Acknowledgements

- [DeepRescore2](https://github.com/bzhanglab/DeepRescore2)
- [AutoRT](https://github.com/bzhanglab/AutoRT)
- [pDeep3](https://github.com/pFindStudio/pDeep3)
- [PhosphoRS](https://github.com/lmsac/phosphoRS-cli)
- [SpectralEntropy](https://github.com/YuanyueLi/SpectralEntropy)

## References

```

```