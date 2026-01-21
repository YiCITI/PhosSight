@echo off
REM ========================================
REM TMT01-TMT17 Complete Workflow (Full SiteQuant)
REM ========================================
REM 
REM Steps:
REM   Step 0: Merge identification results
REM   Step 1: GetIonIntensity - Match MASCI quantification with identification
REM   Step 2: SiteQuant (Full version) - Site quantification
REM   Step 3: AdjustSiteTable - Adjust column names
REM   Step 4: MergeSiteTables - Merge site tables
REM   Step 5: GetUniprotIDGeneName - Add gene names
REM   Step 6: Analysis - Missing value analysis and plotting
REM
REM Configuration:
REM   - 17 TMT plexes: TMT01-TMT17
REM   - 170 TMT channels (17 x 10)
REM   - 153 sample channels (excluding 17 reference channels)
REM   - MASCI results: D:\MASCIresult\QuantificationResults01-17
REM
REM ========================================

setlocal enabledelayedexpansion

echo ========================================
echo TMT01-TMT17 Complete Workflow (Full SiteQuant)
echo ========================================
echo.

set "WORK_DIR=%~dp0"
cd /d "%WORK_DIR%"
echo Working directory: %WORK_DIR%
echo.

REM ===== Configuration =====
set "FASTA_FILE=C:\Users\wk\Desktop\TMTQuantification\TMTQuantification\step5.GetGeneNameBasedOnUniprotID\swiss_prot_human_20190214_target_conts_96Libraries.fasta"

REM ===== Step 0: Merge identification results =====
echo ========================================
echo [Step 0/6] Merging identification results
echo ========================================
echo.

call conda activate R_env
if errorlevel 1 (
    echo ERROR: Failed to activate R_env!
    goto :error
)

echo [0] Merging TMT01-TMT17 PhosSight, PhosphoRS, and DeepRescore2 results...
Rscript merge_identification_results.R
if errorlevel 1 (
    echo ERROR: merge_identification_results failed!
    goto :error
)
echo   Done: PhosSightResults_TMT01_17.txt created
echo   Done: PhosphoRSResults_TMT01_17.txt created
echo   Done: DeepRescore2Results_TMT01_17.txt created
echo.

REM ===== Step 1: GetIonIntensity =====
echo ========================================
echo [Step 1/6] GetIonIntensity (17 plexes)
echo ========================================
echo.

echo [1a] Processing PhosSight (17 plexes)...
Rscript step1_GetIonIntensity_PhosSight_17plex.R
if errorlevel 1 (
    echo ERROR: GetIonIntensity_PhosSight failed!
    goto :error
)
echo   Done: PhosSight_Intensity_full.txt created (170 channels)
echo.

echo [1b] Processing PhosphoRS (17 plexes)...
Rscript step1_GetIonIntensity_PhosphoRS_17plex.R
if errorlevel 1 (
    echo ERROR: GetIonIntensity_PhosphoRS failed!
    goto :error
)
echo   Done: PhosphoRS_Intensity_full.txt created (170 channels)
echo.

echo [1c] Processing DeepRescore2 (17 plexes)...
Rscript step1_GetIonIntensity_DeepRescore2_17plex.R
if errorlevel 1 (
    echo ERROR: GetIonIntensity_DeepRescore2 failed!
    goto :error
)
echo   Done: DeepRescore2_Intensity_full.txt created (170 channels)
echo.

REM ===== Step 2: SiteQuant (Full version) =====
echo ========================================
echo [Step 2/6] SiteQuant (Full version, 17 plexes)
echo ========================================
echo.

call conda activate base
if errorlevel 1 (
    echo ERROR: Failed to activate base environment!
    goto :error
)

echo [2a] Processing PhosSight (Full SiteQuant)...
python step2_sitequant_full_17plex.py -i PhosSight_Intensity_full.txt -f "%FASTA_FILE%" -o PhosSight_full_site_table.tsv -m PhosSight --plexes 1-17
if errorlevel 1 (
    echo ERROR: step2_sitequant_PhosSight failed!
    goto :error
)
echo   Done: PhosSight_full_site_table.tsv created
echo.

echo [2b] Processing PhosphoRS (Full SiteQuant)...
python step2_sitequant_full_17plex.py -i PhosphoRS_Intensity_full.txt -f "%FASTA_FILE%" -o PhosphoRS_full_site_table.tsv -m PhosphoRS --plexes 1-17
if errorlevel 1 (
    echo ERROR: step2_sitequant_PhosphoRS failed!
    goto :error
)
echo   Done: PhosphoRS_full_site_table.tsv created
echo.

echo [2c] Processing DeepRescore2 (Full SiteQuant)...
python step2_sitequant_full_17plex.py -i DeepRescore2_Intensity_full.txt -f "%FASTA_FILE%" -o DeepRescore2_full_site_table.tsv -m DeepRescore2 --plexes 1-17
if errorlevel 1 (
    echo ERROR: step2_sitequant_DeepRescore2 failed!
    goto :error
)
echo   Done: DeepRescore2_full_site_table.tsv created
echo.

REM ===== Step 3: Adjust SiteTable =====
echo ========================================
echo [Step 3/6] Adjusting SiteTable format
echo ========================================
echo.

echo [3a] Processing PhosSight...
python step3_AdjustSiteTable_17plex.py -i PhosSight_full_site_table.tsv -o PhosSight_site_table_Adjust_full.tsv -m PhosSight
if errorlevel 1 (
    echo ERROR: step3_AdjustSiteTable_PhosSight failed!
    goto :error
)
echo   Done: PhosSight_site_table_Adjust_full.tsv created
echo.

echo [3b] Processing PhosphoRS...
python step3_AdjustSiteTable_17plex.py -i PhosphoRS_full_site_table.tsv -o PhosphoRS_site_table_Adjust_full.tsv -m PhosphoRS
if errorlevel 1 (
    echo ERROR: step3_AdjustSiteTable_PhosphoRS failed!
    goto :error
)
echo   Done: PhosphoRS_site_table_Adjust_full.tsv created
echo.

echo [3c] Processing DeepRescore2...
python step3_AdjustSiteTable_17plex.py -i DeepRescore2_full_site_table.tsv -o DeepRescore2_site_table_Adjust_full.tsv -m DeepRescore2
if errorlevel 1 (
    echo ERROR: step3_AdjustSiteTable_DeepRescore2 failed!
    goto :error
)
echo   Done: DeepRescore2_site_table_Adjust_full.tsv created
echo.

REM ===== Step 4: Merge SiteTables =====
echo ========================================
echo [Step 4/6] Merging SiteTables
echo ========================================
echo.

echo [4a] Processing PhosSight...
python step4_MergeSiteTables_17plex.py -i PhosSight_site_table_Adjust_full.tsv -o PhosSight_Merged_full.tsv -m PhosSight
if errorlevel 1 (
    echo ERROR: step4_MergeSiteTables_PhosSight failed!
    goto :error
)
echo   Done: PhosSight_Merged_full.tsv created
echo.

echo [4b] Processing PhosphoRS...
python step4_MergeSiteTables_17plex.py -i PhosphoRS_site_table_Adjust_full.tsv -o PhosphoRS_Merged_full.tsv -m PhosphoRS
if errorlevel 1 (
    echo ERROR: step4_MergeSiteTables_PhosphoRS failed!
    goto :error
)
echo   Done: PhosphoRS_Merged_full.tsv created
echo.

echo [4c] Processing DeepRescore2...
python step4_MergeSiteTables_17plex.py -i DeepRescore2_site_table_Adjust_full.tsv -o DeepRescore2_Merged_full.tsv -m DeepRescore2
if errorlevel 1 (
    echo ERROR: step4_MergeSiteTables_DeepRescore2 failed!
    goto :error
)
echo   Done: DeepRescore2_Merged_full.tsv created
echo.

REM ===== Step 5: Get UniprotID and GeneName =====
echo ========================================
echo [Step 5/6] Getting UniprotID and GeneName
echo ========================================
echo.

echo [5a] Processing PhosSight...
python step5_GetUniprotIDGeneName_17plex.py -i PhosSight_Merged_full.tsv -o PhosSight_UniprotID_GeneName_full.tsv
if errorlevel 1 (
    echo ERROR: step5_GetUniprotIDGeneName_PhosSight failed!
    goto :error
)
echo   Done: PhosSight_UniprotID_GeneName_full.tsv created
echo.

echo [5b] Processing PhosphoRS...
python step5_GetUniprotIDGeneName_17plex.py -i PhosphoRS_Merged_full.tsv -o PhosphoRS_UniprotID_GeneName_full.tsv
if errorlevel 1 (
    echo ERROR: step5_GetUniprotIDGeneName_PhosphoRS failed!
    goto :error
)
echo   Done: PhosphoRS_UniprotID_GeneName_full.tsv created
echo.

echo [5c] Processing DeepRescore2...
python step5_GetUniprotIDGeneName_17plex.py -i DeepRescore2_Merged_full.tsv -o DeepRescore2_UniprotID_GeneName_full.tsv
if errorlevel 1 (
    echo ERROR: step5_GetUniprotIDGeneName_DeepRescore2 failed!
    goto :error
)
echo   Done: DeepRescore2_UniprotID_GeneName_full.tsv created
echo.

REM ===== Step 6: Analysis =====
echo ========================================
echo [Step 6/6] Running analysis
echo ========================================
echo.

call conda activate R_env
if errorlevel 1 (
    echo ERROR: Failed to activate R_env!
    goto :error
)

if not exist analysis mkdir analysis
cd analysis

echo [6a] Calculating missing value cutoff statistics...
Rscript Compare_PhosSight_vs_PhosphoRS_MissingValue_17plex.R
if errorlevel 1 (
    echo WARNING: Compare analysis may have warnings
)
echo   Done: MissingValueCutoff\CountData_PhosSight_vs_PhosphoRS_full.txt created
echo.

echo [6b] Plotting missing value curves...
Rscript Plot_PhosSight_vs_PhosphoRS_17plex.R
if errorlevel 1 (
    echo WARNING: Plotting may have warnings
)
echo   Done: MissingValueCutoff\Plot_PhosSight_vs_PhosphoRS_full.png created
echo.

cd ..

echo ========================================
echo Workflow completed successfully!
echo ========================================
echo.
echo Summary:
echo   - 17 TMT plexes: TMT01-TMT17
echo   - 170 TMT channels total (17 x 10)
echo   - 153 sample channels (excluding 17 reference channels)
echo   - Full SiteQuant with FASTA-based site localization
echo   - Three methods: PhosSight, PhosphoRS, DeepRescore2
echo.
echo Output files:
echo   Step 1: *_Intensity_full.txt
echo   Step 2: *_full_site_table.tsv (Full SiteQuant output)
echo   Step 3: *_site_table_Adjust_full.tsv
echo   Step 4: *_Merged_full.tsv
echo   Step 5: *_UniprotID_GeneName_full.tsv
echo   Step 6: analysis/MissingValueCutoff/*.png (Plots)
echo.

goto :end

:error
echo.
echo ========================================
echo ERROR: Workflow failed!
echo ========================================
exit /b 1

:end
endlocal
