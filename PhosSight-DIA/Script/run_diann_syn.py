import subprocess
import os
import argparse
from pathlib import Path

# TODO: delete these comments
# root_dir = Path("~/PhosSight_analysis").expanduser()
# diann_exe = Path("singularity exec ~/diann-2.2.0/diann-2.2.0.img ~/diann-2.2.0/diann-linux")

###########################################################################
# Batch run DIA-NN with different spectral library filtering parameters
# Filtered by pretrained model
###########################################################################
def batch_run_diann_three_runs(
    output_dir: Path,
    raw_dir: Path,
    spec_lib_dir: Path,
    diann_cmd_prefix: str
):
    condense_ratio_list = [1, 2, 10]

    for cond_factor in condense_ratio_list:
        final_output_dir = output_dir / "for_finetuning_res" / f"{cond_factor}xrep1"
        if not os.path.exists(final_output_dir):
            os.makedirs(final_output_dir)
            
        raw_path = raw_dir / f"PhosphoPooledSynth{cond_factor}xdY_DIA60K10Da30K500_1000Da_160min50cm_19031101.raw"
        spec_lib_parquet_path = spec_lib_dir / "spec_library_original.parquet"
        output_path = final_output_dir / "report.parquet"
        param = f'--f "{raw_path}" --lib "{spec_lib_parquet_path}" --threads 32 --verbose 1 --out "{output_path}" --qvalue 0.01 --matrices  --temp "{final_output_dir}"  --unimod4 --var-mods 1 --var-mod UniMod:21,79.966331,STY --window 7 --mass-acc 8 --mass-acc-ms1 10 --peptidoforms --reanalyse --rt-profiling'
        cmd = f"{diann_cmd_prefix} {param}"
        print(f"Running command: {cmd}")
        result = subprocess.run(cmd, shell=True, capture_output=False, text=True)
        print("Return code:", result.returncode)
        
        
        
###########################################################################
# Run DIA-NN with original spectral library without any filtering
###########################################################################
def batch_run_diann_original(
    output_dir: Path,
    raw_dir: Path,
    spec_lib_dir: Path,
    diann_cmd_prefix: str
):
    final_output_dir = output_dir / "original"
    if not os.path.exists(final_output_dir):
        os.makedirs(final_output_dir)
    
    spec_lib_parquet_path = spec_lib_dir / "spec_library_original.parquet"
    output_path = final_output_dir / "report.parquet"
    
    # Build paths for all raw data files
    raw_files = []
    for concentration in ["1x", "2x", "4x", "10x", "20x"]:
        for replicate in ["1", "2", "3"]:
            raw_file = raw_dir / f"PhosphoPooledSynth{concentration}dY_DIA60K10Da30K500_1000Da_160min50cm_1903110{replicate}.raw"
            raw_files.append(f'--f "{raw_file}"')
    
    raw_files_param = " ".join(raw_files)
    
    param = f'{raw_files_param} --lib "{spec_lib_parquet_path}" --threads 32 --verbose 1 --out "{output_path}" --qvalue 0.01 --matrices  --temp "{final_output_dir}"  --unimod4 --var-mods 1 --var-mod UniMod:21,79.966331,STY --window 7 --mass-acc 8 --mass-acc-ms1 10 --peptidoforms --reanalyse --rt-profiling'
    cmd = f"{diann_cmd_prefix} {param}"
    print(f"Running command: {cmd}")
    result = subprocess.run(cmd, shell=True, capture_output=False, text=True)
    print("Return code:", result.returncode)



###########################################################################
# Batch run DIA-NN with different spectral library filtering parameters
# Filtered by pretrained model
###########################################################################
def batch_run_diann_pretrained_filtered(
    output_dir: Path,
    raw_dir: Path,
    spec_lib_dir: Path,
    diann_cmd_prefix: str
):
    ratio_list = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    ratio_suffix = [f'ratio_{r}' for r in ratio_list]
    suffix = ratio_suffix

    for suf in suffix:
        final_output_dir = output_dir / f"pretrained_{suf}"
        if not os.path.exists(final_output_dir):
            os.makedirs(final_output_dir)
        
        spec_lib_parquet_path = spec_lib_dir / "pretrained_filtered" / f"spec_library_filtered_by_{suf}.parquet"
        output_path = final_output_dir / "report.parquet"
        
        # Build paths for all raw data files
        raw_files = []
        for concentration in ["1x", "2x", "4x", "10x", "20x"]:
            for replicate in ["1", "2", "3"]:
                raw_file = raw_dir / f"PhosphoPooledSynth{concentration}dY_DIA60K10Da30K500_1000Da_160min50cm_1903110{replicate}.raw"
                raw_files.append(f'--f "{raw_file}"')
        
        raw_files_param = " ".join(raw_files)
        
        param = f'{raw_files_param} --lib "{spec_lib_parquet_path}" --threads 32 --verbose 1 --out "{output_path}" --qvalue 0.01 --matrices  --temp "{final_output_dir}"  --unimod4 --var-mods 1 --var-mod UniMod:21,79.966331,STY --window 7 --mass-acc 8 --mass-acc-ms1 10 --peptidoforms --reanalyse --rt-profiling'
        cmd = f"{diann_cmd_prefix} {param}"
        print(f"Running command: {cmd}")
        result = subprocess.run(cmd, shell=True, capture_output=False, text=True)
        print("Return code:", result.returncode)
    
    

###########################################################################
# Batch run DIA-NN with different spectral library filtering parameters
# Filtered by fine-tuned model
###########################################################################
def batch_run_diann_finetuned_filtered(
    output_dir: Path,
    raw_dir: Path,
    spec_lib_dir: Path,
    diann_cmd_prefix: str
):
    ratio_list = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    ratio_suffix = [f'ratio_{r}' for r in ratio_list]
    suffix = ratio_suffix
    for suf in suffix:
        final_output_dir = output_dir / f"finetuned_{suf}"
        if not os.path.exists(final_output_dir):
            os.makedirs(final_output_dir)
        
        spec_lib_parquet_path = spec_lib_dir / "finetuned_filtered" / f"spec_library_filtered_by_{suf}.parquet"
        output_path = final_output_dir / "report.parquet"
        
        # Build paths for all raw data files
        raw_files = []
        for concentration in ["1x", "2x", "4x", "10x", "20x"]:
            for replicate in ["1", "2", "3"]:
                raw_file = raw_dir / f"PhosphoPooledSynth{concentration}dY_DIA60K10Da30K500_1000Da_160min50cm_1903110{replicate}.raw"
                raw_files.append(f'--f "{raw_file}"')
        
        raw_files_param = " ".join(raw_files)
        
        param = f'{raw_files_param} --lib "{spec_lib_parquet_path}" --threads 32 --verbose 1 --out "{output_path}" --qvalue 0.01 --matrices  --temp "{final_output_dir}"  --unimod4 --var-mods 1 --var-mod UniMod:21,79.966331,STY --window 7 --mass-acc 8 --mass-acc-ms1 10 --peptidoforms --reanalyse --rt-profiling'
        cmd = f"{diann_cmd_prefix} {param}"
        print(f"Running command: {cmd}")
        result = subprocess.run(cmd, shell=True, capture_output=False, text=True)
        print("Return code:", result.returncode)
        
        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run DIA-NN batch analysis')
    parser.add_argument('--step', type=str, choices=['three_rep', 'original', 'pretrained', 'finetuned'], 
                        required=True, help='Choose the step to run: 1 (three runs), 2 (pretrained filtering), 3 (fine-tuned filtering)')
    parser.add_argument('--diann-cmd-prefix', type=str, required=True,
                        help='The command prefix to run DIA-NN (inside Singularity container)')
    parser.add_argument('--output_dir', type=str, required=True,
                        help='The output directory for analysis results')
    parser.add_argument('--raw_dir', type=str, required=True,
                        help='The directory for raw data files')
    parser.add_argument('--spec_lib_dir', type=str, required=True,
                        help='The directory for spectral library parquet files')
    
    args = parser.parse_args()
    
    # Convert paths to Path objects and expand ~
    output_dir = Path(args.output_dir).expanduser()
    raw_dir = Path(args.raw_dir).expanduser()
    spec_lib_dir = Path(args.spec_lib_dir).expanduser()
    singularity_diann_cmd_prefix = args.diann_cmd_prefix

    if args.step == 'three_rep':
        print("Running three DIA-NN analyses...")
        batch_run_diann_three_runs(output_dir, raw_dir, spec_lib_dir, singularity_diann_cmd_prefix)
    elif args.step == 'original':
        print("Running DIA-NN analysis with original spectral library...")
        batch_run_diann_original(output_dir, raw_dir, spec_lib_dir, singularity_diann_cmd_prefix)
    elif args.step == 'pretrained':
        print("Running DIA-NN analysis with pretrained model filtering...")
        batch_run_diann_pretrained_filtered(output_dir, raw_dir, spec_lib_dir, singularity_diann_cmd_prefix)
    elif args.step == 'finetuned':
        print("Running DIA-NN analysis with fine-tuned model filtering...")
        batch_run_diann_finetuned_filtered(output_dir, raw_dir, spec_lib_dir, singularity_diann_cmd_prefix)
    
    print("Analysis completed!")
