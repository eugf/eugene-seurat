#!/usr/bin/env bash
#SBATCH --cpus-per-task=1
#SBATCH --mem=100g
#SBATCH --mail-type=BEGIN,END 
#SBATCH --time=3:00:00
#SBATCH --gres=lscratch:25

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# HOW TO RUN:
# - Open a terminal where this script file is and run:
# sbatch CONVERSION-scNanoGPS-merge_matrix_and_matrix_isoform-argparse-sbatch.sh

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# DEFINE VARIABLES FOR ARGPARSE - change these!!!
SCRIPT_PATH='/vf/users/CARDPB/data/snRNA_longread/eugene-seurat/src/CONVERSION-scNanoGPS-merge_matrix_and_matrix_isoform-argparse.py'

# # INPUT PATHS
# INPUT_PATH='/data/CARDPB/data/snRNA_longread/scNanoGPS-neuro/scNanoGPS_res/SH-04-08_singlecell_LR/20240321-1/20240321'
# INPUT_PATH='/data/CARDPB/data/snRNA_longread/scNanoGPS-neuro/scNanoGPS_res/SH-06-25_singlecell_LR/20241105-2/20241105'
# INPUT_PATH='/data/CARDPB/data/snRNA_longread/scNanoGPS-neuro/scNanoGPS_res/SH-07-46_singlecell_LR/20241105-1/20241105'
# INPUT_PATH='/data/CARDPB/data/snRNA_longread/scNanoGPS-neuro/scNanoGPS_res/SH-92-05_singlecell_LR/20240314-1/20240314'
# INPUT_PATH='/data/CARDPB/data/snRNA_longread/scNanoGPS-neuro/scNanoGPS_res/UMARY_4546_FTX_singlecell_LR/20240321-1/20240321'

# OUTPUT PATH
OUTPUT_PATH='/data/CARDPB/data/snRNA_longread/eugene-seurat/output/merged'

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Output script name and the initial command-line arguments to slurm file 
echo "SCRIPT_PATH   : ${SCRIPT_PATH}"
echo "SCRIPT NAME   : $0"
echo "JOB INFO      : $@"

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# IMPORTS
module load python/3.10

# RUN THIS SCRIPT (for converting to a merged matrix + matrix_isoform file)
python ${SCRIPT_PATH} \
    --input_path ${INPUT_PATH} \
    --output_path ${OUTPUT_PATH} || (echo 'Python failed!'; exit 1)