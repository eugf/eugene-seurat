#!/usr/bin/env bash
#SBATCH --cpus-per-task=1
#SBATCH --mem=100g
#SBATCH --mail-type=BEGIN,END 
#SBATCH --time=3:00:00
#SBATCH --gres=lscratch:25

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# HOW TO RUN:
# - Open a terminal where this script file is and run:
# sbatch CONVERSION-scNanoGPS_matrix_and_matrix_isoform_to_barcodes_genes_mtx-argparse-sbatch.sh 

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# DEFINE VARIABLES FOR ARGPARSE -- change these!!!
SCRIPT_PATH='/vf/users/CARDPB/data/snRNA_longread/eugene-seurat/src/CONVERSION-scNanoGPS_matrix_and_matrix_isoform_to_barcodes_genes_mtx-argparse.py'

# INPUT PATHS
INPUT_PATH='<.../scNanoGPS/scNanoGPS_res/...>'

# OUTPUT PATH
OUTPUT_PATH='../output/long_read'

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# LOGGING (to slurm)
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