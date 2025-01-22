#!/usr/bin/env bash
#SBATCH --cpus-per-task=1
#SBATCH --mem=50g
#SBATCH --mail-type=BEGIN,END 
#SBATCH --time=1:00:00
#SBATCH --gres=lscratch:10

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# HOW TO RUN:
# - Open a terminal where this script file is and run:
# sbatch SEURAT-for_short_reads-argparse-sbatch.sh

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# DEFINE VARIABLES FOR ARGPARSE -- change these!!!
SCRIPT_PATH='/vf/users/CARDPB/data/snRNA_longread/eugene-seurat/src/SEURAT-for_short_reads-argparse.R'

# # SHORT READS
# INPUT_DIR='<...output/merged/...>'
INPUT_DIR='/data/CARD_singlecell/Brain_atlas/Cortex_opt/Psomagen102722/4546Cortex/outs/filtered_feature_bc_matrix.h5'

# OUTPUT FOLDER
# NOTE - I'm reusing the same input folder, but you can change it to wherever you want
OUTPUT_DIR='/data/CARDPB/data/snRNA_longread/eugene-seurat/output/short_reads'

# CELLTYPE MARKER TABLES
CELLTYPE_MARKER_DIR='/vf/users/CARD_singlecell/MONOCLE_V3/INPUTS/celltype_marker_tables_split'

# ID -- required for user to input for the short reads!!!
ID='UMARY_4546'

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# LOGGING (to slurm)
# Output script name and the initial command-line arguments to slurm file 
echo "SCRIPT_PATH   : ${SCRIPT_PATH}"
echo "SCRIPT NAME   : $0"
echo "JOB INFO      : $@"

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# IMPORTS 
module load R/4.4

# RUN THIS SCRIPT
Rscript \
    ${SCRIPT_PATH} \
    --input-dir ${INPUT_DIR} \
    --output-dir ${OUTPUT_DIR} \
    --celltype-marker-dir ${CELLTYPE_MARKER_DIR} \
    --id ${ID} || (echo 'R failed!'; exit 1)
