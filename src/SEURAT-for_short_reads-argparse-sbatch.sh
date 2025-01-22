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
INPUT_DIR='<...output/merged/...>'

# OUTPUT FOLDER
# NOTE - I'm reusing the same input folder, but you can change it to wherever you want
OUTPUT_DIR='../output/short_reads'

# CELLTYPE MARKER TABLES
CELLTYPE_MARKER_DIR='/vf/users/CARD_singlecell/MONOCLE_V3/INPUTS/celltype_marker_tables_split'

# ID -- required for user to input for the short reads!!!
ID=<'SAMPLE_ID#'>

# FILTER PARAMETERS
# - These are actually the default values that's already there and is technically optional for me to enter, but I'm leaving them here as an example to remind people to fill this out, you'll need to adjust it later anyway
FILTER_NFEATURE_RNA_LOWER=0
FILTER_NFEATURE_RNA_UPPER=12000
FILTER_PERCENT_MT=5

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
    --id ${ID} \
    --filter_nFeature_RNA_lower ${FILTER_NFEATURE_RNA_LOWER} \
    --filter_nFeature_RNA_upper ${FILTER_NFEATURE_RNA_UPPER} \
    --filter_percent_mt ${FILTER_PERCENT_MT} || (echo 'R failed!'; exit 1)
