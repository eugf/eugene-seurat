#!/usr/bin/env bash
#SBATCH --cpus-per-task=1
#SBATCH --mem=50g
#SBATCH --mail-type=BEGIN,END 
#SBATCH --time=1:00:00
#SBATCH --gres=lscratch:10

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# HOW TO RUN:
# - Open a terminal where this script file is and run:
# sbatch SEURAT-for_scNanoGPS-matrix+isoform+merged-argparse-sbatch.sh

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# # DEFINE VARIABLES FOR ARGPARSE -- change these!!!

# # LONG READS (MATRIX AND MATRIX_ISOFORM FOLDERS)
# INPUT_DIR='<...output/long_reads/...>'

# MERGED FOLDERS
# - (WARNING! disable the matrix and the isoform sections if using the MERGED section, and vice-versa!!!)
INPUT_DIR='<...output/merged/...>'

# OUTPUT FOLDER
# NOTE - I'm reusing the same input folder, but you can change it to wherever you want
OUTPUT_DIR=${INPUT_DIR}

# CELLTYPE MARKER TABLES
CELLTYPE_MARKER_DIR='/vf/users/CARD_singlecell/MONOCLE_V3/INPUTS/celltype_marker_tables_split'

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# LOGGING (to slurm)
# Output script name and the initial command-line arguments to slurm file 
echo "SCRIPT_PATH   : ${SCRIPT_PATH}"
echo "SCRIPT NAME   : $0"
echo "JOB INFO      : $@"

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# IMPORTS 
module load R/4.4

# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# # SCRIPT(S) TO RUN - based on which file(s) you want to use as input

# # MATRIX FILE

# SCRIPT_PATH='/vf/users/CARDPB/data/snRNA_longread/eugene-seurat/src/SEURAT-for_scNanoGPS-matrix-argparse.R'

# # FILTER PARAMETERS
# FILTER_NFEATURE_RNA_LOWER=0
# FILTER_NFEATURE_RNA_UPPER=200
# FILTER_PERCENT_MT=5

# Rscript \
#     ${SCRIPT_PATH} \
#     --input-dir ${INPUT_DIR} \
#     --output-dir ${OUTPUT_DIR} \
#     --celltype-marker-dir ${CELLTYPE_MARKER_DIR} \
#     --filter_nFeature_RNA_lower ${FILTER_NFEATURE_RNA_LOWER} \
#     --filter_nFeature_RNA_upper ${FILTER_NFEATURE_RNA_UPPER} \
#     --filter_percent_mt ${FILTER_PERCENT_MT} || (echo 'R failed!'; exit 1)

# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# # MATRIX ISOFORM FILE

# # SCRIPT_PATH='/vf/users/CARDPB/data/snRNA_longread/eugene-seurat/src/SEURAT-for_scNanoGPS-matrix_isoform-argparse.R'

# # FILTER PARAMETERS
# FILTER_NFEATURE_RNA_LOWER=0
# FILTER_NFEATURE_RNA_UPPER=250
# FILTER_PERCENT_MT=5

# Rscript \
#     ${SCRIPT_PATH} \
#     --input-dir ${INPUT_DIR} \
#     --output-dir ${OUTPUT_DIR} \
#     --celltype-marker-dir ${CELLTYPE_MARKER_DIR} \
#     --filter_nFeature_RNA_lower ${FILTER_NFEATURE_RNA_LOWER} \
#     --filter_nFeature_RNA_upper ${FILTER_NFEATURE_RNA_UPPER} \
#     --filter_percent_mt ${FILTER_PERCENT_MT} || (echo 'R failed!'; exit 1)
# 
# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# MERGED (MATRIX + MATRIX_ISOFORM)
# - (WARNING! disable the matrix and the isoform sections if using the MERGED section, and vice-versa!!!)

SCRIPT_PATH='/vf/users/CARDPB/data/snRNA_longread/eugene-seurat/src/SEURAT-for_scNanoGPS-merged-argparse.R'

# FILTER PARAMETERS
FILTER_NFEATURE_RNA_LOWER=0
FILTER_NFEATURE_RNA_UPPER=2500
FILTER_PERCENT_MT=5

Rscript \
    ${SCRIPT_PATH} \
    --input-dir ${INPUT_DIR} \
    --output-dir ${OUTPUT_DIR} \
    --celltype-marker-dir ${CELLTYPE_MARKER_DIR} \
    --filter_nFeature_RNA_lower ${FILTER_NFEATURE_RNA_LOWER} \
    --filter_nFeature_RNA_upper ${FILTER_NFEATURE_RNA_UPPER} \
    --filter_percent_mt ${FILTER_PERCENT_MT} || (echo 'R failed!'; exit 1)
