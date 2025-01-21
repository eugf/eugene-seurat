#!/usr/bin/env bash
#SBATCH --cpus-per-task=1
#SBATCH --mem=50g
#SBATCH --mail-type=BEGIN,END 
#SBATCH --time=1:00:00
#SBATCH --gres=lscratch:10

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# HOW TO RUN:
# - NOTE: almost all of the user-inputted parameters are saved here as variables you can change and it will send that to the Rscript file
# - Open a terminal where this script file is and run:
# sbatch SEURAT-for_scNanoGPS-matrix+isoform+merged-argparse-sbatch.sh

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# # DEFINE VARIABLES FOR ARGPARSE -- change these!!!

# # LONG READS (MATRIX AND MATRIX_ISOFORM FOLDERS)
# INPUT_DIR='/data/CARDPB/data/snRNA_longread/eugene-seurat/output/long_reads/SH-04-08/20241107'
# INPUT_DIR='/data/CARDPB/data/snRNA_longread/eugene-seurat/output/long_reads/SH-06-25/20241107'
# INPUT_DIR='/data/CARDPB/data/snRNA_longread/eugene-seurat/output/long_reads/SH-07-46/20241107'
# INPUT_DIR='/data/CARDPB/data/snRNA_longread/eugene-seurat/output/long_reads/SH-92-05/20241122'
# INPUT_DIR='/data/CARDPB/data/snRNA_longread/eugene-seurat/output/long_reads/UMARY_4546/20241106'

# # MERGED FOLDERS
# - (WARNING! disable the matrix and the isoform sections then!!!)
# INPUT_DIR='/vf/users/CARDPB/data/snRNA_longread/eugene-seurat/output/merged/SH-04-08/20241219-1'
# INPUT_DIR='/vf/users/CARDPB/data/snRNA_longread/eugene-seurat/output/merged/SH-06-25/20241220'
# INPUT_DIR='/vf/users/CARDPB/data/snRNA_longread/eugene-seurat/output/merged/SH-07-46/20241220'
# INPUT_DIR='/vf/users/CARDPB/data/snRNA_longread/eugene-seurat/output/merged/SH-92-05/20241220'
# INPUT_DIR='/vf/users/CARDPB/data/snRNA_longread/eugene-seurat/output/merged/UMARY_4546/20241219-0'

# NOTE - I'm reusing the same input folder, but you can change it to wherever you want
OUTPUT_DIR=${INPUT_DIR}

CELLTYPE_MARKER_DIR='/vf/users/CARD_singlecell/MONOCLE_V3/INPUTS/celltype_marker_tables_split'

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# LOGGING
# Output script name and the initial command-line arguments to slurm file 
echo "SCRIPT_PATH   : ${SCRIPT_PATH}"
echo "SCRIPT NAME   : $0"
echo "JOB INFO      : $@"

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# IMPORTS 
module load R/4.4

# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# # SCRIPT(S) TO RUN - based on which file you want to use as input

# # MATRIX FILE

# SCRIPT_PATH='/vf/users/CARDPB/data/snRNA_longread/eugene-seurat/src/SEURAT-for_scNanoGPS-matrix-argparse.R'

# Rscript \
#     ${SCRIPT_PATH} \
#     --input-dir ${INPUT_DIR} \
#     --output-dir ${OUTPUT_DIR} \
#     --celltype-marker-dir ${CELLTYPE_MARKER_DIR} || (echo 'R failed!'; exit 1)

# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# # MATRIX ISOFORM FILE

# SCRIPT_PATH='/vf/users/CARDPB/data/snRNA_longread/eugene-seurat/src/SEURAT-for_scNanoGPS-matrix_isoform-argparse.R'

# Rscript \
#     ${SCRIPT_PATH} \
#     --input-dir ${INPUT_DIR} \
#     --output-dir ${OUTPUT_DIR} \
#     --celltype-marker-dir ${CELLTYPE_MARKER_DIR} || (echo 'R failed!'; exit 1)

# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# MERGED (MATRIX + MATRIX_ISOFORM)
# - (WARNING! disable the matrix and the isoform sections then!!!)

SCRIPT_PATH='/vf/users/CARDPB/data/snRNA_longread/eugene-seurat/src/SEURAT-for_scNanoGPS-merged-argparse.R'

Rscript \
    ${SCRIPT_PATH} \
    --input-dir ${INPUT_DIR} \
    --output-dir ${OUTPUT_DIR} \
    --celltype-marker-dir ${CELLTYPE_MARKER_DIR} || (echo 'R failed!'; exit 1)
