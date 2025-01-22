#!/usr/bin/env python
# coding: utf-8

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# # GOAL(S)
# 
# - Convert the outputs from scNanoGPS into a matrix file format that'll be accepted into Seurat
# - Input files needed:
#   - matrix.tsv
#   - matrix_isoform.tsv
# 
# - Output files created:
#   - barcodes.tsv
#   - genes.tsv
#   - matrix.mtx --> in sparse matrix format

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# # HOW TO RUN (for converting to a merged matrix + matrix_isoform file)
# # - Open a terminal where this script file is and run:
# python CONVERSION-scNanoGPS-matrix_and_matrix_isoform_to_barcodes_genes_mtx-argparse.py --input_path <INPUT_PATH> --output_path <OUTPUT_PATH>

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# IMPORTS
import pandas as pd
import os
import numpy as np
from scipy import sparse
from scipy.io import mmwrite
import re
import datetime
import argparse

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ARGUMENT PARSER 
# Create ArgumentParser object to take user inputs
parser = argparse.ArgumentParser(description = 'Process arguments needed to run conversions from scNanoGPS outputs')

# Define arguments for the user to supply
parser.add_argument(
    '--input_path', 
    help        = 'Should be the sample folders inside the /scNanoGPS_res folder ',
    required    = True,
    type        = str,
    )

parser.add_argument(
    '--output_path', 
    help        = 'Base folder for processed files and plots to be saved',
    required    = True,
    type        = str,
    )

# Parse arguments
args = parser.parse_args()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# LOAD DATA
# Input paths - look inside the scNanoGPS output folder (`/scNanoGPS_res`) by patient ID until you reach the parent folder that holds the matrix.tsv + matrix_isoform.tsv files
scNanoGPS_output_dir        = args.input_path
file_path_matrix            = os.path.join(scNanoGPS_output_dir, 'matrix.tsv')
file_path_matrix_isoform    = os.path.join(scNanoGPS_output_dir, 'matrix_isoform.tsv')

# Output dir
output_dir = args.output_path

# PROCESSING
# Extract the patient ID
id = os.path.basename(os.path.dirname(os.path.dirname(scNanoGPS_output_dir)))
separators = ['_FTX_singlecell', '_singlecell']
split_id = re.split('|'.join(separators), id)[0]

# Get today's date
current_date = datetime.date.today()
formatted_date = current_date.strftime("%Y%m%d")

# Define unique output folders based on the patient ID and date
sample_dir = os.path.join(output_dir, split_id, formatted_date)

# CHECK - if the path exists already, adjust name as needed
# Start a counter for repeated folders
counter = 0
if os.path.exists(sample_dir):
    print(f"The directory '{sample_dir}' exists.")
    new_sample_dir = sample_dir + '-' + str(counter)
    
    # Increment on the date folder
    while os.path.exists(new_sample_dir):
        counter += 1
        print('counter =', counter)
        new_sample_dir = sample_dir + '-' + str(counter)
        print('new_sample_dir INSIDE IF BLOCK =', new_sample_dir)
    
    # Reset this variable to call later
    sample_dir = new_sample_dir

else:
    print(f"Sample directory '{sample_dir}' does not exist")

# Make the new directory 
os.makedirs(sample_dir)

# Save full ID to text file (will process this way downstream)
output_id_txt = os.path.join(sample_dir, 'id.txt')
# Open file to write and save ID
with open(output_id_txt, 'w+') as file:
    file.write(id)

# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # MATRIX FILE
# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# LOAD - matrix file into a Pandas DF
df = pd.read_csv(
    file_path_matrix, 
    sep = '\t',     # For TSV file
    skiprows = 0,   # 1st row is metadata, skip it (index = 0)
    header = 1      # Header begins on the 2nd row (index = 1)
    )

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# CONVERSIONS

# ### Create the `genes.tsv` file 
# 
# - Drop the first 2 rows (already done during import)
# - Extract 2 columns ONLY: 
#   - Geneid (ENSEMBLE ID #s) and
#   - gene_name

# Select your desired cols
df_gene_id_and_gene_name = df[['Geneid', 'gene_name']]

# ### EXPORT - `genes.tsv`

# Define output path
file_path_genes_tsv = os.path.join(sample_dir, 'genes.tsv')

# Save file as `genes.tsv`
df_gene_id_and_gene_name.to_csv(
    file_path_genes_tsv, 
    sep = '\t',         # For TSV file
    index = False,      # Remove index col
    header = False,     # Remove header row
    )
print(f'Saved genes.tsv to: {file_path_genes_tsv}')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ### Create the `barcodes.tsv` file
# 
# - Load just the 2nd row of the original `matrix.tsv` file (which is now the header row)
# - Extract all the barcodes, drop the other columns

# RELOAD - matrix file into a Pandas DF
df_header = pd.read_csv(
    file_path_matrix, 
    sep = '\t',     # For TSV file
    skiprows = 1,   # 1st row is metadata, skip it (index = 0)
    nrows = 1,      # Only load the next 1 row (the header row)
    header = None   # Prevents Pandas defaulting to header = True
    )

# Converts to a 1D series
header_series = df_header.squeeze()

# Remove the initial columns from the header, keep the barcodes ONLY
header_series_subset = header_series[7:]

# Change the formatting to match examples and our data
header_series_subset_mod = header_series_subset + '-1'

# ### EXPORT - `barcodes.tsv`

# Define output path
file_path_barcodes_tsv = os.path.join(sample_dir, 'barcodes.tsv')

# Save file as `barcodes.tsv`
header_series_subset_mod.to_csv(
    file_path_barcodes_tsv, 
    sep = '\t',         # For TSV file
    index = False,      # No index col
    header = False,     # Prevents Pandas defaulting to header = True
    )
print(f'Saved barcodes.tsv to: {file_path_barcodes_tsv}')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ### Create the `matrix.mtx` file
# 
# - RELOAD `matrix.tsv` as a DF, 
# - Get rid of cols from `chr` to `gene_name`, just keep `geneid` and then all the counts cols, 
# - then convert to a sparse matrix format (needs matrix library to convert to DCG matrix)

# RELOAD - matrix file into a Pandas DF (but this time get rid of the header)
df_full = pd.read_csv(
    file_path_matrix, 
    sep = '\t',         # For TSV file
    skiprows = [0,1],   # 1st row is metadata, 2nd is header, skip BOTH (index = 0,1),
    header = None,      # Prevents Pandas defaulting to header = True
    )

# Now using the "full" DF, define the corresponding col's we want to drop
cols_to_drop = df_full.columns[00:7]
# APPLY - Drop to those cols in the full DF
df_full_drop = df_full.drop(cols_to_drop, axis = 1)

# ### EXPORT - processed DF to `matrix_dropped.tsv` file 

# Define output path
file_path_matrix_dropped_tsv = os.path.join(sample_dir, 'matrix_dropped.tsv')

# Save file as `matrix_dropped.tsv`
df_full_drop.to_csv(
    file_path_matrix_dropped_tsv, 
    sep = '\t', 
    index = False, 
    header = False,
    )
print(f'Saved matrix_dropped.tsv to: {file_path_matrix_dropped_tsv}')

# ### CONVERSION - function to convert the .TSV file into a sparse matrix .MTX file

def tsv_to_mtx(input_file, output_file):
    
    # Read the TSV file
    df = pd.read_csv(
        input_file, 
        sep = '\t',         # For TSV file
        index_col = None,   # No index col
        header = None,      # Prevents Pandas defaulting to header = True
        )
    
    # CONVERT - DF to numpy array
    data = df.to_numpy()
    
    # CONVERT - numpy array to sparse matrix
    sparse_matrix = sparse.csr_matrix(data)
    
    # EXPORT - save the sparse matrix to a .mtx file
    mmwrite(output_file, sparse_matrix)
    
    # CHECK - saved output file to
    print(f'Saved matrix.mtx file to: {output_file}')

# ### EXPORT - `matrix.mtx` file as a sparse matrix

# Define output path
output_file = os.path.join(sample_dir, 'matrix.mtx')

# Run conversion function
tsv_to_mtx(file_path_matrix_dropped_tsv, output_file)

# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # MATRIX ISOFORM FILE
# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# # LOAD DATA

# DEFINE FILE PATHS
# Output dir - new isoform folder inside the sample directory
output_isoform_dir = os.path.join(sample_dir, 'isoform')
os.makedirs(output_isoform_dir)

# LOAD - matrix isoform file into a Pandas DF
df_isoform = pd.read_csv(
    file_path_matrix_isoform, 
    sep = '\t',     # For TSV file
    )

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ### NOTES
#
# Cols are:
# - `Geneid` (ENS #)
# - MISSING `gene_name`, which should be derived from the isoform's parent?
# - The rest are barcodes
# - NO index col
# 
# Rows:
# - Header row

# # CONVERSION
# 
# ### Create the `genes.tsv` file 
# 
# - Rename col 1 (currently unnamed)
# - Extract 2 columns: 
#   - `Geneid` (Ensemble #s)  
#   - `gene_name` --> MISSING! need to create! but how???
# - Drop the 1st row (header row)

# RENAME - missing header col
df_isoform = df_isoform.rename(columns = {'Unnamed: 0': 'Geneid'})

# CHECK - length before splitting the col names
pre_split = len(df_isoform['Geneid'].unique())

# ### EXPORT - CHECK

# EXPORT - the file before splitting the col names
df_isoform.to_csv(
    os.path.join(output_isoform_dir, 'isoform_before_splitting.tsv'), 
    sep = '\t', 
    index = False, 
    header = True
    )

# NOTES: This col's format:
# 
# - <gene_name>_ENSEMBLE###.#

# Split the `Geneid` column on '_ENST' then add back 'ENST'
# Keep the 0th element to create a new `gene_name` column
df_isoform['gene_name'] = df_isoform['Geneid'].str.split('_ENST').str[0]
# Rename the `Geneid` col using the 1st element of the split and prepending it with 'ENST' to restore the lost characters from the split
df_isoform['Geneid'] = 'ENST' + df_isoform['Geneid'].str.split('_ENST').str[1]

# Select your desired cols
df_isoform_gene_id_and_gene_name = df_isoform[['Geneid', 'gene_name']]

# ### EXPORT - `genes.tsv`

# NOTE: overwriting the old genes.tsv path for matrix.tsv with the matrix_isoform.tsv path
# Define output path
file_path_genes_tsv = os.path.join(output_isoform_dir, 'genes.tsv')

# Save file as `genes.tsv`
df_isoform_gene_id_and_gene_name.to_csv(
    file_path_genes_tsv, 
    sep = '\t',         # For TSV file
    index = False,      # Remove index col???
    header = False,     # Remove header row???
    )
print(f'Saved genes.tsv to: {file_path_genes_tsv}')

# RELOAD - matrix isoform file into a Pandas DF
df_isoform_before_splitting = pd.read_csv(
    os.path.join(output_isoform_dir, 'isoform_before_splitting.tsv'),
    sep = '\t',     # For TSV file
    )

# Select your desired cols
df_isoform_before_splitting_gene_id = df_isoform_before_splitting[['Geneid']]

# OVERWRITE - the genes.tsv file with the original, unchanged, fused column name, but we'll just call it `Geneid`
# Define output path
file_path_genes_tsv = os.path.join(output_isoform_dir, 'genes.tsv')

# EXPORT - Save file as `genes.tsv`
df_isoform_before_splitting_gene_id.to_csv(
    file_path_genes_tsv,
    sep = '\t',         # For TSV file
    index = False,      # Remove index col???
    header = False,     # Remove header row???
    )
print(f'Saved genes.tsv to: {file_path_genes_tsv}')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ### Create the `barcodes.tsv` file
# 
# - Get the 1st row of the `matrix_isoform.tsv` file
# - Drop the 1st col (ENS #s)
# - NOTE: 1st col is empty and there's no other metadata cols

# LOAD - matrix isoform file into a Pandas DF to get the header row only
df_isoform_header = pd.read_csv(
    file_path_matrix_isoform, 
    sep = '\t',         # For TSV file
    nrows = 1,          # Only load the 1st row (the header row)
    header = None,      # Header begins on the 1st row
    )

# Drop the 1st col (the unnamed Geneid col here)
df_isoform_header_skip = df_isoform_header.iloc[:, 1:]

# Convert to a 1D series
isoform_header_series = df_isoform_header_skip.squeeze()

# Change the formatting to match examples and our data
isoform_header_series_mod = isoform_header_series + '-1'

# ### EXPORT - `barcodes.tsv`

# Define output path
file_path_barcodes_tsv = os.path.join(output_isoform_dir, 'barcodes.tsv')

# Save file as `genes.tsv`
isoform_header_series_mod.to_csv(
    file_path_barcodes_tsv, 
    sep = '\t',         # For TSV file
    index = False,      # Remove index col???
    header = False,     # Remove header row???
    )
print(f'Saved barcodes.tsv to: {file_path_barcodes_tsv}')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ### Create the `matrix.mtx` file
# 
# - Import `matrix_isoform.tsv` as a DF, 
# - Keep `geneid` and then all the counts cols, drop the 1st row (header row) tho
# - then convert to a sparse matrix format (needs matrix library to convert to DCG matrix)

# RELOAD - matrix file into a Pandas DF (but this time get rid of the header)
df_isoform_full = pd.read_csv(
    file_path_matrix_isoform, 
    sep = '\t',         # For TSV file
    skiprows = [0],     # 1st row is metadata (index = 0),
    header = None,      # Prevents Pandas defaulting to header = True
    )

# Use the "full" DF
cols_to_drop_isoform = df_isoform_full.columns[00:1]

# APPLY - Drop to those cols in the full DF
df_isoform_full_drop = df_isoform_full.drop(cols_to_drop_isoform, axis = 1)

# ### EXPORT - processed DF to `matrix_isoform_dropped.tsv` file

# Define output path
file_path_matrix_isoform_dropped_tsv = os.path.join(output_isoform_dir, 'matrix_dropped.tsv')

# EXPORT - save as .TSV file
df_isoform_full_drop.to_csv(
    file_path_matrix_isoform_dropped_tsv, 
    sep = '\t',         # For TSV file
    index = False,      # Remove index col
    header = False,     # Remove header row
    )
print(f'Saved matrix_isoform_dropped.tsv to: {file_path_matrix_isoform_dropped_tsv}')

# # CONVERSION - function to convert the TSV file a sparse matrix .MTX file
# 
# ### EXPORT - `matrix.mtx` file as a sparse matrix

# Define output path
output_file = os.path.join(output_isoform_dir, 'matrix.mtx')

# Run conversion function
tsv_to_mtx(file_path_matrix_isoform_dropped_tsv, output_file)
