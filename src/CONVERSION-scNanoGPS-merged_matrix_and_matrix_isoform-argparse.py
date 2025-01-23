#!/usr/bin/env python
# coding: utf-8

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# # GOAL(S)
# - Merge the resulting scNanoGPS matrix and matrix_isoform files together
# 
# - Input files needed:
#   - matrix.tsv
#   - matrix_isoform.tsv
# 
# - Output files created:
#   - merge_matrix_and_matrix_isoform.tsv
# 

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# # HOW TO RUN (for converting to a merged matrix + matrix_isoform file)
# # - Open a terminal where this script file is and run:
# module load python/3.10
# python CONVERSION-scNanoGPS-merged_matrix_and_matrix_isoform-argparse.py --input_path <INPUT_PATH> --output_path <OUTPUT_PATH>

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# IMPORTS
import pandas as pd
import os
import re
import datetime
import numpy as np
from scipy import sparse
from scipy.io import mmwrite
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

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

# OUTPUT PATHS
# Save file as `merged_matrix_and_matrix_isoform.tsv`
file_path_merged_matrix_and_matrix_isoform = os.path.join(sample_dir, 'merged_matrix_and_matrix_isoform.tsv')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# LOAD - matrix file into a Pandas DF
df_matrix = pd.read_csv(
    file_path_matrix, 
    sep = '\t',     # For TSV file
    skiprows = 0,   # 1st row is metadata, skip it (index = 0)
    header = 1      # Header begins on the 2nd row (index = 1)
    )

# ### Steps overview
# - KEEP: Geneid, gene_name, and the rest of the barcodes
# - DROP: in betw cols
# - Then merge on the barcodes

# Drop these cols from the DFs
df_matrix_drop = df_matrix.drop(
    columns = ['Chr', 'Start', 'End', 'Strand', 'Length'], 
    inplace = False
    )

# LOAD - matrix isoform file into a Pandas DF
df_isoform = pd.read_csv(
    file_path_matrix_isoform, 
    sep = '\t',     # For TSV file
    )

# CHECK
df_isoform_shape_initial = df_isoform.shape[1]
# RENAME - missing header col
df_isoform = df_isoform.rename(columns = {'Unnamed: 0': 'Geneid'})

# NOTES: This col's format:
# 
# - <gene_name>_ENSEMBLE###.#
#
# ### (MOVED from lower down barcodes section up here)
# 
# ### Create the `barcodes.tsv` file
# 
# - Get the 2nd row of the original `matrix.tsv` file (which is now the header row)
# - Extract all the barcodes, drop the other columns
# - Drop the first 2 rows (already done during import)
# - Get the barcode column names as a list
# - Extract the columns I want:

# RELOAD - matrix file into a Pandas DF
df_matrix_header = pd.read_csv(
    file_path_matrix, 
    sep = '\t',     # For TSV file
    skiprows = 1,   # 1st row is metadata, skip it (index = 0)
    nrows = 1,      # Only load the 1st row (the header row)
    header = None 
    )

# Convert to a 1D series
header_series = df_matrix_header.squeeze()
# Remove the initial columns from the header
header_series_subset = header_series[7:]
# Remove the initial columns from the header
barcodes_matrix_list = header_series[7:].tolist()

# Copy the 'Geneid" column and paste it as a new column 'gene_name'
df_isoform['gene_name'] = df_isoform['Geneid']

# Create a list of just the column names in the order that you want it in the genes.tsv file
# Then add on the barcodes after that
cols = ['Geneid', 'gene_name']
cols.extend(barcodes_matrix_list)

# Remake the DF in that order
df_isoform = df_isoform[cols]

# CHECK
df_isoform_shape_final = df_isoform.shape[1]
print('df_isoform_shape_initial =', df_isoform_shape_initial)
print('df_isoform_shape_final   =', df_isoform_shape_final)
assert df_isoform_shape_final == df_isoform_shape_initial + 1, f'{df_isoform_shape_final} =/= {df_isoform_shape_initial} + 1'

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ### Apply to the FULL merged DF file
# Merge all
df_merged = pd.merge(df_matrix_drop, df_isoform, how = 'outer')

# # EXPORT - as merged .TSV file
# 
# Next:
# 
# - This will be a new merged .TSV file to run thru the rest of the CONVERSION pipeline
# - run thru the other CONVERSION script
# - run thru Seurat to see if this works, it might be the gene name col that screws it up, even if I specify to use col 1 since the matrix.tsv file needs to use col 2...

# EXPORT
df_merged.to_csv(
    file_path_merged_matrix_and_matrix_isoform, 
    sep = '\t', 
    index = False, 
    )
print(f'Saved merged_matrix_and_matrix_isoform.tsv to: {file_path_merged_matrix_and_matrix_isoform}')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# # CONVERSION
# 
# - from merged matrix + matrix_isoform to:
#   - genes.tsv
#   - barcodes.tsv
#   - matrix.mtx

# RELOAD - the merged matrix file 
df_merged = pd.read_csv(
    file_path_merged_matrix_and_matrix_isoform, 
    sep = '\t',     # For TSV file
    )

# ### Create the `genes.tsv` file 
# 
# - Drop the first 2 rows (already done during import)
# - Extract 2 columns: 
#   - Geneid (the ensemble #s) and 
#   - gene_name

# Select your desired cols
df_gene_id_and_gene_name = df_merged[['Geneid', 'gene_name']]

# ### EXPORT - `genes.tsv`

# Define output path
file_path_genes_tsv = os.path.join(sample_dir, 'genes.tsv')

# Save file as `genes.tsv`
df_gene_id_and_gene_name.to_csv(
    file_path_genes_tsv, 
    sep = '\t', 
    index = False, 
    header = False
    )
print(f'Saved genes.tsv to: {file_path_genes_tsv}')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# # (NOTE: MOVING SOME OF THE BARCODES.TSV CODE UP)

# Change the formatting to match examples and our data
header_series_subset_mod = header_series_subset + '-1'

# Define output path for `barcodes.tsv`
file_path_barcodes_tsv = os.path.join(sample_dir, 'barcodes.tsv')

# EXPORT - file to `barcodes.tsv`
header_series_subset_mod.to_csv(
    file_path_barcodes_tsv, 
    sep = '\t', 
    index = False, 
    header = False
    )
print(f'Saved barcodes.tsv to: {file_path_barcodes_tsv}')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ### Create the `matrix.mtx` file
# 
# - Import `matrix.tsv` as a DF, 
# - Get rid of cols from `chr` to `gene_name`, just keep `geneid` and then all the counts cols, 
# - then convert to a sparse matrix format (needs library matrix to convert to DCG matrix)
# - add dimensions, row, col, 

# RELOAD - merged matrix file w/changes
df_merged_no_header = pd.read_csv(
    file_path_merged_matrix_and_matrix_isoform, 
    sep = '\t',         # For TSV file
    skiprows = [0],     # 1st row is header, skip (index = 0),
    header = None,      # SKIP header this time
    )

# CHECK - using the earlier DF with header names 
cols_to_drop = df_merged.columns[00:2]
# Use the "full" DF
cols_to_drop = df_merged_no_header.columns[00:2]
# APPLY - Drop to those cols in the full DF
df_merged_no_header_drop = df_merged_no_header.drop(cols_to_drop, axis = 1)

# ### EXPORT - processed DF to `matrix_dropped.tsv`
# Define output file
file_path_matrix_dropped_tsv = os.path.join(sample_dir, 'matrix_dropped.tsv')

# EXPORT - the processed DF to a .TSV file
df_merged_no_header_drop.to_csv(
    file_path_matrix_dropped_tsv, 
    sep = '\t', 
    index = False, 
    header = False
    )
print(f'Saved matrix_dropped.tsv to: {file_path_matrix_dropped_tsv}')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ### CONVERSION - function to convert the TSV file to a sparse matrix .MTX file

def tsv_to_mtx(input_file, output_file):
    
    # Read the TSV file
    df = pd.read_csv(
        input_file, 
        sep = '\t',             # For TSV file
        index_col = None,       # Don't make an index col
        header = None,          # Don't make a header row
        # low_memory = False,   # WARNING - will slow it down, but prevents type errors from sampling
        )
    
    # Convert DataFrame to numpy array
    data = df.to_numpy()
    
    # Create sparse matrix
    sparse_matrix = sparse.csr_matrix(data)
    
    # Save the sparse matrix to a .mtx file
    mmwrite(output_file, sparse_matrix)

# ### EXPORT - `matrix.mtx` file as a sparse matrix
# Define output file
output_file = os.path.join(sample_dir, 'matrix.mtx')

# Run conversion function to export to a .MTX file
tsv_to_mtx(file_path_matrix_dropped_tsv, output_file)

# NOTE - in the JNB running in `sinteractive` the mem exceeded 160 GB when the code fails, gets up to 100 GB if it works tho, uses less, like ~50 GB from sbatch
