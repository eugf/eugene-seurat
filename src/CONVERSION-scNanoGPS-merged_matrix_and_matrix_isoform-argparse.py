#!/usr/bin/env python
# coding: utf-8

# 10/4/24
# 
# LAST UPDATED: 1/7/25
# 
# By Eugene Fong
# 
# # #* TODO - add log?
# 
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
#
# # HOW TO RUN (for converting to a merged matrix + matrix_isoform file)
# # - Open a terminal where this script file is and run:
# module load python/3.10
# python CONVERSION-scNanoGPS-merged_matrix_and_matrix_isoform-argparse.py --input_path <INPUT_PATH> --output_path <OUTPUT_PATH>
#
# # IMPORTS

# In[1]:


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
# Create ArgumentParser object
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
    help        = 'Base folder for processed files and plots to be stored',
    required    = True,
    type        = str,
    )

# Parse arguments
args = parser.parse_args()

# CHECK
print('args =', args)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# # LOAD DATA
# 
# See `/scNanoGPS_res` folder for a list of processed sample outputs at:
# 
#     /data/CARDPB/data/snRNA_longread/scNanoGPS-neuro/scNanoGPS_res

# In[2]:


# DEFINE FILE PATHS - CHANGE THIS!!!
# Input paths - look inside the scNanoGPS output folder (`/scNanoGPS_res`) by patient ID until you reach the parent folder that holds the matrix.tsv + matrix_isoform.tsv files
scNanoGPS_output_dir = args.input_path

# From there, create the paths for the matrix and matrix isoform files
file_path_matrix            = os.path.join(scNanoGPS_output_dir, 'matrix.tsv')
file_path_matrix_isoform    = os.path.join(scNanoGPS_output_dir, 'matrix_isoform.tsv')

# CHECK
print('scNanoGPS_output_dir     =', scNanoGPS_output_dir)
print('file_path_matrix         =', file_path_matrix)
print('file_path_matrix_isoform =', file_path_matrix_isoform)


# In[3]:


# # Output dir
output_dir = args.output_path
print('output_dir =', output_dir)


# In[4]:


# Extract the patient ID
id = os.path.basename(os.path.dirname(os.path.dirname(scNanoGPS_output_dir)))
print('id =', id)

# ALT - better way
separators = ['_FTX_singlecell', '_singlecell']
split_id = re.split('|'.join(separators), id)[0]

# CHECK
print('split_id =', split_id)


# In[5]:


# Get today's date
current_date = datetime.date.today()
print('current_date =', current_date)

formatted_date = current_date.strftime("%Y%m%d")
print('formatted_date =', formatted_date)


# In[6]:


# Create unique output folders based on the patient ID and date
sample_dir = os.path.join(output_dir, split_id, formatted_date)
print('sample_dir =', sample_dir)

# Start a counter for repeated folders
counter = 0
   
# CHECK - if the path link exists as text, adjust name as needed
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
    print(f"The directory '{sample_dir}' does not exist")

# Make the new directory 
os.makedirs(sample_dir)


# In[7]:


# Save full ID to text file
output_id_txt = os.path.join(sample_dir, 'id.txt')

# Open file to write and save ID
with open(output_id_txt, 'w+') as file:
    file.write(id)
    # file.write(split_id)

print(f'ID written to: {output_id_txt}')


# In[8]:


# OUTPUT PATHS
# Save file as `merged_matrix_and_matrix_isoform.tsv`
file_path_merged_matrix_and_matrix_isoform = os.path.join(sample_dir, 'merged_matrix_and_matrix_isoform.tsv')

# TEST - subset of 5 rows file
# Save file as `merged_matrix_and_matrix_isoform.tsv`
file_path_merged_matrix_and_matrix_isoform_5 = os.path.join(sample_dir, 'merged_matrix_and_matrix_isoform_5.tsv')


# In[9]:


# LOAD - matrix file into a Pandas DF
df_matrix = pd.read_csv(
    file_path_matrix, 
    sep = '\t',     # For TSV file
    skiprows = 0,   # 1st row is metadata, skip it (index = 0)
    header = 1      # Header begins on the 2nd row (index = 1)
    )

df_matrix


# In[ ]:


# CHECK
print('df_matrix.shape =', df_matrix.shape)


# # PLANNING
# 
# By row/col
# 
# ```
# Row 0   = command (DROP)
# Row 1   = header
# Row 2+  = Transcripts (ENS ID #s) (KEEP)
# 
# ###########################
# 
# Col 0       = Geneid
# Col 1 - 5   = (DROP)
# Col 6       = gene_name (KEEP)
# Col 7+      = barcodes (KEEP)
# ```
# 
# ### Steps overview
# 
# - KEEP: Geneid, gene_name, and the rest of the barcodes
# - DROP: in betw cols
# - Then merge on the barcodes

# In[10]:


# Drop these cols from the DFs
df_matrix_drop = df_matrix.drop(
    columns = ['Chr', 'Start', 'End', 'Strand', 'Length'], 
    inplace = False)
df_matrix_drop


# In[11]:


# LOAD - matrix isoform file into a Pandas DF
df_isoform = pd.read_csv(
    file_path_matrix_isoform, 
    sep = '\t',     # For TSV file
)
df_isoform


# In[12]:


# CHECK
df_isoform_shape_initial = df_isoform.shape[1]
print('df_isoform_shape_initial =', df_isoform_shape_initial)


# In[13]:


# RENAME - missing header col
df_isoform = df_isoform.rename(columns = {'Unnamed: 0': 'Geneid'})
df_isoform.head()


# ### #? Q) Do I need to rebuild the `gene_name` col?
# ### #? A) Maybe? Copy and paste it over entirely
# 
# NOTES: This col's format:
# 
# - <gene_name>_ENSEMBLE###.#

# ### (MOVED from lower down barcodes section up here)
# 
# ### Create the `barcodes.tsv` file
# 
# - Get the 2nd row of the original `matrix.tsv` file (which is now the header row)
# - Extract all the barcodes, drop the other columns

# - Drop the first 2 rows (already done during import)
# - Get the barcode column names as a list
# - Extract the columns I want:

# In[14]:


# RELOAD - matrix file into a Pandas DF
df_matrix_header = pd.read_csv(
    file_path_matrix, 
    sep = '\t',     # For TSV file
    skiprows = 1,   # 1st row is metadata, skip it (index = 0)
    nrows = 1,      # Only load the 1st row (the header row)
    header = None 
    )

df_matrix_header


# In[15]:


# Convert to a 1D series
header_series = df_matrix_header.squeeze()
header_series.head(10)


# In[16]:


# Remove the initial columns from the header
header_series_subset = header_series[7:]
header_series_subset


# In[17]:


# Remove the initial columns from the header
barcodes_matrix_list = header_series[7:].tolist()
print('barcodes_matrix_list[:5] =')
print(barcodes_matrix_list[:5])


# In[18]:


# Copy the 'Geneid" column and paste it as a new column 'gene_name'
df_isoform['gene_name'] = df_isoform['Geneid']

# Create a list of just the column names in the order that you want it the genes.tsv file
cols = ['Geneid', 'gene_name']
cols.extend(barcodes_matrix_list)
# print('cols =', cols)

# Remake the DF in that order
df_isoform = df_isoform[cols]
df_isoform


# In[19]:


# CHECK
df_isoform_shape_final = df_isoform.shape[1]
print('df_isoform_shape_final =', df_isoform_shape_final)

assert df_isoform_shape_final == df_isoform_shape_initial + 1, f'{df_isoform_shape_final} =/= {df_isoform_shape_initial} + 1'


# ### TEST BLOCKS

# In[20]:


# # TEST - subset the DF to only 5 rows
df_matrix_drop_5 = df_matrix_drop.head(5)
df_matrix_drop_5


# In[21]:


# TEST - subset the DF to only 5 rows
df_isoform_5 = df_isoform.head(5)
df_isoform_5


# In[22]:


# TEST - the col dimensions match up w/the matrix dimensions, great!
# Merge all
df_merged_5 = pd.merge(df_matrix_drop_5, df_isoform_5, how = 'outer')
df_merged_5


# In[23]:


# TEST - EXPORT
df_merged_5.to_csv(
    file_path_merged_matrix_and_matrix_isoform_5, 
    sep = '\t', 
    index = False, 
    )

print(f"Saved merged_matrix_and_matrix_isoform_5.tsv to: {file_path_merged_matrix_and_matrix_isoform_5}")


# In[24]:


# TEST - RELOAD - matrix file into a Pandas DF
df_merged_5 = pd.read_csv(
    file_path_merged_matrix_and_matrix_isoform_5, 
    sep = '\t',     # For TSV file
    )

df_merged_5


# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# ### Apply to the FULL merged DF file

# In[25]:


# Merge all
df_merged = pd.merge(df_matrix_drop, df_isoform, how = 'outer')


# In[26]:


# CHECK
df_merged.head()


# In[27]:


# CHECK - see the NaN's in the gene_name col?
df_merged.tail()


# In[28]:


# CHECK - for floats
df_merged.iloc[: , 3].unique()


# In[29]:


# CHECK - for dtypes
# NOTE - floats mostly + 2 objects, uh oh
df_merged.info()


# # EXPORT - as .TSV file
# 
# Next:
# 
# - This will be a new merged .TSV file to run thru the rest of the CONVERSION pipeline
# - run thru the other CONVERSION script
# - run thru Seurat to see if this works, it might be the gene name col that screws it up, even if I specify to use col 1 since the matrix.tsv file needs to use col 2...

# In[30]:


# EXPORT
df_merged.to_csv(
    file_path_merged_matrix_and_matrix_isoform, 
    sep = '\t', 
    index = False, 
    )

print(f'Saved merged_matrix_and_matrix_isoform.tsv to: {file_path_merged_matrix_and_matrix_isoform}')


# # CONVERSION
# 
# - from merged matrix + matrix_isoform to:
#   - genes.tsv
#   - barcodes.tsv
#   - matrix.mtx

# In[31]:


# # Output dir
# output_dir = '/data/CARDPB/data/snRNA_longread/eugene-seurat/output/merged'
# print('output_dir =', output_dir)


# In[32]:


# RELOAD - the merged matrix file 
df_merged = pd.read_csv(
    file_path_merged_matrix_and_matrix_isoform, 
    sep = '\t',     # For TSV file
    )

df_merged

# #! WARNING
# /tmp/ipykernel_3477477/3859502601.py:2: DtypeWarning: Columns (1) have mixed types. Specify dtype option on import or set low_memory=False.
#   df = pd.read_csv(


# # CONVERSION
# 
# ### Create the `genes.tsv` file 
# 
# - Drop the first 2 rows (already done during import)
# - Extract 2 columns: 
#   - Geneid (the ensemble #s) and 
#   - gene_name

# In[34]:


# Select your desired cols
df_gene_id_and_gene_name = df_merged[['Geneid', 'gene_name']]
df_gene_id_and_gene_name.head()


# ### EXPORT - `genes.tsv`

# In[35]:


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


# # (NOTE: MOVING SOME OF THE BARCODES.TSV CODE UP)

# In[36]:


# Change the formatting to match examples and our data
#* NOTE - add the donor ID at a later step in Seurat instead of here
header_series_subset_mod = header_series_subset + '-1'
header_series_subset_mod


# ### EXPORT - `barcodes.tsv`
# 
# #* NOTE - check if the index is still there

# In[37]:


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


# ### Create the `matrix.mtx` file
# 
# - Import `matrix.tsv` as a DF, 
# - Get rid of cols from `chr` to `gene_name`, just keep `geneid` and then all the counts cols, 
# - then convert to a sparse matrix format (needs library matrix to convert to DCG matrix)
# - add dimensions, row, col, 

# In[38]:


# RELOAD - merged matrix file w/changes
df_merged_no_header = pd.read_csv(
    file_path_merged_matrix_and_matrix_isoform, 
    sep = '\t',         # For TSV file
    skiprows = [0],     # 1st row is header, skip (index = 0),
    header = None,      # SKIP header this time
    )

df_merged_no_header


# In[39]:


# CHECK - using the earlier DF with header names 
cols_to_drop = df_merged.columns[00:2]
print('cols_to_drop =', cols_to_drop)

# Use the "full" DF
cols_to_drop = df_merged_no_header.columns[00:2]
print('cols_to_drop =', cols_to_drop)


# In[40]:


# APPLY - Drop to those cols in the full DF
df_merged_no_header_drop = df_merged_no_header.drop(cols_to_drop, axis = 1)
df_merged_no_header_drop.head()


# In[41]:


# CHECK - compare - are there still 2 object datatypes?
df_merged.info()


# In[42]:


# CHECK - dtypes
df_merged_no_header_drop.info()


# ### EXPORT - processed DF to .TSV file
# 
# Save file as `matrix_dropped.tsv`

# In[43]:


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


# ### CONVERSION - function to convert the TSV file to a sparse matrix .MTX file
# 
# NOTES:
# - having non-numerical data types is a problem... what to do about the header and index?
# - prolly drop them for the matrix since its in `barcodes.tsv` and `genes.tsv`

# In[44]:


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

# In[45]:


# Define output file
output_file = os.path.join(sample_dir, 'matrix.mtx')

# Run conversion function to export to a .MTX file
tsv_to_mtx(file_path_matrix_dropped_tsv, output_file)
print(f'Saved matrix.mtx to: {output_file}')

# NOTE - in the JNB running in `sinteractive` the mem exceeded 160 GB when the code fails, gets up to 100 GB if it works tho, uses less, like ~50 GB from sbatch

