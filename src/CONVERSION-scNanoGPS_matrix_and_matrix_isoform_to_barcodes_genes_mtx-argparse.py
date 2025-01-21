#!/usr/bin/env python
# coding: utf-8

# STARTED:      10/4/24
# 
# LAST UPDATED: 12/32/24
# 
# By Eugene Fong
# 
# # #* TODO - add log?
# 
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
# 
# # IMPORTS

# In[1]:


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
# scNanoGPS_output_dir = '/data/CARDPB/data/snRNA_longread/scNanoGPS-neuro/scNanoGPS_res/SH-07-46_singlecell_LR/20241105-1/20241105'
scNanoGPS_output_dir = args.input_path

# From there, append the paths for the matrix files
file_path_matrix            = os.path.join(scNanoGPS_output_dir, 'matrix.tsv')
file_path_matrix_isoform    = os.path.join(scNanoGPS_output_dir, 'matrix_isoform.tsv')

# CHECK
print('scNanoGPS_output_dir     =', scNanoGPS_output_dir)
print('file_path_matrix         =', file_path_matrix)
print('file_path_matrix_isoform =', file_path_matrix_isoform)


# In[3]:


# # Output dir
# output_dir = '/data/CARDPB/data/snRNA_longread/eugene-seurat/output/long_reads'
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
   
# CHECK - if the path exists already, adjust name as needed
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


# In[7]:


# # TEST
# # Turn the make directory checker into a func
# def mkdir_checker(output_dir):
    
#     # Create unique output folders based on the patient ID and date
#     global sample_dir
#     sample_dir = os.path.join(output_dir, split_id, formatted_date)
#     print('sample_dir =', sample_dir)

#     # Start a counter for repeated folders
#     counter = 0
    
#     # CHECK - if the path link exists as text, adjust name as needed
#     if os.path.exists(sample_dir):
#         print(f"The directory '{sample_dir}' exists.")
#         new_sample_dir = sample_dir + '-' + str(counter)
        
#         # Increment on the date folder
#         while os.path.exists(new_sample_dir):
#             counter += 1
#             print('counter =', counter)
#             new_sample_dir = sample_dir + '-' + str(counter)
#             print('new_sample_dir INSIDE IF BLOCK =', new_sample_dir)
        
#         # Reset this variable to call later
#         sample_dir = new_sample_dir

#     else:
#         print(f"The directory '{sample_dir}' does not exist.")

#     # Make the new directory 
#     os.makedirs(sample_dir)

# # RUN - func
# mkdir_checker(output_dir)


# In[8]:


# Save full ID to text file (will process this way downstream)
output_id_txt = os.path.join(sample_dir, 'id.txt')

# Open file to write and save ID
with open(output_id_txt, 'w+') as file:
    file.write(id)

print(f'ID written to: {output_id_txt}')


# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # MATRIX FILE
# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# In[9]:


# LOAD - matrix file into a Pandas DF
df = pd.read_csv(
    file_path_matrix, 
    sep = '\t',     # For TSV file
    skiprows = 0,   # 1st row is metadata, skip it (index = 0)
    header = 1      # Header begins on the 2nd row (index = 1)
    )

df


# In[10]:


# CHECK
print('df.shape =', df.shape)


# In[11]:


# #* TODO - HANDLING METADATA???

# # Assuming your first three columns are: cell_id, sample, isoform_id
# seurat_obj@meta.data[, c("cell_id", "sample")] <- isoform_data[1:2]

# # Add any additional metadata columns
# seurat_obj@meta.data$condition <- isoform_data[3]


# # CONVERSIONS
# 
# ### Create the `genes.tsv` file 
# 
# - Drop the first 2 rows (already done during import)
# - Extract 2 columns ONLY: 
#   - Geneid (ENSEMBLE ID #s) and
#   - gene_name

# In[12]:


# Select your desired cols
df_gene_id_and_gene_name = df[['Geneid', 'gene_name']]
df_gene_id_and_gene_name


# ### EXPORT - `genes.tsv`

# In[13]:


# Define output path
file_path_genes_tsv = os.path.join(sample_dir, 'genes.tsv')

# Save file as `genes.tsv`
df_gene_id_and_gene_name.to_csv(
    file_path_genes_tsv, 
    sep = '\t',         # For TSV file
    index = False,      # Remove index col???
    header = False,     # Remove header row???
    )

print(f'Saved genes.tsv to: {file_path_genes_tsv}')


# ### Create the `barcodes.tsv` file
# 
# - Load just the 2nd row of the original `matrix.tsv` file (which is now the header row)
# - Extract all the barcodes, drop the other columns

# In[14]:


# RELOAD - matrix file into a Pandas DF
df_header = pd.read_csv(
    file_path_matrix, 
    sep = '\t',     # For TSV file
    skiprows = 1,   # 1st row is metadata, skip it (index = 0)
    nrows = 1,      # Only load the next 1 row (the header row)
    header = None   # Prevents Pandas defaulting to header = True
    )

df_header


# In[15]:


# Converts to a 1D series
header_series = df_header.squeeze()
header_series.head(10)


# In[16]:


# Remove the initial columns from the header, keep the barcodes ONLY
header_series_subset = header_series[7:]
header_series_subset


# In[17]:


# Change the formatting to match examples and our data
#* NOTE - add the donor ID at a later step in Seurat instead of here
header_series_subset_mod = header_series_subset + '-1'
header_series_subset_mod


# ### EXPORT - `barcodes.tsv`

# In[18]:


# Define output path
file_path_barcodes_tsv = os.path.join(sample_dir, 'barcodes.tsv')

# Save file as `barcodes.tsv`
header_series_subset_mod.to_csv(
    file_path_barcodes_tsv, 
    sep = '\t',         # For TSV file
    index = False,      # No index col???
    header = False,     # Prevents Pandas defaulting to header = True
    )

print(f'Saved barcodes.tsv to: {file_path_barcodes_tsv}')


# ### Create the `matrix.mtx` file
# 
# - RELOAD `matrix.tsv` as a DF, 
# - Get rid of cols from `chr` to `gene_name`, just keep `geneid` and then all the counts cols, 
# - then convert to a sparse matrix format (needs matrix library to convert to DCG matrix)

# In[19]:


# RELOAD - matrix file into a Pandas DF (but this time get rid of the header)
df_full = pd.read_csv(
    file_path_matrix, 
    sep = '\t',         # For TSV file
    skiprows = [0,1],   # 1st row is metadata, 2nd is header, skip BOTH (index = 0,1),
    header = None,      # Prevents Pandas defaulting to header = True
    )

# CHECK - it shouldn't have words in the "header" but just #'s
df_full.head()


# In[20]:


# CHECK
df_full.tail()


# In[21]:


# CHECK - using the earlier DF with header names 1st
cols_to_drop = df.columns[00:7]
print('cols_to_drop =', cols_to_drop)

# Now using the "full" DF, define the corresponding col's we want to drop
cols_to_drop = df_full.columns[00:7]
print('cols_to_drop =', cols_to_drop)


# In[22]:


# APPLY - Drop to those cols in the full DF
df_full_drop = df_full.drop(cols_to_drop, axis = 1)
df_full_drop.head()


# In[23]:


# CHECK - before and after
df.info()


# In[24]:


# CHECK - dtypes should all be #;s
df_full_drop.info()


# ### EXPORT - CHECK - processed DF to `matrix_dropped.tsv` file 

# In[25]:


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
# 
# NOTES:
# - having non-numerical data types causes an ERROR... what to do about the header and index?
# - prolly drop them for the matrix file since its in `barcodes.tsv` and `genes.tsv`

# In[26]:


def tsv_to_mtx(input_file, output_file):
    
    # Read the TSV file
    df = pd.read_csv(
        input_file, 
        sep = '\t',         # For TSV file
        index_col = None,   # No index col???
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

# In[ ]:


# Define output path
output_file = os.path.join(sample_dir, 'matrix.mtx')

# Run conversion function
tsv_to_mtx(file_path_matrix_dropped_tsv, output_file)


# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # MATRIX ISOFORM FILE
# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# # LOAD DATA

# In[28]:


# DEFINE FILE PATHS
# Output dir - new isoform folder inside the sample directory
output_isoform_dir = os.path.join(sample_dir, 'isoform')

# Make the new directory 
os.makedirs(output_isoform_dir)

# CHECK
print('output_isoform_dir =', output_isoform_dir)


# In[29]:


# LOAD - matrix isoform file into a Pandas DF
df_isoform = pd.read_csv(
    file_path_matrix_isoform, 
    sep = '\t',     # For TSV file
    )

df_isoform


# In[30]:


# CHECK - # of unique rows (might be ALL of them for isoforms)
len(df_isoform['Unnamed: 0'].unique())


# ### NOTES
# 
# Cols are:
# 
# - `Geneid` (ENS #)
# - MISSING `gene_name`, which should be derived from the isoform's parent?
# - The rest are barcodes
# - NO index col
# 
# Rows:
# 
# - Header row

# In[31]:


# # CHECK
# df_isoform['CATTATCTCGCACACA'].unique()

test = df_isoform.iloc[:, 2].unique()
test


# In[32]:


# CHECK
type(test)


# In[33]:


df_isoform.iloc[:,2].info()


# # CONVERSION
# 
# ### Create the `genes.tsv` file 
# 
# - Rename col 1 (currently unnamed)
# - Extract 2 columns: 
#   - `Geneid` (Ensemble #s)  
#   - `gene_name` --> MISSING! need to create! but how???
# - Drop the 1st row (header row)

# In[34]:


# RENAME - missing header col
df_isoform = df_isoform.rename(columns = {'Unnamed: 0': 'Geneid'})
df_isoform.head()


# In[35]:


# CHECK - length before splitting the col names
pre_split = len(df_isoform['Geneid'].unique())
print('pre_split =', pre_split)


# ### EXPORT - CHECK

# In[36]:


# EXPORT - the file before splitting the col names
df_isoform.to_csv(
    os.path.join(output_isoform_dir, 'isoform_before_splitting.tsv'), 
    sep = '\t', 
    index = False, 
    header = True
    )


# ### #? Q) Do I need to rebuild the `gene_name` col?
# 
# ### #? A) will error w/o it
# 
# ### -- OR --
# 
# ### leave it alone, give the col a name, assign Seurat to use col 1 instead
# 
# NOTES: This col's format:
# 
# - <gene_name>_ENSEMBLE###.#
# 
# #? Q) Do I keep the original name? 
# 
# #* A) Tried keeping the original name, but plotting the counts in Seurat was weird
# 
# #* NEXT - try splitting the `gene_name` off col 1, only keep the `Geneid` (Ensemble #)

# In[37]:


# Split the `Geneid` column on '_ENST' then add back 'ENST'

# Keep the 0th element to create a new `gene_name` column
df_isoform['gene_name'] = df_isoform['Geneid'].str.split('_ENST').str[0]

# Rename the `Geneid` col using the 1st element of the split and prepending it with 'ENST' to restore the lost characters from the split
df_isoform['Geneid'] = 'ENST' + df_isoform['Geneid'].str.split('_ENST').str[1]

# CHECK
df_isoform


# In[ ]:


# CHECK
post_split = len(df_isoform['Geneid'].unique())
print('post_split =', post_split)

# HALT SCRIPT IF FAILED
assert pre_split == post_split, "ERROR! values lost!"


# In[39]:


# Select your desired cols
df_isoform_gene_id_and_gene_name = df_isoform[['Geneid', 'gene_name']]
df_isoform_gene_id_and_gene_name


# ### CHECKS

# In[40]:


# CHECK - the # of rows match
num_rows_df_isoform = len(df_isoform)
print('num_rows_df_isoform =', num_rows_df_isoform)

num_rows_df_isoform_gene_id_and_gene_name = len(df_isoform_gene_id_and_gene_name)
print('num_rows_df_isoform_gene_id_and_gene_name =', num_rows_df_isoform_gene_id_and_gene_name)

assert num_rows_df_isoform == num_rows_df_isoform_gene_id_and_gene_name, "WARNING! # of rows in original DF and split col DF do NOT match"


# In[41]:


# CHECK - for dupes in the Geneid column only
duplicates = df_isoform['Geneid'].duplicated(keep=False)  # keep=False marks all duplicates as True
print("Unique values of dupes in Geneid col check = ", duplicates.unique())

# COUNT - T/F values
counts = duplicates.value_counts()
print("Counts of T/F in duplicates:")
print(counts)


# In[42]:


# CHECK - Show duplicate values
duplicate_values = df_isoform['Geneid'][duplicates]
print("Duplicate values in the Geneid column: \n", duplicate_values)


# In[43]:


# EXPORT - as text so i can read the whole thing
duplicate_values.to_csv(
    os.path.join(output_isoform_dir, 'df_isoform_duplicate_values.tsv'), 
    sep = '\t',     # For TSV file
    )


# ### CHECK

# In[44]:


# # LOAD - isoform_var_features.tsv exported from the scNanoGPS Seurat object
# df_isoform_var_features = pd.read_csv(
#     "/home/fonge2/scNanoGPS/eugene-seurat/output/TEST/isoform/isoform_var_features.tsv", 
#     sep = '\t',     # For TSV file
# )
# df_isoform_var_features


# In[45]:


# # CHECK - overlap betw the `genes.tsv` and `isoform_var_features.tsv` files
# # Convert the 'Geneid' columns to sets
# set1 = set(df_isoform['gene_name'])
# set2 = set(df_isoform_var_features['x'])

# # # CHECK - comment out cuz that's way too much text
# # print('set1 = \n', set1)
# # print('set2 = \n', set2)

# # Find common elements
# common_genes = set1.intersection(set2)
# print('common_genes = \n', common_genes)

# # Get the number of common values
# num_common = len(common_genes)
# print(f"Number of common values: {num_common}")


# In[46]:


# # CHECK - extract only the unique values from `isoform_var_features.tsv` --> this shouldn't matter actually

# # CHECK - split the `x` column values on the .period. and only keep the 0th element

# # Keep the 0th element to create a new `gene_name` column
# df_isoform_var_features['x_split'] = df_isoform_var_features['x'].str.split('.').str[0]
# df_isoform_var_features.head(50)


# In[47]:


# # CHECK - overlap betw the `genes.tsv` and `isoform_var_features.tsv` files
# # Convert the 'Geneid' columns to sets
# set1 = set(df_isoform['gene_name'])
# set2 = set(df_isoform_var_features['x_split'])

# # # CHECK
# # print('set1 = \n', set1)
# # print('set2 = \n', set2)

# # Find common elements
# common_genes = set1.intersection(set2)
# print('common_genes = \n', common_genes)

# # Get the number of common values
# num_common = len(common_genes)
# print(f"Number of common values: {num_common}")


# ### EXPORT - `genes.tsv`

# In[48]:


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


# In[49]:


# RELOAD - matrix isoform file into a Pandas DF
df_isoform_before_splitting = pd.read_csv(
    os.path.join(output_isoform_dir, 'isoform_before_splitting.tsv'),
    sep = '\t',     # For TSV file
    )

df_isoform_before_splitting


# In[50]:


# Select your desired cols
df_isoform_before_splitting_gene_id = df_isoform_before_splitting[['Geneid']]
df_isoform_before_splitting_gene_id


# In[ ]:


# OVERWRITE - the genes.tsv file with the original, unchanged, fused column name, but we'll just call it `Geneid`
# Define output path
file_path_genes_tsv = os.path.join(output_isoform_dir, 'genes.tsv')

# EXPORT - Save file as `genes.tsv`
df_isoform_before_splitting_gene_id.to_csv(
    # os.path.join(output_isoform_dir, 'genes.tsv'), 
    file_path_genes_tsv,
    sep = '\t',         # For TSV file
    index = False,      # Remove index col???
    header = False,     # Remove header row???
    )

print(f'Saved genes.tsv to: {file_path_genes_tsv}')


# ### Create the `barcodes.tsv` file
# 
# - Get the 1st row of the `matrix_isoform.tsv` file
# - Drop the 1st col (ENS #s)
# - NOTE: 1st col is empty and there's no other metadata cols

# In[52]:


# LOAD - matrix isoform file into a Pandas DF to get the header row only
df_isoform_header = pd.read_csv(
    file_path_matrix_isoform, 
    sep = '\t',         # For TSV file
    nrows = 1,          # Only load the 1st row (the header row)
    header = None,      # Header begins on the 1st row
)
df_isoform_header


# In[53]:


# Drop the 1st col (the unnamed Geneid col here)
df_isoform_header_skip = df_isoform_header.iloc[:, 1:]
df_isoform_header_skip


# In[54]:


# Convert to a 1D series
isoform_header_series = df_isoform_header_skip.squeeze()
isoform_header_series.head(10)


# In[55]:


# Change the formatting to match examples and our data
#* TODO - add the donor ID at a later step in Seurat instead of here (Cory suggested adding it here tho)
isoform_header_series_mod = isoform_header_series + '-1'
isoform_header_series_mod


# ### EXPORT - `barcodes.tsv`

# In[56]:


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


# ### Create the `matrix.mtx` file
# 
# - Import `matrix_isoform.tsv` as a DF, 
# - Keep `geneid` and then all the counts cols, drop the 1st row (header row) tho
# - then convert to a sparse matrix format (needs matrix library to convert to DCG matrix)

# In[ ]:


# RELOAD - matrix file into a Pandas DF (but this time get rid of the header)
df_isoform_full = pd.read_csv(
    file_path_matrix_isoform, 
    sep = '\t',         # For TSV file
    skiprows = [0],     # 1st row is metadata (index = 0),
    header = None,      # Prevents Pandas defaulting to header = True
    )

df_isoform_full.head()


# In[58]:


# Use the "full" DF
cols_to_drop_isoform = df_isoform_full.columns[00:1]
print('cols_to_drop_isoform =', cols_to_drop_isoform)


# In[59]:


# APPLY - Drop to those cols in the full DF
df_isoform_full_drop = df_isoform_full.drop(cols_to_drop_isoform, axis = 1)
df_isoform_full_drop


# In[60]:


# CHECK - before and after
df_isoform_full.info()


# In[61]:


# CHECK - should have only # dtypes
df_isoform_full_drop.info()


# ### EXPORT - processed DF to `matrix_isoform_dropped.tsv` file

# In[62]:


# Define output path
file_path_matrix_isoform_dropped_tsv = os.path.join(output_isoform_dir, 'matrix_dropped.tsv')

# EXPORT - save as .TSV file
df_isoform_full_drop.to_csv(
    file_path_matrix_isoform_dropped_tsv, 
    sep = '\t',         # For TSV file
    index = False,      # Remove index col???
    header = False,     # Remove header row???
    )

print(f'Saved matrix_isoform_dropped.tsv to: {file_path_matrix_isoform_dropped_tsv}')


# # CONVERSION - function to convert the TSV file a sparse matrix .MTX file
# 
# ### EXPORT - `matrix.mtx` file as a sparse matrix

# In[63]:


# Define output path
output_file = os.path.join(output_isoform_dir, 'matrix.mtx')

#! ERROR - check here for the exponents???
# input_df = pd.read_csv(file_path_matrix_isoform_dropped_tsv)

# Run conversion function
tsv_to_mtx(file_path_matrix_isoform_dropped_tsv, output_file)

# # CHECK
# print(f'Saved matrix.mtx to: {output_file}')

