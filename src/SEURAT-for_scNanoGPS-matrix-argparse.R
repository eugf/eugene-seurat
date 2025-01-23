#!/usr/bin/env Rscript

# # HOW TO RUN:
# # - Open a terminal where this script file is and run:
# Rscript SEURAT-for_scNanoGPS-matrix-argparse.R --input-dir <INPUT_DIR> --output-dir <OUTPUT_DIR> --celltype-marker-dir <CELLTYPE_MARKER_DIR>

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# SEURAT TUTORIAL
# SOURCE: https://satijalab.org/seurat/articles/pbmc3k_tutorial

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# IMPORTS 
require(dplyr)
require(Seurat)
require(patchwork)
require(base)
require(ggplot2)
require(argparse)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# SETUP

# INITIALIZE ARGPARSE
parser <- ArgumentParser()

# DEFINE ARGUMENTS THE USER WILL PROVIDE
parser$add_argument(
  '--input-dir', 
  required    = TRUE,
  type        = 'character'
)

parser$add_argument(
  '--output-dir', 
  required    = TRUE,
  type        = 'character'
)

parser$add_argument(
  '--celltype-marker-dir', 
  required    = TRUE,
  type        = 'character'
)

parser$add_argument(
  '--filter_nFeature_RNA_lower', 
  required    = FALSE,
  type        = 'character',
  default     = 0
)

parser$add_argument(
  '--filter_nFeature_RNA_upper', 
  required    = FALSE,
  type        = 'character',
  default     = 200
)

parser$add_argument(
  '--filter_percent_mt', 
  required    = FALSE,
  type        = 'character',
  default     = 5
)

# PARSE ARGUMENTS
args <- parser$parse_args()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# SETUP

# LOGGING
current_date  <- as.character(Sys.Date())
current_time  <- format(Sys.time(),"%H-%M-%S")

# # DEFINE PATHS
# Input path - takes the results of converting the `/scNanoGPS_res` folder into: barcodes.tsv, genes.tsv, and matrix.mtx
input_dir <- args$input_dir

# Get the Sample ID #
# Shorten ID - Split by underscore, and keep the first 2 elements
id_path <- file.path(input_dir, 'id.txt')
id <- readLines(id_path)
split_ids <- sapply(strsplit(id, "_"), function(x) x[1:2])
split_ids
rejoined_ids <- apply(split_ids, 2, function(x) paste(x, collapse = "_"))
rejoined_ids
id <- rejoined_ids

# # Output path
output_dir <- args$output_dir
output_path <- file.path(output_dir, 'PLOTS-Seurat-matrix')

# CHECK - if the directory exists, if yes, increment the name
create_incremented_dir <- function(base_dir) {
  new_dir <- base_dir
  counter <- 0
  
  # Loop to find an available directory name by incrementing the counter
  while (dir.exists(new_dir)) {
    new_dir <- paste(base_dir, counter, sep ='-')
    counter <- counter + 1
  }
  
  # Create the directory once a unique name is found
  dir.create(new_dir, recursive = TRUE)
  
  return(new_dir)
}

# Create directory with possible increment
output_path <- create_incremented_dir(output_path)
cat('Created directory at: \n', output_path, '\n')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# LOAD DATA

# SETUP THE SEURAT OBJECT
scNanoGPS.data <- Read10X(data.dir = input_dir)

# Initialize the Seurat object with the raw (non-normalized data).
scNanoGPS <- CreateSeuratObject(
  counts  = scNanoGPS.data, 
  project = id)

# Get info
scNanoGPS_initial_info <- print(scNanoGPS)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# STANDARD PRE-PREPROCESSING WORKFLOW
scNanoGPS[["percent.mt"]] <- PercentageFeatureSet(scNanoGPS, pattern = "^MT-")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# MY FUNCTIONS

# FUNCTION - to construct the file name
my_plot_name <- function(plot_name = "") {
  
  # Construct file name and output folder based on user-input and other defined vars
  file_name = paste0('PLOT-', plot_name, '.png')
  
  # Construct the final file path
  file_path = file.path(output_path, file_name)
  cat('file_path = \n', file_path, '\n')
  
  # CHECK - if the output directory exists, if not: create it
  if (!dir.exists(dirname(file_path))) {
    dir.create(dirname(file_path), recursive = TRUE)
  }
  
  return(file_path)
}

# FUNCTION - to save plots
# - Include some sensible defaults for width, height, DPI
# - Take in user-defined name, save last plot
my_plot_save <- function(
    filename  = my_plot_name(),
    width     = 8, 
    height    = 6, 
    dpi       = 300
) 
{  
  # Define the file format and dimensions to open with graphics device
  ggsave(
    filename  = filename,
    width     = width, 
    height    = height, 
    dpi       = dpi
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# # TEMPLATE
# SAVE PLOT
# save_plot_name <- my_plot_name('plot_name')
# my_plot_save(save_plot_name)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# PLOTS
# Visualize QC metrics, and use these to filter cells.
#   - We filter cells that have unique feature counts
#   - We filter cells that have >% mitochondrial counts

# Visualize QC metrics as a violin plot
VlnPlot(
  scNanoGPS, 
  features = c(
    "nFeature_RNA", 
    "nCount_RNA", 
    "percent.mt"
    ), 
  ncol = 3)

# SAVE PLOT
save_plot_violin_counts_features <- my_plot_name('violin_counts_features')
my_plot_save(save_plot_violin_counts_features)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(scNanoGPS, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(scNanoGPS, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# SAVE PLOT
save_plot_feature_scatter <- my_plot_name('feature_scatter')
my_plot_save(
  save_plot_feature_scatter,
  height = 4,
  width = 10
)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# FILTER MANUALLY - based on the violin plot
# filter_nFeature_RNA_lower <- 0
# filter_nFeature_RNA_upper <- 200
# filter_percent.mt         <- 5
filter_nFeature_RNA_lower <- args$filter_nFeature_RNA_lower
filter_nFeature_RNA_upper <- args$filter_nFeature_RNA_upper
filter_percent_mt         <- args$filter_percent_mt

scNanoGPS <- subset(
  scNanoGPS, 
  subset  = nFeature_RNA > filter_nFeature_RNA_lower & 
            nFeature_RNA < filter_nFeature_RNA_upper & 
            percent.mt < filter_percent_mt
  )

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# NORMALIZING THE DATA
# After removing unwanted cells from the dataset, the next step is to normalize the data. 
# By default, we employ a global-scaling normalization method “LogNormalize” 
scNanoGPS <- NormalizeData(scNanoGPS)

# Identification of highly variable features (feature selection)
# We next calculate a subset of features that exhibit high cell-to-cell variation in the dataset (i.e, they are highly expressed in some cells, and lowly expressed in others). 
scNanoGPS <- FindVariableFeatures(scNanoGPS, selection.method = "vst", nfeatures = 10000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(scNanoGPS), 10)

# Plot variable features with and without labels
plot1 <- VariableFeaturePlot(scNanoGPS)

plot2 <- LabelPoints(
  plot = plot1,
  points = top10,
  repel = TRUE,
  xnudge = 0,
  ynudge = 0
)

plot1 + plot2

# SAVE PLOT
save_plot_variable_features <- my_plot_name('variable_features')
my_plot_save(
  save_plot_variable_features,
  width = 12, 
  height = 10
)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# SCALING THE DATA
# Next, we apply a linear transformation (‘scaling’) that is a standard pre-processing step prior to dimensional reduction techniques like PCA. 
all.genes <- rownames(scNanoGPS)
scNanoGPS <- ScaleData(scNanoGPS, features = all.genes)

# In Seurat, we also use the `ScaleData()` function to remove unwanted sources of variation from a single-cell dataset. 
scNanoGPS <- ScaleData(scNanoGPS, vars.to.regress = "percent.mt")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# PERFORM LINEAR DIMENSIONAL REDUCTION
# Next we perform PCA on the scaled data.
scNanoGPS <- RunPCA(scNanoGPS, features = VariableFeatures(object = scNanoGPS))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Get info
cat('Final Seurat object (post filtering) stats: \n') 
scNanoGPS_final_info <- print(scNanoGPS)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Examine and visualize PCA results a few different ways

# PLOT
VizDimLoadings(scNanoGPS, dims = 1:2, reduction = "pca")

# SAVE PLOT
save_plot_viz_dim_loadings <- my_plot_name('viz_dim_loadings')
my_plot_save(
  save_plot_viz_dim_loadings,
  width     = 14
)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# PLOT
DimPlot(scNanoGPS, reduction = "pca") + NoLegend()

# SAVE PLOT
save_plot_dim <- my_plot_name('dim')
my_plot_save(save_plot_dim)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# PLOT
DimHeatmap(scNanoGPS, dims = 1, cells = 500, balanced = TRUE, fast = FALSE)

#* NOTE - Setting `fast = FALSE` is crucial for generating a customizable ggplot object that can be saved. 
#* When fast = TRUE, the function uses the base R image() function, which is faster but doesn't return a plot object.

# SAVE PLOT
save_plot_dim_heatmap <- my_plot_name('dim_heatmap')
my_plot_save(save_plot_dim_heatmap)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# PLOT
DimHeatmap(scNanoGPS, dims = 1:15, cells = 500, balanced = TRUE, fast = FALSE)

# SAVE PLOT
save_plot_dim_heatmap_15 <- my_plot_name('dim_heatmap_15')
my_plot_save(
  save_plot_dim_heatmap_15,
  width = 16,
  height = 20
)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# DETERMINE THE 'DIMENSIONALITY' OF THE DATASET
# PLOT 
ElbowPlot(scNanoGPS)

# SAVE PLOT
save_plot_elbow <- my_plot_name('elbow')
my_plot_save(save_plot_elbow)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# CLUSTER THE CELLS
scNanoGPS <- FindNeighbors(scNanoGPS, dims = 1:10)

scNanoGPS <- FindClusters(scNanoGPS, resolution = 0.2)

# Look at cluster IDs of the first 5 cells
head(Idents(scNanoGPS), 5)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# RUN NON-LINEAR DIMENSIONAL REDUCTION (UMAP/tSNE)
scNanoGPS <- RunUMAP(scNanoGPS, dims = 1:10)

# PLOT
DimPlot(scNanoGPS, reduction = "umap", label = TRUE)

# SAVE PLOT
save_plot_umap <- my_plot_name('umap')
my_plot_save(save_plot_umap)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# FINDING DIFFERENTIALLY EXPRESSED FEATURES (CLUSTER BIOMARKERS)
# Find all markers of cluster 2
cluster2.markers <- FindMarkers(scNanoGPS, ident.1 = 1)
cat('cluster2.markers, n = 5 : \n')
head(cluster2.markers, n = 5)

# Find markers for every cluster compared to all remaining cells, report only the positive ones
scNanoGPS.markers <- FindAllMarkers(scNanoGPS, only.pos = TRUE)

scNanoGPS.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

cluster0.markers <- FindMarkers(scNanoGPS, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

# CHECK
# scNanoGPS@assays[["isoform"]]@meta.data[["var.features"]]
write.table(
  scNanoGPS@assays[["isoform"]]@meta.data[["var.features"]], 
  file = file.path(input_dir, "isoform_var_features.tsv"),
  sep = "\t",
  quote = FALSE
)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Import celltype marker files from Monocle folder
celltype_marker_dir <- args$celltype_marker_dir
celltype_marker_files_list <- list.files(path = celltype_marker_dir, full.names = TRUE)

# Full genes list
#* NOTE - make sure to set the `celltype` to 'all' for the celltype_marker_loop if you use this list
celltype_marker_full_list <- '/vf/users/CARD_singlecell/MONOCLE_V3/INPUTS/celltype_marker_tables-marker_table0124.csv'

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# FUNCTION - loop thru all the celltype marker files for the `genes_list` parameter
celltype_marker_loop <- function(
  celltype_marker_files_list  = celltype_marker_files_list,
  plotting_operation          = plotting_operation,
  celltype                    = '',
  plot_name
) 
{ # Loop through each file and read it from column 1
  for (file in celltype_marker_files_list) {
    
    # Read the CSV file, starting from column 1
      
    # Extract col 1 ('Gene')
    data <- read.csv(file)
    gene_col <- data[['Gene']]
    cat('gene_col = \n')
    cat(gene_col, '\n')
    
      # NOTE: if `celltype = 'all'` is provided as function input, use that as manual override, otherwise extract from single celltype files
      # Extract col 2's header (celltype) ONLY
      if (celltype == 'all') {
        celltype <- 'all'
        
      } else {
        celltype <- names(data)[2]
      }

      # CHECK
      cat('celltype =', celltype, '\n')
    
    # You can now work with the data
    # PLOT
    result <- try({
      eval(plotting_operation) 
    })
    
    # CHECK - if an error occurred
    if (inherits(result, "try-error")) {
      cat("ERROR! With plotting", celltype, ': \n', result, '\n')
      
    } else {
      cat("SUCCESS! Plotted", celltype, ': \n')
      plot(result)
      
      # SAVE PLOT
      save_plot_gene_expression <- my_plot_name(paste0(plot_name, '-', celltype)) 
      my_plot_save(
        save_plot_gene_expression,
        width     = 16,
        height    = 10
      )
    }
  }
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# PLOT
plotting_operation <- expression(VlnPlot(scNanoGPS, features  = gene_col))

#! WARNING, needs to include the defaults even tho I provided it in the function's defaults
celltype_marker_loop(celltype_marker_files_list, plotting_operation, plot_name = 'gene_expression')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# PLOT
plotting_operation <- expression(VlnPlot(scNanoGPS, features = gene_col, slot = "counts", log = TRUE))
celltype_marker_loop(celltype_marker_files_list, plotting_operation, plot_name = 'gene_expression_raw_counts')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# PLOT
plotting_operation <- expression(FeaturePlot(scNanoGPS, features = gene_col))
celltype_marker_loop(celltype_marker_files_list, plotting_operation, plot_name = 'gene_expression_feature')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# `DoHeatmap()` ( https://satijalab.org/seurat/reference/doheatmap ) generates an expression heatmap for given cells and features
scNanoGPS.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10

# PLOT
DoHeatmap(scNanoGPS, features = top10$gene) 

# SAVE PLOT
save_plot_gene_expression_heatmap_top10 <- my_plot_name('gene_expression_heatmap_top10')
my_plot_save(
  save_plot_gene_expression_heatmap_top10,
  height    = 10
)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# PLOT - w/all marker genes
plotting_operation <- expression(DotPlot(scNanoGPS, features = gene_col) + RotatedAxis())
celltype_marker_loop(celltype_marker_full_list, plotting_operation, plot_name = 'gene_expression', celltype = 'all')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# PLOT - w/all marker genes
plotting_operation <- expression(DoHeatmap(scNanoGPS, features = gene_col, size = 3, draw.lines = T) + RotatedAxis())
celltype_marker_loop(celltype_marker_full_list, plotting_operation, plot_name = 'gene_expression_heatmap_celltype_markers', celltype = 'all')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# MANUAL INSPECTION!!!
# - See the UMAP plot for the #'s, see the dot/feature plots for the gene correlations
# - XYLENA: labeling all the clusters might not be necessary until you want to publish/present
# - XYLENA: the dotplot, heatmap, and UMAP should be enough to look at otherwise

# SAVE POINT
# OUTPUT - .RDS PATH
rds_path <- file.path(output_path, 'scNanoGPS.rds')
cat('rds_path = \n', rds_path, '\n')

# EXPORT - .RDS file
saveRDS(scNanoGPS, file = rds_path)

# # RELOAD - .RDS (if needed)
# scNanoGPS <- readRDS(file = rds_path)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# # ASSIGNING CELL TYPE IDENTITY TO CLUSTERS 
# names(new.cluster.ids) <- levels(scNanoGPS)
# scNanoGPS <- RenameIdents(scNanoGPS, new.cluster.ids)
# 
# # PLOT - w/relabeled cluster names
# # DimPlot(scNanoGPS, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
# DimPlot(scNanoGPS, reduction = "umap", label = TRUE, pt.size = 0.5)
#
# # #! ERROR
# # Error in names(new.cluster.ids) <- levels(scNanoGPS) : 
# #   'names' attribute [11] must be the same length as the vector [9]
# 
# # SAVE PLOT
# save_plot_umap_cell_types <- my_plot_name('umap_cell_types')
# my_plot_save(save_plot_umap_cell_types)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# LOGGING 
# - save job metadata, user-provided inputs, and important calculated variables to a plaintext log file

# CREATE - filename + path where you want to save the log
my_log_file <- file.path(output_path, paste0('log-', current_date, '_', current_time, '.txt'))
cat('Creating log file at: \n', my_log_file, '\n')

# OPEN SINK (for logging terminal outputs to a text file)
sink(my_log_file)

# JOB INFO:
script_name <- deparse(substitute(sys.call(-1)))
cat("\nJOB INFO: \n")
cat(paste0("current_datetime  = ", current_date, "_", current_time, " \n"))
cat(paste0("script_name       = ", script_name, " \n"))

# ARGPARSE USER INPUTS
cat('\nUSER-PROVIDED ARGUMENTS: \n')
print(args)

# VARIABLES CALCULATED IN SCRIPT
cat('\nVARIABLES CALCULATED IN SCRIPT: \n')
cat('input_dir = \n', input_dir, '\n')
cat('output_path = \n', output_path, '\n')

cat('\ntop 10 = \n')
paste(top10)

cat('celltype_marker_files_list: \n')
paste(celltype_marker_files_list)
cat('Initial Seurat object created stats: \n')
print(scNanoGPS_initial_info)
cat('Final Seurat object (post filtering) stats: \n')
print(scNanoGPS_final_info)

cat('\nSaved .RDS to: \n', rds_path, '\n')

cat('\nsessionInfo: \n')
sessionInfo()

# CLOSE SINK
sink(NULL)