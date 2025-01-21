# SEURAT TUTORIAL
# SOURCE: https://satijalab.org/seurat/articles/pbmc3k_tutorial

# NOTES:
# - SCRIPT PATH: /vf/users/CARDPB/data/snRNA_longread/eugene-seurat/src/SEURAT-for_short_reads.R
# - CARD folder path: /data/CARD_singlecell/Brain_atlas/NABEC_multiome

# We start by reading in the data. 
# The `Read10X()` function reads in the output of the `cellranger` 
# ( https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger ) 
# pipeline from 10X, returning a unique molecular identified (UMI) count matrix. 
# The values in this matrix represent the number of molecules for each feature (i.e. gene; row) 
# that are detected in each cell (column). Note that more recent versions of cellranger now also output using the 
# `h5 file format` ( https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/h5_matrices ), 
# which can be read in using the `Read10X_h5()` function in Seurat.

# We next use the count matrix to create a `Seurat` object. 
# The object serves as a container that contains both data (like the count matrix) and 
# analysis (like PCA, or clustering results) for a single-cell dataset. 
# For more information, check out our [Seurat object interaction vignette], or our 
# `GitHub Wiki` ( https://github.com/satijalab/seurat/wiki ). 
# For example, in Seurat v5, the count matrix is stored in `pbmc[["RNA"]]$counts`.

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# IMPORTS 
require(dplyr)
require(Seurat)
require(patchwork)
require(base)
require(ggplot2) 

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# SETUP

# #* TODO - NEED TO CHANGE/PUT IN ARGPARSE!!!
# input_dir
# output_dir
# id 

# LOGGING
current_date  <- as.character(Sys.Date())
current_time  <- format(Sys.time(),"%H-%M-%S")

# # DEFINE PATHS
# # Input path - takes a `filtered_feature_bc_matrix.h5`` file from Arc Cellranger output folder
# input_dir <- '/data/CARD_singlecell/Brain_atlas/Cortex_opt/Psomagen102722/4546Cortex/outs/filtered_feature_bc_matrix.h5'
# input_dir <- '/data/CARD_singlecell/Brain_atlas/Cortex_opt/Psomagen102722/SH-92-05/outs/filtered_feature_bc_matrix.h5'
input_dir <- '/data/CARD_singlecell/Brain_atlas/NABEC_multiome/batch5/Multiome/SH-04-08-ARC/outs/filtered_feature_bc_matrix.h5'
cat('input_dir =', input_dir, '\n')

#* TODO - fix the ID section OR make into argparse
id <- 'SH-04-08'
cat('id =', id, '\n')

# # Output path
# output_dir <- input_dir
output_dir <- '/vf/users/CARDPB/data/snRNA_longread/eugene-seurat/output'
output_path <- file.path(output_dir, 'short_reads', id, 'PLOTS-Seurat-short_reads')

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

#* TODO - Cory still said that I'd need to add the patient IDs to the barcode column

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# LOAD DATA

# SETUP THE SEURAT OBJECT
seurat_object.data <- Read10X_h5(input_dir)

# Initialize the Seurat object with the raw (non-normalized data).
# pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
seurat_object <- CreateSeuratObject(counts = seurat_object.data, project = id, min.cells = 3, min.features = 200)

# Get info
seurat_object_initial_info <- print(seurat_object)
cat('Initial Seurat object created stats: \n')
seurat_object_initial_info  
print(seurat_object_initial_info)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# STANDARD PRE-PREPROCESSING WORKFLOW

# The steps below encompass the standard pre-processing workflow for scRNA-seq data in Seurat. 
# These represent the selection and filtration of cells based on QC metrics, data normalization and scaling, 
# and the detection of highly variable features.

# # QC AND SELECTING CELLS FOR FURTHER ANALYSIS
# Seurat allows you to easily explore QC metrics and filter cells based on any user-defined criteria. 
# A few QC metrics commonly used by the community include

# - The number of unique genes detected in each cell.
#     o Low-quality cells or empty droplets will often have very few genes
#     o Cell doublets or multiplets may exhibit an aberrantly high gene count
# - Similarly, the total number of molecules detected within a cell (correlates strongly with unique genes)
# - The percentage of reads that map to the mitochondrial genome
#     o Low-quality / dying cells often exhibit extensive mitochondrial contamination
#     o We calculate mitochondrial QC metrics with the `PercentageFeatureSet()` function, which calculates the percentage of counts originating from a set of features
#     o We use the set of all genes starting with MT- as a set of mitochondrial genes

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^MT-")

# # Where are QC metrics stored in Seurat?
# - The number of unique genes and total molecules are automatically calculated during `CreateSeuratObject()`
#   o You can find them stored in the object meta data
# # Show QC metrics for the first 5 cells
head(seurat_object@meta.data, 5)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# MY FUNCTIONS

# FILE NAMING - Construct the file name
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

# In the example below, we visualize QC metrics, and use these to filter cells.
#   - We filter cells that have unique feature counts > 2,500 or < 200
#   - We filter cells that have >5% mitochondrial counts

# Visualize QC metrics as a violin plot
VlnPlot(seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# SAVE PLOT
save_plot_violin_counts_features <- my_plot_name('violin_counts_features')
my_plot_save(save_plot_violin_counts_features)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# SAVE PLOT
save_plot_feature_scatter <- my_plot_name('feature_scatter')
my_plot_save(
  save_plot_feature_scatter,
  height = 4,
  width = 10
)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# # FILTER MANUALLY - based on the violin plot
# seurat_object <- subset(seurat_object, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
seurat_object <- subset(seurat_object, subset = nFeature_RNA > 0 & nFeature_RNA < 12000 & percent.mt < 5)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# # NORMALIZING THE DATA

# # After removing unwanted cells from the dataset, the next step is to normalize the data. 
# # By default, we employ a global-scaling normalization method “LogNormalize” that normalizes 
# # the feature expression measurements for each cell by the total expression, multiplies this by a scale factor 
# # (10,000 by default), and log-transforms the result. In Seurat v5, Normalized values are stored in `seurat_object[["RNA"]]$data`.
# seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = 10000)

# For clarity, in this previous line of code (and in future commands), we provide the default values 
# for certain parameters in the function call. However, this isn’t required and the same behavior can be achieved with:
seurat_object <- NormalizeData(seurat_object)

# While this method of normalization is standard and widely used in scRNA-seq analysis, 
# global-scaling relies on an assumption that each cell originally contains the same number of RNA molecules. 
# We and others have developed alternative workflows for the single cell preprocessing that do not make these assumptions. 
# For users who are interested, please check out our `SCTransform()` normalization workflow. 
# The method is described in our `paper` ( https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02584-9 ), 
# with a separate vignette using Seurat `here.` ( https://satijalab.org/seurat/articles/sctransform_vignette ). 
# The use of `SCTransform` replaces the need to run `NormalizeData`, `FindVariableFeatures`, or `ScaleData` (described below.)

# Identification of highly variable features (feature selection)
# We next calculate a subset of features that exhibit high cell-to-cell variation in the dataset 
# (i.e, they are highly expressed in some cells, and lowly expressed in others). 
# We and `others`( https://www.nature.com/articles/nmeth.2645 ) have found that focusing on 
# these genes in downstream analysis helps to highlight biological signal in single-cell datasets.

# Our procedure in Seurat is described in detail `here ( https://doi.org/10.1016/j.cell.2019.05.031 )`, 
# and improves on previous versions by directly modeling the mean-variance relationship inherent in single-cell data, 
# and is implemented in the `FindVariableFeatures()` function. By default, we return 2,000 features per dataset. 
# These will be used in downstream analysis, like PCA.
# seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2000)
seurat_object <- FindVariableFeatures(seurat_object, layer = 'data.Gene Expression', selection.method = "vst", nfeatures = 2000)

#! WARNING - MUST SPECIFY LAYER IN SEURAT V5
# `data.` layer is the normalized data, use this!

#* TODO - do i need to specify the layer again?
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seurat_object, layer = 'data.Gene Expression'), 10)
cat('top 10 = \n')
paste(top10)

# Plot variable features with and without labels
plot1 <- VariableFeaturePlot(seurat_object)
# plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2 <- LabelPoints(
  plot = plot1,
  points = top10,
  repel = TRUE,
  xnudge = 0,
  ynudge = 0
)
plot1 + plot2

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# SAVE PLOT
save_plot_variable_features <- my_plot_name('variable_features')
my_plot_save(
  save_plot_variable_features,
  width = 12, 
  height = 10
)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# SCALING THE DATA

# Next, we apply a linear transformation (‘scaling’) that is a standard pre-processing step 
# prior to dimensional reduction techniques like PCA. 
# The `ScaleData()` function:
# - Shifts the expression of each gene, so that the mean expression across cells is 0
# - Scales the expression of each gene, so that the variance across cells is 1
#         o This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate
# - The results of this are stored in `pbmc[["RNA"]]$scale.data`
# - By default, only variable features are scaled.
# - You can specify the `features` argument to scale additional features
all.genes <- rownames(seurat_object)
seurat_object <- ScaleData(seurat_object, features = all.genes)

# # How can I remove unwanted sources of variation
# # In Seurat, we also use the `ScaleData()` function to remove unwanted sources of variation from a single-cell dataset. 
# # For example, we could ‘regress out’ heterogeneity associated with (for example) 
# # `cell cycle stage` ( https://satijalab.org/seurat/articles/cell_cycle_vignette ), or mitochondrial contamination i.e.:
# seurat_object <- ScaleData(seurat_object, vars.to.regress = "percent.mt")

#! ERROR
# Error in LayerData.Assay5(X[[i]], ...) : features are not found

# However, particularly for advanced users who would like to use this functionality, we strongly recommend the use of our new normalization workflow, `SCTransform()`. 
# The method is described in our `paper` ( https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02584-9 ), 
# with a separate vignette using Seurat `here` ( https://satijalab.org/seurat/articles/sctransform_vignette ). 
# As with `ScaleData()`, the function `SCTransform()` also includes a `vars.to.regress` parameter.

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# PERFORM LINEAR DIMENSIONAL REDUCTION

# Next we perform PCA on the scaled data. By default, only the previously determined variable features are used as input, 
# but can be defined using features argument if you wish to choose a different subset 
# (if you do want to use a custom subset of features, make sure you pass these to `ScaleData` first).

# For the first principal components, Seurat outputs a list of genes with the most positive and negative loadings, 
# representing modules of genes that exhibit either correlation (or anti-correlation) across single-cells in the dataset.
seurat_object <- RunPCA(seurat_object, features = VariableFeatures(object = seurat_object))

# Seurat provides several useful ways of visualizing both cells and features that define the PCA, including:
# - `VizDimReduction()`, 
# - `DimPlot()`, and 
# - `DimHeatmap()`

# Examine and visualize PCA results a few different ways
print(seurat_object[["pca"]], dims = 1:5, nfeatures = 5)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Get info
seurat_object_final_info <- print(seurat_object)
cat('Final Seurat object (post filtering) stats: \n')
seurat_object_final_info  
print(seurat_object_final_info)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# PLOT
VizDimLoadings(seurat_object, dims = 1:2, reduction = "pca")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# SAVE PLOT
save_plot_viz_dim_loadings <- my_plot_name('viz_dim_loadings')
my_plot_save(save_plot_viz_dim_loadings)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# PLOT
DimPlot(seurat_object, reduction = "pca") + NoLegend()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# SAVE PLOT
save_plot_dim <- my_plot_name('dim')
my_plot_save(save_plot_dim)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# In particular `DimHeatmap()` allows for easy exploration of the primary sources of heterogeneity in a dataset,
# and can be useful when trying to decide which PCs to include for further downstream analyses. 
# Both cells and features are ordered according to their PCA scores. 
# Setting `cells` to a number plots the ‘extreme’ cells on both ends of the spectrum, 
# which dramatically speeds plotting for large datasets. 
# Though clearly a supervised analysis, we find this to be a valuable tool for exploring correlated feature sets.
# DimHeatmap(seurat_object, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(seurat_object, dims = 1, cells = 500, balanced = TRUE, fast = FALSE)

#* NOTE - Setting `fast = FALSE` is crucial for generating a customizable ggplot object that can be saved. 
#* When fast = TRUE, the function uses the base R image() function, which is faster but doesn't return a plot object.

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# SAVE PLOT
save_plot_dim_heatmap <- my_plot_name('dim_heatmap')
my_plot_save(save_plot_dim_heatmap)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# PLOT
DimHeatmap(seurat_object, dims = 1:15, cells = 500, balanced = TRUE, fast = FALSE)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# SAVE PLOT
save_plot_dim_heatmap_15 <- my_plot_name('dim_heatmap_15')
my_plot_save(
  save_plot_dim_heatmap_15,
  width = 16,
  height = 20
)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# DETERMINE THE 'DIMENSIONALITY' OF THE DATASET

# To overcome the extensive technical noise in any single feature for scRNA-seq data, Seurat clusters cells based on their PCA scores, 
# with each PC essentially representing a ‘metafeature’ that combines information across a correlated feature set. 
# The top principal components therefore represent a robust compression of the dataset. 
# However, how many components should we choose to include? 10? 20? 100?

# In `Macosko et al` ( http://www.cell.com/abstract/S0092-8674(15)00549-8 ),
# we implemented a resampling test inspired by the JackStraw procedure. 
# While still available in Seurat (`see previous vignette` [ https://satijalab.org/seurat/articles/pbmc3k_tutorial ]), 
# this is a slow and computationally expensive procedure, and we is no longer routinely used in single cell analysis.

# An alternative heuristic method generates an ‘Elbow plot’: a ranking of principle components based on the 
# percentage of variance explained by each one (`ElbowPlot()` [ https://satijalab.org/seurat/reference/elbowplot ] function).
# In this example, we can observe an ‘elbow’ around PC9-10, suggesting that the majority of true signal is captured in the first 10 PCs.
ElbowPlot(seurat_object)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# SAVE PLOT
save_plot_elbow <- my_plot_name('elbow')
my_plot_save(save_plot_elbow)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Identifying the true dimensionality of a dataset – can be challenging/uncertain for the user. 
# We therefore suggest these multiple approaches for users. 
# - The first is more supervised, exploring PCs to determine relevant sources of heterogeneity, 
#   and could be used in conjunction with GSEA for example. 
# - The second (`ElbowPlot`). 
# - The third is a heuristic that is commonly used, and can be calculated instantly. 
# In this example, we might have been justified in choosing anything between PC 7-12 as a cutoff.

# We chose 10 here, but encourage users to consider the following:
# - Dendritic cell and NK aficionados may recognize that genes strongly associated with 
#   PCs 12 and 13 define rare immune subsets (i.e. MZB1 is a marker for plasmacytoid DCs). 
# However, these groups are so rare, they are difficult to distinguish from background noise 
# for a dataset of this size without prior knowledge.
# - We encourage users to repeat downstream analyses with a different number of PCs (10, 15, or even 50!). 
#   As you will observe, the results often do not differ dramatically.
# - We advise users to err on the higher side when choosing this parameter. 
#   For example, performing downstream analyses with only 5 PCs does significantly and adversely affect results.

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# CLUSTER THE CELLS

# Seurat applies a graph-based clustering approach, building upon initial strategies in 
# (`Macosko et al` [ http://www.cell.com/abstract/S0092-8674(15)00549-8 ] ). 
# Importantly, the distance metric which drives the clustering analysis (based on previously identified PCs) remains the same.
# However, our approach to partitioning the cellular distance matrix into clusters has dramatically improved. 
# Our approach was heavily inspired by recent manuscripts which applied graph-based clustering approaches to 
# scRNA-seq data [`SNN-Cliq, Xu and Su, Bioinformatics, 2015` ( http://bioinformatics.oxfordjournals.org/content/early/2015/02/10/bioinformatics.btv088.abstract )] 
# and CyTOF data [`PhenoGraph, Levine et al., Cell, 2015` ( http://www.ncbi.nlm.nih.gov/pubmed/26095251 )]. 
# Briefly, these methods embed cells in a graph structure - for example a K-nearest neighbor (KNN) graph, 
# with edges drawn between cells with similar feature expression patterns, 
# and then attempt to partition this graph into highly interconnected ‘quasi-cliques’ or ‘communities’.

# As in PhenoGraph, we first construct a KNN graph based on the euclidean distance in PCA space, 
# and refine the edge weights between any two cells based on the shared overlap in their local neighborhoods (Jaccard similarity). 
# This step is performed using the `FindNeighbors()` ( https://satijalab.org/seurat/reference/findneighbors ) function, 
# and takes as input the previously defined dimensionality of the dataset (first 10 PCs).

# To cluster the cells, we next apply modularity optimization techniques such as the Louvain algorithm (default) 
# or SLM [`SLM, Blondel et al., Journal of Statistical Mechanics` ( http://dx.doi.org/10.1088/1742-5468/2008/10/P10008 )], 
# to iteratively group cells together, with the goal of optimizing the standard modularity function. 
# The `FindClusters()` function implements this procedure, and contains a resolution parameter that sets the ‘granularity’ 
# of the downstream clustering, with increased values leading to a greater number of clusters. 
# We find that setting this parameter between 0.4-1.2 typically returns good results for single-cell datasets of around 3K cells. 
# Optimal resolution often increases for larger datasets. The clusters can be found using the `Idents()` function.

seurat_object <- FindNeighbors(seurat_object, dims = 1:10)
# seurat_object <- FindClusters(seurat_object, resolution = 0.5)
seurat_object <- FindClusters(seurat_object, resolution = 0.2)

# Look at cluster IDs of the first 5 cells
head(Idents(seurat_object), 5)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# RUN NON-LINEAR DIMENSIONAL REDUCTION (UMAP/tSNE)

# Seurat offers several non-linear dimensional reduction techniques, such as tSNE and UMAP, to visualize and explore these datasets. 
# The goal of these algorithms is to learn underlying structure in the dataset, in order to place similar cells together in low-dimensional space. 
# Therefore, cells that are grouped together within graph-based clusters determined above should co-localize on these dimension reduction plots.

# While we and others have routinely found 2D visualization techniques like tSNE and UMAP to be valuable tools for exploring datasets, 
# all visualization techniques have limitations, and cannot fully represent the complexity of the underlying data. 
# In particular, these methods aim to preserve local distances in the dataset 
# (i.e. ensuring that cells with very similar gene expression profiles co-localize), but often do not preserve more global relationships. 
# We encourage users to leverage techniques like UMAP for visualization, but to avoid drawing biological conclusions solely on the basis of visualization techniques.
seurat_object <- RunUMAP(seurat_object, dims = 1:10)

# # note that you can set `label = TRUE` or use the `LabelClusters` function to help label individual clusters
# DimPlot(seurat_object, reduction = "umap")
DimPlot(seurat_object, reduction = "umap", label = TRUE)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# SAVE PLOT
save_plot_umap <- my_plot_name('umap')
my_plot_save(save_plot_umap)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# FINDING DIFFERENTIALLY EXPRESSED FEATURES (CLUSTER BIOMARKERS)

# Seurat can help you find markers that define clusters via differential expression (DE). 
# By default, it identifies positive and negative markers of a single cluster (specified in `ident.1`), compared to all other cells. 
# `FindAllMarkers()` automates this process for all clusters, but you can also test groups of clusters vs. each other, or against all cells.

# In Seurat v5, we use the presto package (as described `here` ( https://www.biorxiv.org/content/10.1101/653253v1 ) 
# and available for installation `here` [https://github.com/immunogenomics/presto ]), 
# to dramatically improve the speed of DE analysis, particularly for large datasets. 
# For users who are not using presto, you can examine the documentation for this function 
# (`?FindMarkers` [https://satijalab.org/seurat/reference/findmarkers ]) 
# to explore the `min.pct` and `logfc.threshold` parameters, 
# which can be increased in order to increase the speed of DE testing.

# #! ERROR
# Error in FindMarkers.StdAssay(object = data.use, latent.vars = latent.vars,  : 
#   data layers are not joined. Please run JoinLayers
# TEST
seurat_object <- JoinLayers(seurat_object)

# # Find all markers of cluster 2
# cluster2.markers <- FindMarkers(seurat_object, ident.1 = 2)
# head(cluster2.markers, n = 5)

# #! ERROR - if no clusters... breaks script!
# Error in WhichCells.Seurat(object = object, idents = ident.1) : 
#   Cannot find the following identities in the object: 2

# Find markers for every cluster compared to all remaining cells, report only the positive ones
seurat_object.markers <- FindAllMarkers(seurat_object, only.pos = TRUE)
seurat_object.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1)

# TEST
cat('Final Seurat object (post filtering) stats: \n')
print(seurat_object)

# #* TODO - comment out?
# # Seurat has several tests for differential expression which can be set with the `test.use` parameter 
# # (see our `DE vignette` ( https://satijalab.org/seurat/articles/de_vignette ) for details). 
# # For example, the ROC test returns the ‘classification power’ for any individual marker (ranging from 0 - random, to 1 - perfect).
# cluster0.markers <- FindMarkers(seurat_object, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

# We include several tools for visualizing marker expression. 
# - `VlnPlot()` (shows expression probability distributions across clusters), and 
# - `FeaturePlot()` (visualizes feature expression on a tSNE or PCA plot) are our most commonly used visualizations. 
# We also suggest exploring:
# - `RidgePlot()`, 
# - `CellScatter()`, and 
# - `DotPlot()` 
# as additional methods to view your dataset.

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# -------------------------------------------

# Import celltype marker files from Monocle folder
celltype_marker_dir <- '/vf/users/CARD_singlecell/MONOCLE_V3/INPUTS/celltype_marker_tables_split'
celltype_marker_files_list <- list.files(path = celltype_marker_dir, full.names = TRUE)
cat('celltype_marker_files_list: \n')
paste(celltype_marker_files_list)

# Full genes list
#* NOTE - make sure to set the `celltype` to 'all' for the celltype_marker_loop if you use this list
celltype_marker_full_list <- '/vf/users/CARD_singlecell/MONOCLE_V3/INPUTS/celltype_marker_tables-marker_table0124.csv'

# -------------------------------------------

# MY  FUNCTION - loop thru all the celltype marker files for the `genes_list` parameter
celltype_marker_loop <- function(
  celltype_marker_files_list  = celltype_marker_files_list,
  plotting_operation          = plotting_operation,
  celltype                    = 'all',
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

      # # NOTE: if `celltype` provided as function input, use that as manual override, otherwise extract from single celltype files
      # # Extract col 2's header (celltype)
      # if (celltype == 'all') {
      # } else {
      #   celltype <- names(data)[2]
      # }
      
      # TEST - it was naming celltype as 'all' but still plotting all the right stuff, shouldn't it just depend on the list I use?
      celltype <- names(data)[2]
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
        #* TODO - fix `celltype` set to `all` instead of each actual cell
        #* 
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

# # Test the gene lists
# VlnPlot(seurat_object, features = genes_list)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# # SAVE PLOT
# save_plot_gene_expression <- my_plot_name('gene_expression')
# my_plot_save(save_plot_gene_expression)

# -------------------------------------------

# RUN
# Define the plotting operation BEFORE each celltype marker plotting call
plotting_operation <- expression(VlnPlot(seurat_object, features  = gene_col))

#! WARNING, needs to include the defaults even tho I provided it in the function's defaults
celltype_marker_loop(celltype_marker_files_list, plotting_operation, plot_name = 'gene_expression')

# -------------------------------------------

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# you can plot raw counts as well
# # VlnPlot(seurat_object, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)
# VlnPlot(seurat_object, features = genes_list, slot = "counts", log = TRUE)

# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# # SAVE PLOT
# save_plot_gene_expression_raw_counts <- my_plot_name('gene_expression_raw_counts')
# my_plot_save(save_plot_gene_expression_raw_counts)

# -------------------------------------------

# RUN
# Define the plotting operation BEFORE each celltype marker plotting call
plotting_operation <- expression(VlnPlot(seurat_object, features = gene_col, slot = "counts", log = TRUE))
celltype_marker_loop(celltype_marker_files_list, plotting_operation, plot_name = 'gene_expression_raw_counts')

# -------------------------------------------

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# FeaturePlot(seurat_object, features = genes_list)

# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# # SAVE PLOT
# save_plot_gene_expression_feature <- my_plot_name('gene_expression_feature')
# my_plot_save(save_plot_gene_expression_feature)

# -------------------------------------------

# RUN
# Define the plotting operation BEFORE each celltype marker plotting call
plotting_operation <- expression(FeaturePlot(seurat_object, features = gene_col))
celltype_marker_loop(celltype_marker_files_list, plotting_operation, plot_name = 'gene_expression_feature')

# -------------------------------------------

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# `DoHeatmap()` ( https://satijalab.org/seurat/reference/doheatmap ) generates an expression heatmap for given cells and features. 
# In this case, we are plotting the top 20 markers (or all markers if less than 20) for each cluster.
seurat_object.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 10) %>%
    ungroup() -> top10

# PLOT - heatmap
# DoHeatmap(seurat_object, features = top10$gene) + NoLegend()
DoHeatmap(seurat_object, features = top10$gene) 

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# SAVE PLOT
save_plot_gene_expression_heatmap_top10 <- my_plot_name('gene_expression_heatmap_top10')
my_plot_save(save_plot_gene_expression_heatmap_top10)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#* NOTE - turn off in argparse
# marker.genes = c(
#   'RORB',
#   'SLC17A6',
#   'SLC17A7',
#   'THEMIS',
#   'GAD1',
#   'GAD2',
#   'PVALB',
#   'SST',
#   'VIP',
#   'APBB1IP',
#   'CD74',
#   'CSF1R',
#   'CX3CR1',
#   'ITGAM',
#   'P2RY12',
#   'PTPRC',
#   'ALDH1L1',
#   'AQP4',
#   'COL5A3',
#   'GFAP',
#   'SLC1A2',
#   'SLC1A3',
#   'CLDN11',
#   'MBP',
#   'MOBP',
#   'OPALIN',
#   'PLP1',
#   'ST18',
#   'LHFPL3',
#   'MEGF11',
#   'PCDH15',
#   'PDGFRA',
#   'VCAN',
#   'CLDN5',
#   'COLEC12',
#   'EPAS1',
#   'VCAM1'
# )

# FROM - SARAH
# # PLOT
# DotPlot(seurat_object, features = marker.genes) + RotatedAxis()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# # SAVE PLOT
# save_plot_gene_expression_dot_plot<- my_plot_name('gene_expression_dot_plot')
# my_plot_save(save_plot_gene_expression_dot_plot)

# -------------------------------------------

#* TODO - fix the naming which adds a '-ExN' like it's the celltype
# RUN - w/all marker genes
plotting_operation <- expression(DotPlot(seurat_object, features = gene_col) + RotatedAxis())
celltype_marker_loop(celltype_marker_full_list, plotting_operation, plot_name = 'gene_expression', celltype = 'all')

# -------------------------------------------

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# # FROM - SARAH
# # DoHeatmap(subset(SH0106_clean), features = marker.genes, size = 3, draw.lines = F)
# DoHeatmap(seurat_object, features = marker.genes, size = 3, draw.lines = T)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# # SAVE PLOT
# save_plot_gene_expression_heatmap_celltype_markers <- my_plot_name('')
# my_plot_save(save_plot_gene_expression_heatmap_celltype_markers)

# -------------------------------------------

#* TODO - fix the naming which adds a '-ExN' like it's the celltype
# RUN - w/all marker genes
plotting_operation <- expression(DoHeatmap(seurat_object, features = gene_col, size = 3, draw.lines = T) + RotatedAxis())
celltype_marker_loop(celltype_marker_full_list, plotting_operation, plot_name = 'gene_expression_heatmap_celltype_markers', celltype = 'all')

# -------------------------------------------

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# NOTE - MANUAL INSPECTION!!!
# - See the UMAP plot for the #'s, see the dot/feature plots for the gene correlations
# - XYLENA: labeling all the clusters might not be necessary until you want to publish/present
# - XYLENA: the dotplot, heatmap, and UMAP should be enough to look at otherwise

# SAVE POINT
# OUTPUT - .RDS PATH
rds_path <- file.path(output_path, 'seurat_object.rds')
cat('rds_path = \n', rds_path, '\n')

# EXPORT - .RDS file
saveRDS(seurat_object, file = rds_path)

# # RELOAD - .RDS (if needed)
# seurat_object <- readRDS(file = rds_path)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# # ASSIGNING CELL TYPE IDENTITY TO CLUSTERS 
# # Fortunately in the case of this dataset, we can use canonical markers to easily match the unbiased clustering to known cell types:
# new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono", "NK", "DC", "Platelet")

# # MY LIST (for SH-04-08):
# 0 = 'Unknown'
# 1 = 'ExN'
# 2 = 'ExN'
# 3 = 'InN'
# 4 = 'InN'
# 5 = 'Unknown'
# 6 = 'Unknown'
# 7 = 'Astro'
# 8 = 'Unknown'
# 9 = 'OPC'
# 10 = 'Unknown'

# new.cluster.ids <- c(
#   'Unknown',
#   'ExN',
#   'ExN',
#   'InN',
#   'InN',
#   'Unknown',
#   'Unknown',
#   'Astro',
#   'Unknown',
#   'OPC',
#   'Unknown'
# )
# 
# names(new.cluster.ids) <- levels(seurat_object)
# seurat_object <- RenameIdents(seurat_object, new.cluster.ids)
# 
# # PLOT - w/relabeled cluster names
# # DimPlot(seurat_object, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
# DimPlot(seurat_object, reduction = "umap", label = TRUE, pt.size = 0.5)
# # #! ERROR
# # Error in names(new.cluster.ids) <- levels(seurat_object) : 
# #   'names' attribute [11] must be the same length as the vector [9]
# 
# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# # SAVE PLOT
# save_plot_umap_cell_types <- my_plot_name('umap_cell_types')
# my_plot_save(save_plot_umap_cell_types)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#* TODO - export some more metadata??? put that and sessionInfo into a text file? or cat to slurm???
# CHECK - Session Info
cat('\nsessionInfo: \n')
sessionInfo()

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
cat('top 10 = \n')
paste(top10)
cat('celltype_marker_files_list: \n')
paste(celltype_marker_files_list)

# Get info
cat('Initial Seurat object created stats: \n')
seurat_object_initial_info  
print(seurat_object_initial_info)

# Get info
cat('Final Seurat object (post filtering) stats: \n')
seurat_object_final_info  
print(seurat_object_final_info)

cat('Saved .RDS to: \n', rds_path, '\n')
cat('\nsessionInfo: \n')
sessionInfo()

# CLOSE SINK
sink(NULL)