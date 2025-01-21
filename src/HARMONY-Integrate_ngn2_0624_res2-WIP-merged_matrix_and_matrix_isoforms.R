#!/usr/bin/env Rscript
# start like this
# Rscript --vanilla Integrate_ngn2_0624_res2.R

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#* TODO - not sure which libraries I'll need
# LOAD LIBRARIES
require(Seurat)
require(dplyr)
require(patchwork)
require(ggplot2)
# require(fields)
require(KernSmooth)
require(Matrix)
require(parallelly)
require(ROCR)
# require(parallel, lib.loc = "/require/Frameworks/R.framework/Versions/4.2/Resources/require")
require(parallel)
require(tidyverse)
require(cowplot)
# require(SeuratObject)
# require(future)
require(future.apply)
# require(harmony, lib.loc='/data/reedx/R/rhel8/4.3/')
require(harmony)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# SETUP

# LOGGING
current_date  <- as.character(Sys.Date())
current_time  <- format(Sys.time(),"%H-%M-%S")

# # DEFINE paths
# # Input path - takes the results of converting the `/scNanoGPS_res` folder into: barcodes.tsv, genes.tsv, and matrix.mtx
# # MATRIX FOLDERS
# input_dir_1 <- '/data/CARDPB/data/snRNA_longread/eugene-seurat/output/SH-04-08/20241107'
# input_dir_2 <- '/data/CARDPB/data/snRNA_longread/eugene-seurat/output/SH-06-25/20241107'
# input_dir_3 <- '/data/CARDPB/data/snRNA_longread/eugene-seurat/output/SH-07-46/20241107'
# input_dir_4 <- '/data/CARDPB/data/snRNA_longread/eugene-seurat/output/SH-92-05/20241122'
# input_dir_5 <- '/data/CARDPB/data/snRNA_longread/eugene-seurat/output/UMARY_4546/20241106'

# # ISOFORM FOLDERS
# input_dir_1_isoform <- file.path(input_dir_1, 'isoform')
# input_dir_2_isoform <- file.path(input_dir_2, 'isoform')
# input_dir_3_isoform <- file.path(input_dir_3, 'isoform')
# input_dir_4_isoform <- file.path(input_dir_4, 'isoform')
# input_dir_5_isoform <- file.path(input_dir_5, 'isoform')

# MERGED FOLDERS
input_dir_1_merged <- '/vf/users/CARDPB/data/snRNA_longread/eugene-seurat/output/merged/SH-04-08'
input_dir_2_merged <- '/vf/users/CARDPB/data/snRNA_longread/eugene-seurat/output/merged/SH-06-25'
input_dir_3_merged <- '/vf/users/CARDPB/data/snRNA_longread/eugene-seurat/output/merged/SH-04-08'
input_dir_4_merged <- '/vf/users/CARDPB/data/snRNA_longread/eugene-seurat/output/merged/SH-07-46'
input_dir_5_merged <- '/vf/users/CARDPB/data/snRNA_longread/eugene-seurat/output/merged/UMARY_4546'

# #* TODO - loop all inputs in list?
# input_dir_list <- c(input_dir_1, input_dir_2, input_dir_3, input_dir_4, input_dir_5)
# cat('input_dir_list = \n')
# paste(input_dir_list)

# Output path - you gotta provide one, also change the name for diff inputs
output_dir <- '/data/CARDPB/data/snRNA_longread/eugene-seurat/output/short_reads'

# ---------------------------------------------------------------------
# # MERGED OUTPUT - uncomment this section to send results to the isoform folder
# output_path <- file.path(output_dir, 'merged', 'PLOTS-Seurat-matrix')
output_path <- file.path(output_dir, 'merged', 'PLOTS-Seurat-merged_matrix_and_matrix_isoform')
# ---------------------------------------------------------------------

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

# # # LOAD DATA
# # ipsc1 <- Read10X_h5("/data/CARD_singlecell/NGN2_diff_Erika/data/sample_iPSC_1/outs/filtered_feature_bc_matrix.h5")
# scNanoGPS.data1 <- Read10X(data.dir = input_dir_1)
# scNanoGPS.data2 <- Read10X(data.dir = input_dir_2)
# scNanoGPS.data3 <- Read10X(data.dir = input_dir_3)
# scNanoGPS.data4 <- Read10X(data.dir = input_dir_4)
# scNanoGPS.data5 <- Read10X(data.dir = input_dir_5)

# scNanoGPS.data1_isoform <- Read10X(data.dir = input_dir_1_isoform, gene.column=1)
# scNanoGPS.data2_isoform <- Read10X(data.dir = input_dir_2_isoform, gene.column=1)
# scNanoGPS.data3_isoform <- Read10X(data.dir = input_dir_3_isoform, gene.column=1)
# scNanoGPS.data4_isoform <- Read10X(data.dir = input_dir_4_isoform, gene.column=1)
# scNanoGPS.data5_isoform <- Read10X(data.dir = input_dir_5_isoform, gene.column=1)

scNanoGPS.data1_merged <- Read10X(data.dir = input_dir_1_merged, gene.column=1)
scNanoGPS.data2_merged <- Read10X(data.dir = input_dir_2_merged, gene.column=1)
scNanoGPS.data3_merged <- Read10X(data.dir = input_dir_3_merged, gene.column=1)
scNanoGPS.data4_merged <- Read10X(data.dir = input_dir_4_merged, gene.column=1)
scNanoGPS.data5_merged <- Read10X(data.dir = input_dir_5_merged, gene.column=1)

# #! ERROR
# Error in Read10X(data.dir = input_dir_4_merged, gene.column = 1) : 
#   Barcode file missing. Expecting barcodes.tsv.gz

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#* TODO - getting the 'id' here would have been nice for naming them

# # # CREATE SEURAT OBJECTS
# # ipsc1 <- CreateSeuratObject(counts = ipsc1, project = "ipsc1", min.cells = 3, min.features = 200)
# scNanoGPS_1 <- CreateSeuratObject(counts = scNanoGPS.data1, project = 'SH-04-08')
# scNanoGPS_2 <- CreateSeuratObject(counts = scNanoGPS.data2, project = 'SH-06-25')
# scNanoGPS_3 <- CreateSeuratObject(counts = scNanoGPS.data3, project = 'SH-07-46')
# scNanoGPS_4 <- CreateSeuratObject(counts = scNanoGPS.data4, project = 'SH-92-05')
# scNanoGPS_5 <- CreateSeuratObject(counts = scNanoGPS.data5, project = 'UMARY_4546')

# scNanoGPS_1_isoform  <- CreateSeuratObject(counts = scNanoGPS.data1_isoform , project = 'SH-04-08_isoform')
# scNanoGPS_2_isoform  <- CreateSeuratObject(counts = scNanoGPS.data2_isoform , project = 'SH-06-25_isoform')
# scNanoGPS_3_isoform  <- CreateSeuratObject(counts = scNanoGPS.data3_isoform , project = 'SH-07-46_isoform')
# scNanoGPS_4_isoform  <- CreateSeuratObject(counts = scNanoGPS.data4_isoform , project = 'SH-92-05_isoform')
# scNanoGPS_5_isoform  <- CreateSeuratObject(counts = scNanoGPS.data5_isoform , project = 'UMARY_4546_isoform')

scNanoGPS_1_merged <- CreateSeuratObject(counts = scNanoGPS.data1_merged , project = 'SH-04-08_merged')
scNanoGPS_2_merged <- CreateSeuratObject(counts = scNanoGPS.data2_merged , project = 'SH-06-25_merged')
scNanoGPS_3_merged <- CreateSeuratObject(counts = scNanoGPS.data3_merged , project = 'SH-07-46_merged')
scNanoGPS_4_merged <- CreateSeuratObject(counts = scNanoGPS.data4_merged , project = 'SH-92-05_merged')
scNanoGPS_5_merged <- CreateSeuratObject(counts = scNanoGPS.data5_merged , project = 'UMARY_4546_merged')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Get info
cat('Initial Seurat object created stats: \n') 
scNanoGPS_1_info_initial <- print(scNanoGPS_1)
scNanoGPS_2_info_initial <- print(scNanoGPS_2)
scNanoGPS_3_info_initial <- print(scNanoGPS_3)
scNanoGPS_4_info_initial <- print(scNanoGPS_4)
scNanoGPS_5_info_initial <- print(scNanoGPS_5)

scNanoGPS_1_isoform_info_initial <- print(scNanoGPS_1_isoform)
scNanoGPS_2_isoform_info_initial <- print(scNanoGPS_2_isoform)
scNanoGPS_3_isoform_info_initial <- print(scNanoGPS_3_isoform)
scNanoGPS_4_isoform_info_initial <- print(scNanoGPS_4_isoform)
scNanoGPS_5_isoform_info_initial <- print(scNanoGPS_5_isoform)

scNanoGPS_1_merged_info_initial <- print(scNanoGPS_1_merged)
scNanoGPS_2_merged_info_initial <- print(scNanoGPS_2_merged)
scNanoGPS_3_merged_info_initial <- print(scNanoGPS_3_merged)
scNanoGPS_4_merged_info_initial <- print(scNanoGPS_4_merged)
scNanoGPS_5_merged_info_initial <- print(scNanoGPS_5_merged)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# # # LABEL EACH DATASET (BY "TYPE" OR IN MY CASE, SAMPLE ID #)
# # ipsc1$type <- "ipsc"
# scNanoGPS_1$type <- 'SH-04-08'
# scNanoGPS_2$type <- 'SH-06-25'
# scNanoGPS_3$type <- 'SH-07-46'
# scNanoGPS_4$type <- 'SH-92-05'
# scNanoGPS_5$type <- 'UMARY_4546'

# scNanoGPS_1_isoform$type <- 'SH-04-08_isoform'
# scNanoGPS_2_isoform$type <- 'SH-06-25_isoform'
# scNanoGPS_3_isoform$type <- 'SH-07-46_isoform'
# scNanoGPS_4_isoform$type <- 'SH-92-05_isoform'
# scNanoGPS_5_isoform$type <- 'UMARY_4546_isoform'

scNanoGPS_1_merged$type <- 'SH-04-08_merged'
scNanoGPS_2_merged$type <- 'SH-06-25_merged'
scNanoGPS_3_merged$type <- 'SH-07-46_merged'
scNanoGPS_4_merged$type <- 'SH-92-05_merged'
scNanoGPS_5_merged$type <- 'UMARY_4546_merged'

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# # # MERGE BY TIMEPOINT
# # ipsc <- merge(ipsc1, y=c(ipsc2, ipsc3, ipsc4), add.cell.ids = c("ipsc1", "ipsc2", "ipsc3", "ipsc4"))
# # nd4 <- merge(nd4.1, y=c(nd4.2, nd4.3, nd4.4), add.cell.ids = c("nd4.1", "nd4.2", "nd4.3", "nd4.4"))
# scNanoGPS <- merge(
# 	scNanoGPS_1,
# 	y = c(scNanoGPS_2, scNanoGPS_3, scNanoGPS_4, scNanoGPS_5, scNanoGPS_1_isoform, scNanoGPS_2_isoform, scNanoGPS_3_isoform, scNanoGPS_4_isoform, scNanoGPS_5_isoform),
# 	add.cell.ids = c('SH-04-08','SH-06-25','SH-07-46','SH-92-05','UMARY_4546', 'SH-04-08_isoform', 'SH-06-25_isoform', 'SH-07-46_isoform', 'SH-92-05_isoform', 'UMARY_4546_isoform') 
# )

#? Q) Should I merge the merged samples here? yeah probably, but take out the individual samples here and merge into (lol) a single merged workflow now

scNanoGPS <- merged(
	scNanoGPS_1_merge,
	y = c( 
		scNanoGPS_2_merged, 
		scNanoGPS_3_merged, 
		scNanoGPS_4_merged, 
		scNanoGPS_5_merged
		),
	add.cell.ids = c(
		'SH-04-08_merged', 
		'SH-06-25_merged', 
		'SH-07-46_merged', 
		'SH-92-05_merged', 
		'UMARY_4546_merged'
		) 
	)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# GET INFO
cat('Initial MERGED Seurat object created stats: \n') 
scNanoGPS_info_initial <- print(scNanoGPS)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# # QC - BY TIMEPOINT
# ipsc[["percent.mt"]] <- PercentageFeatureSet(ipsc, pattern = "^MT-")
# nd4[["percent.mt"]] <- PercentageFeatureSet(nd4, pattern = "^MT-")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# #* TODO - could name with the ID's here if you had IDs or labels saved

# # PLOT - VIOLIN PLOTS
# p1 <- VlnPlot(ipsc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# ggsave("June2024/ipscQC.png", plot = p1, width = 15, height = 5, units = "in", dpi = 300)
VlnPlot(scNanoGPS, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# SAVE PLOT
save_plot_violin_counts_features <- my_plot_name('violin_counts_features')
my_plot_save(save_plot_violin_counts_features)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#? Q) Should I filter each dataset separately? Then merge again?
#? A) Yeah -- can just use 1 filter for 1500

# # # FILTER CELLS - (MANUALLY) BASED ON VIOLIN PLOTS
# # ipsc <- subset(ipsc, subset = nFeature_RNA > 1000 & nFeature_RNA < 10000 & percent.mt < 10)
# # nd4 <- subset(nd4, subset = nFeature_RNA > 1000 & nFeature_RNA < 10000 & percent.mt < 10)
# scNanoGPS_1 <- subset(scNanoGPS_1, subset = nFeature_RNA > 0 & nFeature_RNA < 1500)
# scNanoGPS_2 <- subset(scNanoGPS_2, subset = nFeature_RNA > 0 & nFeature_RNA < 500)
# scNanoGPS_3 <- subset(scNanoGPS_3, subset = nFeature_RNA > 0 & nFeature_RNA < 1000)
# scNanoGPS_4 <- subset(scNanoGPS_4, subset = nFeature_RNA > 0 & nFeature_RNA < 1500)
# scNanoGPS_5 <- subset(scNanoGPS_5, subset = nFeature_RNA > 0 & nFeature_RNA < 500)
scNanoGPS <- subset(scNanoGPS, subset = nFeature_RNA > 0 & nFeature_RNA < 1500)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# # # MERGE TIMEPOINTS
# # ngn2 <- merge(ipsc, y=c(nd4, nd7, nd10, nd14, nd21, nd28, k4), add.cell.ids = c("ipsc", "nd4", "nd7", "nd10", "nd14", "nd21", "nd28", "k4"))
# ngn2 <- NormalizeData(ngn2)
# ngn2 <- FindVariableFeatures(ngn2)
# ngn2 <- ScaleData(ngn2)
# ngn2 <- RunPCA(ngn2)

# # MERGE FILTERED SUBSETS + NORMALIZE + FIND VARIABLE FEATURES + SCALE + PCA
# # NOTE - if not filtering by each sample's subset, then can disable the 2nd merge command
# scNanoGPS <- merge(
# 	scNanoGPS_1,
# 	y = c(scNanoGPS_2, scNanoGPS_3, scNanoGPS_4, scNanoGPS_5),
# 	add.cell.ids = c('SH-04-08','SH-06-25','SH-07-46','SH-92-05','UMARY_4546') )
scNanoGPS <- NormalizeData(scNanoGPS)

scNanoGPS <- FindVariableFeatures(scNanoGPS, nfeatures = 10000)
# #! ERROR
# Error in `.SelectFeatures()`:
#   ! None of the features provided are present in the feature set
#* TODO - manually merge each sample's matrix + matrix_isoform files on the barcodes 

scNanoGPS <- ScaleData(scNanoGPS)
scNanoGPS <- RunPCA(scNanoGPS)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# GET INFO
cat('Final MERGED Seurat object (post filtering) stats: \n')
scNanoGPS_info_final <- print(scNanoGPS)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# PLOT - ELBOW
ElbowPlot(scNanoGPS)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# SAVE PLOT
save_plot_elbow <- my_plot_name('elbow')
my_plot_save(save_plot_elbow)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# #? Q) Before + after?

# #! ERROR - might not be able to run both unintegrated and integrated? or need to specify the layer after?
# Error in `IntegrateLayers()`:
#   ! None of the features provided are found in this assay

# NOTE - integrated and unintegrated plots look very similar, that's fine, batch correction for similar samples should be like that

# # UNINTEGRATED SAMPLES 
# scNanoGPS <- FindNeighbors(scNanoGPS, dims = 1:10, reduction = "pca")
# scNanoGPS <- FindClusters(scNanoGPS, resolutions = 0.2, cluster.name = "unintegrated_clusters")
# scNanoGPS <- RunUMAP(scNanoGPS, dims = 1:10, reduction = "pca", reduction.name = "umap.unintegrated")
# 
# # PLOT - UNINTEGRATED SAMPLES
# # Octopus
# p1 <- DimPlot(scNanoGPS, reduction = "umap.unintegrated", group.by = c("type", "seurat_clusters"))
# plot(p1)
# #ggsave("June2024/unintegrated_ngn2.png", plot = p1, width = 15, height = 5, units = "in", dpi = 300)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# INTEGRATE USING HARMONY
scNanoGPS <- IntegrateLayers(object = scNanoGPS, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony", verbose = FALSE)
scNanoGPS <- FindNeighbors(scNanoGPS, reduction = "harmony", dims = 1:10)
scNanoGPS <- FindClusters(scNanoGPS, resolution = 0.2)
scNanoGPS <- RunUMAP(scNanoGPS, dims = 1:10, reduction = "harmony")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# # PLOT
# p1 <- DimPlot(ngn2, reduction="umap", group.by=c("type", "seurat_clusters"), label=TRUE)
# ggsave("June2024/harmony_ngn2all.png", plot = p1, width = 15, height = 5, units = "in", dpi = 300)
DimPlot(scNanoGPS, reduction = "umap", group.by = c("type", "seurat_clusters"), label = TRUE)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# SAVE PLOT
save_plot_dim_groupby_clusters <- my_plot_name('dim_groupby_clusters')
my_plot_save(save_plot_dim_groupby_clusters)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# p1 <- DimPlot(ngn2, reduction="umap", split.by="type", label = TRUE)
# ggsave("June2024/TypeH_ngn2all.png", plot = p1, width = 15, height = 5, units = "in", dpi = 300)
DimPlot(scNanoGPS, reduction="umap", split.by="type", label = TRUE)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# SAVE PLOT
save_plot_dim_splitby_type <- my_plot_name('dim_splitby_type')
my_plot_save(save_plot_dim_splitby_type)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# # TEST - SHAHROZE - REQUIRED!!!
# ?JoinLayers
scNanoGPS <- JoinLayers(scNanoGPS)

# XYLENA - SKIP!
cat('\nMarkers \n')
# FindAllMarkers is looking for DEG
scNanoGPS.markers <- FindAllMarkers(object = scNanoGPS)
print(scNanoGPS.markers)
 
# #* WARNING 
# Calculating cluster 0
# Calculating cluster 1
# Calculating cluster 2
# Calculating cluster 3
# Calculating cluster 4
# Calculating cluster 5
# Calculating cluster 6
# Calculating cluster 7
# Calculating cluster 8
# Warning: No DE genes identified
# Warning: The following tests were not performed: 
#   Warning: When testing 0 versus all:
#   data layers are not joined. Please run JoinLayers
# Warning: When testing 1 versus all:
#   data layers are not joined. Please run JoinLayers
# Warning: When testing 2 versus all:
#   data layers are not joined. Please run JoinLayers
# Warning: When testing 3 versus all:
#   data layers are not joined. Please run JoinLayers
# Warning: When testing 4 versus all:
#   data layers are not joined. Please run JoinLayers
# Warning: When testing 5 versus all:
#   data layers are not joined. Please run JoinLayers
# Warning: When testing 6 versus all:
#   data layers are not joined. Please run JoinLayers
# Warning: When testing 7 versus all:
#   data layers are not joined. Please run JoinLayers
# Warning: When testing 8 versus all:
#   data layers are not joined. Please run JoinLayers

# # EXPORT - to .CSV file
# ngn2.markers <- FindAllMarkers(object = ngn2)
# write.csv(ngn2.markers, file = "June2024/ngn2_markers_all.csv")
marker_table_path <- file.path(output_path, 'scNanoGPS_markers_all.csv')
write.csv(scNanoGPS.markers, file = marker_table_path)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# # SAVE - .RDS FILE
# saveRDS(ngn2, file = "June2024/ngn2_062024.rds")

# SAVE POINT
# OUTPUT - .RDS PATH
rds_path <- file.path(output_path, 'scNanoGPS.rds')
cat('rds_path = \n', rds_path, '\n')

# EXPORT - .RDS file
saveRDS(scNanoGPS, file = rds_path)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# # TEST
# # Identify the 10 most highly variable genes
# top10 <- head(VariableFeatures(scNanoGPS), 10)
# cat('top 10 = \n')
# paste(top10)

# # TEST 
# scNanoGPS.markers %>%
#   group_by(.data) %>%
#   dplyr::filter(avg_log2FC > 1)
# 
# scNanoGPS.markers %>%
#   group_by(.data) %>%
#   dplyr::filter(avg_log2FC > 1) %>%
#   slice_head(n = 10) %>%
#   ungroup() -> top10
#
# #! ERROR
# Error in `group_by()`:
#   ! Must group by variables found in `.data`.
# ✖ Column `.data` is not found.

# #* TODO - change to top 10 for isoforms
# # PLOT - expression of all markers
# marker.genes <- c("RBFOX3", "TUBB3", "MAPT", "MAP2", "SYN1", "DLG4", "CUX1", "CUX2", "BCL11B", "FOXP2", "TBR1", "POU3F3", "POU3F2",
# 	"SLC17A7", "SLC17A6", "GRIN2B", "GRIN1", "GRIN2B", "RORB", "THEMIS", "GAD1", "SST", "VIP", "PVALB", "ISL1", "SLC5A7", "DCN", "PRPH",
# 	"POU4F1", "LHX9", "GPM6A", "CD99", "SOX2", "POU5F1", "NES", "NANOG")

marker.genes <- c(
	'ORB',
	'LC17A6',
	'LC17A7',
	'HEMIS',
	'AD1',
	'AD2',
	'VALB',
	'ST',
	'IP',
	'PBB1IP',
	'D74',
	'SF1R',
	'X3CR1',
	'TGAM',
	'2RY12',
	'TPRC',
	'LDH1L1',
	'QP4',
	'OL5A3',
	'FAP',
	'LC1A2',
	'LC1A3',
	'LDN11',
	'BP',
	'OBP',
	'PALIN',
	'LP1',
	'T18',
	'HFPL3',
	'EGF11',
	'CDH15',
	'DGFRA',
	'CAN',
	'LDN5',
	'OLEC12',
	'PAS1',
	'CAM1'
	)

# p1 <- DotPlot(object=ngn2, features = marker.genes) + RotatedAxis()
# ggsave("June2024/ngn2Dotplotall.png", plot = p1, width = 12, height = 5, units = "in", dpi = 300)	
DotPlot(object = scNanoGPS, features = marker.genes) + RotatedAxis()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# SAVE PLOT
save_plot_dotplot_marker_genes <- my_plot_name('dotplot_marker_genes')
my_plot_save(save_plot_dotplot_marker_genes)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# # PLOT - FEATURE PLOTS - FOR CELL TYPE SPECIFIC GENE EXPRESSION
# p1 <- FeaturePlot(ngn2, features = c("RBFOX3", "TUBB3", "MAPT", "MAP2"))
# ggsave("June2024/all_neurons.png", plot = p1, width = 10, height = 5, units = "in", dpi = 300)	
# p1 <- FeaturePlot(ngn2, features = c("CUX1", "CUX2", "POU3F2", "BCL11B"))
# ggsave("June2024/all_CorticL.png", plot = p1, width = 10, height = 5, units = "in", dpi = 300)	

# Put our celltype markers list
marker_genes 		<- list()
marker_genes$Astro 	<- c('ALDH1L1', 'AQP4', 'COL5A3', 'GFAP', 'SLC1A2', 'SLC1A3')
marker_genes$ExN 	<- c('RORB', 'SLC17A6', 'SLC17A7', 'THEMIS')
marker_genes$InN 	<- c('GAD1', 'GAD2', 'PVALB', 'SST', 'VIP')
marker_genes$MG 	<- c('APBB1IP', 'CD74', 'CSF1R', 'CX3CR1', 'ITGAM', 'P2RY12', 'PTPRC')
marker_genes$Oligo 	<- c('CLDN11', 'MBP', 'MOBP', 'OPALIN', 'PLP1', 'ST18')
marker_genes$OPC 	<- c('LHFPL3', 'MEGF11', 'PCDH15', 'PDGFRA', 'VCAN')
marker_genes$VC 	<- c('CLDN5', 'COLEC12', 'EPAS1', 'VCAM1') 

# PLOT - feature plot with each of the cell types (separately)
FeaturePlot(scNanoGPS, features = marker_genes$Astro)
# SAVE PLOT
save_plot_feature_Astro <- my_plot_name('feature-Astro')
my_plot_save(save_plot_feature_Astro)

# FeaturePlot(scNanoGPS, features = marker_genes$Exn)
# # SAVE PLOT
# save_plot_feature_ExN <- my_plot_name('feature-ExN')
# my_plot_save(save_plot_feature_ExN)

# # #! ERROR
# # Error in `FeaturePlot()`:
# #   ! None of the requested features were found:  in slot  data

FeaturePlot(scNanoGPS, features = marker_genes$InN)
# SAVE PLOT
save_plot_feature_InN <- my_plot_name('feature-InN')
my_plot_save(save_plot_feature_InN)

FeaturePlot(scNanoGPS, features = marker_genes$MG)
# SAVE PLOT
save_plot_feature_MG <- my_plot_name('feature-MG')
my_plot_save(save_plot_feature_MG)

FeaturePlot(scNanoGPS, features = marker_genes$Oligo)
# SAVE PLOT
save_plot_feature_Oligo <- my_plot_name('feature-Oligo')
my_plot_save(save_plot_feature_Oligo)

FeaturePlot(scNanoGPS, features = marker_genes$OPC)
# SAVE PLOT
save_plot_feature_OPC <- my_plot_name('feature-OPC')
my_plot_save(save_plot_feature_OPC)

FeaturePlot(scNanoGPS, features = marker_genes$VC)
# SAVE PLOT
save_plot_feature_VC <- my_plot_name('feature-VC')
my_plot_save(save_plot_feature_VC)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# # PLOT - HEATMAP
# p1 <- DoHeatmap(object=ngn2, features = marker.genes, size = 3, draw.lines = F)
# ggsave("June2024/ngn2heatmapall.png", plot = p1, width = 12, height = 5, units = "in", dpi = 300)
DoHeatmap(object = scNanoGPS, features = marker.genes, size = 3, draw.lines = T)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# SAVE PLOT
save_plot_heatmap_marker_genes <- my_plot_name('heatmap-marker_genes')
my_plot_save(save_plot_heatmap_marker_genes)

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

# VARIABLES CALCULATED IN SCRIPT
cat('\nVARIABLES CALCULATED IN SCRIPT: \n')

cat('\nInput directories: \n')
cat('input_dir_1 = \n', input_dir_1, '\n') 
cat('input_dir_2 = \n', input_dir_2, '\n') 
cat('input_dir_3 = \n', input_dir_3, '\n') 
cat('input_dir_4 = \n', input_dir_4, '\n') 
cat('input_dir_5 = \n', input_dir_5, '\n') 

cat('input_dir_1_isoform  = \n', input_dir_1_isoform , '\n') 
cat('input_dir_2_isoform  = \n', input_dir_2_isoform , '\n') 
cat('input_dir_3_isoform  = \n', input_dir_3_isoform , '\n') 
cat('input_dir_4_isoform  = \n', input_dir_4_isoform , '\n') 
cat('input_dir_5_isoform  = \n', input_dir_5_isoform , '\n') 

cat('\nOutput directories: \n')
cat('output_path = \n', output_path, '\n')

cat('\nInitial Seurat object created stats: \n') 
scNanoGPS_1_info_initial
scNanoGPS_2_info_initial
scNanoGPS_3_info_initial
scNanoGPS_4_info_initial
scNanoGPS_5_info_initial

scNanoGPS_1_isoform_info_initial
scNanoGPS_2_isoform_info_initial
scNanoGPS_3_isoform_info_initial
scNanoGPS_4_isoform_info_initial
scNanoGPS_5_isoform_info_initial

cat('\nInitial MERGED Seurat object created stats: \n') 
scNanoGPS_info_initial

cat('\nFinal MERGED Seurat object (post filtering) stats: \n')
scNanoGPS_info_final <- print(scNanoGPS)

cat('\nMarkers \n')
print(scNanoGPS.markers)

cat('\nmarker_genes lists: \n')
print(marker_genes)

# Output files
cat('marker_table_path = \n', marker_table_path)
cat('rds_path = \n', rds_path, '\n')

cat('\nsessionInfo: \n')
sessionInfo()

# CLOSE SINK
sink(NULL)