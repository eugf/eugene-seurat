#!/usr/bin/env Rscript
# start like this
# Rscript --vanilla Integrate_ngn2_0624_res2.R

#load libraries
library(Seurat, lib.loc='/data/reedx/R/rhel8/4.3/')
library(dplyr)
library(patchwork)
library(ggplot2)
library(fields)
library(KernSmooth)
library(Matrix)
library(parallelly)
library(ROCR)
library(parallel, lib.loc = "/Library/Frameworks/R.framework/Versions/4.2/Resources/library")
library(tidyverse)
library(cowplot)
library(SeuratObject)
library(future)
library(future.apply)
library(harmony, lib.loc='/data/reedx/R/rhel8/4.3/')

#load data
ipsc1 <- Read10X_h5("/data/CARD_singlecell/NGN2_diff_Erika/data/sample_iPSC_1/outs/filtered_feature_bc_matrix.h5")
ipsc2 <- Read10X_h5("/data/CARD_singlecell/NGN2_diff_Erika/data/sample_iPSC_2/outs/filtered_feature_bc_matrix.h5")
ipsc3 <- Read10X_h5("/data/CARD_singlecell/NGN2_diff_Erika/data/sample_iPSC_3/outs/filtered_feature_bc_matrix.h5")
ipsc4 <- Read10X_h5("/data/CARD_singlecell/NGN2_diff_Erika/data/sample_iPSC_4/outs/filtered_feature_bc_matrix.h5")
d4.1 <- Read10X_h5("/data/CARD_singlecell/NGN2_diff_Erika/data/sample_Nd4_1/outs/filtered_feature_bc_matrix.h5")
d4.2 <- Read10X_h5("/data/CARD_singlecell/NGN2_diff_Erika/data/sample_Nd4_2/outs/filtered_feature_bc_matrix.h5")
d4.3 <- Read10X_h5("/data/CARD_singlecell/NGN2_diff_Erika/data/sample_Nd4_3/outs/filtered_feature_bc_matrix.h5")
d4.4 <- Read10X_h5("/data/CARD_singlecell/NGN2_diff_Erika/data/sample_Nd4_4/outs/filtered_feature_bc_matrix.h5")
d7.1 <- Read10X_h5("/data/CARD_singlecell/NGN2_diff_Erika/data/sample_Nd7_1/outs/filtered_feature_bc_matrix.h5")
d7.2 <- Read10X_h5("/data/CARD_singlecell/NGN2_diff_Erika/data/sample_Nd7_2/outs/filtered_feature_bc_matrix.h5")
d7.3 <- Read10X_h5("/data/CARD_singlecell/NGN2_diff_Erika/data/sample_Nd7_3/outs/filtered_feature_bc_matrix.h5")
d7.4 <- Read10X_h5("/data/CARD_singlecell/NGN2_diff_Erika/data/sample_Nd7_4/outs/filtered_feature_bc_matrix.h5")
d10.1 <- Read10X_h5("/data/CARD_singlecell/NGN2_diff_Erika/data/sample_Nd10_1/outs/filtered_feature_bc_matrix.h5")
d10.2 <- Read10X_h5("/data/CARD_singlecell/NGN2_diff_Erika/data/sample_Nd10_2/outs/filtered_feature_bc_matrix.h5")
d10.3 <- Read10X_h5("/data/CARD_singlecell/NGN2_diff_Erika/data/sample_Nd10_3/outs/filtered_feature_bc_matrix.h5")
d10.4 <- Read10X_h5("/data/CARD_singlecell/NGN2_diff_Erika/data/sample_Nd10_4/outs/filtered_feature_bc_matrix.h5")
d14.1 <- Read10X_h5("/data/CARD_singlecell/NGN2_diff_Erika/data/sample_Nd14_1/outs/filtered_feature_bc_matrix.h5")
d14.2 <- Read10X_h5("/data/CARD_singlecell/NGN2_diff_Erika/data/sample_Nd14_2/outs/filtered_feature_bc_matrix.h5")
d14.3 <- Read10X_h5("/data/CARD_singlecell/NGN2_diff_Erika/data/sample_Nd14_3/outs/filtered_feature_bc_matrix.h5")
d14.4 <- Read10X_h5("/data/CARD_singlecell/NGN2_diff_Erika/data/sample_Nd14_4/outs/filtered_feature_bc_matrix.h5")
d21.1 <- Read10X_h5("/data/CARD_singlecell/NGN2_diff_Erika/data/sample_Nd21_1/outs/filtered_feature_bc_matrix.h5")
d21.2 <- Read10X_h5("/data/CARD_singlecell/NGN2_diff_Erika/data/sample_Nd21_2/outs/filtered_feature_bc_matrix.h5")
d21.3 <- Read10X_h5("/data/CARD_singlecell/NGN2_diff_Erika/data/sample_Nd21_3/outs/filtered_feature_bc_matrix.h5")
d21.4 <- Read10X_h5("/data/CARD_singlecell/NGN2_diff_Erika/data/sample_Nd21_4/outs/filtered_feature_bc_matrix.h5")
d28.1 <- Read10X_h5("/data/CARD_singlecell/NGN2_diff_Erika/data/Nd28-1n/outs/filtered_feature_bc_matrix.h5")
d28.2 <- Read10X_h5("/data/CARD_singlecell/NGN2_diff_Erika/data/Nd28-2n/outs/filtered_feature_bc_matrix.h5")
d28.3 <- Read10X_h5("/data/CARD_singlecell/NGN2_diff_Erika/data/Nd28-3n/outs/filtered_feature_bc_matrix.h5")
d28.4 <- Read10X_h5("/data/CARD_singlecell/NGN2_diff_Erika/data/Nd28-4n/outs/filtered_feature_bc_matrix.h5")
k4.1 <- Read10X_h5("/data/CARD_singlecell/NGN2_diff_Erika/data/K4-1n/outs/filtered_feature_bc_matrix.h5")
k4.2 <- Read10X_h5("/data/CARD_singlecell/NGN2_diff_Erika/data/K4-2n/outs/filtered_feature_bc_matrix.h5")
k4.3 <- Read10X_h5("/data/CARD_singlecell/NGN2_diff_Erika/data/K4-3n/outs/filtered_feature_bc_matrix.h5")
k4.4 <- Read10X_h5("/data/CARD_singlecell/NGN2_diff_Erika/data/K4-4n/outs/filtered_feature_bc_matrix.h5")

#Create Seurat Objects
ipsc1 <- CreateSeuratObject(counts = ipsc1, project = "ipsc1", min.cells = 3, min.features = 200)
ipsc2 <- CreateSeuratObject(counts = ipsc2, project = "ipsc2", min.cells = 3, min.features = 200)
ipsc3 <- CreateSeuratObject(counts = ipsc3, project = "ipsc3", min.cells = 3, min.features = 200)
ipsc4 <- CreateSeuratObject(counts = ipsc4, project = "ipsc4", min.cells = 3, min.features = 200)
nd4.1 <- CreateSeuratObject(counts = d4.1, project = "nd4.1", min.cells = 3, min.features = 200)
nd4.2 <- CreateSeuratObject(counts = d4.2, project = "nd4.2", min.cells = 3, min.features = 200)
nd4.3 <- CreateSeuratObject(counts = d4.3, project = "nd4.3", min.cells = 3, min.features = 200)
nd4.4 <- CreateSeuratObject(counts = d4.4, project = "nd4.4", min.cells = 3, min.features = 200)
nd7.1 <- CreateSeuratObject(counts = d7.1, project = "nd7.1", min.cells = 3, min.features = 200)
nd7.2 <- CreateSeuratObject(counts = d7.2, project = "nd7.2", min.cells = 3, min.features = 200)
nd7.3 <- CreateSeuratObject(counts = d7.3, project = "nd7.3", min.cells = 3, min.features = 200)
nd7.4 <- CreateSeuratObject(counts = d7.4, project = "nd7.4", min.cells = 3, min.features = 200)
nd10.1 <- CreateSeuratObject(counts = d10.1, project = "nd10.1", min.cells = 3, min.features = 200)
nd10.2 <- CreateSeuratObject(counts = d10.2, project = "nd10.2", min.cells = 3, min.features = 200)
nd10.3 <- CreateSeuratObject(counts = d10.3, project = "nd10.3", min.cells = 3, min.features = 200)
nd10.4 <- CreateSeuratObject(counts = d10.4, project = "nd10.4", min.cells = 3, min.features = 200)
nd14.1 <- CreateSeuratObject(counts = d14.1, project = "nd14.1", min.cells = 3, min.features = 200)
nd14.2 <- CreateSeuratObject(counts = d14.2, project = "nd14.2", min.cells = 3, min.features = 200)
nd14.3 <- CreateSeuratObject(counts = d14.3, project = "nd14.3", min.cells = 3, min.features = 200)
nd14.4 <- CreateSeuratObject(counts = d14.4, project = "nd14.4", min.cells = 3, min.features = 200)
nd21.1 <- CreateSeuratObject(counts = d21.1, project = "nd21.1", min.cells = 3, min.features = 200)
nd21.2 <- CreateSeuratObject(counts = d21.2, project = "nd21.2", min.cells = 3, min.features = 200)
nd21.3 <- CreateSeuratObject(counts = d21.3, project = "nd21.3", min.cells = 3, min.features = 200)
nd21.4 <- CreateSeuratObject(counts = d21.4, project = "nd21.4", min.cells = 3, min.features = 200)
nd28.1 <- CreateSeuratObject(counts = d28.1, project = "nd28.1", min.cells = 3, min.features = 200)
nd28.2 <- CreateSeuratObject(counts = d28.2, project = "nd28.2", min.cells = 3, min.features = 200)
nd28.3 <- CreateSeuratObject(counts = d28.3, project = "nd28.3", min.cells = 3, min.features = 200)
nd28.4 <- CreateSeuratObject(counts = d28.4, project = "nd28.4", min.cells = 3, min.features = 200)
k4.1 <- CreateSeuratObject(counts = k4.1, project = "k4.1", min.cells = 3, min.features = 200)
k4.2 <- CreateSeuratObject(counts = k4.2, project = "k4.2", min.cells = 3, min.features = 200)
k4.3 <- CreateSeuratObject(counts = k4.3, project = "k4.3", min.cells = 3, min.features = 200)
k4.4 <- CreateSeuratObject(counts = k4.4, project = "k4.4", min.cells = 3, min.features = 200)


#label each dataset
ipsc1$type <- "ipsc"
ipsc2$type <- "ipsc"
ipsc3$type <- "ipsc"
ipsc4$type <- "ipsc"
nd4.1$type <- "nd4"
nd4.2$type <- "nd4"
nd4.3$type <- "nd4"
nd4.4$type <- "nd4"
nd7.1$type <- "nd7"
nd7.2$type <- "nd7"
nd7.3$type <- "nd7"
nd7.4$type <- "nd7"
nd10.1$type <- "nd10"
nd10.2$type <- "nd10"
nd10.3$type <- "nd10"
nd10.4$type <- "nd10"
nd14.1$type <- "nd14"
nd14.2$type <- "nd14"
nd14.3$type <- "nd14"
nd14.4$type <- "nd14"
nd21.1$type <- "nd21"
nd21.2$type <- "nd21"
nd21.3$type <- "nd21"
nd21.4$type <- "nd21"
nd28.1$type <- "nd28"
nd28.2$type <- "nd28"
nd28.3$type <- "nd28"
nd28.4$type <- "nd28"
k4.1$type <- "k4"
k4.2$type <- "k4"
k4.3$type <- "k4"
k4.4$type <- "k4"


#merge by timepoint
ipsc <- merge(ipsc1, y=c(ipsc2, ipsc3, ipsc4), add.cell.ids = c("ipsc1", "ipsc2", "ipsc3", "ipsc4"))
nd4 <- merge(nd4.1, y=c(nd4.2, nd4.3, nd4.4), add.cell.ids = c("nd4.1", "nd4.2", "nd4.3", "nd4.4"))
nd7 <- merge(nd7.1, y=c(nd7.2, nd7.3, nd7.4), add.cell.ids = c("nd7.1", "nd7.2", "nd7.3", "nd7.4"))
nd10 <- merge(nd10.1, y=c(nd10.2, nd10.3, nd10.4), add.cell.ids = c("nd10.1", "nd10.2", "nd10.3", "nd10.4"))
nd14 <- merge(nd14.1, y=c(nd14.2, nd14.3, nd14.4), add.cell.ids = c("nd14.1", "nd14.2", "nd14.3", "nd14.4"))
nd21 <- merge(nd21.1, y=c(nd21.2, nd21.3, nd21.4), add.cell.ids = c("nd21.1", "nd21.2", "nd21.3", "nd21.4"))
nd28 <- merge(nd28.1, y=c(nd28.2, nd28.3, nd28.4), add.cell.ids = c("nd28.1", "nd28.2", "nd28.3", "nd28.4"))
k4 <- merge(k4.1, y=c(k4.2, k4.3, k4.4), add.cell.ids = c("k4.1", "k4.2", "k4.3", "k4.4"))

#QC by timepoint
ipsc[["percent.mt"]] <- PercentageFeatureSet(ipsc, pattern = "^MT-")
nd4[["percent.mt"]] <- PercentageFeatureSet(nd4, pattern = "^MT-")
nd7[["percent.mt"]] <- PercentageFeatureSet(nd7, pattern = "^MT-")
nd10[["percent.mt"]] <- PercentageFeatureSet(nd10, pattern = "^MT-")
nd14[["percent.mt"]] <- PercentageFeatureSet(nd14, pattern = "^MT-")
nd21[["percent.mt"]] <- PercentageFeatureSet(nd21, pattern = "^MT-")
nd28[["percent.mt"]] <- PercentageFeatureSet(nd28, pattern = "^MT-")
k4[["percent.mt"]] <- PercentageFeatureSet(k4, pattern = "^MT-")

p1 <- VlnPlot(ipsc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave("June2024/ipscQC.png", plot = p1, width = 15, height = 5, units = "in", dpi = 300)
p1 <- VlnPlot(nd4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave("June2024/nd4QC.png", plot = p1, width = 15, height = 5, units = "in", dpi = 300)
p1 <- VlnPlot(nd7, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave("June2024/nd7QC.png", plot = p1, width = 15, height = 5, units = "in", dpi = 300)
p1 <- VlnPlot(nd10, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave("June2024/nd10QC.png", plot = p1, width = 15, height = 5, units = "in", dpi = 300)
p1 <- VlnPlot(nd14, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave("June2024/nd14QC.png", plot = p1, width = 15, height = 5, units = "in", dpi = 300)
p1 <- VlnPlot(nd21, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave("June2024/nd21QC.png", plot = p1, width = 15, height = 5, units = "in", dpi = 300)
p1 <- VlnPlot(nd28, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave("June2024/nd28QC.png", plot = p1, width = 15, height = 5, units = "in", dpi = 300)
p1 <- VlnPlot(k4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave("June2024/k4QC.png", plot = p1, width = 15, height = 5, units = "in", dpi = 300)


#filter cells
ipsc <- subset(ipsc, subset = nFeature_RNA > 1000 & nFeature_RNA < 10000 & percent.mt < 10)
nd4 <- subset(nd4, subset = nFeature_RNA > 1000 & nFeature_RNA < 10000 & percent.mt < 10)
nd7 <- subset(nd7, subset = nFeature_RNA > 1000 & nFeature_RNA < 10000 & percent.mt < 10)
nd10 <- subset(nd10, subset = nFeature_RNA > 1000 & nFeature_RNA < 10000 & percent.mt < 10)
nd14 <- subset(nd14, subset = nFeature_RNA > 1000 & nFeature_RNA < 10000 & percent.mt < 10)
nd21 <- subset(nd21, subset = nFeature_RNA > 1000 & nFeature_RNA < 10000 & percent.mt < 10)
nd28 <- subset(nd28, subset = nFeature_RNA > 1000 & nFeature_RNA < 10000 & percent.mt < 10)
k4 <- subset(k4, subset = nFeature_RNA > 1000 & nFeature_RNA < 10000 & percent.mt < 10)

#merge all timepoints
ngn2 <- merge(ipsc, y=c(nd4, nd7, nd10, nd14, nd21, nd28, k4), add.cell.ids = c("ipsc", "nd4", "nd7", "nd10", "nd14", "nd21", "nd28", "k4"))
ngn2 <- NormalizeData(ngn2)
ngn2 <- FindVariableFeatures(ngn2)
ngn2 <- ScaleData(ngn2)
ngn2 <- RunPCA(ngn2)

#Plot unintegrated samples
#ngn2 <- FindNeighbors(ngn2, dims=1:30, reduction="pca")
#ngn2 <- FindClusters(ngn2, resolutions =0.3, cluster.name = "unintegrated_clusters")
#ngn2 <- RunUMAP(ngn2, dims=1:30, reduction ="pca", reduction.name = "umap.unintegrated")
#p1 <- DimPlot(ngn2, reduction="umap.unintegrated", group.by=c("type", "seurat_clusters"))
#ggsave("June2024/unintegrated_ngn2.png", plot = p1, width = 15, height = 5, units = "in", dpi = 300)

#Integrate using harmony
ngn2 <- IntegrateLayers(object=ngn2, method = HarmonyIntegration, orig.reduction= "pca", new.reduction = "harmony", verbose= FALSE)
ngn2 <- FindNeighbors(ngn2, reduction = "harmony", dims = 1:30)
ngn2 <- FindClusters(ngn2, resolution = 0.2)
ngn2 <- RunUMAP(ngn2, dims=1:30, reduction = "harmony")
p1 <- DimPlot(ngn2, reduction="umap", group.by=c("type", "seurat_clusters"), label=TRUE)
ggsave("June2024/harmony_ngn2all.png", plot = p1, width = 15, height = 5, units = "in", dpi = 300)

p1 <- DimPlot(ngn2, reduction="umap", split.by="type", label = TRUE)
ggsave("June2024/TypeH_ngn2all.png", plot = p1, width = 15, height = 5, units = "in", dpi = 300)

saveRDS(ngn2, file = "June2024/ngn2_062024.rds")

p1 <- FeaturePlot(ngn2, features = c("NANOG", "DCN", "MAP2", "PRPH"))
ggsave("June2024/ngn_all1.png", plot = p1, width = 10, height = 5, units = "in", dpi = 300)	
p1 <- FeaturePlot(ngn2, features = c("POU4F1", "LHX9", "GPM6A", "CD99"))
ggsave("June2024/ngn_all2.png", plot = p1, width = 10, height = 5, units = "in", dpi = 300)

ngn2.markers <- FindAllMarkers(object = ngn2)
write.csv(ngn2.markers, file = "June2024/ngn2_markers_all.csv")

#Plot expression of all markers
marker.genes <- c("RBFOX3", "TUBB3", "MAPT", "MAP2", "SYN1", "DLG4", "CUX1", "CUX2", "BCL11B", "FOXP2", "TBR1", "POU3F3", "POU3F2",
	"SLC17A7", "SLC17A6", "GRIN2B", "GRIN1", "GRIN2B", "RORB", "THEMIS", "GAD1", "SST", "VIP", "PVALB", "ISL1", "SLC5A7", "DCN", "PRPH",
	"POU4F1", "LHX9", "GPM6A", "CD99", "SOX2", "POU5F1", "NES", "NANOG")

p1 <- DotPlot(object=ngn2, features = marker.genes) + RotatedAxis()
ggsave("June2024/ngn2Dotplotall.png", plot = p1, width = 12, height = 5, units = "in", dpi = 300)	

#FeaturePlots for cell type specific gene expression
p1 <- FeaturePlot(ngn2, features = c("RBFOX3", "TUBB3", "MAPT", "MAP2"))
ggsave("June2024/all_neurons.png", plot = p1, width = 10, height = 5, units = "in", dpi = 300)	
p1 <- FeaturePlot(ngn2, features = c("CUX1", "CUX2", "POU3F2", "BCL11B"))
ggsave("June2024/all_CorticL.png", plot = p1, width = 10, height = 5, units = "in", dpi = 300)	
p1 <- FeaturePlot(ngn2, features = c("FOXP2", "TBR1", "POU3F3", "POU3F2"))
ggsave("June2024/all_CorticL2.png", plot = p1, width = 10, height = 5, units = "in", dpi = 300)	
p1 <- FeaturePlot(ngn2, features = c("SLC17A6", "SLC17A7", "GRIN1", "GRIN2B"))
ggsave("June2024/all_Glut.png", plot = p1, width = 10, height = 5, units = "in", dpi = 300)	
p1 <- FeaturePlot(ngn2, features = c("SOX2", "POU5F1", "NES", "HES1"))
ggsave("June2024/all_Pluri.png", plot = p1, width = 10, height = 5, units = "in", dpi = 300)	

p1 <- FeaturePlot(ngn2, features = c("GRIN2B", "RORB", "THEMIS"))
ggsave("June2024/all_Glut2.png", plot = p1, width = 10, height = 5, units = "in", dpi = 300)	
p1 <- FeaturePlot(ngn2, features = c("GAD1", "SST", "VIP", "PVALB"))
ggsave("June2024/all_GABA.png", plot = p1, width = 10, height = 5, units = "in", dpi = 300)
p1 <- FeaturePlot(ngn2, features = c("PLP1", "MOBP", "ST18", "OPALIN"))
ggsave("June2024/all_Oligo.png", plot = p1, width = 10, height = 5, units = "in", dpi = 300)	
p1 <- FeaturePlot(ngn2, features = c("PDGFRA", "VCAN", "LHFPL3", "MYT1"))
ggsave("June2024/all_OPC.png", plot = p1, width = 10, height = 5, units = "in", dpi = 300)	
p1 <- FeaturePlot(ngn2, features = c("P2RY12", "CSF1R", "APBB1IP", "CD74"))
ggsave("June2024/all_MGL.png", plot = p1, width = 10, height = 5, units = "in", dpi = 300)	
p1 <- FeaturePlot(ngn2, features = c("ALDH1L1", "COL5A3", "GFAP", "AQP4"))
ggsave("June2024/all_Astro.png", plot = p1, width = 10, height = 5, units = "in", dpi = 300)	
p1 <- FeaturePlot(ngn2, features = c("CLDN5", "EPAS1", "VCAM1"))
ggsave("June2024/all_VC.png", plot = p1, width = 10, height = 5, units = "in", dpi = 300)	
p1 <- FeaturePlot(ngn2, features = c("SYN1", "DLG4", "SLC5A7", "ISL1"))
ggsave("June2024/all_Mature_chol.png", plot = p1, width = 10, height = 5, units = "in", dpi = 300)

p1 <- DoHeatmap(object=ngn2, features = marker.genes, size = 3, draw.lines = F)
ggsave("June2024/ngn2heatmapall.png", plot = p1, width = 12, height = 5, units = "in", dpi = 300)	

marker.genes2 <- c("SYN1", "DLG4", "SLC5A7", "ISL1", "NANOG", "DCN", "MAP2", "PRPH", "POU4F1", "LHX9", 
	"GPM6A", "CD99", "RBFOX3", "TUBB3", "MAPT", "MAP2", "CUX1", "CUX2", "POU3F2", "BCL11B", 
	"FOXP2", "TBR1", "POU3F3", "POU3F2", "SLC17A6", "SLC17A7", "GRIN1", "GRIN2B", "SOX2", "POU5F1", "NES", "HES1", 
	"GRIN2B", "RORB", "THEMIS", "GAD1", "SST", "VIP", "PVALB")

p1 <- DoHeatmap(object=ngn2, features = marker.genes2, size = 3, draw.lines = F)
ggsave("plots/ngn2heatmapELF.png", plot = p1, width = 10, height = 5, units = "in", dpi = 300)	

p1 <- DotPlot(object=ngn2, features = marker.genes2) + RotatedAxis()
ggsave("plots/ngn2DotplotELF.png", plot = p1, width = 10, height = 5, units = "in", dpi = 300)	















