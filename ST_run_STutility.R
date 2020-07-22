
# Title: Create Seurat object from ST counts with STutility
# Author: Andrew Ji, Khavari Lab
# 07-01-20
# For more details on STutility package and tutorials, please see: https://ludvigla.github.io/STUtility_web_site/index.html
# For illustrative purposes, this script will create a Seurat object from ST pipeline outputs
# Please download processed example p2st.Rds to further explore if desired

# Install STutility package
devtools::install_github(
  "jbergenstrahle/STUtility"
  ,build_opts=c("--no-resave-data", "--no-manual"), force = T
)

# Load necessary packages
library(knitr)
knitr::opts_chunk$set(echo = TRUE, global.par=TRUE)
library(STutility)
library(Seurat)
library(dplyr)
library(magrittr)

# Load in sequencing data
samples <- list.files(pattern = "stdata.tsv", path = "~/ST_output/", full.names = T, recursive = T)
spotfiles <- list.files(pattern = "selection", path = "~/ST_output/", full.names = T, recursive = T)
imgs <- list.files(pattern = ".jpg", path = "~/ST_output/", full.names = T, recursive = T)
infoTable <- data.frame(samples, spotfiles, imgs, stringsAsFactors = F)
infoTable <- infoTable[1:3, ]
infoTable

# Make Seurat object
p2st <- InputFromTable(infotable = infoTable, 
                         transpose = T, 
                         min.gene.count = 10, 
                         min.gene.spots = 2,
                         min.spot.count = 200)

# Open interactive annotation session (if desired)
p2st <- ManualAnnotation(object=p2st)

# Load in images, mask images
p2st <- LoadImages(object = p2st, verbose = T, xdim = 800)
p2st <- MaskImages(object = p2st, verbose = T)

# Normalize expression across spots
p2st <- SCTransform(p2st, vars.to.regress = c("sample", "nFeature_RNA"), return.only.var.genes = F)

# Run ICA
p2st <- RunICA(p2st) # Can also try RunPCA
ProjectDim(p2st, reduction = "ica", dims = 1:20)

# Cluster
keep.dims <- 1:20
p2st <- FindNeighbors(object = p2st, dims = keep.dims, verbose = FALSE, reduction = "ica") #ICA
#p2st <- FindNeighbors(object = p2st, dims = keep.dims, verbose = FALSE, reduction = "pca") #PCA

# Try clustering at different resolution
p2st <- FindClusters(object = p2st, verbose = FALSE, resolution = 0.6)
p2st <- FindClusters(object = p2st, verbose = FALSE, resolution = 0.8)
p2st <- FindClusters(object = p2st, verbose = FALSE, resolution = 1)
p2st <- FindClusters(object = p2st, verbose = FALSE, resolution = 1.5)

# UMAP
p2st <- RunUMAP(object = p2st, dims = keep.dims, verbose = FALSE, reduction = "ica") # change to reduction = "pca" for PCA after running RunPCA

# Plot spots only in UMAP space
DimPlot(p2st, group.by = "SCT_snn_res.0.8", reduction = "umap",label=T)

# Find markers 
res = "SCT_snn_res.0.8"
Idents(p2st) = res
p2st.markers <- FindAllMarkers(object = p2st, only.pos = T)
write.table(p2st.markers,file="p2st_res0.8_seurat_markers.csv", sep=",", row.names = T, col.names = T)

p2st_top10 = p2st.markers %>% group_by(cluster) %>% top_n(10, avg_logFC) 

# Plot spatial feature plot 
ST.FeaturePlot(p2st,features = c("VWF","PECAM1"))

# Plot spatial overlay feature plot
FeatureOverlay(p2st, features = "SCT_snn_res.0.8", sample.index = 2)

# Save object
saveRDS(p2st,file="p2st.Rds")



