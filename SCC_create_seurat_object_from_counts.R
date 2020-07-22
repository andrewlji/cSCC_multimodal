# Title: Create Seurat object from counts table and run dimensionality reduction, clustering
# Note: Data presented in Ji et al., Cell 2020, Figure 1 was initially created with Seurat v2 followed by UMAP from Seurat v3
# This script is for Seurat v3 and may generate UMAP that looks slightly different from that of Figure 1
# Author: Andrew Ji, Khavari Lab
# 06-30-2020

library(reticulate)
library(Seurat)
library(dplyr)
library(cowplot)


# Load GEO counts table and patient_metadata table
patient_metadata = read.table("patient_metadata_new.txt.gz",sep="\t",row.names=1,header=T,stringsAsFactors=F)
counts = read.table("cSCC_counts.txt.gz",sep="\t",row.names=1,header=T,stringsAsFactors=F)

# Filter out top 2 rows of metadata in counts table
counts_filt = counts[3:32740,]

# Create Seurat Object
scc_all = CreateSeuratObject(counts = counts_filt, project = "test", meta.data = patient_metadata)

# Set mitochondrial gene %
scc_all[["percent.mt"]] <- PercentageFeatureSet(scc_all, pattern = "^MT-")

# Normalize data and scale data
scc_all <- NormalizeData(scc_all)
scc_all <- FindVariableFeatures(scc_all, selection.method = "vst", nfeatures = 2000)

# Scale data (this step can take awhile; for speed, consider using parallelization with future package, details: https://satijalab.org/seurat/v3.0/future_vignette.html)
all.genes <- rownames(scc_all)
scc_all <- ScaleData(scc_all, features = all.genes, vars.to.regress = c("nCount_RNA","percent.mt"))

# Run PCA
scc_all <- RunPCA(scc_all, features = VariableFeatures(object = scc_all))

# Cluster
scc_all <- FindNeighbors(scc_all, dims = 1:15, verbose = FALSE)
scc_all <- FindClusters(scc_all, resolution = 0.1, verbose = FALSE)

# Run UMAP
scc_all <- RunUMAP(scc_all, dims = 1:15, verbose = TRUE)



# Plot by meta.data feature
png("scc_all_level1_celltype.png", units = "in", height = 5, width = 6.5, res=300)
p=DimPlot(object = scc_all,group.by = "level1_celltype")
print(p)
dev.off()

png("scc_all_level2_celltype.png", units = "in", height = 5, width = 8, res=300)
p=DimPlot(object = scc_all,group.by = "level2_celltype")
print(p)
dev.off()

png("scc_all_level3_celltype.png", units = "in", height = 5, width = 8, res=300)
p=DimPlot(object = scc_all,group.by = "level3_celltype")
print(p)
dev.off()

png("scc_all_patient.png", units = "in", height = 5, width = 6, res=300)
p=DimPlot(object = scc_all,group.by = "patient")
print(p)
dev.off()

png("scc_all_tumnorm.png", units = "in", height = 5, width = 6.5, res=300)
p=DimPlot(object = scc_all,group.by = "tum.norm")
print(p)
dev.off()

png("scc_all_res.0.1.png", units = "in", height = 5, width = 6, res=300)
p=DimPlot(object = scc_all,label=T)
p = p + NoAxes()
print(p)
dev.off()


# Save scc_all
saveRDS(scc_all,file="scc_all.Rds")




