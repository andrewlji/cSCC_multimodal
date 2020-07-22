# cSCC_multimodal
Analysis scripts for cSCC scRNA-seq and ST analysis from Ji et al., Cell 2020.

SCC_create_seurat_object_from_counts.R
- Creates Seurat object from scRNA-seq counts and metadata tables uploaded in GEO GSE144240

SCC_ligand_receptor_scRNA_analysis.R
- Workflow for ligand-receptor analysis in scRNA-seq data to generate p-values based average expression values

ST_run_STutility.R
- Create Seurat object using STutility package (for illustrative purposes)
- For ease, you can download the fully analyzed Seurat object from the script here: https://andrewji.s3-us-west-1.amazonaws.com/p2st.Rds

ST_leading_edge_proximity_analysis.R
- Workflow for ligand-receptor proximity at leading edge based on ST expression values (use example p2st.Rds Seurat object from our patient 2 ST data)

ST_calculate_nearest_neighbor.R
- Tabulates cluster identities of nearest neighbors for spots in ST data and compares to randomized data

Useful files:

LR_barcode_celltype.txt
- Barcodes and cell type annotations for data used to generate scRNA-seq ligand-receptor analysis shown in Figure 6 of Ji et al., Cell 2020. See "SCC_ligand_receptor_scRNA_analysis.R" for additional information.

ligand_receptor_list.csv
- Compiled list of ligand-receptor pairs used for scRNA-seq and ST ligand-receptor analysis pulled from Ramilowski et al., Nature Communications 2015.
