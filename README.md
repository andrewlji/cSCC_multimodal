# cSCC_multimodal
Analysis scripts for cSCC scRNA-seq and ST analysis

SCC_create_seurat_object_from_counts.R
- Creates Seurat object from counts and metadata tables uploaded in GEO GSE144240

SCC_ligand_receptor_scRNA_analysis.R
- Workflow for ligand-receptor analysis in scRNA-seq data to generate p-values based average expression values

ST_run_STutility.R
- Create Seurat object using STutility package (for illustrative purposes)

ST_leading_edge_proximity_analysis.R
- Workflow for ligand-receptor proximity at leading edge based on ST expression values (use example p2st.Rds Seurat object from our patient 2 ST data)

ST_calculate_nearest_neighbor.R
- Tabulates cluster identities of nearest neighbors for spots in ST data and compares to randomized data
