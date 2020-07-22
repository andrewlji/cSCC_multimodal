# Title: Calculate average ligand-receptor pair expression and p-values across normalized scRNA-seq data
# Author: Andrew Ji, Khavari Lab
# 07-01-20

#####
# Overview
# Ligand-receptor analysis is based on CellPhoneDB method described in Vento-Tormo et al., Nature 2018
# Workflow can be divided into 3 parts:
# Part 1: Calculate true average expression of ligand-receptor pairs across cell type pairs
# Part 2: Calculate null distribution of average expression across randomized data
# Part 3: Calculate p-value by comparing true average to null distribution

####
# Note: data from Ji et al. Cell 2020 shown in Figures 6 were generated with Seurat v2 normalization and on a subset of cells annotated earlier than final cell type annotations
# Therefore, running this script may result in slightly different results than what is shown in Figure 6
# To recreate Figure 6 analysis, use the barcodes and annotations in the file "LR_barcode_celltype.txt" to subset the expression data
# Please reach out directly for any additional questions/clarifications
####

#############
# Part 1: Calculate true average expression of ligand-receptor pairs across cell type pairs
library(reticulate)
library(Seurat)
library(dplyr)
library(cowplot)
require(scales)
library(gtools)

# Read in Seurat object
scc_all = readRDS("scc_all.Rds")

# Read in ligand-receptor table
# ligand_receptor_list.csv pulled from Ramilowski et al., Nature Comms. 2016 supplemental data
ligrec_table = read.table("ligand_receptor_list.csv",sep=",",header=T,stringsAsFactors = F)

# Subset tumor data from seven patients with TSK cells 
p7 = c("P1","P2","P5","P6","P7","P9","P10")
tum_lr = subset(scc_all, patient %in% p7 & tum.norm=="Tumor")


# Subset object further to cell types of interest
lr_types = c("Endothelial Cell","Fibroblast","Melanocyte","LC","PDC","CLEC9A","CD1C","Mac","ASDC","Tcell","MDSC",
	"B Cell","NK","Tumor_KC_Basal","Tumor_KC_Cyc","Tumor_KC_Diff","TSK")
tum_lr = subset(tum_lr, level2_celltype %in% lr_types)


# Save object for ease of loading later
saveRDS(tum_lr,file="tum_lr.Rds")


# Make list of cell type pair permutations
ct_pairs_perm = permutations(n=length(lr_types),r=2,v=lr_types,repeats.allowed=T)

# Convert to vector of all cell type pairs
ct_pair_name = c()
for (j in 1:nrow(ct_pairs_perm))
{
	pair = ct_pairs_perm[j,]	
	ct_pair_name = c(ct_pair_name,paste(pair[1],pair[2],sep="_"))
}

# Calculate average expression of each ligand-receptor pair across cell type pairs ()
Idents(tum_lr) = "level2_celltype"
tum_lr_data = as.matrix(tum_lr[["RNA"]]@data)
all_pair_avg = c()
ligrec_name = c()
for (i in 1:nrow(ligrec_table))
{
	lig = ligrec_table[i,2]
	rec = ligrec_table[i,4]
	if (!is.na(match(lig,rownames(tum_lr))) & !is.na(match(rec,rownames(tum_lr))))
	{
		all_avg = c()
		for (j in 1:nrow(ct_pairs_perm))
		{
			pair = ct_pairs_perm[j,]
			cells1 = WhichCells(tum_lr, idents = pair[1])
			cells2 = WhichCells(tum_lr, idents = pair[2])
			pctlig1 = length(which(tum_lr_data[lig,cells1]>0))/length(cells1)
			pctlig2 = length(which(tum_lr_data[lig,cells2]>0))/length(cells2)
			pctrec1 = length(which(tum_lr_data[rec,cells1]>0))/length(cells1)
			pctrec2 = length(which(tum_lr_data[rec,cells2]>0))/length(cells2)
			
			# If at least one gene is expressed in at least 10% of cells, calculate avg
			if ((max(pctlig1,pctlig2,pctrec1,pctrec2))>=0.1)
			{
				avg_pair = mean(c(mean(tum_lr_data[lig,cells1]),mean(tum_lr_data[rec,cells2])))
				all_avg = c(all_avg,avg_pair)
			}

			# If neither lig nor rec gene is expressed in at least 10% of cells, set average to 0
			if ((max(pctlig1,pctlig2,pctrec1,pctrec2))<0.1)
			{
				all_avg = c(all_avg,0)
			}
		
		}
		all_pair_avg = rbind(all_pair_avg,all_avg)
		print(dim(all_pair_avg))
		ligrec_name = c(ligrec_name,paste(lig,rec,sep="_"))	# Lig_rec names for rownames
		print(paste(lig,rec,sep="_"))

	}
	print(length(all_avg))
	print(i)
	print(length(ligrec_name))
	rownames(all_pair_avg) = ligrec_name
	colnames(all_pair_avg) = ct_pair_name

}
colnames(all_pair_avg) = ct_pair_name
rownames(all_pair_avg) = ligrec_name

# Save true average expression values (matrix is 2550 ligand-receptor pairs x 289 cell type pairs)
write.table(all_pair_avg,file="tum_lr_ligrec_avg_pair_exp.txt",sep="\t",row.names=T,col.names=T)


##########################
# Part 2: Calculate null distribution of average expression across randomized data
# Steps:
# 2A) Create permuted gene expression matrix of ligands and receptors across single-cell dataset
# 2B) For each permuted matrix, calculate average expression of each ligand-receptor pair across each cell-type pair (this is null distribution)
# Note: the permuted matrix is a large object; ensure that your computing nodes have sufficient memory to handle it

# PART 2A
# Permute ligand-receptor genes by cells matrix 1000 times
ligrec_genes = unique(c(ligrec_table[,2],ligrec_table[,4]))
ligrec_genes_match = ligrec_genes[which(!is.na(match(ligrec_genes,rownames(tum_lr_data))))]
all_mat = as.matrix(tum_lr_data[ligrec_genes_match,])
cell_id = seq(1,length(ncol(tum_lr)),1)
id_perm = replicate(1000, sample(cell_id,length(cell_id)))

perm_mat_list = list()
for (i in 1:ncol(id_perm))
{
	perm_mat = all_mat[,id_perm[,i]]
	perm_mat_list[[i]] = perm_mat
}

# perm_mat_list = list of 1000 randomized matrices of ligand-receptor gene expression across single cells
save(perm_mat_list,file="perm_mat_list_tum_lr.Robj")

rm(perm_mat_list)
rm(tum_lr)

## Saved list of cell names per cell type
cell_list = list()
for (i in 1:length(lr_types))
{
	type = lr_types[i]
	cells = WhichCells(tum_lr, idents = type)
	cell_list[[i]] = cells
}
names(cell_list) = lr_types
save(cell_list,file="tum_lr_cell_list.Robj")
write.table(all_mat,file="tum_lr_ligrec_exp_mat.txt",sep="\t",row.names=T,col.names=T)

###################
# PART 2B
# Calc null distribution with perm_mat_list
library(reticulate)
library(Seurat)
library(dplyr)
library(cowplot)
require(scales)
library(gtools)


tum_lr_mat = read.table(file="tum_lr_ligrec_exp_mat.txt",sep="\t",row.names=1,header=T,stringsAsFactors=F)
load(file="tum_lr_cell_list.Robj")

setwd("/scratch/PI/khavari/users/andrewji/gene_sig/")
ligrec_table = read.table("ligand_receptor_list.csv",sep=",",header=T,stringsAsFactors = F)

setwd("/scratch/PI/khavari/users/andrewji/10xdata/merge12ptsnew/tum")
load(file="perm_mat_list_tum_lr.Robj")


lr_types = c("Endothelial Cell","Fibroblast","Melanocyte","LC","PDC","CLEC9A","CD1C","Mac","ASDC","Tcell","MDSC",
	"B Cell","NK","Tumor_KC_Basal","Tumor_KC_Cyc","Tumor_KC_Diff","TSK")

Idents(tum_lr) = "level2_celltype"

setwd("/scratch/PI/khavari/users/andrewji/10xdata/merge12ptsnew/tum")

ct_pairs_perm = permutations(n=length(lr_types),r=2,v=lr_types,repeats.allowed=T)


ct_pair_name = c()
for (j in 1:nrow(ct_pairs_perm))
{
	pair = ct_pairs_perm[j,]	
	ct_pair_name = c(ct_pair_name,paste(pair[1],pair[2],sep="_"))
}

# Null distibution = average lig-rec pair expression from random data x number of times (1000 in this case)
# 1000 average pair expression values per ligand-receptor pair
# result should be 2550 x 1000 matrix for each cell type pair (289 cell type pairs) 
# list of length 2550 where each item in list is 1000 x 289 matrix
null_d = list()
ligrec_name = c()
for (i in 1:nrow(ligrec_table))
{
	lig = ligrec_table[i,2]
	rec = ligrec_table[i,4]
	print(i)
	if (!is.na(match(lig,rownames(tum_lr_mat))) & !is.na(match(rec,rownames(tum_lr_mat))))
	{
		all_avg = c()
		for (j in 1:nrow(ct_pairs_perm))
		{
			pair = ct_pairs_perm[j,]
			match1 = match(pair[1],lr_types)
			match2 = match(pair[2],lr_types)
			cells1 = cell_list[[match1]]
			cells2 = cell_list[[match2]]
			cell1_index = match(cells1,colnames(tum_lr_mat))
			cell2_index = match(cells2,colnames(tum_lr_mat))
			lig_index = match(lig,rownames(tum_lr_mat))
			rec_index = match(rec,rownames(tum_lr_mat))
			pctlig1 = length(which(tum_lr_mat[lig,cells1]>0))/length(cells1)
			pctlig2 = length(which(tum_lr_mat[lig,cells2]>0))/length(cells2)
			pctrec1 = length(which(tum_lr_mat[rec,cells1]>0))/length(cells1)
			pctrec2 = length(which(tum_lr_mat[rec,cells2]>0))/length(cells2)
			
			if ((max(pctlig1,pctlig2,pctrec1,pctrec2))>=0.1)
			{
				ligc1 = lapply(perm_mat_list,function(x) x[lig_index,cell1_index])
				ligc1_avg = lapply(ligc1,mean)
				recc2 = lapply(perm_mat_list,function(x) x[rec_index,cell2_index])
				recc2_avg = lapply(recc2,mean)
				avg_pair = (unlist(ligc1_avg)+unlist(recc2_avg))/2
				all_avg = cbind(all_avg,avg_pair)
			}
			if ((max(pctlig1,pctlig2,pctrec1,pctrec2))<0.1)
			{
				all_avg = cbind(all_avg,rep(0,1000))
			}

		}
		null_d[[i]] = all_avg
		print(dim(all_avg))
		ligrec_name = c(ligrec_name,paste(lig,rec,sep="_"))	# Lig_rec names for pairs that pass threshold

	}

	print(dim(all_avg))
	print(i)
	print(ligrec_name[i])
	print(length(ligrec_name))

}
#names(null_d) = ligrec_name
names(null_d) = ligrec_table[,1]
save(null_d,file="tum_lr_null_d_list.Robj")

rm(null_d)
rm(tum_lr)
rm(perm_mat_list)


######################

# Part 3: Calculate p-value by comparing true average to null distribution

# Read in TRUE average expression and null distribution
all_pair_avg = read.table(file="tum_lr_ligrec_avg_pair_exp.txt",sep="\t",row.names=1,header=T,stringsAsFactors=F) # 2550 x 298 matrix
load(file="tum_lr_null_d_list.Robj")

# Some entries may be null in tum_lr_null_d_list.Robj; filter out NULL entries; this will match length of list to number of L-R pairs (nrow) of all_pair_avg
x = null_d
x[sapply(x, is.null)] <- NULL
names(x) = rownames(all_pair_avg)
null_d_filt = x 

# Loop through all averages and compare to the 1000 permutated null distribution
# Want a matrix of 289 (celltype pairs) x 2550 (L-R pairs) with p-vals, which is proportion of null_d avgs that are equal to or exceed actual avg
pval = c()
for (i in 1:nrow(all_pair_avg))
{
	lr_null = null_d_filt[[i]] # lr_null is 1000 x 289 matrix of all randomized averages for all 289 cell type pairs given L-R pair i
	all_prop = c()
	for (j in 1:ncol(all_pair_avg))
	{
		real_avg = all_pair_avg[i,j] # Real avg for L-R pair i and CT pair j
		ct_null_avg = lr_null[,j] # Null d for L-R pair i and CT pair j, should be vector of 1000 values
		prop = length(which(ct_null_avg>=real_avg))/1000
		all_prop = c(all_prop,prop)
	}
	pval = rbind(pval,all_prop) # all_prop is 289 length vector with all p-values of each CT pair for that ligand-receptor pair
	print(dim(pval))
	
}
print(dim(pval))
rownames(pval) = rownames(all_pair_avg)
colnames(pval) = colnames(all_pair_avg)

write.table(pval,file="tum_lr_ligrec_pval.txt",sep="\t",row.names=T,col.names=T) # pval is 2550 x 289 pvalues (same dim as all_pair_avg)

rm(null_d)











