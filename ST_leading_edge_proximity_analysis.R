# Title: Calculate average ligand-receptor pair expression and p-values across leading edge ST data
# Author: Andrew Ji, Khavari Lab
# 07-01-20
# Very similar to SCC_ligand_receptor_scRNA_analysis.R script

#####
# Overview
# Ligand-receptor analysis is based on CellPhoneDB method described in Vento-Tormo et al., Nature 2018
# Instead of cell type pairs, however, expression is averaged over "neighborhoods" or "sliding windows" of spots consisting of a central spot and its 4 nearest neighbors (fewer NN if on the edge of tissue)
# Workflow can be divided into 3 parts:
# Part 1: Calculate true average expression of ligand-receptor pairs across neighborhoods of spot and 4 nearest neighbors (fewer NN if on the edge of tissue)
# Part 2: Calculate null distribution of average expression across neighborhoods in randomized data
# Part 3: Calculate p-value by comparing true average of L-R expression across leading edge to the null distribution


######
# Part 1: Calculate average expression of each ligand-receptor pair in "neighborhoods"
library(reticulate)
library(Seurat)
library(dplyr)

p2st = readRDS("p2st.Rds")

ligrec_table = read.table("ligand_receptor_list.csv",sep=",",header=T,stringsAsFactors = F)
ligrec_genes = unique(c(ligrec_table[,2],ligrec_table[,4]))
ligrec_genes_match = ligrec_genes[ligrec_genes %in% rownames(p2st)]
keep_ligrec_table = c()
for (i in 1:nrow(ligrec_table))
{
  lig = ligrec_table[i,2]
  rec = ligrec_table[i,4]
  if (lig %in% rownames(p2st) & rec %in% rownames(p2st)) 
  {keep_ligrec_table = rbind(keep_ligrec_table,ligrec_table[i,])}
}

# Get sample spots
Idents(p2st) = "sample"
s1_spots = WhichCells(p2st, idents = 1)
s2_spots = WhichCells(p2st, idents = 2)
s3_spots = WhichCells(p2st, idents = 3)
all_spots_list = list(s1_spots,s2_spots,s3_spots)


# Get lists of all nearest cluster identities for each spot across all 3 sections
all_lr_mean = c()
for (i in 1:nrow(keep_ligrec_table))
{
  lig = keep_ligrec_table[i,2]
  rec = keep_ligrec_table[i,4]
  lr_mean_spot = c()
  for (j in 1:3) {
    spots_to_test = all_spots_list[[j]]

    for (k in 1:length(spots_to_test))
    {
      spot = spots_to_test[k]
      coor = strsplit(spots_to_test[k],"_")[[1]][1]
      x1 = as.numeric(strsplit(coor,"x")[[1]][1])
      y1 = as.numeric(strsplit(coor,"x")[[1]][2])
      range_x = c(x1-1,x1,x1+1)
      range_y = c(y1-1,y1,y1+1)
      all_combo = expand.grid(range_x,range_y)
      nearest_n = apply(all_combo,1,function(x) paste(x[1],x[2],sep="x"))
      nearest_n_add = paste(nearest_n,"_",j,sep="")
      nearest_n_match = nearest_n_add[which(!is.na(match(nearest_n_add,all_spots_list[[j]])))]

      lig_exp = p2st[["SCT"]]@data[lig,nearest_n_match]
      rec_exp = p2st[["SCT"]]@data[rec,nearest_n_match]
      mean_lig_exp = mean(lig_exp)
      mean_rec_exp = mean(rec_exp)
      mean_lr = (mean_lig_exp+mean_rec_exp)/2
      lr_mean_spot = c(lr_mean_spot, mean_lr)

    }
  }
  all_lr_mean = rbind(all_lr_mean,lr_mean_spot)
  print(i)
}
rownames(all_lr_mean) = keep_ligrec_table[,1]
colnames(all_lr_mean) = unlist(all_spots_list)

write.table(all_lr_mean,file="p2st_spots_lr_mean_nn.txt",sep="\t",row.names = T,col.names = T)

#rm(p2st)

#############

# Part 2: Calculate null distribution of average expression across neighborhoods in randomized data
ligrec_table = read.table("ligand_receptor_list.csv",sep=",",header=T,stringsAsFactors = F)
ligrec_genes = unique(c(ligrec_table[,2],ligrec_table[,4]))
ligrec_genes_match = ligrec_genes[ligrec_genes %in% rownames(p2st)]
keep_ligrec_table = c()
for (i in 1:nrow(ligrec_table))
{
  lig = ligrec_table[i,2]
  rec = ligrec_table[i,4]
  if (lig %in% rownames(p2st) & rec %in% rownames(p2st)) 
  {keep_ligrec_table = rbind(keep_ligrec_table,ligrec_table[i,])}
}

unique_lr_genes = unique(c(keep_ligrec_table[,2],keep_ligrec_table[,4]))

# Get sample spots
Idents(p2st) = "sample"
s1_spots = WhichCells(p2st, idents = 1)
s2_spots = WhichCells(p2st, idents = 2)
s3_spots = WhichCells(p2st, idents = 3)
all_spots_list = list(s1_spots,s2_spots,s3_spots)

# Generate randomized data 
random_s1 = replicate(1000,sample(s1_spots,length(s1_spots)))
random_s2 = replicate(1000,sample(s2_spots,length(s2_spots)))
random_s3 = replicate(1000,sample(s3_spots,length(s3_spots)))

# Generate list of 1000 randomized spots while maintaining sample structure
shuffled_lists = list()
for (i in 1:1000) {
  com_list = list(random_s1[,i],random_s2[,i],random_s3[,i])
  shuffled_lists[[i]] = com_list
}

# Get 1000 expression matrices using randomized spot structure
exp_mat = as.matrix(p2st[["SCT"]]@data[unique_lr_genes,])
exp_mat_list = list()
for (i in 1:1000) {
  spot_id = shuffled_lists[[i]]
  exp_mat_shuffle = exp_mat[,unlist(spot_id)]
  exp_mat_list[[i]] = exp_mat_shuffle
}

# Get true nearest neighbor list
true_nn_list = list()
for (j in 1:3) {
    spots_to_test = all_spots_list[[j]]
    true_nn_list[[j]] = list()
    for (k in 1:length(spots_to_test))
    {
      spot = spots_to_test[k]
      coor = strsplit(spots_to_test[k],"_")[[1]][1]
      x1 = as.numeric(strsplit(coor,"x")[[1]][1])
      y1 = as.numeric(strsplit(coor,"x")[[1]][2])
      range_x = c(x1-1,x1,x1+1)
      range_y = c(y1-1,y1,y1+1)
      all_combo = expand.grid(range_x,range_y)
      nearest_n = apply(all_combo,1,function(x) paste(x[1],x[2],sep="x"))
      nearest_n_add = paste(nearest_n,"_",j,sep="")
      #nearest_n_filt = nearest_n_add[! nearest_n_add %in% spot]
      nearest_n_match = nearest_n_add[which(!is.na(match(nearest_n_add,all_spots_list[[j]])))]
      nearest_n_id = match(nearest_n_match,colnames(p2st))
      true_nn_list[[j]][[k]] = nearest_n_id
    }
}

# Get list of NN for all spots in relation to entire object
true_nn_list_all = unlist(true_nn_list, recursive=F)


# Get list of null distribution: 1540 L-R pairs of 1000 x 1949 spot neighborhood 
lr_spot_nulld = list()
for (k in 1:nrow(keep_ligrec_table)) {
  lig = keep_ligrec_table[k,2]
  rec = keep_ligrec_table[k,4]
  lig_index = match(lig,rownames(exp_mat))
  rec_index = match(rec,rownames(exp_mat))
  all_lr_mean = c()
  for (j in 1:ncol(p2st))
  {
    nearest_n_match = true_nn_list_all[[j]]
    lig_exp = lapply(exp_mat_list,function(x) x[lig_index,nearest_n_match]) # list of 1000 vectors of lig exp
    mean_lig_exp = lapply(lig_exp,mean)
    rec_exp = lapply(exp_mat_list,function(x) x[rec_index,nearest_n_match])
    mean_rec_exp = lapply(rec_exp,mean)
    mean_lr = (unlist(mean_lig_exp)+unlist(mean_rec_exp))/2
    all_lr_mean = cbind(all_lr_mean,mean_lr)
  }
  lr_spot_nulld[[k]] = all_lr_mean
  print(dim(all_lr_mean))
  print(k)

}

saveRDS(lr_spot_nulld,file="lr_spot_nulld.Rds")
#rm(p2st)
#rm(lr_spot_nulld)

################################

# Part 3: Calculate p-value by comparing true average of L-R expression across leading edge to the null distribution
library(reticulate)
library(Seurat)
library(dplyr)

lr_spot_nulld = readRDS(file="lr_spot_nulld.Rds")


# Calculate p-values from null distribution

# Load in true average data
all_pair_avg = read.table(file="p2st_spots_lr_mean_nn.txt",sep="\t",row.names=1,header=T,stringsAsFactors=F) # 1540 x 1949 matrix

# Calculate p-val of leading edge L-R average (isolate leading edge spots and calculate average over null distibution)

# Isolate all leading edge spots (structures such as the leading edge can be annotated in several ways, including with STutility package for 2K ST data, 10X Loupe Browser, or manually in a program such as Adobe Illustrator)
# In this case, we used STutility
Idents(p2st) = "labels"
le_spots = WhichCells(p2st, ident = "LeadEdge")

le_spot_idx = match(le_spots,colnames(p2st))
all_pair_avg_le = all_pair_avg[,le_spot_idx] # Subsetted matrix of ligand-receptor averages across leading edge spots
all_pair_avg_le_lr = apply(all_pair_avg_le,1,mean) # True average across all leading edge spots of each ligand-receptor pair
pval_le = c()
ligrec_name = c()
all_prop = c()
for (i in 1:nrow(all_pair_avg_le))
{
  lr_null = lr_spot_nulld[[i]] # lr_null is 1000 x 1949 matrix 
  lr_null_le = lr_null[,le_spot_idx] #1000 x length(le_spots_idx) matrix
  real_avg = all_pair_avg_le_lr[i] 
  ct_null_avg = apply(lr_null_le,1,mean) # Null d for L-R pair i over le_spots, should be 1000 values
  prop = length(which(ct_null_avg>=real_avg))/1000
  all_prop = c(all_prop,prop)
  print(i)
}
names(all_prop) = rownames(all_pair_avg_le)

# Save p-values in table
write.table(all_prop,file="p2st_lr_nn_pval_le.txt",sep="\t",row.names=T,col.names=T) # pval is named vector of pvalues (same dim as all_pair_avg)

# Merge expression and p-values into one table
le_avg_pair_pval = cbind(all_pair_avg_le_lr,all_prop)
le_avg_pair_pval_order = le_avg_pair_pval[order(le_avg_pair_pval[,1],decreasing=T),]
colnames(le_avg_pair_pval) = c("Average_Pair_Exp","P-val")
write.table(le_avg_pair_pval,file="p2st_lr_le_avg_exp_pval_table.txt",sep="\t",row.names=T,col.names=T)





