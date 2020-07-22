# Title: Nearest neighbor analysis in ST
# Author: Andrew Ji, Khavari Lab
# 07-01-20
# This script will calculate number of nearest neighbors for every spot in each cluster and compare it to randomized data 
# Workflow can be broken into 3 steps if needed
# 1. Get true nearest neighbor identities
# 2. Get randomized identities of nearest neighbors for null distribution
# 3. Compare true nearest neighbor identity counts to null distribution

library(Seurat)
library(dplyr)
library(magrittr)

# Load object (can be downloaded from ***lab website***)
p2st = readRDS("p2st.Rds")


# Set cluster resolution, separate spots by sample
res = 0.8
Idents(p2st) = paste("SCT_snn_res.",res,sep="")
all_clus = levels(Idents(p2st))
Idents(p2st) = "sample"
s1_spots = WhichCells(p2st, idents = 1)
s2_spots = WhichCells(p2st, idents = 2)
s3_spots = WhichCells(p2st, idents = 3)

##### STEP 1: Get true identities of NN
# Get lists of all nearest cluster identities for each spot across all 3 sections (TRUE DATA)
all_spots_list = list(s1_spots,s2_spots,s3_spots)
all_nn_spots = list()
for (j in 1:3) {
  spots_to_test = all_spots_list[[j]]
  clus = as.vector(p2st@meta.data[spots_to_test,paste("SCT_snn_res.",res,sep="")])
  spot_nn = list()
  for (i in 1:length(spots_to_test))
  {
    spot = spots_to_test[i]
    coor = strsplit(spots_to_test[i],"_")[[1]][1]
    x1 = as.numeric(strsplit(coor,"x")[[1]][1])
    y1 = as.numeric(strsplit(coor,"x")[[1]][2])
    range_x = c(x1-1,x1,x1+1)
    range_y = c(y1-1,y1,y1+1)
    all_combo = expand.grid(range_x,range_y)
    nearest_n = apply(all_combo,1,function(x) paste(x[1],x[2],sep="x"))
    nearest_n_add = paste(nearest_n,"_",j,sep="")
    nearest_n_filt = nearest_n_add[! nearest_n_add %in% spot]
    nearest_n_match = nearest_n_filt[which(!is.na(match(nearest_n_filt,all_spots_list[[j]])))]
    spot_nn[[i]] = clus[match(nearest_n_match,spots_to_test)]
  }
  names(spot_nn) = clus
  all_nn_spots[[j]] = spot_nn
}

# For all 3 sections, tabulate all nearest neighbor clusters
# Result is n_spots x n_spots matrix with counts of nearest neighbors for each cluster
nn_mat = matrix(data=0,nrow=length(all_clus),ncol=length(all_clus))
for (k in 1:3) {
  spot_nn = all_nn_spots[[k]]
  for (i in 1:length(spot_nn))
  {
    clus_nn = all_nn_spots[[k]][[i]]
    #clus_nn = spot_nn[[i]]
    clus_name = as.numeric(names(spot_nn)[i])
    for (j in 1:length(clus_nn))
    {
      nn_index = clus_nn[j]
      x = match(nn_index,all_clus)
      nn_mat[clus_name+1,x] = nn_mat[clus_name+1,x]+1 
    }
  }
}
colnames(nn_mat) = all_clus
rownames(nn_mat) = all_clus

##########################
####### STEP 2: Randomize data for a null distribution
# 1000 permutations, result is a 1000-length list of nearest neighbor lists (like the true list calculated above, but now from randomized data)
spot_clus = as.vector(p2st@meta.data[,paste("SCT_snn_res.",res,sep="")])
clus_perm = replicate(1000, sample(spot_clus,length(spot_clus)))
p2_nulld = list()
for (k in 1:ncol(clus_perm))
{
  clus_scr = clus_perm[,k]
  all_nn_null_mat = list()
  for (j in 1:3) {
    spot_nn_null = list()
    spots_to_test = all_spots_list[[j]]
    clus = clus_scr[match(spots_to_test,colnames(p2st))]
    spot_nn = list()
    for (i in 1:length(spots_to_test))
    {
      spot = spots_to_test[i]
      coor = strsplit(spots_to_test[i],"_")[[1]][1]
      x1 = as.numeric(strsplit(coor,"x")[[1]][1])
      y1 = as.numeric(strsplit(coor,"x")[[1]][2])
      range_x = c(x1-1,x1,x1+1)
      range_y = c(y1-1,y1,y1+1)
      all_combo = expand.grid(range_x,range_y)
      nearest_n = apply(all_combo,1,function(x) paste(x[1],x[2],sep="x"))
      nearest_n_add = paste(nearest_n,"_",j,sep="")
      nearest_n_filt = nearest_n_add[! nearest_n_add %in% spot]
      nearest_n_match = nearest_n_filt[which(!is.na(match(nearest_n_filt,all_spots_list[[j]])))]
      spot_nn[[i]] = clus_scr[match(nearest_n_match,colnames(p2st))]
    }
    names(spot_nn) = clus
    all_nn_null_mat[[j]] = spot_nn
  }
  p2_nulld[[k]] = all_nn_null_mat
  print(k)
}

# Change the lists into n_spots x n_spots matrices
all_mat_null = list()
for (h in 1:1000) {
  nn_mat_null = matrix(data=0,nrow=length(all_clus),ncol=length(all_clus))
  for (k in 1:3) {
    spot_nn = p2_nulld[[h]][[k]]
    for (i in 1:length(spot_nn))
    {
      clus_nn = p2_nulld[[h]][[k]][[i]]
      #clus_nn = spot_nn[[i]]
      clus_name = as.numeric(names(spot_nn)[i])
      for (j in 1:length(clus_nn))
      {
        nn_index = clus_nn[j]
        x = match(nn_index,all_clus)
        nn_mat_null[clus_name+1,x] = nn_mat_null[clus_name+1,x]+1 
        
      }
      
    }
  }
  colnames(nn_mat_null) = all_clus
  rownames(nn_mat_null) = all_clus
  all_mat_null[[h]] = nn_mat_null
  print(h)
}

# change list into array for comparison to TRUE values
p2_nulld_array = array(as.numeric(unlist(all_mat_null)), dim=c(length(all_clus), length(all_clus), 1000))

#### STEP 3: Compare random data to true data
#### Count number of random values >= observed data = P-value
pval_mat = matrix(data=NA,nrow=length(all_clus),ncol=length(all_clus))
for (i in 1:nrow(nn_mat))
{
  for (j in 1:ncol(nn_mat))
  {
    null_d = p2_nulld_array[i,j,]
    print(range(null_d))
    prop = length(which(null_d>nn_mat[i,j]))/1000
    pval_mat[i,j] = prop
  }
  
}
colnames(pval_mat) = all_clus
rownames(pval_mat) = all_clus

write.table(pval_mat,file="p2st_res0.8_pval_mat.txt",sep="\t",row.names = T,col.names = T)
write.table(nn_mat,file="p2st_res0.8_nn_mat.txt",sep="\t",row.names = T,col.names = T)

nn_mat = read.table("p2st_res0.8_nn_mat.txt",sep="\t",row.names=1,header=T)
pval_mat = read.table("p2st_res0.8_pval_mat.txt",sep="\t",row.names=1,header=T)


##### PLOT RESULTS
# Function to get lower half of symmetrical matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}

# Make data frame for ggplots plotting
lower_tri <- get_lower_tri(nn_mat)
lower_tri_p = get_lower_tri(pval_mat)
sig_vec_lt = rep(NA,length(as.vector(lower_tri_p)))
sig_vec_lt[which(lower_tri_p<0.05)] = "*"

rownames(lower_tri) = c("C0","C1","C2","C3","C4","C5","C6","C7","C8","C9","C10","C11")
colnames(lower_tri) = c("C0","C1","C2","C3","C4","C5","C6","C7","C8","C9","C10","C11")
clus2 = c()
for (i in 1:nrow(lower_tri))
{
  clus2 = c(clus2,rownames(lower_tri)[i:nrow(lower_tri)])
}

library(reshape2)
melted_lower_tri <- melt(lower_tri, na.rm = TRUE, value.name = "num_n")
melted_lower_tri_p = melt(lower_tri_p, na.rm=T, value.name = "p_val")
melted_lower_tri_cut = melted_lower_tri
melted_lower_tri_cut$num_n[which(melted_lower_tri_cut$num_n>100)]=100
sig_vec_lt = rep(NA,nrow(melted_lower_tri_p))
sig_vec_lt[which(melted_lower_tri_p$p_val<0.05)] = "*"

melted_lower_tri_cut = cbind(melted_lower_tri_cut,melted_lower_tri_p$p_val,sig_vec_lt,clus2)

p1 = ggplot(melted_lower_tri_cut, aes(x=variable, y=clus2, fill=num_n)) + scale_x_discrete(limits=colnames(lower_tri)) +
  scale_y_discrete(limits=rev(colnames(lower_tri))) + 
  geom_tile() +
  scale_fill_gradient2(low = "white", high = "red", 
                       limit = c(0,100), space = "Lab", 
                       name="Number\nNeighbors") +
  geom_text(aes(variable, clus2, label = sig_vec_lt), color = "black", size = 4) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    #panel.grid.major = element_blank(),
    #panel.border = element_blank(),
    axis.line = element_line(colour = "black"),
    panel.background = element_blank(),
    axis.ticks = element_blank())
print(p1)

ggsave(plot=p1,height=5,width=6.5,dpi=300, filename="p2st_res0.8_nn_heatmap_sig_lower_tri.pdf", useDingbats=FALSE)

