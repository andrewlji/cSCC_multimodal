

########################
########################

# Process MIBI annotated table to perform analysis in Figure 5 of Ji et al. 2020 
# Includes extra analysis not contained in paper figures

########################
########################

# Load annotated MIBI table of segmented cells with positions and labels 

mat = read.table("MIBI_celltype_annotations.csv.gz",stringsAsFactors=F,sep=",",header=T,row.names=1)

########################
########################

# Load modules

require(gplots)
require(RColorBrewer)
require(raster)
require(tiff)
require(scales)
require(Rtsne)
require(ggplot2)
require(grid)
require(gridExtra)
require(viridis)

#######################
#######################

# Make cell label scatter plots
# Overlay cell type annotation for each FOV

unique(mat$cell_type_label)
cell_types = c("Keratinocyte","LeadingEdge","CyclingK","HairFollicle","Fibroblast","Neutrophil","CD4","B","Macrophage","CD8","Tregs","DC","NK","Endothelial")

par(mfrow=c(3,6),mar=c(2,2,2,2))

my_palette <- colorRampPalette(brewer.pal(6,"Set2"))(n =length(cell_types))

for (idx1 in c(1:17)){
fov_sel=idx1
point_sel = mat$PointNum==fov_sel
#pdf(paste("Image_FOV_",fov_sel,".pdf",sep=""),width=5,height=4)
#par(mar=c(4, 4, 4, 8), xpd=TRUE)
plot(mat$x_cent[point_sel],2000-mat$y_cent[point_sel],cex=0.1,pch=19,col="black",xlab="X (pixels)",ylab="Y (pixels)",main=paste("FOV ",fov_sel)) # plot Y as 2000-Y

for (idx2 in c(1:length(cell_types))){
	sel_type = mat$cell_type_label == cell_types[idx2]
	if (cell_types[idx2]=="Fibroblast"){
		sel_pch = 19
		sel_cex = 0.7
		sel_color = my_palette[idx2]}
	else {
		sel_pch=19
		sel_cex=0.8
		sel_color  = my_palette[idx2]}
	
	points(mat$x_cent[point_sel&sel_type],2000-mat$y_cent[point_sel&sel_type],cex=sel_cex,pch=sel_pch,col=sel_color)}

legend(2150,2000, cell_types,fill=my_palette,horiz=F,inset=0,cex=0.8)
#dev.off()
}

plot(mat$x_cent[point_sel],2000-mat$y_cent[point_sel],cex=0.1,pch=19,col="white",xlab="X (pixels)",ylab="Y (pixels)",main="") # plot Y as 2000-Y
legend(200,2000, cell_types,fill=my_palette,horiz=F,inset=0,cex=0.9)


########################
########################

# Make cell label scatter plots
# Overlay tumor and non-tumor classes for each FOV

unique(mat$cell_type_label)
cell_types = c("Keratinocyte","LeadingEdge","CyclingK","HairFollicle","Fibroblast","Neutrophil","CD4","B","Macrophage","CD8","Tregs","DC","NK","Endothelial")


sel_tumor = mat$cell_type_label %in% c("Keratinocyte","LeadingEdge","CyclingK","HairFollicle")
sel_nontumor = mat$cell_type_label %in% c("Fibroblast","Neutrophil","CD4","B","Macrophage","CD8","Tregs","DC","NK","Endothelial")

sel_cex = 0.2
sel_pch = 19

par(mfrow=c(3,6),mar=c(2,2,2,2))

for (idx1 in c(1:17)){
fov_sel=idx1
point_sel = mat$PointNum==fov_sel

plot(mat$x_cent[point_sel],2000-mat$y_cent[point_sel],cex=0.1,pch=19,col="black",xlab="X (pixels)",ylab="Y (pixels)",main=paste("FOV ",fov_sel)) # plot Y as 2000-Y

	points(mat$x_cent[point_sel&sel_tumor],2000-mat$y_cent[point_sel& sel_tumor],cex=sel_cex,pch=sel_pch,col=alpha("royalblue1",0.5))
		points(mat$x_cent[point_sel& sel_nontumor],2000-mat$y_cent[point_sel& sel_nontumor],cex=sel_cex,pch=sel_pch,col=alpha("red3",0.5))

}

plot(mat$x_cent[point_sel],2000-mat$y_cent[point_sel],cex=0.1,pch=19,col="white",xlab="X (pixels)",ylab="Y (pixels)",main="") # plot Y as 2000-Y
legend(200,2000, c("Tumor","Non-tumor"),fill=c("royalblue1","red3"),horiz=F,inset=0,cex=0.9)

########################
########################

# Make cell label scatter plots
# Overlay a single cell type annotation in each plot, for one FOV (adjust "sel_fov" variable)
# Label epithelial cells blue

sel_fov = 13 # Select a single FOV (1-17)

unique(mat$cell_type_label)
cell_types = unique(mat$cell_type_label)
cell_types = c("Keratinocyte","LeadingEdge","CyclingK","HairFollicle","Fibroblast","Neutrophil","CD4","B","Macrophage","CD8","Tregs","DC","NK","Endothelial")

my_palette <- colorRampPalette(brewer.pal(6,"Set2"))(n =length(cell_types))


par(mfrow=c(3,5),mar=c(2,2,2,2))

point_sel = mat$PointNum==sel_fov

for (idx2 in c(1:length(cell_types))){
	sel_pch=19
	sel_cex=1
	plot(mat$x_cent[point_sel],2000-mat$y_cent[point_sel],cex=0.1,pch=19,col="white",xlab="X (pixels)",ylab="Y (pixels)",main=paste("FOV ", sel_fov,cell_types[idx2])) # plot Y as 2000-Y

	sel_kc = mat$cell_type_label %in% c("Keratinocyte","LeadingEdge","CyclingK","HairFollicle")	
	points(mat$x_cent[point_sel& sel_kc],2000-mat$y_cent[point_sel& sel_kc],cex=sel_cex,pch=sel_pch,col=alpha("dodgerblue",0.5))

	sel_type = mat$cell_type_label == cell_types[idx2]

	sel_color  = my_palette[idx2]
	points(mat$x_cent[point_sel&sel_type],2000-mat$y_cent[point_sel&sel_type],cex=sel_cex,pch=sel_pch,col=sel_color)
	points(mat$x_cent[point_sel&sel_type],2000-mat$y_cent[point_sel&sel_type],cex=sel_cex,pch=1,col="black")
}

plot(mat$x_cent[point_sel],2000-mat$y_cent[point_sel],cex=0.1,pch=19,col="white",xlab="X (pixels)",ylab="Y (pixels)",main="") # plot Y as 2000-Y
legend(200,2000, cell_types,fill=my_palette,horiz=F,inset=0,cex=0.7)


########################
########################

# Generate heat maps of marker correlation across cells and expression in each cell type


# Heat map: correlation of marker signals across cells
#mat_markers = mat[,8:39]
mat_markers = mat[,c(8,10,9,11,13,14,15,17,18,20,21,25,26,27,28,29,30,31,32,33,35,36,38,39)] # limit to subset of "high confidence" markers exhibiting consistent expression with 10x scRNA
mat_cor = cor(mat_markers,method="pearson")

colors = c(seq(0,1,length=22))
my_palette <- rev(colorRampPalette(brewer.pal(6,"RdBu"))(n = 21))

hm2_call = heatmap.2(mat_cor,col=my_palette,breaks=colors,density.info="none",trace="none",Rowv=T,Colv=T,dendrogram="both",symm=T,labRow=colnames(mat_markers),labCol=colnames(mat_markers),margins=c(5,5),scale="none",cexRow=0.6,cexCol=0.5,rowsep=c(0:50),colsep=c(0:50),sepcolor="black",sepwidth=c(0.0001,0.0001))


# Heat map: average marker expression for each *annotated* cell type 
#cell_types = unique(mat$cell_type_label)
cell_types = c("Keratinocyte","LeadingEdge","CyclingK","Fibroblast","Neutrophil","CD4","B","Macrophage","CD8","Tregs","Endothelial")

cell_type_mat = matrix(0,nrow=length(cell_types),ncol=ncol(mat_markers))

for (idx in c(1:length(cell_types))){
	print(idx)
	row_match = mat$cell_type_label==cell_types[idx] 
	print(length(which(row_match)))
	cell_type_mat[idx,] = colMeans(mat_markers[row_match,])
}


# Heat map of 0-1 values
colors = c(seq(0,1,length=22))
my_palette <- rev(colorRampPalette(brewer.pal(6,"RdBu"))(n = 21))
my_palette <- colorRampPalette(brewer.pal(6,"PuRd"))(n = 21)

hm2_call = heatmap.2(cell_type_mat,col=my_palette,breaks=colors,density.info="none",trace="none",Rowv=F,Colv=T,dendrogram="col",symm=F,labRow=cell_types,labCol=colnames(mat_markers),margins=c(6,6),scale="none",cexRow=0.7,cexCol=0.7,rowsep=c(0:20),colsep=c(0:32),sepcolor="black",sepwidth=c(0.0001,0.0001))


# Heat map of z-score values across cell types
colors = c(seq(-3,3,length=22))
my_palette <- rev(colorRampPalette(brewer.pal(6,"RdBu"))(n = 21))

hm2_call = heatmap.2(cell_type_mat,col=my_palette,breaks=colors,density.info="none",trace="none",Rowv=F,Colv=F,dendrogram="none",symm=F,labRow=cell_types,labCol=colnames(mat_markers),margins=c(6,6),scale="col",cexRow=0.6,cexCol=0.5,rowsep=c(0:20),colsep=c(0:32),sepcolor="black",sepwidth=c(0.0001,0.0001))


########################
########################

# Explore non-tumor cell type abundance across FOVs
#   Generate stacked bar plot, heatmap for correlation of proportions across FOVs, 
#	scatter plot of abundance for two cell types across FOVs, and heat map of proportions in each FOV


# Stacked bar plot of cell types in each point and each patient 
#cell_types = c("Keratinocyte","LeadingEdge","CyclingK","HairFollicle","Fibroblast","Neutrophil","CD4","DC","B","Macrophage","CD8","Tregs","NK","Endothelial")
cell_types = c("Fibroblast","Neutrophil","CD4","B","Macrophage","CD8","Tregs","Endothelial")
fov_names = c("1_T21","2_T21","3_T26","4_T26","5_T27","6_T27","7_T23","8_T23","9_T26","10_T26","11_T21","12_T27","13_T1","14_T1","15_T1","16_T10","17_T10")

point_ids = unique(mat$PointNum)

# Rows = cell types, columns = points
cell_type_count = matrix(0,nrow=length(cell_types),ncol=length(point_ids))

# Count number of each cell type in each point
for (idx in c(1:length(point_ids))){
	print(idx)
	cell_types_in_point = mat$cell_type_label[mat$PointNum==point_ids[idx]]
	for (idx2 in c(1:length(cell_types))){
		ct_count = length(which(cell_types_in_point==cell_types[idx2]))
		cell_type_count[idx2,idx] = ct_count
}}


# Color for each cell type
my_palette <- colorRampPalette(brewer.pal(6,"Set2"))(n=length(cell_types))

rownames(cell_type_count) = cell_types

# Normalize points (colums) to 1
cell_type_count_norm = apply(cell_type_count,2,function(x) x/sum(x))

# Reorder columns by row (cell type) values of interest
col_order = order(cell_type_count_norm[1,])
#col_order = hm2_call$colInd # Use this to rearrange columns by clustering with later heatmap.2 call 

# Reorder rows (cell types) by total abundance of each cell type
row_order = order(apply(cell_type_count_norm,1,sum))

par(mar=c(2, 4, 4, 8), xpd=TRUE)
bp = barplot(cell_type_count_norm[row_order,col_order],col=my_palette,axisnames=F,cex.names=0.6)
text(bp,par("usr")[3],labels= fov_names[col_order],srt=0.2,adj=c(0.6,1.1),xpd=T,cex=0.5)
legend(21,1,rownames(cell_type_count_norm)[row_order],fill=my_palette,horiz=F,inset=0,cex=0.6)


# Bar plot of total non-tumor cell type counts
cell_count_total = apply(cell_type_count,2,sum)
bp = barplot(cell_count_total[col_order],col="black")
text(bp,par("usr")[3],labels=c(1:17)[col_order],srt=0,adj=c(1.1,1.1),xpd=T,cex=0.9)


# Heat map: correlation of cell type counts across FOVs
cell_types_cor = cor(t(cell_type_count_norm),method="pearson")

colors = c(seq(-1,1,length=22))
my_palette <- rev(colorRampPalette(brewer.pal(6,"RdBu"))(n = 21))

hm2_call = heatmap.2(cell_types_cor,col=my_palette,breaks=colors,density.info="none",trace="none",Rowv=T,Colv=T,dendrogram="both",symm=T,labRow=rownames(cell_type_count),labCol=rownames(cell_type_count),margins=c(6,6),scale="none",cexRow=0.6,cexCol=0.5,rowsep=c(0:12),colsep=c(0:12),sepcolor="black",sepwidth=c(0.0001,0.0001))


# Scatter plot of abundance of any two cell types
plot(cell_type_count_norm[6,], cell_type_count_norm[7,],pch=19,xlab="CD8 proportion",ylab="Treg proportion",xlim=c(0,0.3),ylim=c(0,0.3))
abline(0,1,lty=2,col="darkred",lwd=3)
cor(cell_type_count_norm[6,],cell_type_count_norm[7,])

# Heat map of non-tumor cell type proportions by FOV
colors = c(seq(0,0.75,length=20))
my_palette <- plasma(19)
fov_names = c("1_T21","2_T21","3_T26","4_T26","5_T27","6_T27","7_T23","8_T23","9_T26","10_T26","11_T21","12_T27","13_T1","14_T1","15_T1","16_T10","17_T10")

hm2_call = heatmap.2(cell_type_count_norm,col=my_palette,breaks=colors,density.info="none",trace="none",Rowv=T,Colv=T,dendrogram="both",symm=F,labRow=rownames(cell_type_count),labCol=fov_names,margins=c(6,6),scale="none",cexRow=0.6,cexCol=0.8,rowsep=c(0:12),colsep=c(0:18),sepcolor="black",sepwidth=c(0.0001,0.0001))



########################
########################

# ** NOTE this section is NOT a final analysis included in the paper. **
#   Overall we did not find that assessing variable marker expression based on proximity yielded robust results, however this framework may be refined or applied more succesffuly to other datasets. 

# Explore the relationship between cell type proximity and marker expression
# Are features expressed differently for cell type 1 close vs far from cell type 2?

cell_type1 = "CD8" 
cell_type2 = "Tregs"

label2_mean_dists = c()
mat_label1_all = matrix(,nrow=0,ncol=ncol(mat)) # make mat subset for cells included in analysis 

for (pt_idx in c(1:17)){
	
mat_point = mat[mat$PointNum==pt_idx,]	
print(pt_idx)
print(nrow(mat_point))
print("")

mat_label1 = mat_point[mat_point$cell_type_label== cell_type1,]
mat_label2 = mat_point[mat_point$cell_type_label== cell_type2,]

if (nrow(mat_label1)==0 | nrow(mat_label2)==0){
	next
}

mat_label1_all = rbind(mat_label1_all,mat_label1)

for (idx1 in c(1:nrow(mat_label1))){
	coord1 = c(mat_label1$x_cent[idx1], mat_label1 $y_cent[idx1])
	coord2_list = cbind(mat_label2$x_cent,mat_label2$y_cent)
	label2_dists = pointDistance(coord1, coord2_list,lonlat=F)
	label2_mean_dists = c(label2_mean_dists,mean(label2_dists))
	}
}

## End of loop

# Check mean/min distance distribution 
plot(density(label2_mean_dists),"CD8 distance to Treg",xlab="Distance (pixels)")
abline(v=label2_mean_dists[order(label2_mean_dists)][101],lty=2,col="red")
abline(v=label2_mean_dists[order(label2_mean_dists,decreasing=T)][101],lty=2,col="red")

hist(label2_mean_dists,n=60,main="CD8 distance to Treg")

# Compare feature expression in proximal vs distal cell type 1 
label2_mean_dists_z = (label2_mean_dists-mean(label2_mean_dists))/sd(label2_mean_dists)
sel1 = label2_mean_dists_z<(-2)
sel2 = label2_mean_dists_z>2
length(which(sel1))
length(which(sel2))

# What is feature expression distribution of proximal vs distal populations? 
par(mfrow=c(4,8),mar=c(2,2,2,2),cex.main=0.7)
for (idx in c(1:32)){
	box_list = list(mat_label1_all[sel1,idx+7], mat_label1_all[sel2,idx+7])
	ksp = round(ks.test(box_list[[1]],box_list[[2]])$p.value,4)
boxplot(box_list,outline=F,names=c("Prox","Dist"),col=c("cornflowerblue","firebrick2"),main=paste(colnames(mat)[idx+7],ksp,sep="_"))
	print(colnames(mat)[idx+7])
	print(round(cor(label2_mean_dists, mat_label1_all[idx+7]),4))
	print("")
}

# What proportion of cells are "high" for a feature in proximal vs distal populations?
par(mfrow=c(4,8),mar=c(2,2,2,2),cex.main=0.7)
for (idx in c(1:32)){
	high1 = length(which(sel1 & mat_label1_all[,idx+7]>0.5))/length(which(sel1))
	high2 = length(which(sel2 & mat_label1_all[,idx+7]>0.5))/length(which(sel2))
	barplot(c(high1,high2),main=colnames(mat_label1_all)[idx+7])
}

# Do high vs low expressing cells have distinct label2 mean distance distribution? 
par(mfrow=c(4,8),mar=c(2,2,2,2),cex.main=0.7,cex.axis=0.6)
for (idx in c(1:32)){
	sel_high = mat_label1_all[,idx+7]>0.8
	sel_low = mat_label1_all[,idx+7]<0.2
	box_list = list(label2_mean_dists[sel_low], label2_mean_dists[sel_high])
	ksp=NA
boxplot(box_list,outline=T,names=c("Low","High"),col=c("cornflowerblue","firebrick2"),main=paste(colnames(mat)[idx+7],ksp,sep="_"))
}


# Setup for ggplot2 violin plots
mat_label1_all$Treg_prox = rep("NA",nrow(mat_label1_all))
mat_label1_all$Treg_prox[sel1] = "Prox"
mat_label1_all$Treg_prox[sel2] = "Dist"
mat_label1_all$Treg_prox = as.factor(mat_label1_all$Treg_prox)
p = list()
for (idx in c(1:32)){
feature = colnames(mat_label1_all)[idx+7]
p[[idx]] = ggplot(mat_label1_all, aes_string(x="Treg_prox",y=feature),main=feature) + geom_violin(trim=F) + scale_x_discrete(limits=c("Prox","Dist")) + geom_boxplot(width=0.1,color=c("red","blue"))
}
do.call(grid.arrange,p)


# Heat map of median or mean feature expression in each proximity selection group
mat_cell = mat_label1_all[,c(7:39)]
mat_sel_avg = cbind(apply(mat_cell,2,function(x) median(x[sel1])),apply(mat_cell,2,function(x) median(x[sel2])))
colnames(mat_sel_avg) = c("Prox","Dist")

colors = c(seq(0,0.5,length=22))
my_palette <- colorRampPalette(brewer.pal(6,"YlGnBu"))(n = 21)

hm2_call = heatmap.2(mat_sel_avg,col=my_palette,breaks=colors,density.info="none",trace="none",Rowv=T,Colv=F,dendrogram="row",symm=F,labRow= rownames(mat_sel_avg),labCol=colnames(mat_sel_avg),margins=c(3,5),scale="none",cexRow=0.6,cexCol=0.8,rowsep=c(0:50),colsep=c(0:3),sepcolor="black",sepwidth=c(0.0001,0.0001))


#######################
#######################

# Are particular cell types co-localized in stroma? 
# Calculate distance distribution for between two cell types for real and non-tumor sampled positions
# Null distribution from 100x permutations of label set 1 

non_tumor_types = c("Fibroblast","Neutrophil","CD4","B","Macrophage","Endothelial","CD8","Tregs","DC")
label1="CD8"
label2="Tregs"
num_perm = 100 # how many permutations for FDR?

# Find distribution of distances to nearest and mean type Y

par(mfrow=c(3,6),mar=c(2,2,2,2))

for (pt_idx in c(1:17)){

print(pt_idx)
mat_point = mat[mat$PointNum==pt_idx,]	

mat_label1 = mat_point[mat_point$cell_type_label==label1,]
mat_label2 = mat_point[mat_point$cell_type_label==label2,]
print(nrow(mat_label1))
print(nrow(mat_label2))

if (nrow(mat_label1)<20 | nrow(mat_label2)<20){
	print("")
	print("")
	label2_mean_dists = rep(0,10)
	label2_mean_dists_shuf = rep(0,10)
	boxplot(list(label2_mean_dists, label2_mean_dists_shuf),col=c("cornflowerblue","firebrick3"),ylab="Mean distance",main=paste(pt_idx,"NA","NA",sep="_"),names=c(nrow(mat_label1),nrow(mat_label2)))
	next
}

# Shuffle (select shuffled set of non-tumor cells for label1)
shuffle_sel = (mat_point$cell_type_label %in% non_tumor_types) & !(mat_point$cell_type_label %in% c(label1,label2))
mat_point_nont = mat_point[shuffle_sel,]

coord2_dist_mat = matrix(,nrow=nrow(mat_label1),ncol=nrow(mat_label2))
coord1 = cbind(mat_label1$x_cent,mat_label1$y_cent) # every row is a cell in label1
coord2_shuf_dist_array = array(dim=c(nrow(mat_label1),nrow(mat_label2), num_perm))

for (idx2 in c(1:nrow(mat_label2))){
	coord2 = c(mat_label2$x_cent[idx2], mat_label2$y_cent[idx2]) 
	coord2_dist_mat[,idx2] = pointDistance(coord1,coord2,lonlat=F)
}

for (idx3 in c(1:num_perm)){
	mat_label1_shuf = mat_point_nont[sample(c(1:nrow(mat_point_nont)),nrow(mat_label1)),]
	coord1_shuf = cbind(mat_label1_shuf$x_cent, mat_label1_shuf$y_cent)
	for (idx2 in c(1:nrow(mat_label2))){	
		coord2 = c(mat_label2$x_cent[idx2], mat_label2$y_cent[idx2]) 
		coord2_shuf_dist_array[,idx2,idx3] = pointDistance(coord1_shuf,coord2,lonlat=F)
	}
}

label2_mean_dists = rowMeans(coord2_dist_mat)
label2_mean_dists_median = median(label2_mean_dists) # observed median of mean dists
label2_shuf_medians = rep(0,num_perm) # store shuffled medians of mean dists
# Loop through all shuffle label2 mean distance sets and score as higher or lower median vs observed
fdr_count = c(0,0) # higher, lower shuffle median vs. observed median
for (idx4 in c(1: num_perm)){
	label2_mean_dists_shuf = rowMeans(coord2_shuf_dist_array[,,idx4])
	label2_shuf_medians[idx4] = median(label2_mean_dists_shuf)
	if (label2_shuf_medians[idx4]>label2_mean_dists_median){
		fdr_count[1] = fdr_count[1]+1
	}
	else{
		fdr_count[2] = fdr_count[2]+1
	}
}

delta = round(log2(median(label2_mean_dists))-log2(median(label2_shuf_medians)),3)
fdr = min(fdr_count)
par(cex.main=0.7)
#boxplot(list(label2_mean_dists, label2_shuf_medians),col=c("cornflowerblue","firebrick3"),ylab="Mean distance",main=paste(pt_idx,fdr,delta,sep="_"),names=c(nrow(mat_label1),nrow(mat_label2)))
print("")

d = density(label2_mean_dists)
plot(d,lwd=1,main=paste(pt_idx,fdr,delta,sep="_"))
polygon(d,col=alpha("cornflowerblue",0.2))
#abline(v=mean(label2_shuf_medians),col="firebrick2",lty=2,lwd=2)
d2 = density(label2_shuf_medians)
lines(d2,lwd=1)
polygon(d2,col=alpha("firebrick2",0.2))
}

### ** End of loop *** ###


# Heatmap summary of co-localization calls for CD8, CD4, NK, Macrophage
# 1 = co-localized, 0 = NA, -1 = NOT co-localized, -2 = excluded
# flags from FDR cutoff: 1k permutations, FDR < 0.01 -> FDR < 10/1000
colocal_cd8 = c(0,0,1,-1,-1,-1,-1,0,0,-1,0,-1,-1,1,1,1,1)
colocal_cd4 = c(0,0,1,1,-1,-1,-1,1,0,-1,0,-1,-1,1,0,-1,-1)
colocal_mac = c(0,0,1,-1,-1,-1,1,1,0,1,0,-1,-1,-1,-1,-1,1)
colocal_nk = c(0,0,-1,1,-2,-1,0,-2,0,-1,0,-1,-2,-2,1,-2,-2)

colocal_flags = rbind(colocal_cd8, colocal_cd4,colocal_mac,colocal_nk)

colors = c(seq(-2,1,length=5))
my_palette <- rev(colorRampPalette(brewer.pal(6,"RdBu"))(n = 4))
fov_names = c("1_T21","2_T21","3_T26","4_T26","5_T27","6_T27","7_T23","8_T23","9_T26","10_T26","11_T21","12_T27","13_T1","14_T1","15_T1","16_T10","17_T10")

reorder = order(colocal_flags[1,],decreasing=T)

hm2_call = heatmap.2(colocal_flags[,reorder],col=my_palette,breaks=colors,density.info="none",trace="none",Rowv=F,Colv=F,dendrogram="none",symm=F,labRow= c("CD8","CD4","Mac","NK"),labCol=fov_names[reorder],margins=c(3,3),scale="none",cexRow=0.6,cexCol=0.6,rowsep=c(0:5),colsep=c(0:18),sepcolor="black",sepwidth=c(0.0001,0.0001))


#######################
#######################


# What are the features of FOVs exhibiting co-localization of cell types?

# confused on 7 and 10 (inconsistent KS p.val, but generally not significant co-localizing)
fov_set1 = c(3,14,15,16,17)  # CD8 co-localized with Treg (n=5)
fov_set2 = c(4,5,6,7,10,12,13) # CD8 NOT co-localized with Treg (n=7)
fov_set3 = c(1,2,8,9,11) # Not called (n=5)

cell_type_count
cbind(c(1:9),cell_types)
cell_type_count_log = log10(cell_type_count+1)

cell_type_prop = apply(cell_type_count,2,function(x) x/sum(x))

ratio1 = cell_type_count[6,]/(cell_type_count[7,]+1) # 6 = CD8, 7 = Tregs
ratio1 = cell_type_prop[6,]

box_list = list(ratio1[fov_set1],ratio1[fov_set2],ratio1[fov_set3])

# Box plot of ratios divided by co-localized vs not vs NA
par(cex.axis=0.7)
boxplot(box_list,ylab=c("CD8:Treg ratio"),names=c("Yes","No","NA"),outline=T)
stripchart(box_list,vertical=T, add=T, method="jitter",pch=20,cex=2,col=c("cornflowerblue","firebrick3","darkgray"))
abline(h=1,lty=2,lwd=3,col="forestgreen")


# Bar plot of ratios divided by co-localized vs not vs NA
vals = unlist(box_list)
col = c(rep("darkorange2",5),rep("Forestgreen",7),rep("gray",5))
names = c(fov_set1,fov_set2,fov_set3)
reorder_notNA = order(vals[1:12],decreasing=T)
vals_reorder = c(vals[reorder_notNA],vals[13:17])
col_reorder = c(col[reorder_notNA],col[13:17])
names_reorder = c(names[reorder_notNA], names[13:17])
barplot(vals_reorder,col= col_reorder,names= names_reorder,ylab="CD8 proportion")


plot(cell_type_count_log[6,fov_set1], cell_type_count_log[7,fov_set1],col="blue",xlim=c(0,3),ylim=c(0,3),xlab=c("CD8 count in FOV"),ylab=c("Treg count in FOV"),pch=19,main="Cell type counts and co-localizatoin per FOV")
points(cell_type_count_log[6,fov_set2], cell_type_count_log[7,fov_set2],col="red",pch=19)
points(cell_type_count_log[6,fov_set3], cell_type_count_log[7,fov_set3],col="darkgray",pch=19)
abline(0,1,lty=2,col="black")
legend("bottomright",c("Co-localized","NOT co-localized","Not called"),pch=19,col=c("blue","red","darkgray"),cex=0.6)


plot(cell_type_prop[6,fov_set1], cell_type_prop[7,fov_set1],col="blue",xlim=c(0,0.3),ylim=c(0,0.3),xlab=c("CD8 prop in FOV"),ylab=c("Treg prop in FOV"),pch=19,main="Cell type counts and co-localizatoin per FOV")
points(cell_type_prop[6,fov_set2], cell_type_prop[7,fov_set2],col="red",pch=19)
points(cell_type_prop[6,fov_set3], cell_type_prop[7,fov_set3],col="darkgray",pch=19)
abline(0,1,lty=2,col="black")
legend("bottomright",c("Co-localized","NOT co-localized","Not called"),pch=19,col=c("blue","red","darkgray"),cex=0.6)


######################
######################


# Are non-tumor cells preferentially located in parenchyma, tumor-stroma border, or stroma?
# For each non-tumor cell, determine the proportion of nearest neighbors that are epithelial (tumor) cells
# Stratify non-tumor cells by proportion of tumor cell neighbors: high (infiltrated into tumor parenchyma), medium (tumor-stroma border), and low (stroma)


tumor_cell_types = c("Keratinocyte","LeadingEdge","CyclingK")
#tumor_cell_types = c("Keratinocyte","CyclingK","Epithelial") # T1 tiled mat

mat_high = matrix(,nrow=0,ncol=ncol(mat))
mat_med = matrix(,nrow=0,ncol=ncol(mat))
mat_low = matrix(,nrow=0,ncol=ncol(mat))

for (pt_idx in c(2,4,5,6,7,8,10,12,13,14,15)){
	print(pt_idx)
	mat_point = mat[mat$PointNum==pt_idx,]
	sel_nontumor = !(mat_point$cell_type_label %in% tumor_cell_types)  
	kc_neighbor_prop = c()
	for (idx1 in which(sel_nontumor)){
		coord1 = c(mat_point$x_cent[idx1], mat_point$y_cent[idx1])
		coord2 = cbind(mat_point$x_cent,mat_point$y_cent)
		dists = pointDistance(coord1,coord2,lonlat=F)
		neighbors = order(dists)[1:30]
		kc_neighbor_prop = c(kc_neighbor_prop,length(which(mat_point$cell_type_label[neighbors]%in% tumor_cell_types)))
	}
	mat_high = rbind(mat_high,mat_point[sel_nontumor,][kc_neighbor_prop>19,])
	mat_med = rbind(mat_med,mat_point[sel_nontumor,][(kc_neighbor_prop>4) & (kc_neighbor_prop<14),])
	mat_low = rbind(mat_low,mat_point[sel_nontumor,][kc_neighbor_prop<1,])
}


### END OF LOOP


# FOV scatter images of parenchyma vs stroma non-tumor cells
par(mfrow=c(3,6),mar=c(2,2,2,2))
for (pt_idx in c(1:17)){
mat_point = mat[mat$PointNum==pt_idx,]
plot(mat_point$x_cent,2000-mat_point$y_cent,cex=0.1,pch=19,col="black",xlab="X (pixels)",ylab="Y (pixels)",main=pt_idx)
points(mat_high$x_cent[mat_high $PointNum==pt_idx],2000-mat_high$y_cent[mat_high $PointNum==pt_idx],cex=0.7,pch=19,col="red")
points(mat_med$x_cent[mat_med $PointNum==pt_idx],2000-mat_med$y_cent[mat_med $PointNum==pt_idx],cex=0.7,pch=19,col="orange")
points(mat_low$x_cent[mat_low $PointNum==pt_idx],2000-mat_low$y_cent[mat_low $PointNum==pt_idx],cex=0.7,pch=19,col="blue")

}

plot(mat_point$x_cent,2000-mat_point$y_cent,cex=0.1,pch=19,col="white",xlab="X (pixels)",ylab="Y (pixels)")
legend("topleft",c("Parenchyma","Border","Stroma"),fill=c("red","orange","blue"),cex=0.6)

cell_types = c("Fibroblast","Neutrophil","CD4","B","Macrophage","CD8","Tregs","Endothelial")


# 1.) Analyze each FOV separately to determine bias in loalization for each cell type
# Calculate proportion of each non-tumor cell type in high/med/low in each FOV

mat_enrich_p = matrix(,nrow=11,ncol=length(cell_types)) # rows = FOVs, cols = cell types
mat_enrich_b = matrix(,nrow=11,ncol=length(cell_types))
mat_enrich_s = matrix(,nrow=11,ncol=length(cell_types))
mat_prop_p = matrix(,nrow=11,ncol=length(cell_types)) # rows = FOVs, cols = cell types
mat_prop_b = matrix(,nrow=11,ncol=length(cell_types))
mat_prop_s = matrix(,nrow=11,ncol=length(cell_types))

idx_step = 1
for (fov_idx in c(2,4,5,6,7,8,10,12,13,14,15)){
cell_type_count = matrix(0,nrow=length(cell_types),ncol=3)
fov_sel1 = mat_high$PointNum==fov_idx
fov_sel2 = mat_med$PointNum==fov_idx
fov_sel3 = mat_low$PointNum==fov_idx

for (idx1 in c(1:length(cell_types))){
	cell_type_count[idx1,1] = length(which(mat_high$cell_type_label==cell_types[idx1] & fov_sel1))
	cell_type_count[idx1,2] = length(which(mat_med$cell_type_label==cell_types[idx1] & fov_sel2))
	cell_type_count[idx1,3] = length(which(mat_low$cell_type_label==cell_types[idx1] & fov_sel3))
}

cell_type_prop = apply(cell_type_count,2,function(x) x/(sum(x)+1))
mat_prop_p[idx_step,] = cell_type_prop[,1]
mat_prop_b[idx_step,] = cell_type_prop[,2]
mat_prop_s[idx_step,] = cell_type_prop[,3]

cell_type_enrich = apply(cell_type_prop,1,function(x) x-mean(x))
mat_enrich_p[idx_step,] = cell_type_enrich[1,]
mat_enrich_b[idx_step,] = cell_type_enrich[2,]
mat_enrich_s[idx_step,] = cell_type_enrich[3,]

idx_step = idx_step+1
}


par(mfrow=c(3,1),mar=c(2,2,2,2),cex.axis=0.7)

reorder = order(apply(mat_enrich_p,2,function(x) median(x)),decreasing=T)
colors = my_palette <- colorRampPalette(brewer.pal(6,"Set2"))(n =length(cell_types))

boxplot(mat_enrich_p[,reorder],names=cell_types[reorder],cex.names=0.2,col= colors,outline=F,main="Parenchyma enrichment")
abline(h=0,lty=2,col="firebrick2")

boxplot(mat_enrich_b[,reorder],names=cell_types[reorder],cex.names=0.2,col= colors,outline=F,main="Border enrichment")
abline(h=0,lty=2,col="firebrick2")

boxplot(mat_enrich_s[,reorder],names=cell_types[reorder],cex.names=0.2,col= colors,outline=F,main="Stroma enrichment")
abline(h=0,lty=2,col="firebrick2")


# 2.) Analyze all FOVs together to determine bias in loalization for each cell type
cell_type_count = matrix(0,nrow=length(cell_types),ncol=3)
for (idx1 in c(1:length(cell_types))){
	cell_type_count[idx1,1] = length(which(mat_high$cell_type_label==cell_types[idx1]))
	cell_type_count[idx1,2] = length(which(mat_med$cell_type_label==cell_types[idx1]))
	cell_type_count[idx1,3] = length(which(mat_low$cell_type_label==cell_types[idx1]))
}

cell_type_prop = apply(cell_type_count,2,function(x) x/sum(x))


# Heatmap of cell type proportions in parenchyma and stroma
colors = c(seq(0,0.4,length=22))
my_palette <- plasma(21)
hm2_call = heatmap.2(t(cell_type_prop),col=my_palette,breaks=colors,density.info="none",trace="none",Rowv=F,Colv=T,dendrogram="col",symm=F,labRow= c("Parenchyma","Border","Stroma"),labCol=cell_types,margins=c(4,6),scale="none",cexRow=0.6,cexCol=0.6,colsep=c(0:10),rowsep=c(0:4),sepcolor="black",sepwidth=c(0.001,0.001),cellnote=round(t(cell_type_prop),3))


colors = c(seq(-1.15,1.15,length=11))
my_palette <- rev(colorRampPalette(brewer.pal(6,"RdBu"))(n = 10))
hm2_call = heatmap.2(t(cell_type_prop),col=my_palette,breaks=colors,density.info="none",trace="none",Rowv=F,Colv=T,dendrogram="none",symm=F,labRow= c("Parenchyma","Border","Stroma"),labCol=cell_types,margins=c(4,6),scale="col",cexRow=0.6,cexCol=0.6,colsep=c(0:10),rowsep=c(0:4),sepcolor="black",sepwidth=c(0.001,0.001),cellnote=round(t(cell_type_prop),3),notecol="black")

######################
######################
