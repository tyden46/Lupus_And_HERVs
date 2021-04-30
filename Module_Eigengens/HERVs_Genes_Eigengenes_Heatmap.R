#### Written by Tyson Dawson

#### This script takes the summary_eigengene.tsv file from cemitool
#### and removes the least informative modules (in this case)
#### M14 and the Not.Correlated module and plots them in a heatmap.
#### In this case, we have gene+HERV modules from 99 lupus patients
#### from PRJNA294187
library(vroom)
setwd("C:\\Users\\tyson\\OneDrive\\Desktop\\Lupus\\PRJNA294187")
#Read files
eigenGenes=vroom("C:/Users/tyson/OneDrive/Desktop/Lupus/PRJNA294187/cemitool_results/summary_eigengene.tsv")
library(ComplexHeatmap)
library(dendextend)
#Convert Counts to matrix
mat=as.matrix(eigenGenes[,2:length(colnames(eigenGenes))])
#Remove uninformative modules
mat=mat[-which(eigenGenes$modules=="M14"),]
eigenGenes=eigenGenes[-which(eigenGenes$modules=="M14"),]
mat=mat[-which(eigenGenes$modules=="Not.Correlated"),]
eigenGenes=eigenGenes[-which(eigenGenes$modules=="Not.Correlated"),]
#Set rownames
rownames(mat)=eigenGenes$modules
#Make row dendrogram, color two most prominent branches
row_dend = as.dendrogram(hclust(dist(mat)))
row_dend = color_branches(row_dend, k = 2)
#Save
pdf(file="HERV_Genes_Eigengenes-Heatmap.pdf",width=10,height=10)
Heatmap(mat,
        cluster_rows = row_dend)
dev.off()
