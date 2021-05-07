### Author: Tyson Dawson
This Script takes the counts from PRJNA294187 and conducts differential expression analysis
library(DESeq2)
library(gplots)
library(RColorBrewer)
library(genefilter)
library(biomaRt)
library(Biobase)
library(limma)
library(edgeR)
library(HTSFilter)
library(ggplot2)
library(BiocParallel)
library(GEOquery)
library(stringr)
library(pcaExplorer)
library(airway)
library(dplyr)
library(RCurl)

### Gather all counts files
setwd("C:/Users/tyson/OneDrive/Desktop/Lupus/PRJNA294187/")
myTable=read.table("PRJNA294187-counts.txt",header=TRUE)
counts.data <- mutate_all(myTable[,c(7:length(colnames(myTable)))], function(x) as.numeric(as.character(x)))
rownames(counts.data)=myTable$Geneid
#counts.data=counts.data[,-c(46,51,75)]
### Cross reference the counts data with sample meta data
filenames <- download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE72nnn/GSE72509/matrix/GSE72509_series_matrix.txt.gz", destfile="metadata.txt")
counts.data <- counts.data[rowSums(counts.data) > 1, ] #Filter out those features that are only zeros
gseFile=getGEOfile(GEO="GSE72509", destdir = ".")
gse=getGEO(filename=gseFile)
DISEASE=c()
index=1
listOfColNames=c()
for(x in names(gse@gsms)){
  if(index==1){
    for(q in as.character(names(gse@gsms[[x]]@header))){
      listOfColNames=append(listOfColNames,q)
      if(length(gse@gsms[[x]]@header[[q]])>1){
        for(y in 2:length(gse@gsms[[x]]@header[[q]])){
          listOfColNames=append(listOfColNames,paste(q,y,sep="__"))
        }
      }
    }
    mydf=as.data.frame(t(listOfColNames))
    colnames(mydf)=mydf[1,]
  }
  if(str_count(gse@gsms[[x]]@header$characteristics_ch1[1], "lupus")>0){
    DISEASE=append(DISEASE, "lupus")
  }else{
    DISEASE=append(DISEASE, "control")
  }
  #mydf=as.data.frame(rbind(mydf, as.character(gse@gsms[[x]]@header)))
  myCount=1
  for(y in gse@gsms[[x]]@header){
    i=as.character(names(gse@gsms[[x]]@header))[myCount]
    if(myCount==1){
      newDataFrame=data.frame()
    }
    for(w in 1:length(y)){
      thisi=i
      if(w>1){
        thisi=paste(i,w,sep="__")
      }
      newDataFrame[1,thisi]=y[w]
    }
    last2=substr(i,nchar(i)-2,nchar(i)-1)
    myCount=myCount+1
  }
  index=index+1
  mydf=suppressMessages(full_join(mydf, newDataFrame))
}
mydf=mydf[2:length(row.names(mydf)),]
#write.csv(mydf, "Metadata.csv")
metastuff=data.frame(DISEASE)
metastuff$DISEASE <- as.factor(DISEASE)
metastuff$DISEASE=relevel(metastuff$DISEASE, ref = "control")
ismStatus=data.frame(mydf$characteristics_ch1__4)
ismStatus$mydf.characteristics_ch1__4=as.factor(mydf$characteristics_ch1__4)
ismStatus$mydf.characteristics_ch1__4=relevel(ismStatus$mydf.characteristics_ch1__4, ref="ism: control")
counts.data[is.na(counts.data)] <- 0
dds <- DESeqDataSetFromMatrix(countData = counts.data,
                              colData = ismStatus,
                              design = ~ mydf.characteristics_ch1__4)

### Filter raw reads with HTSFilter to eliminate uninformative and repetitive transcripts by Jaccard Index
idx <- as.character(dds$DISEASE) #Make a vector of WHO variables for HTS filter
pdf("Roche-HTS-Filter_Plot.pdf", width = 16, height = 9)
filtered <- HTSFilter(dds,
                      conds = idx,
                      pAdjustMethod = "BH",
                      s.len = 100,
                      s.min = 1,
                      s.max = 200,
                      normalization = "DESeq",
                      plot = T)
dev.off()
filtered <- filtered$filteredData #Get the data frame of the filtered transcript reads
dim(filtered) #11756 29 remaining features

dds <- DESeq(filtered) #Normalizes via DESeq2 method; already done with HTSFilter, but results will be the same
pdf(file="Dispersion_Estimates.pdf",height = 9,width = 16)
plotDispEsts( dds, ylim = c(1e-8, 1e1) ) #Check dispersion estimates 
dev.off()
rownames(dds) <- substr(rownames(dds),1,15)
#annotation <- gene.list[match(rownames(dds), gene.list$ensembl),]
#all(rownames(dds) == annotation$ensembl)
#mcols(dds) <- cbind(mcols(dds),annotation)
### Differential Analyses
## Compare the SLE versus CTL
res<- results(dds,
              contrast=c("DISEASE","lupus","control"),
              alpha = 0.05,
              pAdjustMethod = "BH",
              independentFiltering=FALSE,  #Since HTSFilter was used, no need for DESeq to filter
              #BPPARAM = SnowParam(4),
              #parallel = T
)
res <- res[order(res$padj),]
index <- which(res$padj<0.2)
res.sig <- res[index,]
res.sig <- as.data.frame(res.sig)
res.sig <- res.sig[ order(res.sig$padj), ] 


gene.list <- as.data.frame(rownames(res))
colnames(gene.list) <- "ensembl"
gene.list <- substr(gene.list$ensembl,1,15)
gene.list <- as.data.frame(gene.list)
colnames(gene.list) <- c("ensembl")

ensembl <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl",
                      mirror="www")
genemap <- getBM(attributes = c("ensembl_gene_id",
                                "entrezgene_id",
                                "hgnc_symbol",
                                "chromosome_name",
                                "start_position",
                                "end_position"),
                 filters = "ensembl_gene_id",
                 values = gene.list,#$ensembl,
                 mart = ensembl)
idx <- match(gene.list$ensembl, genemap$ensembl_gene_id)
gene.list$entrez <- genemap$entrezgene_id[idx]
gene.list$symbol <- genemap$hgnc_symbol[idx]
gene.list$chromosome <- genemap$chromosome_name[idx]
gene.list$start <- genemap$start_position[idx]
gene.list$end <- genemap$end_position[idx]
res=res[order(row.names(res)),]
gene.list=gene.list[order(gene.list$ensembl),]
res$entrez=gene.list$entrez
res$symbol=gene.list$symbol
res$chromosome=gene.list$chromosome
res$start=gene.list$start
res$end=gene.list$end

write.csv(as.data.frame(res), file = "Res.csv")

# ### Check MA plots of differential data
pdf(file="EA-SLE-CTL_MA_plots.pdf",height = 9,width = 16)
DESeq2::plotMA( dds, ylim = c(-10, 10), main = "All genes")
dev.off() 
# ### Data Transformations
# Log2 transformation
log2 <- normTransform(dds)
log2.counts <- assay(log2)
log2.counts <- as.data.frame(log2.counts)
write.csv(log2.counts, "Log2Counts.csv")
## Alternatively:
# glom.log2.counts <- log2(1+counts(ddssva.glom,normalized=TRUE))
# Regularized log transformation
rld <- rlogTransformation(dds, blind=FALSE)
rld.counts <- assay(rld)
rld.counts <- as.data.frame(rld.counts)
# Variance stablizing transformation
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
vsd.counts <- assay(vsd)
vsd.counts <- as.data.frame(vsd.counts) 
# ### Compare Data Transformation Methods
pdf(file="Transformation_Comparison_1.pdf",height=9,width=16)
par(mfrow=c(1,3))
plot(assay(log2)[,1:4],col="#00000020",pch=20,cex=0.3,main="Log2 Transformed")
plot(assay(rld)[,1:4],col="#00000020",pch=20,cex=0.3,main="Regularized Logarithm (rLog)")
plot(assay(vsd)[,1:4],col="#00000020",pch=20,cex=0.3,main="Variance Stabilization Tranformation (VST)")
dev.off() 
library(gplots)
library(RColorBrewer)
library(genefilter)
library(ggplot2)
library(gplots)
library(ggrepel)
## Read in the data
genes <- res.sig
genes$symbol=row.names(res.sig)
pdf(file="LabelledDifferentiallyExpressed.pdf",height = 12,width = 12) 
sigUp <- which(genes$log2FoldChange>0 & genes$padj<4*(10^-14))
sigDown <- which(genes$log2FoldChange<0 & genes$padj<4*(10^-14))
nonSigUp <- which(genes$log2FoldChange>0 & genes$padj>4*(10^-14))
nonSigDown <- which(genes$log2FoldChange<0 & genes$padj>4*(10^-14))
genes$Significance <- "FDR<4e-14, Upregulated"
genes[ sigUp, "Significance"] <- "FDR<4e-14, Upregulated"
genes[ sigDown, "Significance"] <- "FDR<4e-14, Downregulated"
genes[ nonSigUp, "Significance"] <- "Less Significant, Upregulated"
genes[ nonSigDown, "Significance"] <- "Less Significant, Downregulated" 
ggplot(genes, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(aes(color = Significance)) +
  scale_color_manual(values = c("green1", "red", "#b8d9b4",  "#d9b4b4")) +
  xlab(expression('Log'[2]*'Fold Change in Expression Value')) +
  ylab(expression('-Log'[10]*'P-Value')) +
  ggtitle("Differentially Expressed HERVs in Alzheimers") +
  theme(text = element_text(size = 40),
        legend.position = "right", plot.title = element_text(hjust = 0.5, size=40),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.key.height=unit(2, "cm")) +
  geom_text_repel(
    data = subset(genes, padj < 1.5*(10^-3)),
    aes(label = symbol),
    size = 4,
    box.padding = unit(0.4, "lines"),
    point.padding = unit(0.3, "lines")
  )
ggsave("VolcanoPlot.png", plot = last_plot(),
       height = 10, width=15, device = "png", dpi = 500, limitsize = FALSE)
write.csv(res.sig, "AlzheimersHervs-PRJEB28518.csv", quote=FALSE) 
