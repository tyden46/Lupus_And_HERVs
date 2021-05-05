#### Author: Tyson Dawson
#### This script takes as input a dataframe with genes and HERVs that are correlated. It fetches the positions of these loci
#### and plots them on a circos plot.

#Load Libraries
library(BioCircos)
library(Rgb)
library(biomaRt)
library(rentrez)
library(stringr)
setwd("C:/Users/tyson/OneDrive/Desktop/Lupus/PRJNA294187")
top100Links=read.table("top100Hervs.tsv",skip=1)
gtf=read.gtf("transcripts.gtf")
col1Chroms=rep("A",100)
col1Start=rep("A",100)
col1End=rep("A",100)
col2Chroms=rep("A",100)
col2Start=rep("A",100)
col2End=rep("A",100)
linksFrame=as.data.frame(cbind(col1Chroms,col1Start,col1End,
                               col2Chroms,col2Start,col2End))

#Col 1 Hervs
col1Hervs=which(!str_detect(top100Links$V1, "ENSG"))
locCol1Hervs=which(gtf$gene_id %in% top100Links$V1[col1Hervs])
linksFrame$col1Chroms[col1Hervs]=gtf$seqname[match(top100Links$V1[col1Hervs],gtf$gene_id)]
linksFrame$col1Start[col1Hervs]=gtf$start[match(top100Links$V1[col1Hervs],gtf$gene_id)]
linksFrame$col1End[col1Hervs]=gtf$end[match(top100Links$V1[col1Hervs],gtf$gene_id)]

#Col 2 HERVs
col2Hervs=which(!str_detect(top100Links$V2, "ENSG"))
locCol2Hervs=which(gtf$gene_id %in% top100Links$V2[col2Hervs])
linksFrame$col2Chroms[col2Hervs]=gtf$seqname[match(top100Links$V2[col2Hervs],gtf$gene_id)]
linksFrame$col2Start[col2Hervs]=gtf$start[match(top100Links$V2[col2Hervs],gtf$gene_id)]
linksFrame$col2End[col2Hervs]=gtf$end[match(top100Links$V2[col2Hervs],gtf$gene_id)]

#Col 1 Genes
gene.list <- as.data.frame(top100Links$V1[which(str_detect(top100Links$V1, "ENSG"))])
colnames(gene.list) <- "ensembl"
gene.list <- substr(gene.list$ensembl,1,15)
gene.list <- as.data.frame(gene.list)
colnames(gene.list) <- c("ensembl")
ensembl <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl",
                      mirror="www")
genemap <- getBM(attributes = c("ensembl_gene_id",
                                "chromosome_name",
                                "hgnc_symbol",
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

linksFrame$col1Chroms[which(str_detect(top100Links$V1, "ENSG"))]=gene.list$chromosome
linksFrame$col1Start[which(str_detect(top100Links$V1, "ENSG"))]=gene.list$start
linksFrame$col1End[which(str_detect(top100Links$V1, "ENSG"))]=gene.list$end


#Col 2 Genes
gene.list <- as.data.frame(top100Links$V2[which(str_detect(top100Links$V2, "ENSG"))])
colnames(gene.list) <- "ensembl"
gene.list <- substr(gene.list$ensembl,1,15)
gene.list <- as.data.frame(gene.list)
colnames(gene.list) <- c("ensembl")
ensembl <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl",
                      mirror="www")
genemap <- getBM(attributes = c("ensembl_gene_id",
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

linksFrame$col2Chroms[which(str_detect(top100Links$V2, "ENSG"))]=gene.list$chromosome
linksFrame$col2Start[which(str_detect(top100Links$V2, "ENSG"))]=gene.list$start
linksFrame$col2End[which(str_detect(top100Links$V2, "ENSG"))]=gene.list$end

linksFrame$col1Chroms=str_remove_all(linksFrame$col1Chroms,"chr")
linksFrame$col2Chroms=str_remove_all(linksFrame$col2Chroms,"chr")



# Make Circos
links_chromosomes_1 = linksFrame$col1Chroms # Chromosomes on which the links should start
links_chromosomes_2 = linksFrame$col2Chroms # Chromosomes on which the links should end

links_pos_1 = as.numeric(linksFrame$col1Start)
links_pos_2 = as.numeric(linksFrame$col2Start)

tracklist = BioCircosBackgroundTrack("myBackgroundTrack", minRadius = 0, maxRadius = 1,
                                     borderSize = 0, fillColors = "#EEFFEE")  

tracklist = tracklist + BioCircosLinkTrack('myLinkTrack', links_chromosomes_1, links_pos_1,
                                           links_pos_1, links_chromosomes_2, links_pos_2, links_pos_2,
                                           maxRadius = 1, width="0.02em")

pdf(file="HervsGenesCircos.pdf")
BioCircos(tracklist, genomeFillColor = "Spectral",
          chrPad = 0.02, displayGenomeBorder = FALSE, yChr =  TRUE,
          genomeTicksDisplay = FALSE,  genomeLabelTextSize = "8pt", genomeLabelDy = 0)
dev.off()
