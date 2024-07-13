## This script loads, analyzes, and plots RNAseq differential expression

library(readxl) #loads library for reading excel files
library(tidyverse) #loads tidyverse library
library(DESeq2) #loads DESeq2
library(pheatmap)
library(dplyr)


excelsheet = "G:/My Drive/My Papers/Barbara_2024/BarbaraRNAseq_AH240526.xlsx"

######### Analysis of all data together

#import dataset for all genotypes
genecounts <- read_excel(excelsheet, sheet =1, col_names = TRUE) #loads data from excel
genecounts <- as.data.frame(genecounts) #turns tibble output from above into data.frame
genecounts <- genecounts[!duplicated(genecounts[,1]),] #need to remove duplicate rows
rownames(genecounts) <- genecounts[,1] #makes first column into rownames
genecounts$Gene_Symbol <- NULL #removes the Gene_name column (rownames have this info)
genecounts <- genecounts %>% mutate_if(is.numeric,round) #round values so they can be imported into deseq2. I don't love this solution- I'd rather have raw counts

#make DESeqDataSet, needs matrix of conditions to organize comparisons
condition <- read_excel(excelsheet, sheet =2, col_names = TRUE) #loads data from excel
condition <-as.data.frame(condition)
rownames(condition) <- condition[,1] #makes first column into rownames
condition$Library <- NULL #removes the Gene_name column (rownames have this info)

dds <- DESeqDataSetFromMatrix(countData = genecounts, colData = condition, design = ~Cell_Type+Genotype)

rlogcounts = rlog(dds,blind=FALSE) #DESeq2 normalization that log2 transforms and scales for library size
rlogcountstable = assay(rlogcounts)

#PCA
pcaData<-plotPCA(rlogcounts, intgroup=c('Genotype','Cell_Type'),returnData=TRUE) #plotting using plotPCA
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1,PC2, color=Genotype, shape=Cell_Type)) + 
  geom_point(size=4) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

## plot variable genes
countVar <- apply(rlogcountstable, 1, var) #get variance of each row
highVar <- order(countVar, decreasing=TRUE)[1:500] #find top 150 most variable genes
hmDat <- rlogcountstable[highVar,] #pull top 250 most variable genes

pheatmap(hmDat, cluster_rows=TRUE, show_rownames = FALSE, cluster_cols=TRUE, scale = "row")


############ Analysis of just CNs
genecounts <- read_excel(excelsheet, sheet =3, col_names = TRUE) #loads data from excel
genecounts <- as.data.frame(genecounts) #turns tibble output from above into data.frame
genecounts <- genecounts[!duplicated(genecounts[,1]),] #need to remove duplicate rows
rownames(genecounts) <- genecounts[,1] #makes first column into rownames
genecounts$Gene_Symbol <- NULL #removes the Gene_name column (rownames have this info)
genecounts <- genecounts %>% mutate_if(is.numeric,round) #round values so they can be imported into deseq2. I don't love this solution- I'd rather have raw counts

#make DESeqDataSet, needs matrix of conditions to organize comparisons
condition <- read_excel(excelsheet, sheet =4, col_names = TRUE) #loads data from excel
condition <-as.data.frame(condition)
rownames(condition) <- condition[,1] #makes first column into rownames
condition$Library <- NULL #removes the Gene_name column (rownames have this info)

dds <- DESeqDataSetFromMatrix(countData = genecounts, colData = condition, design = ~Batch+Genotype)

rlogcounts = rlog(dds,blind=FALSE) #DESeq2 normalization that log2 transforms and scales for library size
rlogcountstable = assay(rlogcounts)

#PCA
pcaData<-plotPCA(rlogcounts, intgroup=c('Genotype','Batch'),returnData=TRUE) #plotting using plotPCA
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1,PC2, color=Genotype, shape=Batch)) + 
  geom_point(size=4) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

## plot variable genes
countVar <- apply(rlogcountstable, 1, var) #get variance of each row
highVar <- order(countVar, decreasing=TRUE)[1:500] #find top 150 most variable genes
hmDat <- rlogcountstable[highVar,] #pull top 250 most variable genes

pheatmap(hmDat, cluster_rows=TRUE, show_rownames = FALSE, cluster_cols=TRUE, scale = "row")


############ Analysis of just MNs
genecounts <- read_excel(excelsheet, sheet =5, col_names = TRUE) #loads data from excel
genecounts <- as.data.frame(genecounts) #turns tibble output from above into data.frame
genecounts <- genecounts[!duplicated(genecounts[,1]),] #need to remove duplicate rows
rownames(genecounts) <- genecounts[,1] #makes first column into rownames
genecounts$Gene_Symbol <- NULL #removes the Gene_name column (rownames have this info)
genecounts <- genecounts %>% mutate_if(is.numeric,round) #round values so they can be imported into deseq2. I don't love this solution- I'd rather have raw counts

#make DESeqDataSet, needs matrix of conditions to organize comparisons
condition <- read_excel(excelsheet, sheet =6, col_names = TRUE) #loads data from excel
condition <-as.data.frame(condition)
rownames(condition) <- condition[,1] #makes first column into rownames
condition$Library <- NULL #removes the Gene_name column (rownames have this info)

dds <- DESeqDataSetFromMatrix(countData = genecounts, colData = condition, design = ~Genotype)

rlogcounts = rlog(dds,blind=FALSE) #DESeq2 normalization that log2 transforms and scales for library size
rlogcountstable = assay(rlogcounts)

#PCA
pcaData<-plotPCA(rlogcounts, intgroup=c('Genotype'),returnData=TRUE) #plotting using plotPCA
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1,PC2, color=Genotype)) + 
  geom_point(size=4) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

## plot variable genes
countVar <- apply(rlogcountstable, 1, var) #get variance of each row
highVar <- order(countVar, decreasing=TRUE)[1:500] #find top 150 most variable genes
hmDat <- rlogcountstable[highVar,] #pull top 250 most variable genes

pheatmap(hmDat, cluster_rows=TRUE, show_rownames = FALSE, cluster_cols=TRUE, scale = "row")

########### GSEA on Barbara's LMNs

library(fgsea)
library(data.table)
library(ggplot2)
library(tidyverse)
library(readxl) #loads library for reading excel files


#import pathways of interest
pathways.barbara<- gmtPathways("G:/My Drive/My Papers/Barbara_2024/GO/Barbara_genelist1.gmt") #these are relevant GOs from above

#import fold change for Barbara LMN
excelsheet = "G:/My Drive/My Papers/Barbara_2024/BarbaraRNAseq_AH240526.xlsx"
genes <- read_excel(excelsheet, sheet =7, col_names = TRUE) #loads data from excel
genes<-as.data.frame(genes)

subgroup<- select(genes,All_Genes,Barbara_log2FoldChange) #here I've changed it to look at Barbara's LMNs column
values = subgroup$Barbara_log2FoldChange #for some reason, NA values aren't recognized as numeric. The next few lines circumvent this
values<-as.numeric(values)
subgroup$Barbara_log2FoldChange<-values
subgroup<-na.omit(subgroup)
ranks <-deframe(subgroup)

fgseaRes <- fgsea(pathways=pathways.barbara, stats=ranks, eps = 0.0) #analysis of Barbara's pathways

fgseaout <- apply(fgseaRes,2,as.character)
write.csv(fgseaout,file='G:/My Drive/My Papers/Barbara_2024/GO/Barbara_LMN_fgsea.csv')

#make GSEA plots
plotEnrichment(pathways.barbara[["Axonal_transport"]],ranks) + labs(title="Axonal_transport") #just change the strings to plot different groups (see below)

plotEnrichment(pathways.barbara[["DNA Damage from Christine"]],ranks) + labs(title="DNA Damage Sensor") #just change the strings to plot different groups (see below)

plotEnrichment(pathways.barbara[["DNA_Damage_response"]],ranks) + labs(title="DNA Damage Response") #just change the strings to plot different groups (see below)

plotEnrichment(pathways.barbara[["positive_cell_cycle_regulation"]],ranks) + labs(title="Positive Cell Cycle Regulation") #just change the strings to plot different groups (see below)

plotEnrichment(pathways.barbara[["Random"]],ranks) + labs(title="Random 500 Genes") #just change the strings to plot different groups (see below)

plotEnrichment(pathways.barbara[["Type1_interferon_GO0060337"]],ranks) + labs(title="Type1_interferon") #just change the strings to plot different groups (see below)

plotEnrichment(pathways.barbara[["Innate_Immune_GO0045087"]],ranks) + labs(title="Innate_Immune") #just change the strings to plot different groups (see below)

############ GSEA on Barbara's CSMNs
#import pathways of interest
pathways.barbara<- gmtPathways("G:/My Drive/My Papers/Barbara_2024/GO/Barbara_genelist1.gmt") #these are relevant GOs from above

#import fold change for Barbara LMN
excelsheet = "G:/My Drive/My Papers/Barbara_2024/BarbaraRNAseq_AH240526.xlsx"
genes <- read_excel(excelsheet, sheet =9, col_names = TRUE) #loads data from excel
genes<-as.data.frame(genes)

subgroup<- select(genes,All_Genes,Barbara_log2FoldChange) #here I've changed it to look at Barbara's LMNs column
values = subgroup$Barbara_log2FoldChange #for some reason, NA values aren't recognized as numeric. The next few lines circumvent this
values<-as.numeric(values)
subgroup$Barbara_log2FoldChange<-values
subgroup<-na.omit(subgroup)
ranks <-deframe(subgroup)

fgseaRes <- fgsea(pathways=pathways.barbara, stats=ranks, eps = 0.0) #analysis of Barbara's pathways

fgseaout <- apply(fgseaRes,2,as.character)
write.csv(fgseaout,file='G:/My Drive/My Papers/Barbara_2024/GO/Barbara_CSMN_fgsea.csv')

#make GSEA plots
plotEnrichment(pathways.barbara[["Axonal_transport GO0098930"]],ranks) + labs(title="Axonal_transport") #just change the strings to plot different groups (see below)

plotEnrichment(pathways.barbara[["DNA Damage sensor go0140612"]],ranks) + labs(title="DNA Damage Sensor") #just change the strings to plot different groups (see below)

plotEnrichment(pathways.barbara[["DNA_Damage_response GO0006974"]],ranks) + labs(title="DNA Damage Response") #just change the strings to plot different groups (see below)

plotEnrichment(pathways.barbara[["positive_cell_cycle_regulation GO0090068"]],ranks) + labs(title="Positive Cell Cycle Regulation") #just change the strings to plot different groups (see below)

plotEnrichment(pathways.barbara[["Random"]],ranks) + labs(title="Random 500 Genes") #just change the strings to plot different groups (see below)

