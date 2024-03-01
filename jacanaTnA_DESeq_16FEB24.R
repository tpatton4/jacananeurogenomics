#DESeq2 on all TnA samples

library(DESeq2)
#read in gene count data
setwd("C:/Users/tmp15/Downloads/R/jacana")
counts_genes <- as.matrix(read.csv("jacanaTnA_counts_8FEB24.csv",row.names="Gene.stable.ID"))
counts <- subset(counts_genes,select=-c(Gene.name)) #remove gene names column for DESeq steps
counts <- as.matrix(counts)
mode(counts) <- "integer"

treat <- read.csv("jacana_allsamples_sexStageInfo.csv",header=TRUE,row.names=1)
#remove POA info
treat <- treat[!(treat$Brain.Region=="POA"),]
treat <- subset(treat,select=-c(Sex,Stage,Brain.Region))

all(rownames(treat) %in% colnames(counts))
#true - all sample IDs are in 
all(rownames(treat) == colnames(counts))
#true - all sample IDs are in the correct order

#make a DESeq object
dds <- DESeqDataSetFromMatrix(countData=counts, colData=treat, design= ~ SexStage)
dds
yes = rowSums(counts > 10) >= 5
dds2 <- dds[yes,]
#make F the reference
dds2$Sex <-relevel(dds$SexStage, ref="Female")
dds2<-DESeq(dds2)

#Make comparisons between groups
resFvMP <- results(dds2, contrast=c("SexStage", "Female", "MaleP"))
resFvMC <- results(dds2, contrast=c("SexStage", "Female", "MaleC"))
resMCvMP <- results(dds2, contrast=c("SexStage", "MaleC", "MaleP"))

summary(resFvMP)
summary(resFvMC)
summary(resMCvMP)

####write results to files####
allgenesFvMP <- merge(as.matrix(resFvMP), counts_genes, by="row.names") #add gene names and counts
allgenesFvMC <- merge(as.matrix(resFvMC), counts_genes, by="row.names")
allgenesMCvMP <- merge(as.matrix(resMCvMP), counts_genes, by="row.names")
#optionally, write a file with DE data for all genes, not just significantly differentially expressed
#write.csv(allgenesFvMP, "jacanaTnAallgenes_FvMP_DESeq2_16FEB24.csv")
#write.csv(allgenesFvMC, "jacanaTnAallgenes_FvMC_DESeq2_16FEB24.csv")
#write.csv(allgenesMCvMP, "jacanaTnAallgenes_MCvMP_DESeq2_16FEB24.csv")

#write a file with signifiantly differentially expressed genes, as well as their gene names and counts
siggenesFvMP <- as.matrix(subset(resFvMP,padj <0.05)) #adjusted pval less than 0.05
siggenesFvMC <- as.matrix(subset(resFvMC,padj <0.05))
siggenesMCvMP <- as.matrix(subset(resMCvMP,padj <0.05))
sigcountsFvMP <- merge(siggenesFvMP, counts_genes, by="row.names", all=F) #merge the matrices by gene IDs, all=F so we don't have a bunch of NA lines for non-significant genes
sigcountsFvMC <- merge(siggenesFvMC, counts_genes, by="row.names", all=F)
sigcountsMCvMP <- merge(siggenesMCvMP, counts_genes, by="row.names", all=F)
write.csv(sigcountsFvMP, "jacanaTnA_FvMP_DESeq2Out_16FEB24.csv")
write.csv(sigcountsFvMC, "jacanaTnA_FvMC_DESeq2Out_16FEB24.csv")
write.csv(sigcountsMCvMP, "jacanaTnA_MCvMP_DESeq2Out_16FEB24.csv")


####PCA####
vsd <- vst(dds2)
pcaData <- plotPCA(vsd, intgroup="SexStage", returnData=TRUE)
plotPCA(vsd, intgroup="SexStage")
#identify genes that contribute the most to PCs
install.packages("factoextra")
library(factoextra)
vsd_mat <- assay(vsd)
pca <- prcomp(t(vsd_mat), scale=TRUE)
fviz_eig(pca) #scree plot
var_coord_func <- function(loadings, comp.sdev) 
  
####sample heatmap####
install.packages("pheatmap")
library(pheatmap)
rlogData <- rlog(dds2, blind=TRUE) #blind=TRUE for unsupervised clustering
#Selecting the top 50 genes by variance
#topVarGenes <- head(order(rowVars(assay(rlogData)), decreasing = TRUE), 50)
variances <- rowVars(assay(rlogData), useNames=TRUE)
orderedVariances <- order(variances, decreasing = TRUE)
#topVarGenes <- head(orderedVariances, 50)
#mat <- assay(rlogData)[topVarGenes, ]
mat <- assay(rlogData)
rlog_cor <- cor(mat)
head(rlog_cor)
pheatmap(rlog_cor, annotation= treat)

####Volcano Plots####
#FvMP
with(allgenesFvMP, plot(log2FoldChange, -log10(padj), pch=19, main="Females vs Parenting Males TnA", cex=0.7, xlab=bquote(~Log[2]~fold~change), ylab=bquote(~-log[10]~P~value), xlim=c(-3,13), ylim=c(-1,300)))

#Add lines for cut-offs: logFC>0.5 and p-value cut-off at padj<0.05
abline(h=-log10(0.05), col="black", lty=3, lwd=1)
abline(v=-0.5, col="black", lty=3, lwd=1)
abline(v=0.5, col="black", lty=3, lwd=1)

#color significant genes based on whether they were up- or down-regulated
with(subset(allgenesFvMP, padj<0.05 & log2FoldChange< -0.5), points(log2FoldChange, -log10(padj), pch=19, col="skyblue3", cex=0.7))
with(subset(allgenesFvMP, padj<0.05 & log2FoldChange>0.5), points(log2FoldChange, -log10(padj), pch=19, col="tomato3", cex=0.7))

#Add genes names to the significant genes using the code below. Try adjusting how many gene names are shown.
cutoff=sort(allgenesFvMP$padj)[10] #selects the top 10 smallest p values
sign.genes=which(allgenesFvMP$padj <= cutoff)
text(x=allgenesFvMP$log2FoldChange[sign.genes] , y=-log10(allgenesFvMP$padj[sign.genes]), label=allgenesFvMP$Gene.name[sign.genes], cex=0.5)

#FvMC
with(allgenesFvMC, plot(log2FoldChange, -log10(padj), pch=19, main="Females vs Courting Males TnA", cex=0.7, xlab=bquote(~Log[2]~fold~change), ylab=bquote(~-log[10]~P~value), xlim=c(-3,15), ylim=c(-1,300)))

#Add lines for cut-offs: logFC>0.5 and p-value cut-off at padj<0.05
abline(h=-log10(0.05), col="black", lty=3, lwd=1)
abline(v=-0.5, col="black", lty=3, lwd=1)
abline(v=0.5, col="black", lty=3, lwd=1)

#color significant genes based on whether they were up- or down-regulated
with(subset(allgenesFvMC, padj<0.05 & log2FoldChange< -0.5), points(log2FoldChange, -log10(padj), pch=19, col="skyblue3", cex=0.7))
with(subset(allgenesFvMC, padj<0.05 & log2FoldChange>0.5), points(log2FoldChange, -log10(padj), pch=19, col="tomato3", cex=0.7))

#Add genes names to the significant genes using the code below. Try adjusting how many gene names are shown.
cutoff=sort(allgenesFvMC$padj)[10] #selects the top 10 smallest p values
sign.genes=which(allgenesFvMC$padj <= cutoff)
text(x=allgenesFvMC$log2FoldChange[sign.genes] , y=-log10(allgenesFvMC$padj[sign.genes]), label=allgenesFvMC$Gene.name[sign.genes], cex=0.5)

#MCvMP
with(allgenesMCvMP, plot(log2FoldChange, -log10(padj), pch=19, main="Courting vs Parenting Males TnA", cex=0.7, xlab=bquote(~Log[2]~fold~change), ylab=bquote(~-log[10]~P~value)))

#Add lines for cut-offs: logFC>0.5 and p-value cut-off at padj<0.05
abline(h=-log10(0.05), col="black", lty=3, lwd=1)
abline(v=-0.5, col="black", lty=3, lwd=1)
abline(v=0.5, col="black", lty=3, lwd=1)

#color significant genes based on whether they were up- or down-regulated
with(subset(allgenesMCvMP, padj<0.05 & log2FoldChange< -0.5), points(log2FoldChange, -log10(padj), pch=19, col="skyblue3", cex=0.7))
with(subset(allgenesMCvMP, padj<0.05 & log2FoldChange>0.5), points(log2FoldChange, -log10(padj), pch=19, col="tomato3", cex=0.7))

#Add genes names to the significant genes using the code below. Try adjusting how many gene names are shown.
cutoff=sort(allgenesMCvMP$padj)[10] #selects the top 10 smallest p values
sign.genes=which(allgenesMCvMP$padj <= cutoff)
text(x=allgenesMCvMP$log2FoldChange[sign.genes] , y=-log10(allgenesMCvMP$padj[sign.genes]), label=allgenesMCvMP$Gene.name[sign.genes], cex=0.5)