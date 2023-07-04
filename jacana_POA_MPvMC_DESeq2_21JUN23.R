##### Courting vs parenting males
##### POA

library(DESeq2)
#sessionInfo()
#read in gene count data
counts <- as.matrix(read.csv("counts_with_egenes_POA.csv",row.names="Gene.stable.ID"))
counts <- subset(counts,select=-c(X,Gene.name,Annotation.Names))

stage1 <- read.csv("jacana_allsamples_sexStageInfo.csv",header=TRUE,row.names=1)
#remove TNA info and Females
stage <- stage1[!(stage1$Brain.Region=="TnA" | stage1$SexStage=="Female"),]
stage <- subset(stage,select=-c(Sex,Stage,Brain.Region))

#remove those parenting Male individuals from counts
countskeep <- intersect(colnames(counts),rownames(stage)) 
counts <- counts[,countskeep,drop=FALSE]
counts <- as.matrix(counts)
mode(counts) <- "integer"

all(rownames(stage) %in% colnames(counts))
#true - all sample IDs are in 
#for dESeq, needed to get rid of countdat and treatment info for one of the groups to do two at a time
#so just doing Parenting and Courting males
all(rownames(stage) == colnames(counts))
#DESeq dataset
dds <- DESeqDataSetFromMatrix(countData=counts, colData=stage, design=~ SexStage)
dds
yes = rowSums(counts > 10) >= 5
dds2 <- dds[yes,]
#make MP the reference
dds2$Sex <-relevel(dds$SexStage, ref="MaleP")
dds2<-DESeq(dds2)
res <- results(dds2, alpha=0.05, contrast=c("SexStage", "MaleP", "MaleC"))
summary(res)

#write a file with DE data for all genes
write.csv(res, "jacana_POA_MPvMC_DESeq2allGenes.csv")
#add the gene names
annot <- read.csv("counts_with_egenes_POA.csv")
rownames(annot) <- annot[,2]
allgenes <- read.csv("jacana_POA_MPvMC_DESeq2allGenes.csv", row.names="X")
genematrix <- merge(allgenes, annot, by="row.names", all=FALSE)
head(genematrix)
write.csv(genematrix, "jacana_POA_MPvMC_DESeq2allGenes_annot_23JUN21.csv")

#looking at just the significant genes
res_sig <- subset(res,padj <0.05)
head(res_sig)
write.csv(res_sig, "jacana_POA_MPvMC_siggenes_13JUN23.csv")
#log2fold scatterplot
DESeq2::plotMA(res,main="Log2-Fold Change vs. Base Mean", ylim=c(-3,3))
##PCA
#extract transformed values
vsd <- vst(dds2)
pcaData <- plotPCA(vsd, intgroup="Sex", returnData=TRUE)
plotPCA(vsd, intgroup="Sex")


#writing a file with the counts, gene IDs, and gene names for significant genes
allcounts <- as.matrix(read.csv("counts_with_egenes_POA.csv", row.names="Gene.stable.ID"))
siggenes <- as.matrix(res_sig) #turn the significant gene results into a matrix,make sure row names are the Gene IDs
head(siggenes)
sigcounts <- merge(siggenes, allcounts, by="row.names", all=F) #merge the matrices by gene IDs, all=F so we don't have a bunch of NA lines for non-significant genes
head(sigcounts)


deg<-read.csv("jacana_POA_MPvMC_siggenes_13JUN23.csv", na.strings=c("","NA"))
rownames(deg) <- deg[,1]
annot <- read.csv("counts_with_egenes_POA.csv")
rownames(annot) <- annot[,2]
genemat <- merge(deg, annot, by="row.names", all=FALSE)
deg <- genemat
#deg<- read.csv("~/Dropbox/Rosvall_Postdoc/RNAseq/Heat Stress/hot_vs_control.DESeq2_Results.csv")
head(deg)
#with(deg, plot(logFC, -log10(padj), pch=19, main="Mouse DEG", cex=0.7, xlab=bquote(~Log[2]~fold~change), ylab=bquote(~-log[10]~P~value)))
with(deg, plot(log2FoldChange, -log10(padj), pch=19, main="Parenting vs Courting Males POA", cex=0.7, xlab=bquote(~Log[2]~fold~change), ylab=bquote(~-log[10]~P~value)))

write.csv(deg, "jacana_sigcounts_POA_MPvMC_13JUN21.csv")

#Add lines for cut-offs: logFC>0.5 and p-value cut-off at padj<0.05
abline(h=-log10(0.05), col="black", lty=3, lwd=1)
abline(v=-0.5, col="black", lty=3, lwd=1)
abline(v=0.5, col="black", lty=3, lwd=1)

#color significant genes based on whether they were up- or down-regulated
with(subset(deg, padj<0.05 & log2FoldChange< -0.5), points(log2FoldChange, -log10(padj), pch=19, col="skyblue3", cex=0.7))
with(subset(deg, padj<0.05 & log2FoldChange>0.5), points(log2FoldChange, -log10(padj), pch=19, col="tomato3", cex=0.7))

#1) Adjust the log fold change cutoff to 2 and p-value cutoff to 0.1.

#2) Add genes names to the significant genes using the code below. Try adjusting how many gene names are shown.
cutoff=sort(deg$padj)[10] #selects the top 10 smallest p values

sign.genes=which(deg$padj <= cutoff)
#sign.gene=which(deg$log2FoldChange>0.5 | deg$log2FoldChange< -0.5) & which(deg$padj < 0.05)
#text(x=deg$log2FoldChange[sign.gene] , y=-log10(deg$padj[sign.gene]), label=deg$Gene.name[sign.gene], cex=0.9)
text(x=deg$log2FoldChange[sign.genes] , y=-log10(deg$padj[sign.genes]), label=deg$Gene.name[sign.genes], cex=0.5)
