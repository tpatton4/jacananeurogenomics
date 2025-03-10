library(DESeq2)
#read in gene count data
setwd("C:/Users/tmp15/Downloads/R/jacana")
counts_genes <- as.matrix(read.csv("jacana_POATnAcounts_14FEB24.csv",row.names="Gene.stable.ID"))
counts_POA <- read.csv("jacana_POATnAcounts_14FEB24.csv",row.names="Gene.stable.ID")
counts_TNA <- read.csv("jacana_POATnAcounts_14FEB24.csv",row.names="Gene.stable.ID")
POA <- grep("POA", names(counts_POA), value = TRUE)
TNA <- grep("TnA", names(counts_TNA), value = TRUE)
counts_POA <- as.matrix(counts_POA[,POA])
counts_TNA <- as.matrix(counts_TNA[,TNA])
counts <- subset(counts_genes,select=-c(Gene.name)) #remove gene names column for DESeq steps
counts <- as.matrix(counts)
mode(counts) <- "integer"

treat <- read.csv("jacana_allsamples_sexStageInfo.csv",header=TRUE,row.names=1)
treat_POA <- treat[1:24,]
treat_TNA <- treat[25:48,]
treat_POA$group <- paste(treat_POA$Brain.Region,treat_POA$SexStage,sep="_")
treat_TNA$group <- paste(treat_TNA$Brain.Region,treat_TNA$SexStage,sep="_")
#treat <- subset(treat,select=-c(Sex,Stage))

all(rownames(treat_POA) %in% colnames(counts_POA))
#true - all sample IDs are in 
all(rownames(treat_POA) == colnames(counts_POA))
#true - all sample IDs are in the correct order
all(rownames(treat_TNA) %in% colnames(counts_TNA))
#true - all sample IDs are in 
all(rownames(treat_TNA) == colnames(counts_TNA))

#make a DESeq objectfor POA
POAdds <- DESeqDataSetFromMatrix(countData=counts_POA, colData=treat_POA, design= ~ group)
POAdds
yes = rowSums(counts_POA > 10) >= 5
POAdds2 <- POAdds[yes,]
#make F the reference
#dds2$Sex <-relevel(dds$SexStage, ref="Female")
POAdds2<-DESeq(POAdds2)

#make a DESeq object for TnA
TNAdds <- DESeqDataSetFromMatrix(countData=counts_TNA, colData=treat_TNA, design= ~ group)
TNAdds
yes = rowSums(counts_TNA > 10) >= 5
TNAdds2 <- TNAdds[yes,]
#make F the reference
#dds2$Sex <-relevel(dds$SexStage, ref="Female")
TNAdds2<-DESeq(TNAdds2)

#Make comparisons between groups
resPOAFvMP <- results(POAdds2, contrast=c("group", "POA_Female", "POA_MaleP"))
resPOAFvMC <- results(POAdds2, contrast=c("group", "POA_Female", "POA_MaleC"))
resPOAMCvMP <- results(POAdds2, contrast=c("group", "POA_MaleC", "POA_MaleP"))
resTnAFvMP <- results(TNAdds2, contrast=c("group", "TnA_Female", "TnA_MaleP"))
resTnAFvMC <- results(TNAdds2, contrast=c("group", "TnA_Female", "TnA_MaleC"))
resTnAMCvMP <- results(TNAdds2, contrast=c("group", "TnA_MaleC", "TnA_MaleP"))

summary(resPOAFvMP)
summary(resPOAFvMC)
summary(resPOAMCvMP)
summary(resTnAFvMP)
summary(resTnAFvMC)
summary(resTnAMCvMP)

####write results to files####
allgenesPOAFvMP <- merge(as.matrix(resPOAFvMP), counts_POA, by="row.names") #add gene names and counts
allgenesPOAFvMC <- merge(as.matrix(resPOAFvMC), counts_POA, by="row.names")
allgenesPOAMCvMP <- merge(as.matrix(resPOAMCvMP), counts_POA, by="row.names")
allgenesTnAFvMP <- merge(as.matrix(resTnAFvMP), counts_TNA, by="row.names") #add gene names and counts
allgenesTnAFvMC <- merge(as.matrix(resTnAFvMC), counts_TNA, by="row.names")
allgenesTnAMCvMP <- merge(as.matrix(resTnAMCvMP), counts_TNA, by="row.names")
#optionally, write a file with DE data for all genes, not just significantly differentially expressed
#write.csv(allgenesFvMP, "jacanaPOAallgenes_FvMP_DESeq2_16FEB24.csv")
#write.csv(allgenesFvMC, "jacanaPOAallgenes_FvMC_DESeq2_16FEB24.csv")
#write.csv(allgenesMCvMP, "jacanaPOAallgenes_MCvMP_DESeq2_16FEB24.csv")

#write a file with signifiantly differentially expressed genes, as well as their gene names and counts
#adjusted pval less than 0.05
siggenesPOAFvMP <- as.matrix(subset(resPOAFvMP,padj <0.05)) 
siggenesPOAFvMC <- as.matrix(subset(resPOAFvMC,padj <0.05))
siggenesPOAMCvMP <- as.matrix(subset(resPOAMCvMP,padj <0.05))
siggenesTnAFvMP <- as.matrix(subset(resTnAFvMP,padj <0.05)) 
siggenesTnAFvMC <- as.matrix(subset(resTnAFvMC,padj <0.05))
siggenesTnAMCvMP <- as.matrix(subset(resTnAMCvMP,padj <0.05))

sigcountsPOAFvMP <- merge(siggenesPOAFvMP, counts_POA, by="row.names", all=F) #merge the matrices by gene IDs, all=F so we don't have a bunch of NA lines for non-significant genes
sigcountsPOAFvMC <- merge(siggenesPOAFvMC, counts_POA, by="row.names", all=F)
sigcountsPOAMCvMP <- merge(siggenesPOAMCvMP, counts_POA, by="row.names", all=F)
sigcountsTnAFvMP <- merge(siggenesTnAFvMP, counts_TNA, by="row.names", all=F)
sigcountsTnAFvMC <- merge(siggenesTnAFvMC, counts_TNA, by="row.names", all=F)
sigcountsTnAMCvMP <- merge(siggenesTnAMCvMP, counts_TNA, by="row.names", all=F)

write.csv(sigcountsPOAFvMP, "jacanaPOA_FvMP_DESeq2Out_6MAR25.csv")
write.csv(sigcountsPOAFvMC, "jacanaPOA_FvMC_DESeq2Out_6MAR25.csv")
write.csv(sigcountsPOAMCvMP, "jacanaPOA_MCvMP_DESeq2Out_6MAR25.csv")
write.csv(sigcountsTnAFvMP, "jacanaTnA_FvMP_DESeq2Out_6MAR25.csv")
write.csv(sigcountsTnAFvMC, "jacanaTnA_FvMC_DESeq2Out_6MAR25.csv")
write.csv(sigcountsTnAMCvMP, "jacanaTnA_MCvMP_DESeq2Out_6MAR25.csv")

