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

# add a column for Z linked genes
Z_POA_FvMP <- read.csv("jacana_POA_FvMP_Zlinked.csv")
sigcounts_POAFvMP_Z <- sigcountsPOAFvMP |> 
  mutate(Z_linked = sigcountsPOAFvMP$Gene.name %in% Z_POA_FvMP$Name)
write.csv(sigcounts_POAFvMP_Z, "sigcounts_POAFvMP_Z_19MAR2025.csv")
Z_POA_FvMC <- read.csv("jacana_POA_FvMC_Zlinked.csv")
sigcounts_POAFvMC_Z <- sigcountsPOAFvMC |> 
  mutate(Z_linked = sigcountsPOAFvMC$Gene.name %in% Z_POA_FvMC$Name)
write.csv(sigcounts_POAFvMC_Z, "sigcounts_POAFvMC_Z_19MAR2025.csv")
Z_TnA_FvMP <- read.csv("jacana_TnA_FvMP_Zlinked.csv")
sigcounts_TnAFvMP_Z <- sigcountsTnAFvMP |> 
  mutate(Z_linked = sigcountsTnAFvMP$Gene.name %in% Z_TnA_FvMP$Name)
write.csv(sigcounts_TnAFvMP_Z, "sigcounts_TnAFvMP_Z_19MAR2025.csv")
Z_TnA_FvMC <- read.csv("jacana_TnA_FvMC_Zlinked.csv")
sigcounts_TnAFvMC_Z <- sigcountsTnAFvMC |> 
  mutate(Z_linked = sigcountsTnAFvMC$Gene.name %in% Z_TnA_FvMC$Name)
write.csv(sigcounts_TnAFvMC_Z, "sigcounts_TnAFvMC_Z_19MAR2025.csv")

####Volcano Plots####
#POA FvMP
z_POAFvMP <-  read.csv("jacana_POA_FvMP_Zlinked.csv")
z_POAFvMP <- z_POAFvMP$Name

with(allgenesPOAFvMP, plot(log2FoldChange, -log10(padj), pch=19, main="Females vs Parenting Males POA", cex=0.7, xlab=bquote(~Log[2]~fold~change), ylab=bquote(~-log[10]~P~value), xlim=c(-3,13), ylim=c(-1,300)))

#Add lines for cut-offs: logFC>0.5 and p-value cut-off at padj<0.05
abline(h=-log10(0.05), col="black", lty=3, lwd=1)
abline(v=-0.5, col="black", lty=3, lwd=1)
abline(v=0.5, col="black", lty=3, lwd=1)

#color genes on the Z in grey
with(subset(allgenesPOAFvMP, Gene.name %in% z_POAFvMP), points(log2FoldChange, -log10(padj), pch=19, col="grey", cex=0.7))

#color significant genes based on whether they were up- or down-regulated
with(subset(allgenesPOAFvMP, padj<0.05 & log2FoldChange< -0.5 & !Gene.name %in% z_POAFvMP), points(log2FoldChange, -log10(padj), pch=19, col="blue", cex=0.7))
with(subset(allgenesPOAFvMP, padj<0.05 & log2FoldChange>0.5), points(log2FoldChange, -log10(padj), pch=19, col="red", cex=0.7))


#Add genes names to the significant genes using the code below. Try adjusting how many gene names are shown.
#cutoff=sort(allgenesPOAFvMP$padj)[10] #selects the top 10 smallest p values
#sign.genes=which(allgenesPOAFvMP$padj <= cutoff)
#text(x=allgenesPOAFvMP$log2FoldChange[sign.genes] , y=-log10(allgenesPOAFvMP$padj[sign.genes]), label=allgenesPOAFvMP$Gene.name[sign.genes], cex=0.5)

#FvMC POA
with(allgenesPOAFvMC, plot(log2FoldChange, -log10(padj), pch=19, main="Females vs Courting Males POA", cex=0.7, xlab=bquote(~Log[2]~fold~change), ylab=bquote(~-log[10]~P~value), xlim=c(-3,13), ylim=c(-1,300)))

#Add lines for cut-offs: logFC>0.5 and p-value cut-off at padj<0.05
abline(h=-log10(0.05), col="black", lty=3, lwd=1)
abline(v=-0.5, col="black", lty=3, lwd=1)
abline(v=0.5, col="black", lty=3, lwd=1) 

#color significant genes based on whether they were up- or down-regulated
with(subset(allgenesPOAFvMC, padj<0.05 & log2FoldChange< -0.5), points(log2FoldChange, -log10(padj), pch=19, col="lightblue", cex=0.7))
with(subset(allgenesPOAFvMC, padj<0.05 & log2FoldChange>0.5), points(log2FoldChange, -log10(padj), pch=19, col="red", cex=0.7))

#Add genes names to the significant genes using the code below. Try adjusting how many gene names are shown.
#cutoff=sort(allgenesFvMC$padj)[10] #selects the top 10 smallest p values
#sign.genes=which(allgenesFvMC$padj <= cutoff)
#text(x=allgenesFvMC$log2FoldChange[sign.genes] , y=-log10(allgenesFvMC$padj[sign.genes]), label=allgenesFvMC$Gene.name[sign.genes], cex=0.5)

#MCvMP POA
with(allgenesPOAMCvMP, plot(log2FoldChange, -log10(padj), pch=19, main="Courting vs Parenting Males POA", cex=0.7, xlab=bquote(~Log[2]~fold~change), ylab=bquote(~-log[10]~P~value), xlim=c(-3,3), ylim=c(0,5)))

#Add lines for cut-offs: logFC>0.5 and p-value cut-off at padj<0.05
abline(h=-log10(0.05), col="black", lty=3, lwd=1)
abline(v=-0.5, col="black", lty=3, lwd=1)
abline(v=0.5, col="black", lty=3, lwd=1)

#color significant genes based on whether they were up- or down-regulated
with(subset(allgenesPOAMCvMP, padj<0.05 & log2FoldChange< -0.5), points(log2FoldChange, -log10(padj), pch=19, col="blue", cex=0.7))
with(subset(allgenesPOAMCvMP, padj<0.05 & log2FoldChange>0.5), points(log2FoldChange, -log10(padj), pch=19, col="lightblue", cex=0.7))

#Add genes names to the significant genes using the code below. Try adjusting how many gene names are shown.
#cutoff=sort(allgenesMCvMP$padj)[10] #selects the top 10 smallest p values
#sign.genes=which(allgenesMCvMP$padj <= cutoff)
#text(x=allgenesMCvMP$log2FoldChange[sign.genes] , y=-log10(allgenesMCvMP$padj[sign.genes]), label=allgenesMCvMP$Gene.name[sign.genes], cex=0.5)

# F v MP TnA
with(allgenesTnAFvMP, plot(log2FoldChange, -log10(padj), pch=19, main="Females vs Parenting Males TnA", cex=0.7, xlab=bquote(~Log[2]~fold~change), ylab=bquote(~-log[10]~P~value), xlim=c(-3,13), ylim=c(-1,300)))

#Add lines for cut-offs: logFC>0.5 and p-value cut-off at padj<0.05
abline(h=-log10(0.05), col="black", lty=3, lwd=1)
abline(v=-0.5, col="black", lty=3, lwd=1)
abline(v=0.5, col="black", lty=3, lwd=1)

#color significant genes based on whether they were up- or down-regulated
with(subset(allgenesTnAFvMP, padj<0.05 & log2FoldChange< -0.5), points(log2FoldChange, -log10(padj), pch=19, col="blue", cex=0.7))
with(subset(allgenesTnAFvMP, padj<0.05 & log2FoldChange>0.5), points(log2FoldChange, -log10(padj), pch=19, col="tomato", cex=0.7))

#Add genes names to the significant genes using the code below. Try adjusting how many gene names are shown.
#cutoff=sort(allgenesFvMP$padj)[10] #selects the top 10 smallest p values
#sign.genes=which(allgenesFvMP$padj <= cutoff)
#text(x=allgenesFvMP$log2FoldChange[sign.genes] , y=-log10(allgenesFvMP$padj[sign.genes]), label=allgenesFvMP$Gene.name[sign.genes], cex=0.5)

#FvMC TnA
with(allgenesTnAFvMC, plot(log2FoldChange, -log10(padj), pch=19, main="Females vs Courting Males TnA", cex=0.7, xlab=bquote(~Log[2]~fold~change), ylab=bquote(~-log[10]~P~value), xlim=c(-3,13), ylim=c(-1,300)))

#Add lines for cut-offs: logFC>0.5 and p-value cut-off at padj<0.05
abline(h=-log10(0.05), col="black", lty=3, lwd=1)
abline(v=-0.5, col="black", lty=3, lwd=1)
abline(v=0.5, col="black", lty=3, lwd=1) 

#color significant genes based on whether they were up- or down-regulated
with(subset(allgenesTnAFvMC, padj<0.05 & log2FoldChange< -0.5), points(log2FoldChange, -log10(padj), pch=19, col="lightblue", cex=0.7))
with(subset(allgenesTnAFvMC, padj<0.05 & log2FoldChange>0.5), points(log2FoldChange, -log10(padj), pch=19, col="tomato", cex=0.7))

#Add genes names to the significant genes using the code below. Try adjusting how many gene names are shown.
cutoff=sort(allgenesFvMC$padj)[10] #selects the top 10 smallest p values
sign.genes=which(allgenesFvMC$padj <= cutoff)
text(x=allgenesFvMC$log2FoldChange[sign.genes] , y=-log10(allgenesFvMC$padj[sign.genes]), label=allgenesFvMC$Gene.name[sign.genes], cex=0.5)

#MCvMP TnA
with(allgenesTnAMCvMP, plot(log2FoldChange, -log10(padj), pch=19, main="Courting vs Parenting Males TnA", cex=0.7, xlab=bquote(~Log[2]~fold~change), ylab=bquote(~-log[10]~P~value), xlim=c(-3,3), ylim=c(0,5)))

#Add lines for cut-offs: logFC>0.5 and p-value cut-off at padj<0.05
abline(h=-log10(0.05), col="black", lty=3, lwd=1)
abline(v=-0.5, col="black", lty=3, lwd=1)
abline(v=0.5, col="black", lty=3, lwd=1)

#color significant genes based on whether they were up- or down-regulated
with(subset(allgenesTnAMCvMP, padj<0.05 & log2FoldChange< -0.5), points(log2FoldChange, -log10(padj), pch=19, col="blue", cex=0.7))
with(subset(allgenesTnAMCvMP, padj<0.05 & log2FoldChange>0.5), points(log2FoldChange, -log10(padj), pch=19, col="lightblue", cex=0.7))

#Add genes names to the significant genes using the code below. Try adjusting how many gene names are shown.
#cutoff=sort(allgenesMCvMP$padj)[10] #selects the top 10 smallest p values
#sign.genes=which(allgenesMCvMP$padj <= cutoff)
#text(x=allgenesMCvMP$log2FoldChange[sign.genes] , y=-log10(allgenesMCvMP$padj[sign.genes]), label=allgenesMCvMP$Gene.name[sign.genes], cex=0.5)




#### Boxplots ####
library(ggplot2)
library(patchwork)
gene <- "PRLR"
genecounts <- as.data.frame(read.csv("jacana_POATnAcounts_14FEB24.csv", row.names="Gene.stable.ID"))
P <- genecounts[genecounts$Gene.name == gene, ]
P <- t(as.matrix(P, row.names=1))
P <- P[2:49,]#get rid of gene name label
P
treat <- read.csv("jacana_allsamples_sexStageInfo.csv",header=TRUE,row.names=1)

genemat <- merge(P, treat, by="row.names", all=TRUE)
genemat <- as.data.frame(genemat)
genemat$group <- paste(genemat$Brain.Region,genemat$SexStage,sep="_")
genemat$x = as.numeric(as.character(genemat$x))#make numbers into numbers lol
POA <- genemat[grepl("POA$",genemat[,1]), ]
TnA <- genemat[grepl("TnA$",genemat[,1]), ]

#ggplot(POA, aes(x=SexStage, y=x, fill=SexStage))+geom_boxplot(color="black", alpha=0.5)+geom_jitter(aes(color=SexStage),position = position_jitter(width = 0.2), shape = 16, size = 2, alpha = 1)+scale_color_manual(values=c("Female"="red","MaleP"="blue", "MaleC"="lightblue"))+scale_fill_manual(values = c("Female" = "red", MaleP="blue","MaleC" = "lightblue"))+labs(title=paste(gene,"Expression vs Sex"),subtitle="POA", x="Sex",y=paste(gene,"Expression"))+theme_minimal()+theme(plot.title=element_text(hjust=0.5),plot.subtitle = element_text(hjust=0.5),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(color="grey",linewidth=0.5))

POA$SexStage <- factor(POA$SexStage, levels = c("Female", "MaleP", "MaleC"))
TnA$SexStage <- factor(POA$SexStage, levels = c("Female", "MaleP", "MaleC"))

Y_lim <- 300

TnAplot <- ggplot(TnA, aes(x=SexStage, y=x, fill=SexStage)) +
  geom_boxplot(color="black", alpha=0.5) +
  geom_jitter(aes(color=SexStage), position = position_jitter(width = 0.2), shape = 16, size = 2, alpha = 1) +
  scale_color_manual(values=c("Female"="red", "MaleP"="blue", "MaleC"="lightblue")) +
  scale_fill_manual(values = c("Female" = "red", "MaleP" = "blue", "MaleC" = "lightblue")) +
  labs(subtitle="TnA", x="Sex", y=paste(gene, "Expression")) +
  coord_cartesian(ylim = c(0, Y_lim)) + # Ensure same y-axis limit
  theme_minimal() +
  theme(plot.title = element_text(hjust=0.5),
        plot.subtitle = element_text(hjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color="grey", linewidth=0.5))
POAplot <- ggplot(POA, aes(x=SexStage, y=x, fill=SexStage)) +
  geom_boxplot(color="black", alpha=0.5) +
  geom_jitter(aes(color=SexStage), position = position_jitter(width = 0.2), shape = 16, size = 2, alpha = 1) +
  scale_color_manual(values=c("Female"="red", "MaleP"="blue", "MaleC"="lightblue")) +
  scale_fill_manual(values = c("Female" = "red", "MaleP" = "blue", "MaleC" = "lightblue")) +
  labs(subtitle="POA", x="Sex", y=paste(gene, "Expression")) +
  coord_cartesian(ylim = c(0, Y_lim)) + # Ensure same y-axis limit
  theme_minimal() +
  theme(plot.title = element_text(hjust=0.5),
        plot.subtitle = element_text(hjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color="grey", linewidth=0.5))
plotp <- POAplot + TnAplot + plot_layout(guides="collect")
plotp

#library(rvg)
#library(officer)
editable_graph <- dml(ggobj = plotp)
doc <- read_pptx("jacana_suppboxplots_14FEB25.pptx")
doc <- add_slide(doc)
doc <- ph_with(x = doc, editable_graph,
               location = ph_location_type(type = "body") )
print(doc, target = "jacana_suppboxplots_14FEB25.pptx")


