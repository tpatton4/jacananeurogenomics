# WGCNA for jacana TnA samples 20FEB2024

# Set working directory
#setwd("C:/Users/tmp15/Downloads/R/jacana")
setwd("~/Dropbox/Rosvall_Postdoc/Neurogenomics of SRR - Jacanas/WGCNA")

# Install packages
# if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
# BiocManager::install(c("impute", "preprocessCore", "GO.db", "AnnotationDbi"))
# BiocManager::install(c("genefilter"))
# BiocManager::install(c("apeglm","pheatmap","ashr","goseq","biomaRt","clusterProfiler","enrichplot","ggupset"))
# BiocManager::install(c("ensembldb"))

# Load the WGCNA package
library(WGCNA);
library("genefilter")
options(stringsAsFactors = FALSE);
enableWGCNAThreads()

#If you've already run through the code once, load the workspace back in (won't have to re-run sft or TOM)
#load("TnA_WGCNA_20FEB24.RData")

# Data formatting  (Part 1)

#Read in gene counts - can use raw or normalized counts
TNA1 = read.csv("jacana_counts_genenames_TnA_8FEB24.csv") # Sex is binary, females = 0, males = 1; # Male breeding stage, parenting = 0, courting = 1
TNA <- subset(TNA1, select=-c(Gene.name, GSF2960.NOF6.TnA)) #remove gene names and outlier NOF6
head(TNA) #just counts and Gene.stable.ID

# Take a quick look at what is in the data set:
dim(TNA); # 14588 24
names(TNA); # ENSEMBL ID is the first column

# Reformat - get rid of unneccessary columns to only include count data
TNA_normcounts <- data.frame(TNA, row.names=TNA[,1])
TNA_normcounts = TNA_normcounts[-1]
dim(TNA_normcounts) # 14588  23
names(TNA_normcounts)
#remove low counts
use = rowSums(TNA_normcounts > 10) >=5 # 5 individuals must have at least a count of 10
countMatrixFiltered = TNA_normcounts[ use, ]
dim(countMatrixFiltered) # 11818    23
names(countMatrixFiltered)

#now we normalize using variance stabilizing transformation
library(DESeq2)
library(ggplot2)
library(dplyr)

#need numerical trait data to stay, but first 4 columns should be factors
traitData <-  read.csv("jacana_TnA_WGCNATraitData.csv")
traitData <- data.frame(traitData, row.names=traitData[,1])
traitData <- traitData[-1] #remove names since they're colnames now
traitData <- traitData[row.names(traitData) != "GSF2960.NOF6.TnA", ] #remove F6 outlier
traitData[] <- data.frame(
  lapply(seq_along(traitData), function(i) {
    if (i<=4) {
     as.factor(traitData[[i]])
   } else{
     as.numeric(as.character(traitData[[i]]))
   }
  }),
    stringsAsFactors=FALSE)

TNAdds <- DESeqDataSetFromMatrix(
  countData=countMatrixFiltered, 
  colData= traitData,
  design= ~ Group)

TNAdds
TNAvsd <- vst(TNAdds, blind = FALSE) #blind = FALSE bc we will use experimental design info to estimate dispersion trend
TNAvst_counts <- assay(TNAvsd)

# Check how well normalization worked
rld <- rlog(TNAdds, blind = FALSE)
dds <- estimateSizeFactors(TNAdds)

df <- bind_rows(
  as_data_frame(log2(counts(dds, normalized=TRUE)[, 1:2]+1)) %>%
    mutate(transformation = "log2(x + 1)"),
  as_data_frame(assay(TNAvsd)[, 1:2]) %>% mutate(transformation = "vst"),
  as_data_frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"))

colnames(df)[1:2] <- c("x", "y")  

lvls <- c("log2(x + 1)", "vst", "rlog")
df$transformation <- factor(df$transformation, levels=lvls)

ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation)  


# Start with normalized counts (vst) and log-transform them using log2(x+1).
#Whether one uses RPKM, FPKM, or simply normalized counts doesn't make a whole lot of difference for WGCNA analysis as long as all samples were processed the same way.
log2countMatrixFiltered = log2(TNAvst_counts+1)

#Remove the genes with the lowest variance (lower 25 percentile).
log2countMatrixFiltered = as.matrix(log2countMatrixFiltered)
rv = rowVars(log2countMatrixFiltered)
q25 = quantile(rowVars(log2countMatrixFiltered), .25, na.rm=TRUE)
Filtered = log2countMatrixFiltered[rv > q25, ]
dim(Filtered) # Left with 8863  genes
#write.csv(Filtered,file="jacana_TnA_logreduced25_20FEB24.csv")

#### Read gene annotations from gene counts file - so that we can match gene names with Ensembl IDs ####
annot = subset(TNA1,select=c(Gene.stable.ID,Gene.name))
#write.csv(annot, file="jacana_TnA_geneAnnotations.csv")

#### Read in trait data #####
traitData = read.csv("jacana_TnA_WGCNATraitData.csv");
traitData <- traitData[traitData$Individual != "GSF2960.NOF6.TnA", ] #get rid of NOF6, an outlier
dim(traitData) # 23 26
names(traitData)
#view(traitData)

# remove columns that hold information we do not need.
library(tidyverse)

# For heatmap in supplement
#datTrait = subset(traitData, select=-c(Group, Individual_Region, Testosterone, Facial.Shield.Length, Facial.Shield.Width, Bill.length, Bill.width, Wing.length))

# For heatmap in main text
datTrait = subset(traitData, select=-c(Group,log.T., WingSpread, Flyovers,Avg.Testis.Mass, Ovary.mass, Julian.Collected, Brain.Time, Euthanize.Time, Individual_Region, Testosterone, Bill.length, Bill.width, Facial.Shield.Length, Facial.Shield.Width, Wing.length))

dim(datTrait) # 23 9
names(datTrait) #these are the traits that will show up in the heatmap
datTrait[,1] #individuals for which we have trait data

#### One morphology variable
sizemorphtest <- dplyr::select(datTrait, Mass_min, Tarsus.length, Avg.Spur)
scaled_traits <- scale(sizemorphtest)
pca_results <- prcomp(scaled_traits, center=TRUE, scale.=TRUE)
  # error for missing data in facial shield width
pca_results$rotation #looks like everything except facial shield loads together
pc1 <- pca_results$x[,1] #pc1 scores (first column in scores matrix)
datTrait$PCAtrait <- pc1 *(-1) # reverse negative to positive to reflect larger body size in females
print(datTrait)

# Remove extra morphology traits now summaried in PC1
datTrait <- subset(datTrait, select=-c(Mass_min, Tarsus.length, Avg.Spur)) #PCAtrait now reflects these measurements, which loaded together on PC1
#biplot(pca_results, scale=0)
datTrait <- subset(datTrait, select=c(Individual, Sex, Male.Breeding.Stage, Distance, Swoops, Vocalizations, PCAtrait)) #PCAtrait now reflects these measurements, which loaded together on PC1

# Form a data frame analogous to expression data that will hold the traits.
datTraits <- datTrait[,-1] #new df with individual ID as rowname
rownames(datTraits) <- datTrait[,1]
datExpr<-Filtered
dim(datExpr) # 8863  23

# Note that each row corresponds to a gene and column to a sample or auxiliary information.
# We now remove the auxiliary data and transpose the expression data for further analysis.
datExpr0 = t(datExpr)
#names(datExpr0) = datExpr$X   #Use the first column name after '$'
#rownames(datExpr0) = names(datExpr)[-1]

gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK # TRUE

sampleTree = hclust(dist(datExpr0), method = "average");
sizeGrWindow(12,9)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
# NOF6 is an outlier, if you havent already removed it. Go back and take her out


# Re-cluster samples
sampleTree2 = hclust(dist(datExpr0), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
head(datTraits)
datTraits.matrix <- data.matrix(datTraits, rownames.force=NA)

traitColors = numbers2colors(datTraits.matrix, signed = TRUE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")


####Soft connectivity#####
# Choose a set of soft-thresholding powers

powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft <- pickSoftThreshold(t(countMatrixFiltered), powerVector=powers, verbose=5)

# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2),mar=c(4.5,4.5,0.5,0.5),oma=c(0,0,0,0));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n");
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.80,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

sft$fitIndices

#We now calculate the adjacencies:
#power of 18 for TNA? 12 doesn't look much different, so using that for now
adjacencymatrix = adjacency(datExpr0, power=12, corFnc="bicor", 
                            corOptions = list(maxPOutliers =0.1),type="signed hybrid")

#To minimize effects of noise and spurious associations, we transform the adjacency into Topological Overlap Matrix,
#and calculate the corresponding dissimilarity:
TOM = TOMsimilarity(adjacencymatrix, TOMType="signed") #this takes a minute to run
dissTOM = 1-TOM

#We now use hierarchical clustering to produce a hierarchical clustering tree (dendrogram) of genes.
# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");

# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);
# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2,cutHeight = 0.99, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)
# Convert numeric labels into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)

# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

#The Dynamic Tree Cut may identify modules whose expression profiles are very similar. It may be prudent to merge such modules since their genes are highly co-expressed.

# Calculate eigengenes
MEList = moduleEigengenes(datExpr0, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
MEDissThres = 0.25 # Corresponds to similarity of.75 to merge
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Yes, looks like we should merge 

# Call an automatic merging function
merge = mergeCloseModules(datExpr0, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;
sizeGrWindow(12, 9)
#pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
# Save module colors and labels for use in subsequent parts
#save(MEs, moduleLabels, moduleColors, geneTree, file = "jacana_TnA_normcounts_modules.RData")

#identify modules that are significantly associated with the measured clinical traits.
#Since we already have a summary profile (eigengene) for each module, we simply correlate eigengenes with external traits and look for the most significant associations:
# Define numbers of genes and samples
nGenes =  ncol(datExpr0);
nSamples = nrow(datExpr0);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr0, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
#write.csv(MEs,file="jacana_TnA_ModuleEigengenes_4NOV25.csv")

moduleTraitCor = bicor(MEs, datTraits.matrix, use = "p",robustY = FALSE,maxPOutliers =0.1);

#moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));

#color palette for heat map
palette <- colorRampPalette(c("deepskyblue3","white","coral"))(50)
# color changes for heat map
library(RColorBrewer)
dim(moduleTraitCor)
dim(datTraits)
length(colnames(datTraits))
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = palette,
               textMatrix = NULL, # can change this to NULL for no text
               setStdMargins = FALSE,
               cex.text = 0.5,cex.lab=0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))


table(moduleColors)
#black         blue        brown         cyan  greenyellow         grey midnightblue         pink       purple 
#182          678          703           72          110         4379           46          139          501 
#salmon    turquoise       yellow 
#74         1670          309 


#identify modules that are significantly associated with the measured traits
#we already have a summary profile (eigengene) for each module, we simply correlate eigengenes with external traits and look for the most significant associations:




# Female-specific heatmap
traitData <-  read.csv("jacana_TnA_WGCNATraitData_noNOF6_females.csv")

# For heatmap in main text
datTraits = subset(traitData, select=-c(Distance,PCAtrait,Vocalizations,Individual, Male.Breeding.Stage,Sex,Group, WingSpread, Flyovers,Avg.Testis.Mass, Ovary.mass, Julian.Collected, Brain.Time, Euthanize.Time, Individual_Region, Testosterone, Bill.length, Bill.width, Facial.Shield.Length, Facial.Shield.Width, Wing.length,Mass_min, Tarsus.length, Avg.Spur, Avg.Testis.Mass))
head(datTraits)

# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr0, moduleColors)$eigengenes
MEs = orderMEs(MEs0)

moduleTraitCor = bicor(MEs, datTraits, use = "p",robustY = FALSE,maxPOutliers =0.1);
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
#sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)

# Display the correlation values within a heatmap plot
palette <- colorRampPalette(c("deepskyblue3","white","coral"))(50)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = palette,
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1))



# Male-specific heatmap
traitData <-  read.csv("jacana_TnA_WGCNATraitData_males.csv")

# For heatmap in main text
datTraits = subset(traitData, select=-c(Flyovers,PCAtrait,Individual,Sex,Group, Avg.Testis.Mass, log.T., Ovary.mass, Julian.Collected, Brain.Time, Euthanize.Time, Individual_Region, Testosterone, Bill.length, Bill.width, Facial.Shield.Length, Facial.Shield.Width, Wing.length,Mass_min, Tarsus.length, Avg.Spur))
head(datTraits)

# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr0, moduleColors)$eigengenes
MEs = orderMEs(MEs0)

moduleTraitCor = bicor(MEs, datTraits, use = "p",robustY = FALSE,maxPOutliers =0.1);
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
#sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)

# Display the correlation values within a heatmap plot
palette <- colorRampPalette(c("deepskyblue3","white","coral"))(50)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = palette,
               textMatrix = NULL,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1))







#focus on one trait - Distance
# Define variable weight containing the weight column of datTrait
Distance = as.data.frame(datTraits$Distance);
names(Distance) = "Distance"

# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr0, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(bicor(datExpr0, Distance, use = "p",robustY = FALSE,maxPOutliers =0.1));

#geneTraitSignificance
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(Distance), sep="");
names(GSPvalue) = paste("p.GS.", names(Distance), sep="");

#link trait with module
module = "pink"
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for Distance",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
#GS and MM are highly correlated, illustrating that genes highly significantly
#associated with a trait are often also the most important (central) elements of modules associated with the trait


# # Now for cyan and PCAtrait
PCAtrait = as.data.frame(datTraits$PCAtrait);
names(PCAtrait) = "PCAtrait"
# 
# # names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr0, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(bicor(datExpr0, PCAtrait, use = "p",robustY = FALSE,maxPOutliers =0.1));
# 
# #geneTraitSignificance
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(PCAtrait), sep="");
names(GSPvalue) = paste("p.GS.", names(PCAtrait), sep="");
# 
# #link trait with module
module = "cyan"
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                    abs(geneTraitSignificance[moduleGenes, 1]),
                    xlab = paste("Module Membership in", module, "module"),
                    ylab = "Gene significance for Morphology PC1",
                    main = paste("Module membership vs. gene significance\n"),
                    cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

# Summary output of network analysis results 
head(annot)
dim(annot) # 14588     2
#annot<-annot[,c(1,18,19)]
names(annot)

# Need to add column names back to dataframe
names(datExpr0) = row.names(datExpr)
probes = row.names(datExpr)
probes2annot = match(probes, annot$Gene.stable.ID)

# The following is the number or probes without annotation:
sum(is.na(probes2annot))
# Should return 0

#We now create a data frame holding the following information for all probes: gene ID, gene symbol, module color, gene significance for Distance, and module membership and p-values in all modules. 
#Themodules will be ordered by their significance for weight, with the most significant ones to the left.


# Create the starting data frame
geneInfo0 = data.frame(geneID = probes,
                       genename = annot$Gene.name[probes2annot],
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)
# Order modules by their significance for trait of interest
modOrder = order(-abs(cor(MEs, PCAtrait, use = "p")))

# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}

# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
#geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.Distance))

geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$PCAtrait))
geneInfo = geneInfo0[geneOrder, ]

write.csv(geneInfo, file = "jacana_TnA_WGCNAmodule_4Nov25_PCAtrait.csv")
