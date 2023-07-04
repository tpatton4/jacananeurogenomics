# WGCNA for ZEFI gonad hot vs control

# Set working directory
setwd("C:/Users/tmp15/Downloads/R/jacana")

# Load the WGCNA package
library(WGCNA);
library("genefilter")
options(stringsAsFactors = FALSE);
enableWGCNAThreads()

#If you've already run through the code once, load the workspace back in (won't have to re-run sft or TOM)
load("TnA_WGCNA_21MAR2023.RData")

# Data formatting  (Part 1)
#Read in gene counts - can use raw or normalized counts
TNA = read.csv("counts_with_egenes_TnA.csv")
head(TNA)
TNA <- subset(TNA, select=-c(Gene.name,Annotation.Names,X,GSF2960.NOF6.TnA))
head(TNA) #just counts and Gene.stable.ID
# Take a quick look at what is in the data set:
dim(TNA); # 14588 24
names(TNA); # ENSEMBL ID is the first column

# Reformat - get rid of unnecessary columns to only include count data
PTR_normcounts <- data.frame(TNA, row.names=TNA[,1]) #got rid of [,-1] after PTR because that erases NOF1
head(PTR_normcounts)
dim(PTR_normcounts)
PTR_normcounts = PTR_normcounts[-1]
head(PTR_normcounts)
dim(PTR_normcounts) # 14588  23
names(PTR_normcounts)
head(PTR_normcounts)
#We suggest removing features whose counts are consistently low (for example, removing all 
#features that have a count of less than say 10 in more than 90% of the samples) 
#because such low-expressed features tend to reflect noise and correlations based on 
#counts that are mostly zero aren't really meaningful. 
use = rowSums(PTR_normcounts > 10) >=20 # 90% (20 individuals) must have at least a normalized count of 10
countMatrixFiltered = PTR_normcounts[ use, ]
dim(countMatrixFiltered) # 11267    24
names(countMatrixFiltered)
head(countMatrixFiltered)
#We then recommend a variance-stabilizing transformation. Start with normalized counts 
#(or RPKM/FPKM data) and log-transform them using log2(x+1).
#Whether one uses RPKM, FPKM, or simply normalized counts doesn't make a whole lot of difference for WGCNA analysis as long as all samples were processed the same way.
log2countMatrixFiltered = log2(countMatrixFiltered+1)

#Remove the genes with the lowest variance (lower 25 percentile).
log2countMatrixFiltered = as.matrix(log2countMatrixFiltered)
head(log2countMatrixFiltered)
rv = rowVars(log2countMatrixFiltered)
q25 = quantile(rowVars(log2countMatrixFiltered), .25, na.rm=TRUE)
Filtered = log2countMatrixFiltered[rv > q25, ]
dim(Filtered) # Left with 8455  genes
head(Filtered)
write.csv(Filtered,file="jacana_TnAall_normCounts_logreduced25_threshold4.csv")

# Read gene annotation file - so that we can match gene names with Ensembl IDs
TNA = read.csv("counts_with_egenes_TnA.csv")
annot = subset(TNA,select=c(Gene.stable.ID,Gene.name))
head(annot)

########################################### Read in trait data ##########################################
traitData = read.csv("jacana_TnA_WGCNATraitData.csv");
dim(traitData) # 24 48
names(traitData)

# remove columns that hold information we do not need.
library(tidyverse)
datTrait = subset(traitData, select=-c(Individual_Region, Testosterone, Swoops, WingSpread, Flyovers, Facial.Shield.Length, Bill.length, Bill.width, Wing.length))
datTrait
dim(datTrait) # 24 `25
names(datTrait)
head(datTrait,n=10)

# Form a data frame analogous to expression data that will hold the traits.
#datTraits <- datTrait[,-1] # Get rid of first column... did this do anything?
#datTraits <- datTrait[,-1] # Get rid of first column... did this do anything?
#dim(datTrait)
rownames(datTrait) <- datTrait[,1] #make indiv IDs the row names
head(datTrait)
datTrait <- subset(datTrait,select=-c(Individual)) #Ids are row names, so we can remove them from the matrix
head(datTrait)
datTrait=datTrait[-9,] #remove F6 (she's in row 9)
dim(datTrait) #23 24
head(datTrait,n=10)

datExpr<-read.csv("jacana_TnAall_normCounts_logreduced25_threshold4.csv")
dim(datExpr) # 8455   24

# Note that each row corresponds to a gene and column to a sample or auxiliary information.
# We now remove the auxiliary data and transpose the expression data for further analysis.
#datExpr0 = t(datExpr[, -1]) # this gets rid of Ensembl names...
datExpr0 = t(datExpr[,-1])
names(datExpr0) = datExpr$X                    #Use the first column name after '$'
rownames(datExpr0) = names(datExpr)[-1]

gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK # TRUE

sampleTree = hclust(dist(datExpr0), method = "average");
sizeGrWindow(12,9)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
# No obvious outliers now that F6 is gone

# Re-cluster samples
sampleTree2 = hclust(dist(datExpr0), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
head(datTrait)
datTraits1 <- datTrait
datTrait1 <- data.matrix(datTrait,rownames.force=TRUE)
head(datTrait)
head(datTrait1)

datTraits <- datTrait1


######################################################################
######################################################################
################ SIGNED vs UNSIGNED ##################################
######################################################################
######################################################################

traitColors = numbers2colors(datTraits, signed = TRUE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTrait),
                    main = "Sample dendrogram and trait heatmap")


############## Soft connectivity #######

# Choose a set of soft-thresholding powers
# R^2>0.80
# mean connectivity should be high so that the network contains enough information (e.g. for module detection)
#Scale-free topology fit needs to reach values above 0.8 for reasonable powers 
#(less than 15 for unsigned or signed hybrid networks, and less than 30 for signed) 
#and the mean connectivity remains relatively high (in the hundreds or above).
#If this doesn't happen, chances are that the data exhibit a strong driver that 
#makes a subset of the samples globally different from the rest. 
#The difference causes high correlation among large groups of genes which 
#invalidates the assumption of the scale-free topology approximation.
#If you can't meet these criteria,
# default choice of B: for an unsigned/signedhybrid; B=9, for a signed network B=12. for <20 samples (conservative)
#unsigned = positive or negative correlations equal a connection
#signed = negatively correlated genes are considered unconnected
#signed hybrid puts negative correlations at 0 for cosmetic reasons - not that different from signed
#the resulting matrix is always non-negative #s so best to use signed hybrid to avoid confusion
powers = c(c(1:10), seq(from = 12, to=20, by=2))
#sft = pickSoftThreshold(datExpr, networkType="signed hybrid",corOptions = list(maxPOutliers =0.1),
# powerVector = powers, verbose = 5,corFnc= "bicor")
#### Why doesn't this work with datExpr0, but does work with datExpr?
### Idk what that sft was supposed to do, but the one below works and is from Horvath tutorial
sft <- pickSoftThreshold(t(countMatrixFiltered), powerVector=powers, verbose=5)

# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2),mar=c(4.5,4.5,0.5,0.5),oma=c(0,0,0,0));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.80,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

#We now calculate the adjacencies:
#power of 9 looks best now
adjacencymatrix = adjacency(datExpr0, power=9, corFnc="bicor", 
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
#black        blue       brown       green greenyellow     grey 
#314        1174         676         546         141        1201
#magenta    pink      purple         red         tan   turquoise
#237         285         188         545         114        2449
#yellow 
#585 

# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

#The Dynamic Tree Cut may identify modules whose expression profiles are very similar. It may be prudent to merge
#such modules since their genes are highly co-expressed.

### How do you know if your modules need to be merged?

# Calculate eigengenes
MEList = moduleEigengenes(datExpr0, colors = dynamicColors)
moduleColors=MEList$validColors
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
# Yes, looks like we should merge a few modules (black + green)

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
#save(MEs, moduleLabels, moduleColors, geneTree, file = "jacana_POA_normcounts_modules.RData")

#identify modules that are significantly associated with the measured clinical traits.
#Since we already have a summary profile (eigengene) for each module, we simply correlate eigengenes with external
#traits and look for the most significant associations:
# Define numbers of genes and samples
nGenes = ncol(datExpr0);
nSamples = nrow(datExpr0);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr0, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
write.csv(MEs,file="jacana_TnA_ModuleEigengenes_noF6signed.csv")

moduleTraitCor = bicor(MEs, datTraits, use = "p",robustY = FALSE,maxPOutliers =0.1);

#moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));

# color changes for heat map
library(RColorBrewer)
dim(moduleTraitCor)
#write.csv(moduleTraitCor,file="jacana_TnA_ModuleTraitCor_noF6signed.csv")
dim(datTraits)
length(colnames(datTraits))
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = colorRampPalette(brewer.pal(8, "Spectral"))(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,cex.lab=0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
table(moduleColors)
#black        blue       brown greenyellow        grey     magenta 
#1048        1174        3670         141        1201         237 
#pink         tan      yellow 
#285         114         585 
#write.csv(MEs, file = "f")

#moduleColors

#focus on one trait
# Define variable weight containing the weight column of datTrait
Mass_min = as.data.frame(datTraits1$Mass_min);
names(Mass_min) = "Mass_min"

# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr0, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(bicor(datExpr0, Mass_min, use = "p",robustY = FALSE,maxPOutliers =0.1));

geneTraitSignificance
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(Mass_min), sep="");
names(GSPvalue) = paste("p.GS.", names(Mass_min), sep="");

#link trait with module
module = "blue"
  column = match(module, modNames);
  moduleGenes = moduleColors==module;
  sizeGrWindow(7, 7);
  par(mfrow = c(1,1));
  verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                     abs(geneTraitSignificance[moduleGenes, 1]),
                     xlab = paste("Module Membership in", module, "module"),
                     ylab = "Gene significance for Mass_min",
                     main = paste("Module membership vs. gene significance\n"),
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
  #GS and MM are highly correlated, illustrating that genes highly significantly
  #associated with a trait are often also the most important (central) elements of modules associated with the trait
  
  #### What genes are in each module??
  #black.ID <- colnames(datExpr)[moduleColors=="black"] ## Couldn't get this to work
  #write.csv(black.ID, file = "ZF_gonad_normcounts_black_Ensembl.ID.csv")
  
  # Summary output of network analysis results 
  #annot = read.csv(file = "ZF_PTR_annGeneNames.csv")
  head(annot)
  dim(annot) # 14588     2
  #annot<-annot[,c(1,18,19)]
  names(annot)
  
  # Need to add column names back to dataframe
  names(datExpr0) = datExpr$X
  probes = datExpr$X
  probes2annot = match(probes, annot$Gene.stable.ID)
  
  # The following is the number or probes without annotation:
  sum(is.na(probes2annot))
  # Should return 0
  
  #We now create a data frame holding the following information for all probes: gene ID, gene symbol, Locus Link ID
  #(Entrez code), module color, gene significance for weight, and module membership and p-values in all modules. The
  #modules will be ordered by their significance for weight, with the most significant ones to the left.
  
  
  # Create the starting data frame
  geneInfo0 = data.frame(geneID = probes,
                         genename = annot$Gene.name[probes2annot],
                         moduleColor = moduleColors,
                         geneTraitSignificance,
                         GSPvalue)
  # Order modules by their significance for trait of interest
  modOrder = order(-abs(cor(MEs, Mass_min, use = "p")))
  #modOrder = order(-abs(cor(MEs, Neutral.Posture, use = "p")))
  #modOrder = order(-abs(cor(MEs, Piloerect, use = "p")))
  
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
  geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.Mass_min))
  #geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.Neutral.Posture))
  #geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.Piloerect))
  
  geneInfo = geneInfo0[geneOrder, ]
  
  write.csv(geneInfo, file = "jacana_TnA_MassMin_NoF6_signed_nomerge.csv")

  
  
save.image(file="TnA_WGCNA_21MAR2023.RData")
  