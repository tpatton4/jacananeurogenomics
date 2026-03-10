#tissue specificity Tau for gene expression

library(readr)
library(DESeq2)
library(dplyr)

data <- read.csv("jacana_counts_genenames_8FEB24.csv")
#NOF12, NOM2, NOM3
gene_names <- data[[3]] #extract gene names
gene_IDs <- data[[1]]
count_data <- data %>%
  select(contains(c("NOF12","NOM2","NOM3"))) #keep only indivs with all tissue types

#convert to matrix form
count_matrix <- as.matrix(count_data)
rownames(count_matrix) <- gene_IDs

#reorder col names by tissue type
tissues <- c("Brain","Blood","Gonad","POA","TnA")
tissue <- sub("^([^.]+)\\.([^.]+)\\.([^.]+)\\.bam$", "\\3", colnames(count_matrix)) #pull tissue (3rd dot separated chunk) from colname
ord <- order(match(tissue, tissues))
count_matrix2 <- count_matrix[, ord]

colData <- data.frame(
  row.names = colnames(count_matrix2),
  condition = factor(tissue[ord], levels = tissues)
)

dds <- DESeqDataSetFromMatrix(countData = count_matrix2, colData = colData, design = ~ condition)


dds <- DESeq(dds)
normalized_counts <- counts(dds, normalized = TRUE)

#avg expression per tissue
mean_by_tissue <- sapply(tissues, function(t) {
  rowMeans(normalized_counts[, colData$condition == t, drop = FALSE])
})

#make tau function
tau_ts <- function(x) {
  if (all(x == 0)) return(NA_real_)  # undefined if gene not expressed
  x <- x / max(x)
  sum(1 - x) / (length(x) - 1)
}

#calcualte tau for all genes
tau_results <- apply(mean_by_tissue, 1, tau_ts)
tau_df <- data.frame(
  gene = rownames(mean_by_tissue),
  tau  = tau_results
)

#tau is represented for each gene (like unevennes in expression for each gene), so to summarize across tissues, need to do something else

#distribution of tau across all genes (used chatGPT here- not sure who does this. Maybe ask Kim if this is how she represented Tau?)
hist(tau_results,
     breaks = 50,
     col = "grey",
     border = "white",
     xlab = "Tissue specificity (τ)",
     main = "Genome-wide tissue specificity for\nblood, brain, gonad, POA, TnA")
#fraction of tissue-specific genes?
mean(tau_results > 0.8, na.rm = TRUE)
#abt 32% of genes show strong tissue specificity
#so upset plots show that most genes are expressed in all tissues, but tau tells us that those genes are expressed unevenly across tissues (biased expression)

#### Now we want to separate autosomes and Z genes ####
zgenes <- read.csv("ZgeneIDs_Dec25.csv")
mean_Z <- mean(
  tau_df$tau[tau_df$gene %in% zgenes$x],
  na.rm = TRUE
)

mean_autosome <- mean(
  tau_df$tau[!tau_df$gene %in% zgenes$x],
  na.rm = TRUE
)

Z_tau <- tau_df[tau_df$gene %in% zgenes$x, ]
mean(Z_tau$tau > 0.8, na.rm = TRUE) #29.2% of Z genes have tau > 0.8

autosome_tau <- tau_df[!tau_df$gene %in% zgenes$x, ]
mean(autosome_tau$tau > 0.8, na.rm = TRUE) #32.1% of autosomal genes have tau > 0.8
