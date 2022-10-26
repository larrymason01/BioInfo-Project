#Directories
counts_dir <- "./data/GSE127500/GSE127500_Counts_TE.txt"
metadata_dir <- "./data/GSE127500/SraRunTable.txt"

counts <- read.table(counts_dir, sep = "\t", row.names = 1, header=TRUE)
metadata <- read.csv(metadata_dir, sep=",")

#clean metadata file
metadata <- subset(metadata, LibrarySelection != "cDNA")
counts <- counts[,order(colnames(counts))]
metadata <- metadata[order(metadata$Genotype),]
metadata <- metadata[c(25:28,1:24),]
rownames(metadata) <- colnames(counts)
if (!all(rownames(metadata) == colnames(counts))) {
  stop("Rownames != colnames")
}

#dds and filtering
library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = metadata,
                              design = ~ Genotype)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

#Plot PCA
vsd <- vst(dds, blind=TRUE)
library(apeglm)
library(umap)
library(magrittr)
library(ggplot2)
#diff analysis
dds_diff <- DESeq(dds)
dds_diff_results <- results(dds_diff, contrast=c("Genotype", "Control", "7B"))
results_df <- dds_diff_results %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Gene") %>%
  dplyr::mutate(threshold = padj < 0.05) %>%
  dplyr::arrange(pvalue)
  #dplyr::arrange(dplyr::desc(pvalue))

N = 100
topN <- head(results_df, N)
topN_names <- topN$Gene
rowdata <- rowData(dds_diff)

#w/o DESeq
newdatanls <- rowdata[c("baseMean", "baseVar")]
newdatanls <- newdatanls[topN_names,]
newdata <- data.frame(log(newdatanls[,1]), log(newdatanls[,2]))
rownames(newdata) <- topN_names
colnames(newdata) <- colnames(newdatanls)

#from DESeq
# newdata <- topN[c("log2FoldChange", "pvalue")]
# rownames(newdata) <- topN$Gene

#test
# rns <- rownames(newdata)
# cns <- colnames(newdata)
# newdata <- data.frame(replicate(2,sample(0:100,5000,rep=TRUE)))
# rownames(newdata) <- rns
# colnames(newdata) <- cns

cols <- c("#210b1a", "#18125c", "#f54242", "#f5bc42", "#66f542", "#4287f5", "#c842f5")

#KMEANS
K <- 8
km2 <- kmeans(newdata, centers=K, iter.max=1000)
centers <- km2$centers[order(km2$centers[,1],decreasing=FALSE),]
rownames(centers) <- 1:8
km2 <- kmeans(newdata, centers=centers, iter.max=1000)
library(factoextra)
fviz_cluster(km2, data = newdata,
             palette = cols[1:K], 
             geom = "point",
             ellipse.type = "convex", 
             ggtheme = theme_bw()
)

#PAM
library(cluster)
pamCluster = pam(newdata, K)
fviz_cluster(pamCluster, data = newdata,
             palette = cols[1:K],
             geom = "point",
             ellipse.type = "convex",
             ggtheme = theme_bw()
)


###CHI SQUARED TESTING###
chisquaredtest <- function(cluster_vector)
{
  top_counts <- log(counts[topN_names,] + 1)
  chi_table <- data.frame(matrix(0, nrow=K, ncol=28))
  colnames(chi_table) <- colnames(top_counts)
  for (rn in rownames(top_counts)) {
    clusternum <- cluster_vector[rn]
    chi_table[clusternum,] <- chi_table[clusternum,] + top_counts[rn,]
  }
  
  grouped_chi_table <- data.frame(matrix(0, nrow=K, ncol=8))
  colnames(grouped_chi_table) <- list("Control", "1B", "2B", "3B", "4B", "5B", "6B", "7B")
  groups_counted <- data.frame(matrix(0, nrow=1, ncol=8))
  colnames(groups_counted) <- colnames(grouped_chi_table)
  for (cn in colnames(chi_table)) {
    group <- metadata[cn, "Genotype"]
    grouped_chi_table[,group] <- grouped_chi_table[,group]+ chi_table[,cn]
    groups_counted[, group] <- groups_counted[, group] + 1
  }
  for (rn in rownames(grouped_chi_table)) {
    grouped_chi_table[rn,] <- grouped_chi_table[rn,] / groups_counted
  }
  return(grouped_chi_table)
}

grouped_chi_table <- chisquaredtest(pamCluster$clustering)
chisq.test(grouped_chi_table)
grouped_chi_table
