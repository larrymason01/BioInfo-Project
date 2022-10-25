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

N = 5000
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
#newdata <- topN[c("log2FoldChange", "pvalue")]
#rownames(newdata) <- topN$Gene

cols <- c("#210b1a", "#18125c", "#f54242", "#f5bc42", "#66f542", "#4287f5", "#c842f5")

K <- 8
km2 <- kmeans(newdata, K, iter.max=1000)
library(factoextra)
fviz_cluster(km2, data = newdata,
             palette = cols[1:K], 
             geom = "point",
             ellipse.type = "convex", 
             ggtheme = theme_bw()
)
