#Directories
counts_dir <- "./data/GSE127500/GSE127500_Counts_TE.txt"
metadata_dir <- "./data/GSE127500/SraRunTable.txt"

counts <- read.table(counts_dir, sep = "\t", row.names = 1, header=TRUE)
metadata <- read.csv(metadata_dir, sep=",")

###Q1###
#log-scale data
counts_log <- log(counts + 1)

#Calculate ranges and plot density function
ranges <- apply(X = counts_log, MARGIN = 1, FUN = range)
diff <- ranges[2, ] - ranges[1, ]
png(filename="./plots/Density_plot.png")
plot(density(diff))
dev.off()

###Q2###
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
png(filename="./plots/PCA.png")
plotPCA(vsd, intgroup=c("Genotype"))
dev.off()

#Plot UMAP
library(umap)
library(magrittr)
library(ggplot2)

#prepare vsd-scaled data
normalized_vsd <- assay(vsd) %>% t()
vsd_umap <- umap(normalized_vsd)
vsd_umap_df <- data.frame(vsd_umap$layout) %>%
  tibble::rownames_to_column("sample")

#group sample name based on count
vsd_umap_df$count <- rep('control', nrow(vsd_umap_df))
for (c in 1:7)
{
  vsd_umap_df$count[grep(paste('^X', c, 'B', sep=""), vsd_umap_df$sample)] <- c
}

png(filename="./plots/UMAP.png")
ggplot(
  vsd_umap_df,
  aes(
    x = X1,
    y = X2,
    color=count
  )
) +
  geom_point(size = 4)
dev.off()


###Q3###
library(apeglm)
#differential analysis, for now we are using control / 3B
dds_diff <- DESeq(dds)
dds_diff_results <- results(dds_diff, contrast=c("Genotype", "Control", "3B"))

results_df <- dds_diff_results %>%
  as.data.frame() %>%
  # the gene names are row names -- let's make them a column for easy display
  tibble::rownames_to_column("Gene") %>%
  # add a column for significance threshold results
  dplyr::mutate(threshold = padj < 0.05) %>%
  # sort by statistic -- the highest values will be genes with
  # higher expression in RPL10 mutated samples
  dplyr::arrange(dplyr::desc(log2FoldChange))
library(readr)
readr::write_tsv(
  results_df,
  "./plots/difftable.tsv"
)
