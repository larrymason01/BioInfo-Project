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

##Plot PCA
vsd <- vst(dds, blind=TRUE)
png(filename="./plots/PCA.png")
plotPCA(vsd, intgroup=c("Genotype"))
dev.off()

##Plot UMAP
library(umap)
library(magrittr)
library(ggplot2)

#prepare vsd-scaled data
normalized_vsd <- assay(vsd) %>% t()
vsd_umap <- umap(normalized_vsd)
vsd_umap_df <- data.frame(vsd_umap$layout) %>%
  tibble::rownames_to_column("sample")

#group sample name based on count; create new count column expressing this data
vsd_umap_df$count <- rep('control', nrow(vsd_umap_df))
for (c in 1:7)
{
  vsd_umap_df$count[grep(paste('^X', c, 'B', sep=""), vsd_umap_df$sample)] <- c
}
#plotting
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
#differential analysis, for now we are using control / 1B
dds_diff <- DESeq(dds)
dds_diff_results <- results(dds_diff, contrast=c("Genotype", "Control", "1B"))

#clean up table
results_df <- dds_diff_results %>%
  as.data.frame() %>%
  # the gene names are row names -- let's make them a column for easy display
  tibble::rownames_to_column("Gene") %>%
  # add a column for significance threshold results
  dplyr::mutate(threshold = padj < 0.05) %>%
  # sort by pval
  dplyr::arrange(padj)

#save as tsv
library(readr)
readr::write_tsv(
  results_df,
  "./plots/difftable.tsv"
)

#volcano plot
library(EnhancedVolcano)
volcano_plot <- EnhancedVolcano::EnhancedVolcano(
  results_df,
  lab = results_df$Gene,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.01, # Loosen the cutoff since we supplied corrected p-values
  labSize = 3.0
)
png(filename="./plots/Volcano.png")
volcano_plot
dev.off()


###heatmap###
top <- head(results_df, 30)
crows <- as.matrix(counts[top$Gene,])
library("circlize")
library("ComplexHeatmap")
col_fun = colorRamp2(c(0, 1000, 10000), c("green", "white", "red"))
column_ha <- HeatmapAnnotation(Chromosome_count = rep(0:7, times=c(4, 4, 4, 4, 4, 2, 3, 3)))
hm <- Heatmap(crows, name="Counts", col=col_fun, top_annotation = column_ha, row_names_gp = grid::gpar(fontsize = 10))
png(filename="./plots/Heatmap.png")
hm
dev.off()



