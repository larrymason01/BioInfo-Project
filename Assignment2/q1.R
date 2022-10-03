counts_dir <- "./data/GSE127500/GSE127500_Counts_TE.txt"
counts <- read.table(counts_dir, header = TRUE, sep = "\t")

###Q1###
##Names rows based on gene names
#row.names(counts) <- counts$GeneId
counts <- subset(counts, select = -c(1))

##counts matrix already appears to be logarithmic
countsLog <- log(counts + 1)

#Calculate ranges and (WIP) plot density function
ranges <- apply(X = countsLog, MARGIN = 1, FUN = range)
diff <- ranges[2, ] - ranges[1, ]
plot(density(diff))


###Q2###
metadata_dir <- "./data/GSE127500/SraRunTable.txt"
metadata <- read.csv(metadata_dir, sep=",")

##clean metadata
metadata <- subset(metadata, LibrarySelection != "cDNA")
counts <- counts[,order(colnames(counts))]
metadata <- metadata[order(metadata$Genotype),]
metadata <- metadata[c(25:28,1:24),]
rownames(metadata) <- colnames(counts)

##Deseq2
library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = round(counts),
                              colData = metadata,
                              design = ~ Genotype)
dds <- DESeq(dds)
res <- results(dds, name="Genotype_comparison")
# or to shrink log fold changes association with condition: 
res <- lfcShrink(dds, coef="Genotype_2B_vs_1B", type="apeglm")


vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
plotPCA(vsd, intgroup=c("Genotype"))

##UMAP
library(umap)
library(ggplot2)
library(magrittr)
set.seed(12345)
normalized_counts <- assay(vsd) %>%
  t()

umap_results <- umap::umap(normalized_counts)
# umap_plot_df <- data.frame(umap_results$layout) %>%
#   tibble::rownames_to_column("idname") %>%
#   dplyr::inner_join(metadata, by = "idname")

umap_plot_df <- data.frame(umap_results$layout)
ggplot(
  umap_plot_df,
  aes(
    x = X1,
    y = X2
  )
) +
  geom_point() # Plot individual points to make a scatterplot

