counts_dir <- "./data/GSE127500/GSE127500_Counts_TE.txt"
counts <- read.table(counts_dir, header = TRUE, sep = "\t")

##Names rows based on gene names
row.names(counts) <- counts$GeneId
counts <- subset(counts, select = -c(1))

##counts matrix already appears to be logarithmic
countsLog <- log(counts + 1)

#Calculate ranges and (WIP) plot density function
ranges <- apply(X = countsLog, MARGIN = 1, FUN = range)
diff <- ranges[2, ] - ranges[1, ]
plot(density(diff))


###Q2
metadata_dir <- "./data/GSE127500/SraRunTable.txt"
metadata <- read.csv(metadata_dir, sep=",")

##Cut cDNA portion of metadata
metadata <- subset(metadata, LibrarySelection != "cDNA")


counts <- counts[,order(colnames(counts))]
metadata <- metadata[order(metadata$Genotype),]
metadata <- metadata[c(25:28,1:24),]
rownames(metadata) <- colnames(counts)
# colnames(counts) <- sub("_EB[0-9]+", "", colnames(counts))
# colnames(counts) <- sub("X", "", colnames(counts))
# colnames(counts) <- sub("Ct_D[0-9]", "Control", colnames(counts))

library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = round(counts),
                              colData = metadata,
                              design = ~ Genotype)
dds

dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients
res <- results(dds, name="Genotype_2B_vs_1B")
# or to shrink log fold changes association with condition: 
BiocManager::install("apeglm")
res <- lfcShrink(dds, coef="Genotype_2B_vs_1B", type="apeglm")

vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 3)

plotPCA(vsd, intgroup=c("Genotype"))

library(M3C)
data(pollen)
tsne(pollen$data,colvec=c('gold'))


