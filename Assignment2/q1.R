counts_dir <- "./data/GSE127500/GSE127500_Counts_TE.txt"
counts <- read.table(counts_dir, header = TRUE, sep = "\t")

##Names rows based on gene names
row.names(counts) <- counts$GeneId
counts <- subset(counts, select = -c(1))

##counts matrix already appears to be logarithmic
counts <- log(counts + 1)

#Calculate ranges and (WIP) plot density function
ranges <- apply(X = counts, MARGIN = 1, FUN = range)
diff <- ranges[2, ] - ranges[1, ]
plot(density(diff))


###Q2
metadata_dir <- "./data/GSE127500/SraRunTable.txt"
metadata <- read.csv(metadata_dir, sep=",")
metadata <- subset(metadata, LibrarySelection != "cDNA")
#metadata <- metadata[,c("Genotype")]

counts <- counts[,order(colnames(counts))]
colnames(counts) <- sub("_EB[0-9]+", "", colnames(counts))
colnames(counts) <- sub("X", "", colnames(counts))
colnames(counts) <- sub("Ct_D[0-9]", "Control", colnames(counts))

rownames(metadata) <- metadata$Run
metadata <- metadata[order(metadata$Genotype),]
metadata <- metadata[c(25:28,1:24),]


# rownames(metadata) <- colnames(counts)
#library("DESeq2")
# dds <- DESeqDataSetFromMatrix(countData = counts,
#                               colData = metadata,
#                               design = ~ Genotype)
# dds
