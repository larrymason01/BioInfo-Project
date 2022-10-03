counts_dir <- "./data/GSE127500/GSE127500_Counts_gene.txt"
counts <- read.table(counts_dir, header = TRUE, sep = "\t")

##Names rows based on gene names
row.names(counts) <- counts$GeneId
counts <- subset(counts, select = -c(1))

##counts matrix already appears to be logarithmic
#counts <- log(counts)

#Calculate ranges and (WIP) plot density function
ranges <- apply(X = counts, MARGIN = 1, FUN = range)
diff <- ranges[2, ] - ranges[1, ]
plot(density(log(diff)))

