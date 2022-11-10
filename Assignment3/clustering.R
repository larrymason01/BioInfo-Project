setwd("Bioinfo-Project/Assignment3")
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
# newdata <- topN[c("log2FoldChange", "pvalue")]
# rownames(newdata) <- topN$Gene

#test
# rns <- rownames(newdata)
# cns <- colnames(newdata)
# newdata <- data.frame(replicate(2,sample(0:100,N,rep=TRUE)))
# rownames(newdata) <- rns
# colnames(newdata) <- cns

cols <- c("#210b1a", "#18125c", "#f54242", "#f5bc42", "#66f542", "#4287f5", "#c842f5")

#KMEANS
K <- 8
km2 <- kmeans(newdata, centers=K, iter.max=1000)
centers <- km2$centers[order(km2$centers[,1],decreasing=TRUE),]
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

#Hierarchical
dist_mat <- dist(newdata, method = 'euclidean')
hcl <- hclust(dist_mat, method = 'average')
hierarchicalCl <- cutree(hcl, k = 8)
cl <- as.data.frame(hierarchicalCl)

cldata <- merge(cl, counts, by='row.names')
names(cldata)[1] <- "Gene"
names(cldata)[2] <- "Cluster #"

logdata <- log(cldata[, 3:30] + 1)

logdata["Gene"] <- cldata[,1]
logdata["ClusterNum"] <- cldata[,2]
logdata <- logdata[order(logdata$ClusterNum),]

testld <- logdata[,1:28]

testld <- as.matrix(testld)

library("ComplexHeatmap")
row_ha <- rowAnnotation(foo1= logdata$ClusterNum)
hm <- Heatmap(testld, cluster_rows = FALSE, row_names_gp = gpar(fontsize=4), right_annotation = row_ha)
hm

library(mclust)
res <- Mclust(newdata)
res1 <- Mclust(newdata, x = res$BIC)
GMMres <- res1$classification
plot(res1, what = "classification")
resDR <- MclustDR(res)
plot(resDR, what = "contour")
# #ConsensusClustering
# library(ConsensusClusterPlus)
# library(data.table)
# t_newdata <- t(as.matrix(newdata))
# ccp <- ConsensusClusterPlus(t_newdata, maxK = 8)
# ccp_res <- ccp[[8]]$consensusClass

# Alluvial diagram

# For creating the 10k input in the gmm diagram, a similar process was used
# to create the diagrams for other clustering methods

cl <- as.data.frame(res1$classification)
freq10k_gmm<- table(cl['res1$classification'])
freq10k_gmm <- as.data.frame(freq10k_gmm)
freq10k_gmm$Freq <- as.numeric(freq10k_gmm$Freq) / 10000

freq10k_gmm <- cbind(genecount = 10000, freq10k_gmm)


# As the workspace in R saves variables, the above code was altered each time 
# to create each input based on the number of genes. 
alluvtable <- rbind(freq10_gmm, freq100_gmm, freq1k_gmm, freq10k_gmm)
alluvtable <- cbind(map = TRUE, alluvtable)

ggplot(as.data.frame(alluvtable),
       aes(y = Freq, axis1 = genecount, axis2 = res1.classification)) +
  geom_alluvium(aes(fill = map), width = 1/12) +
  geom_stratum(width = 1/12, fill = "black", color = "grey") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Gene Count", "Cluster"), expand = c(0.05, 0.05)) +
  scale_fill_brewer(type = "qual", palette = "Set1") +
  ggtitle("Cluster Size Comparison by Percentage (GMM Clustering)")


#heatmap
cldata <- merge(GMMres, counts, by='row.names')
names(cldata)[1] <- "Gene"
names(cldata)[2] <- "Cluster #"
logdata <- log(cldata[, 3:30] + 1)
logdata["Gene"] <- cldata[,1]
logdata["ClusterNum"] <- cldata[,2]
logdata <- logdata[order(logdata$ClusterNum),]
testld <- logdata[,1:28]
testld <- as.matrix(testld)
library("ComplexHeatmap")
library("circlize")
col_fun = colorRamp2(c(0, 8), c("white", "green"))
col_fun2 = colorRamp2(c(0, 8), c("white", "blue"))
row_ha <- rowAnnotation(col=list(Cluster_Num = col_fun2), Cluster_Num= logdata$ClusterNum)
column_ha <- HeatmapAnnotation(col=list(Chromosome_count = col_fun), Chromosome_count = rep(0:7, times=c(4, 4, 4, 4, 4, 2, 3, 3)))
hm <- Heatmap(testld, name = "Log Scaled Counts", row_names_gp = gpar(fontsize=4), top_annotation = column_ha, right_annotation = row_ha, row_labels = rep("", N) )
hm



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

grouped_chi_table_pam <- chisquaredtest(pamCluster$clustering)
grouped_chi_table_km <- chisquaredtest(km2$cluster)
grouped_chi_table_hierarchical <- chisquaredtest(hierarchicalCl)
grouped_chi_table_GMM <- chisquaredtest(GMMres)

#temp
# cc <- read.csv(file = 'data/consensus_clustering_vect.csv')
# cc_vector <-  cc$ccp_res
# names(cc_vector) <- cc[,1]

sum_pam <- rowSums(grouped_chi_table_pam)
sum_km <- rowSums(grouped_chi_table_km)
sum_hierarchical <- rowSums(grouped_chi_table_hierarchical)
sum_GMM <- rowSums(grouped_chi_table_GMM)

pam_km <- data.frame(sum_pam, sum_km)
pam_hierarchical <- data.frame(sum_pam, sum_hierarchical)
pam_GMM <- data.frame(sum_pam, sum_GMM)
km_hierarchical <- data.frame(sum_km, sum_hierarchical)
km_GMM <- data.frame(sum_km, sum_GMM)
GMM_hierarchical <- data.frame(sum_hierarchical, sum_GMM)

chisq.test(pam_km)
chisq.test(pam_hierarchical)
chisq.test(km_hierarchical)
p <- list(chisq.test(pam_km)$p.val, chisq.test(pam_hierarchical)$p.val, chisq.test(km_hierarchical)$p.val)
p.adjust(p, method = p.adjust.methods, n = length(p))
