library(future)
plan(multiprocess)

library("ComplexHeatmap")
set.seed(123)

ddsVST <- vst(dds)
ddsVST <- assay(ddsVST)

ddsVST <- as.data.frame(ddsVST)
ddsVST$Gene <- rownames(ddsVST)
head(ddsVST)

resDataFrame <- as.data.frame(res)

sigGenes <- rownames(resDataFrame[resDataFrame$padj <= .05 & abs(resDataFrame$log2FoldChange) > 3,])
ddsVST <- ddsVST[ddsVST$Gene %in% sigGenes,]

Heatmap(sigGenes)
