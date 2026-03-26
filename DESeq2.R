##############################
# 1. Install & load packages
##############################
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("DESeq2","EnhancedVolcano","pheatmap"))

library(DESeq2)
library(data.table)
library(tibble)
library(ggplot2)
library(pheatmap)
library(EnhancedVolcano)

##############################
# 2. Load data
##############################
countdata <- fread("FC.txt")

##############################
# 3. Fix column names
##############################
colnames(countdata) <- c("Gene","C1","C2","C3","T1","T2","T3")

##############################
# 4. Set rownames
##############################
count_matrix <- countdata %>%
  column_to_rownames(var = "Gene")

count_matrix <- as.matrix(count_matrix)

##############################
# 5. Sample metadata
##############################
condition <- factor(c("Control","Control","Control",
                      "Treatment","Treatment","Treatment"))

colData <- data.frame(row.names = colnames(count_matrix),
                      condition)

##############################
# 6. Create DESeq dataset
##############################
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = colData,
                              design = ~ condition)

##############################
# 7. Filter low counts
##############################
dds <- dds[rowSums(counts(dds)) > 10, ]

##############################
# 8. Run DESeq2
##############################
dds <- DESeq(dds)

##############################
# 9. Results
##############################
res <- results(dds)
resOrdered <- res[order(res$padj),]

write.csv(as.data.frame(resOrdered),
          file = "DESeq2_results.csv")

##############################
# 10. Significant DEGs
##############################
resSig <- subset(resOrdered,
                 padj < 0.05 & abs(log2FoldChange) > 1)

write.csv(as.data.frame(resSig),
          file = "Significant_DEGs.csv")

##############################
# 11. Transform data
##############################
vsd <- vst(dds, blind=FALSE)

##############################
# 12. PCA Plot (Publication)
##############################
pcaData <- plotPCA(vsd, intgroup="condition", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=5) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw() +
  ggtitle("PCA Plot")

##############################
# 13. Volcano Plot (Advanced)
##############################
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 2.5,
                labSize = 3,
                col = c("grey30","forestgreen","royalblue","red2"),
                title = "Volcano Plot")

##############################
# 14. Heatmap (Top 50 genes)
##############################
topGenes <- rownames(resOrdered)[1:50]
mat <- assay(vsd)[topGenes, ]

annotation_col <- data.frame(condition)
rownames(annotation_col) <- colnames(mat)

pheatmap(mat,
         scale="row",
         annotation_col = annotation_col,
         show_rownames=FALSE,
         fontsize_col=12)

##############################
# 15. Sample distance heatmap
##############################
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)

pheatmap(sampleDistMatrix,
         main="Sample Distance")

##############################
# 16. MA plot
##############################
plotMA(res, ylim=c(-5,5))

##############################
# 17. DEG summary barplot
##############################
res_df <- as.data.frame(res)

res_df$group <- "NS"
res_df$group[res_df$padj < 0.05 & res_df$log2FoldChange > 1] <- "Up"
res_df$group[res_df$padj < 0.05 & res_df$log2FoldChange < -1] <- "Down"

barplot(table(res_df$group),
        col=c("blue","grey","red"),
        main="DEG Summary")

##############################
# 18. Dispersion plot
##############################
plotDispEsts(dds)