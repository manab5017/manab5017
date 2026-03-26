##############################
# Load libraries
##############################
library(edgeR)
library(data.table)
library(dplyr)
library(tibble)

##############################
# 1. Read file
##############################
countdata <- fread("FC.txt")

##############################
# 2. Fix column names
##############################
colnames(countdata) <- c("Gene","C1","C2","C3","T1","T2","T3")

##############################
# 3. Set gene names as rownames
##############################
count_matrix <- countdata %>%
  column_to_rownames(var = "Gene")

##############################
# 4. Convert to matrix
##############################
count_matrix <- as.matrix(count_matrix)

##############################
# 5. Define sample groups
##############################
sample_info <- factor(c("Control","Control","Control",
                        "Treated","Treated","Treated"))

##############################
# 6. Create DGEList
##############################
y <- DGEList(counts = count_matrix, group = sample_info)

##############################
# 7. Filter low counts
##############################
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes = FALSE]

##############################
# 8. Normalize
##############################
y <- calcNormFactors(y, method = "TMM")

##############################
# 9. Estimate dispersion
##############################
y <- estimateDisp(y)

##############################
# 10. Differential expression
##############################
et <- exactTest(y)

##############################
# 11. Get DEGs
##############################
top_degs <- topTags(et, n = Inf)

##############################
# 12. Save results
##############################
write.csv(top_degs$table, "Expression.csv")

##############################
# 13. Summary
##############################
summary(decideTests(et, lfc = 1, p.value = 0.05))

##############################
# 14. Plot
##############################
plotMD(et)
abline(h = c(-1, 1), col = "green")