############################################################
# ==========================================================
#        BULK RNA-seq DIFFERENTIAL EXPRESSION ANALYSIS
# ==========================================================
#
# Dataset:
# Lung cancer bulk RNA-seq
#
# Biological comparison:
# CONTROL vs CANCER
#
# Samples:
# Cancer1, Cancer2, Cancer3
# Control1, Control2, Control3
#
# Main objective:
# Identify genes whose expression differs between
# Control and Cancer samples.
#
# Statistical method:
# DESeq2
#
# Significance criteria:
# Adjusted P-value < 0.05
# |log2FC| >= 1
#
# Visualization:
# ggplot2 + pheatmap
#
# Output:
# Publication-quality PDF
# 600 DPI PNG
#
############################################################


############################################################
# ==========================================================
# 1. INSTALL REQUIRED PACKAGES
# ==========================================================
#
# WHY?
#
# RNA-seq analysis requires specialized Bioconductor and
# CRAN packages.
#
# DESeq2:
# Statistical differential expression analysis.
#
# ggplot2:
# Publication-quality graphical visualization.
#
# ggrepel:
# Prevents overlapping gene labels.
#
# pheatmap:
# Heatmap visualization.
#
############################################################

if (!requireNamespace(
  "BiocManager",
  quietly = TRUE
)) {
  
  install.packages("BiocManager")
}


cran_packages <- c(
  
  "ggplot2",
  "dplyr",
  "tidyr",
  "tibble",
  "ggrepel",
  "pheatmap",
  "RColorBrewer",
  "scales"
)


for (pkg in cran_packages) {
  
  if (!requireNamespace(
    pkg,
    quietly = TRUE
  )) {
    
    install.packages(pkg)
  }
}


if (!requireNamespace(
  "DESeq2",
  quietly = TRUE
)) {
  
  BiocManager::install(
    
    "DESeq2",
    
    ask = FALSE,
    
    update = FALSE
  )
}


############################################################
# ==========================================================
# 2. LOAD LIBRARIES
# ==========================================================
############################################################

library(DESeq2)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)
library(scales)


############################################################
# ==========================================================
# 3. CREATE OUTPUT DIRECTORIES
# ==========================================================
#
# WHY?
#
# Keeping tables, QC results and figures separately makes
# the analysis reproducible and easier to prepare for
# publication.
#
############################################################

dir.create(
  "DESeq2_results",
  showWarnings = FALSE
)

dir.create(
  "DESeq2_results/plots",
  showWarnings = FALSE
)

dir.create(
  "DESeq2_results/tables",
  showWarnings = FALSE
)

dir.create(
  "DESeq2_results/qc",
  showWarnings = FALSE
)


############################################################
# ==========================================================
# 4. DEFINE PUBLICATION COLOR PALETTE
# ==========================================================
#
# WHY?
#
# A consistent color scheme makes all figures visually
# coherent.
#
# Control:
# Blue
#
# Cancer:
# Red
#
# These colors are used consistently throughout the
# analysis.
#
############################################################

condition_colors <- c(
  
  "Control" = "#2C7FB8",
  
  "Cancer" = "#D7301F"
)


deg_colors <- c(
  
  "Higher in Control" = "#2C7FB8",
  
  "Higher in Cancer" = "#D7301F",
  
  "Not significant" = "#BDBDBD"
)


############################################################
# ==========================================================
# 5. PUBLICATION GRAPH THEME
# ==========================================================
#
# WHY?
#
# Default ggplot figures are not optimized for journal
# presentation.
#
# This theme provides:
#
# - Clear typography
# - Strong axis labels
# - Minimal background
# - Consistent spacing
# - Professional appearance
#
############################################################

publication_theme <- theme_classic(
  
  base_size = 14,
  
  base_family = "sans"
  
) +
  
  theme(
    
    plot.title =
      element_text(
        
        size = 17,
        
        face = "bold",
        
        hjust = 0.5
      ),
    
    plot.subtitle =
      element_text(
        
        size = 12,
        
        hjust = 0.5
      ),
    
    axis.title =
      element_text(
        
        size = 14,
        
        face = "bold"
      ),
    
    axis.text =
      element_text(
        
        size = 11
      ),
    
    legend.title =
      element_text(
        
        size = 12,
        
        face = "bold"
      ),
    
    legend.text =
      element_text(
        
        size = 11
      ),
    
    legend.position =
      "right",
    
    panel.border =
      element_rect(
        
        fill = NA,
        
        linewidth = 0.6
      ),
    
    plot.margin =
      margin(
        10,
        15,
        10,
        10
      )
  )


############################################################
# ==========================================================
# 6. READ FEATURECOUNTS FILE
# ==========================================================
#
# WHY?
#
# featureCounts produces gene-level read counts.
#
# DESeq2 requires a matrix where:
#
# rows    = genes
# columns = samples
# values  = raw integer read counts
#
############################################################

counts <- read.delim(
  
  "gene_count.text",
  
  comment.char = "#",
  
  check.names = FALSE,
  
  stringsAsFactors = FALSE
)


############################################################
# ==========================================================
# 7. CHECK INPUT DATA
# ==========================================================
#
# WHY?
#
# Before statistical analysis, we need to make sure the
# featureCounts file has been read correctly.
#
############################################################

cat(
  "\n========================================\n"
)

cat(
  "INPUT DATA QC\n"
)

cat(
  "========================================\n"
)

cat(
  "Number of genes:",
  nrow(counts),
  "\n"
)

cat(
  "Number of columns:",
  ncol(counts),
  "\n\n"
)

print(
  colnames(counts)
)


############################################################
# ==========================================================
# 8. CHECK FEATURECOUNTS STRUCTURE
# ==========================================================
############################################################

if (
  ncol(counts) < 12
) {
  
  stop(
    
    paste(
      
      "The count file does not contain the expected",
      "6 sample columns.",
      "Please check gene_count.text."
    )
  )
}


############################################################
# ==========================================================
# 9. EXTRACT GENE IDs
# ==========================================================
#
# WHY?
#
# Ensembl IDs sometimes contain version numbers:
#
# ENSG00000123456.5
#
# Removing ".5" makes downstream annotation easier.
#
############################################################

gene_ids <- counts$Geneid


gene_ids <- sub(
  
  "\\..*$",
  
  "",
  
  gene_ids
)


############################################################
# ==========================================================
# 10. EXTRACT COUNT MATRIX
# ==========================================================
#
# featureCounts structure:
#
# Geneid
# Chr
# Start
# End
# Strand
# Length
# Sample counts...
#
# Therefore columns 7 onward are the count matrix.
#
############################################################

count_matrix <- counts[
  
  ,
  
  7:ncol(counts)
]


############################################################
# ==========================================================
# 11. ASSIGN SAMPLE NAMES
# ==========================================================
############################################################

colnames(count_matrix) <- c(
  
  "Cancer1",
  "Cancer2",
  "Cancer3",
  
  "Control1",
  "Control2",
  "Control3"
)


############################################################
# ==========================================================
# 12. SET GENE IDs AS ROW NAMES
# ==========================================================
############################################################

rownames(
  count_matrix
) <- gene_ids


############################################################
# ==========================================================
# 13. CONVERT TO INTEGER MATRIX
# ==========================================================
#
# WHY?
#
# DESeq2 expects raw integer count data.
#
############################################################

count_matrix <- as.matrix(
  count_matrix
)


storage.mode(
  count_matrix
) <- "integer"


############################################################
# ==========================================================
# 14. CHECK DUPLICATE GENES
# ==========================================================
#
# WHY?
#
# A gene should appear only once in the final count matrix.
#
############################################################

duplicate_count <-
  
  sum(
    
    duplicated(
      rownames(count_matrix)
    )
  )


cat(
  
  "\nDuplicated gene IDs:",
  
  duplicate_count,
  
  "\n"
)


if (
  duplicate_count > 0
) {
  
  count_matrix <- rowsum(
    
    count_matrix,
    
    group =
      rownames(
        count_matrix
      )
  )
}


############################################################
# ==========================================================
# 15. REMOVE ZERO-COUNT GENES
# ==========================================================
#
# WHY?
#
# Genes with zero reads in every sample contain no
# information for differential expression testing.
#
############################################################

count_matrix <- count_matrix[
  
  rowSums(
    count_matrix
  ) > 0,
  
]


############################################################
# ==========================================================
# 16. SAMPLE INFORMATION
# ==========================================================
#
# WHY?
#
# DESeq2 needs to know which samples belong to which
# biological condition.
#
# CONTROL is explicitly defined as the reference condition.
#
############################################################

sample_info <- data.frame(
  
  row.names =
    colnames(count_matrix),
  
  condition = factor(
    
    c(
      
      "Cancer",
      "Cancer",
      "Cancer",
      
      "Control",
      "Control",
      "Control"
    ),
    
    levels = c(
      
      "Control",
      "Cancer"
    )
  )
)


print(
  sample_info
)


############################################################
# ==========================================================
# 17. CREATE DESeq2 OBJECT
# ==========================================================
#
# WHY?
#
# This combines:
#
# - raw read counts
# - sample information
# - experimental design
#
# design = ~ condition
#
# tells DESeq2 to test expression differences according
# to the Control/Cancer condition.
#
############################################################

dds <- DESeqDataSetFromMatrix(
  
  countData =
    count_matrix,
  
  colData =
    sample_info,
  
  design =
    ~ condition
)


############################################################
# ==========================================================
# 18. LOW-COUNT FILTERING
# ==========================================================
#
# WHY?
#
# Genes with extremely low counts provide little reliable
# statistical information and increase the multiple-testing
# burden.
#
# Here we keep genes with:
#
# >= 10 reads in at least 3 samples.
#
############################################################

keep <- rowSums(
  
  counts(dds) >= 10
  
) >= 3


cat(
  
  "\nGenes before filtering:",
  
  nrow(dds),
  
  "\n"
)


cat(
  
  "Genes after filtering:",
  
  sum(keep),
  
  "\n"
)


cat(
  
  "Genes removed:",
  
  sum(!keep),
  
  "\n"
)


############################################################
# ==========================================================
# 19. APPLY FILTER
# ==========================================================
############################################################

dds <- dds[
  keep,
]


############################################################
# ==========================================================
# 20. SAFETY CHECK
# ==========================================================
############################################################

if (
  nrow(dds) < 100
) {
  
  warning(
    
    paste(
      
      "Only",
      
      nrow(dds),
      
      "genes remain after filtering."
    )
  )
}


############################################################
# ==========================================================
# 21. LIBRARY SIZE
# ==========================================================
#
# WHY?
#
# Library size shows the total number of reads available
# for each sample.
#
# Large differences between samples can indicate sequencing
# depth differences that should be considered during QC.
#
############################################################

library_size <- colSums(
  count_matrix
)


library_df <- data.frame(
  
  Sample =
    names(library_size),
  
  LibrarySize =
    as.numeric(
      library_size
    )
)


library_df$Condition <- factor(
  
  ifelse(
    
    grepl(
      "^Cancer",
      library_df$Sample
    ),
    
    "Cancer",
    
    "Control"
  ),
  
  levels =
    c(
      "Control",
      "Cancer"
    )
)


############################################################
# 22. SAVE LIBRARY SIZE
############################################################

write.csv(
  
  library_df,
  
  "DESeq2_results/qc/library_size.csv",
  
  row.names = FALSE
)


############################################################
# ==========================================================
# 23. LIBRARY SIZE PLOT
# ==========================================================
############################################################

p1 <- ggplot(
  
  library_df,
  
  aes(
    
    x = Sample,
    
    y = LibrarySize,
    
    fill = Condition
  )
) +
  
  geom_col(
    
    width = 0.65,
    
    color = "white",
    
    linewidth = 0.4
  ) +
  
  scale_fill_manual(
    
    values =
      condition_colors
  ) +
  
  scale_y_continuous(
    
    labels =
      comma,
    
    expand =
      expansion(
        mult = c(
          0,
          0.05
        )
      )
  ) +
  
  publication_theme +
  
  labs(
    
    title =
      "Sequencing library size",
    
    subtitle =
      "Total mapped reads per sample",
    
    x = NULL,
    
    y =
      "Total reads",
    
    fill =
      "Condition"
  )


ggsave(
  
  "DESeq2_results/plots/Figure1_library_size.png",
  
  p1,
  
  width = 8,
  
  height = 6,
  
  dpi = 600
)


ggsave(
  
  "DESeq2_results/plots/Figure1_library_size.pdf",
  
  p1,
  
  width = 8,
  
  height = 6
)


############################################################
# ==========================================================
# 24. RUN DESeq2
# ==========================================================
#
# WHY?
#
# DESeq2 models count data using a negative binomial model
# and estimates:
#
# - size factors
# - dispersion
# - statistical significance
#
############################################################

dds <- DESeq(
  dds
)


############################################################
# ==========================================================
# 25. SAVE DESeq2 OBJECT
# ==========================================================
############################################################

saveRDS(
  
  dds,
  
  "DESeq2_results/dds.rds"
)


############################################################
# ==========================================================
# 26. NORMALIZED COUNTS
# ==========================================================
#
# WHY?
#
# Raw read counts cannot be directly compared between
# samples because sequencing depth differs.
#
# DESeq2 therefore calculates normalized counts.
#
############################################################

normalized_counts <- counts(
  
  dds,
  
  normalized = TRUE
)


############################################################
# SAVE NORMALIZED COUNTS
############################################################

write.csv(
  
  as.data.frame(
    normalized_counts
  ),
  
  "DESeq2_results/tables/normalized_counts.csv"
)


############################################################
# ==========================================================
# 27. DIFFERENTIAL EXPRESSION
# ==========================================================
#
# COMPARISON:
#
# CONTROL vs CANCER
#
# This is the central statistical comparison.
#
# Positive log2FC:
# Higher expression in CONTROL.
#
# Negative log2FC:
# Higher expression in CANCER.
#
############################################################

res <- results(
  
  dds,
  
  contrast = c(
    
    "condition",
    
    "Control",
    
    "Cancer"
  ),
  
  alpha = 0.05
)


############################################################
# ==========================================================
# 28. LOG2FC SHRINKAGE
# ==========================================================
#
# WHY?
#
# Very large fold changes from low-count genes can be noisy.
#
# Shrinkage provides more stable effect-size estimates.
#
############################################################

res_shrunk <- tryCatch(
  
  {
    
    lfcShrink(
      
      dds,
      
      contrast = c(
        
        "condition",
        
        "Control",
        
        "Cancer"
      ),
      
      type = "normal"
    )
    
  },
  
  error = function(e) {
    
    message(
      
      "LFC shrinkage unavailable; original results retained."
    )
    
    res
  }
)


############################################################
# ==========================================================
# 29. CONVERT RESULTS TO DATA FRAME
# ==========================================================
############################################################

res_df <- as.data.frame(
  res
)


res_df$GeneID <-
  rownames(res_df)


res_df <- res_df[
  
  ,
  
  c(
    
    "GeneID",
    
    setdiff(
      
      colnames(res_df),
      
      "GeneID"
    )
  )
]


############################################################
# ==========================================================
# 30. DEFINE DEG CATEGORIES
# ==========================================================
#
# Significance criteria:
#
# padj < 0.05
# |log2FC| >= 1
#
############################################################

res_df <- res_df %>%
  
  mutate(
    
    significance = case_when(
      
      !is.na(padj) &
        
        padj < 0.05 &
        
        log2FoldChange >= 1
      
      ~ "Higher in Control",
      
      
      !is.na(padj) &
        
        padj < 0.05 &
        
        log2FoldChange <= -1
      
      ~ "Higher in Cancer",
      
      
      TRUE
      
      ~ "Not significant"
    )
  )


############################################################
# ==========================================================
# 31. SAVE COMPLETE RESULTS
# ==========================================================
############################################################

write.csv(
  
  res_df,
  
  "DESeq2_results/tables/Control_vs_Cancer_all_results.csv",
  
  row.names = FALSE
)


############################################################
# ==========================================================
# 32. EXTRACT SIGNIFICANT GENES
# ==========================================================
############################################################

sig_genes <- res_df %>%
  
  filter(
    
    !is.na(padj),
    
    padj < 0.05,
    
    abs(log2FoldChange) >= 1
  )


control_genes <- res_df %>%
  
  filter(
    
    !is.na(padj),
    
    padj < 0.05,
    
    log2FoldChange >= 1
  ) %>%
  
  arrange(
    
    desc(log2FoldChange)
  )


cancer_genes <- res_df %>%
  
  filter(
    
    !is.na(padj),
    
    padj < 0.05,
    
    log2FoldChange <= -1
  ) %>%
  
  arrange(
    
    log2FoldChange
  )


############################################################
# 33. SAVE DEG TABLES
############################################################

write.csv(
  
  sig_genes,
  
  "DESeq2_results/tables/Control_vs_Cancer_significant_DEGs.csv",
  
  row.names = FALSE
)


write.csv(
  
  control_genes,
  
  "DESeq2_results/tables/genes_higher_in_Control.csv",
  
  row.names = FALSE
)


write.csv(
  
  cancer_genes,
  
  "DESeq2_results/tables/genes_higher_in_Cancer.csv",
  
  row.names = FALSE
)


############################################################
# ==========================================================
# 34. RESULTS SUMMARY
# ==========================================================
############################################################

cat(
  "\n========================================\n"
)

cat(
  "CONTROL VS CANCER RESULTS\n"
)

cat(
  "========================================\n"
)

cat(
  "Total genes:",
  nrow(res_df),
  "\n"
)

cat(
  "Higher in Control:",
  nrow(control_genes),
  "\n"
)

cat(
  "Higher in Cancer:",
  nrow(cancer_genes),
  "\n"
)

cat(
  "Total significant DEGs:",
  nrow(sig_genes),
  "\n"
)


############################################################
# ==========================================================
# 35. VST TRANSFORMATION
# ==========================================================
#
# WHY?
#
# DESeq2 count data are not ideal for PCA, clustering and
# correlation visualization.
#
# VST reduces the dependence of variance on expression level.
#
# IMPORTANT:
#
# varianceStabilizingTransformation() is used directly
# instead of vst() to avoid the nsub error encountered
# previously.
#
############################################################

vsd <- varianceStabilizingTransformation(
  
  dds,
  
  blind = FALSE
)


############################################################
# 36. EXTRACT VST MATRIX
############################################################

vsd_matrix <- assay(
  vsd
)


write.csv(
  
  vsd_matrix,
  
  "DESeq2_results/tables/VST_expression_matrix.csv"
)


############################################################
# ==========================================================
# 37. NORMALIZED EXPRESSION BOXPLOT
# ==========================================================
#
# WHY?
#
# This checks whether normalized expression distributions
# are broadly comparable across samples.
#
############################################################

norm_long <- as.data.frame(
  
  normalized_counts
  
) %>%
  
  rownames_to_column(
    "GeneID"
  ) %>%
  
  pivot_longer(
    
    cols = -GeneID,
    
    names_to = "Sample",
    
    values_to = "Expression"
  ) %>%
  
  left_join(
    
    sample_info %>%
      
      rownames_to_column(
        "Sample"
      ),
    
    by = "Sample"
  )


p2 <- ggplot(
  
  norm_long,
  
  aes(
    
    x = Sample,
    
    y =
      log2(
        Expression + 1
      ),
    
    fill = condition
  )
) +
  
  geom_boxplot(
    
    width = 0.65,
    
    outlier.size = 0.25,
    
    linewidth = 0.4,
    
    color = "grey20"
  ) +
  
  scale_fill_manual(
    
    values =
      condition_colors
  ) +
  
  publication_theme +
  
  labs(
    
    title =
      "Normalized expression distribution",
    
    subtitle =
      "Expression distributions after library-size normalization",
    
    x = NULL,
    
    y =
      "log2(normalized counts + 1)",
    
    fill =
      "Condition"
  )


ggsave(
  
  "DESeq2_results/plots/Figure2_normalized_expression.png",
  
  p2,
  
  width = 8,
  
  height = 6,
  
  dpi = 600
)


ggsave(
  
  "DESeq2_results/plots/Figure2_normalized_expression.pdf",
  
  p2,
  
  width = 8,
  
  height = 6
)


############################################################
# ==========================================================
# 38. EXPRESSION DENSITY
# ==========================================================
#
# WHY?
#
# Density curves show the overall distribution of gene
# expression across samples.
#
# Similar distributions suggest successful normalization.
#
############################################################

p3 <- ggplot(
  
  norm_long,
  
  aes(
    
    x =
      log2(
        Expression + 1
      ),
    
    color = Sample,
    
    group = Sample
  )
) +
  
  geom_density(
    
    linewidth = 0.9,
    
    alpha = 0.8
  ) +
  
  publication_theme +
  
  labs(
    
    title =
      "Normalized expression density",
    
    subtitle =
      "Distribution of gene expression across samples",
    
    x =
      "log2(normalized counts + 1)",
    
    y =
      "Density",
    
    color =
      "Sample"
  )


ggsave(
  
  "DESeq2_results/plots/Figure3_expression_density.png",
  
  p3,
  
  width = 9,
  
  height = 6,
  
  dpi = 600
)


ggsave(
  
  "DESeq2_results/plots/Figure3_expression_density.pdf",
  
  p3,
  
  width = 9,
  
  height = 6
)


############################################################
# ==========================================================
# 39. PCA
# ==========================================================
#
# WHY?
#
# PCA summarizes the largest sources of variation in the
# dataset.
#
# It helps determine:
#
# - Whether Control and Cancer separate
# - Whether biological replicates cluster
# - Whether any sample behaves as an outlier
#
############################################################

pca_data <- plotPCA(
  
  vsd,
  
  intgroup =
    "condition",
  
  returnData =
    TRUE
)


percent_var <- round(
  
  100 *
    attr(
      
      pca_data,
      
      "percentVar"
    ),
  
  1
)


p4 <- ggplot(
  
  pca_data,
  
  aes(
    
    x = PC1,
    
    y = PC2,
    
    color = condition,
    
    label = name
  )
) +
  
  geom_point(
    
    size = 5,
    
    alpha = 0.9
  ) +
  
  geom_text_repel(
    
    size = 4.2,
    
    box.padding = 0.6,
    
    point.padding = 0.3
  ) +
  
  scale_color_manual(
    
    values =
      condition_colors
  ) +
  
  publication_theme +
  
  labs(
    
    title =
      "Principal component analysis",
    
    subtitle =
      "Sample-level expression variation",
    
    x =
      paste0(
        "PC1 (",
        percent_var[1],
        "%)"
      ),
    
    y =
      paste0(
        "PC2 (",
        percent_var[2],
        "%)"
      ),
    
    color =
      "Condition"
  )


ggsave(
  
  "DESeq2_results/plots/Figure4_PCA.png",
  
  p4,
  
  width = 8,
  
  height = 7,
  
  dpi = 600
)


ggsave(
  
  "DESeq2_results/plots/Figure4_PCA.pdf",
  
  p4,
  
  width = 8,
  
  height = 7
)


############################################################
# SAVE PCA DATA
############################################################

write.csv(
  
  pca_data,
  
  "DESeq2_results/qc/PCA_coordinates.csv",
  
  row.names = FALSE
)


############################################################
# ==========================================================
# 40. SAMPLE CORRELATION
# ==========================================================
#
# WHY?
#
# Biological replicates should generally show relatively
# high expression correlation.
#
# This helps identify unusual or potentially problematic
# samples.
#
############################################################

cor_matrix <- cor(
  
  vsd_matrix,
  
  method = "pearson"
)


write.csv(
  
  cor_matrix,
  
  "DESeq2_results/qc/sample_correlation_matrix.csv"
)


############################################################
# CORRELATION HEATMAP
############################################################

png(
  
  "DESeq2_results/plots/Figure5_sample_correlation.png",
  
  width = 2600,
  
  height = 2400,
  
  res = 300
)

pheatmap(
  
  cor_matrix,
  
  color =
    colorRampPalette(
      
      c(
        
        "#F7FBFF",
        
        "#6BAED6",
        
        "#08306B"
      )
      
    )(100),
  
  border_color =
    "white",
  
  display_numbers =
    TRUE,
  
  number_format =
    "%.2f",
  
  fontsize =
    11,
  
  fontsize_number =
    9,
  
  main =
    "Pearson correlation between samples"
)

dev.off()


pdf(
  
  "DESeq2_results/plots/Figure5_sample_correlation.pdf",
  
  width = 8,
  
  height = 7
)

pheatmap(
  
  cor_matrix,
  
  color =
    colorRampPalette(
      
      c(
        
        "#F7FBFF",
        
        "#6BAED6",
        
        "#08306B"
      )
      
    )(100),
  
  border_color =
    "white",
  
  display_numbers =
    TRUE,
  
  number_format =
    "%.2f",
  
  fontsize =
    11,
  
  fontsize_number =
    9,
  
  main =
    "Pearson correlation between samples"
)

dev.off()


############################################################
# ==========================================================
# 41. SAMPLE DISTANCE
# ==========================================================
#
# WHY?
#
# Sample-distance analysis measures how different samples
# are from one another globally.
#
# Similar samples should have smaller distances.
#
############################################################

sample_distance <- dist(
  
  t(
    vsd_matrix
  )
)


sample_distance_matrix <-
  
  as.matrix(
    sample_distance
  )


write.csv(
  
  sample_distance_matrix,
  
  "DESeq2_results/qc/sample_distance_matrix.csv"
)


############################################################
# SAMPLE DISTANCE HEATMAP
############################################################

png(
  
  "DESeq2_results/plots/Figure6_sample_distance.png",
  
  width = 2600,
  
  height = 2400,
  
  res = 300
)

pheatmap(
  
  sample_distance_matrix,
  
  color =
    colorRampPalette(
      
      c(
        
        "#F7FCF5",
        
        "#74C476",
        
        "#00441B"
      )
      
    )(100),
  
  border_color =
    "white",
  
  display_numbers =
    TRUE,
  
  number_format =
    "%.2f",
  
  fontsize =
    11,
  
  fontsize_number =
    9,
  
  main =
    "Sample-to-sample distance"
)

dev.off()


############################################################
# ==========================================================
# 42. HIERARCHICAL CLUSTERING
# ==========================================================
#
# WHY?
#
# Hierarchical clustering provides an additional view of
# sample relationships.
#
# Biological replicates should generally cluster together
# if the biological signal is strong and technical variation
# is limited.
#
############################################################

cluster <- hclust(
  
  sample_distance,
  
  method =
    "complete"
)


png(
  
  "DESeq2_results/plots/Figure7_hierarchical_clustering.png",
  
  width = 2600,
  
  height = 2000,
  
  res = 300
)

plot(
  
  cluster,
  
  main =
    "Hierarchical clustering of samples",
  
  xlab = "",
  
  sub = "",
  
  cex = 1.3,
  
  hang = -1
)

dev.off()


############################################################
# ==========================================================
# 43. MA PLOT
# ==========================================================
#
# WHY?
#
# The MA plot displays:
#
# X-axis:
# Mean expression
#
# Y-axis:
# Log2 fold change
#
# It helps evaluate the magnitude of differential expression
# across the expression range.
#
############################################################

png(
  
  "DESeq2_results/plots/Figure8_MA_plot.png",
  
  width = 2600,
  
  height = 2200,
  
  res = 300
)

plotMA(
  
  res,
  
  alpha =
    0.05,
  
  main =
    "MA plot: Control vs Cancer"
)

abline(
  
  h = 0,
  
  lty = 2,
  
  linewidth = 1
)

dev.off()


############################################################
# ==========================================================
# 44. VOLCANO PLOT
# ==========================================================
#
# WHY?
#
# Volcano plots simultaneously show:
#
# X-axis:
# Effect size (log2FC)
#
# Y-axis:
# Statistical significance
#
# This makes it easy to identify genes with both large
# expression changes and strong statistical evidence.
#
############################################################

volcano_df <- res_df %>%
  
  filter(
    
    !is.na(padj),
    
    !is.na(log2FoldChange),
    
    is.finite(padj),
    
    is.finite(log2FoldChange),
    
    padj > 0
  ) %>%
  
  mutate(
    
    negLog10Padj =
      -log10(padj)
  )


############################################################
# TOP GENES TO LABEL
############################################################

top_labels <- volcano_df %>%
  
  filter(
    
    significance !=
      "Not significant"
  ) %>%
  
  arrange(
    
    padj
  ) %>%
  
  head(15)


############################################################
# VOLCANO PLOT
############################################################

p9 <- ggplot(
  
  volcano_df,
  
  aes(
    
    x = log2FoldChange,
    
    y = negLog10Padj,
    
    color = significance
  )
) +
  
  geom_point(
    
    size = 1.8,
    
    alpha = 0.65
  ) +
  
  geom_vline(
    
    xintercept =
      c(
        -1,
        1
      ),
    
    linetype =
      "dashed",
    
    linewidth =
      0.5,
    
    color =
      "grey40"
  ) +
  
  geom_hline(
    
    yintercept =
      -log10(0.05),
    
    linetype =
      "dashed",
    
    linewidth =
      0.5,
    
    color =
      "grey40"
  ) +
  
  geom_point(
    
    data =
      top_labels,
    
    size = 2.8,
    
    alpha = 1
  ) +
  
  geom_text_repel(
    
    data =
      top_labels,
    
    aes(
      label = GeneID
    ),
    
    size = 3.4,
    
    box.padding =
      0.6,
    
    point.padding =
      0.25,
    
    max.overlaps =
      15,
    
    show.legend =
      FALSE
  ) +
  
  scale_color_manual(
    
    values =
      deg_colors
  ) +
  
  publication_theme +
  
  labs(
    
    title =
      "Differential gene expression",
    
    subtitle =
      "Control vs Cancer",
    
    x =
      "log2 fold change (Control / Cancer)",
    
    y =
      expression(
        -log[10](
          "adjusted P-value"
        )
      ),
    
    color =
      NULL
  )


ggsave(
  
  "DESeq2_results/plots/Figure9_volcano_Control_vs_Cancer.png",
  
  p9,
  
  width = 9,
  
  height = 7,
  
  dpi = 600
)


ggsave(
  
  "DESeq2_results/plots/Figure9_volcano_Control_vs_Cancer.pdf",
  
  p9,
  
  width = 9,
  
  height = 7
)


############################################################
# ==========================================================
# 45. TOP 50 DEG HEATMAP
# ==========================================================
#
# WHY?
#
# The heatmap shows the expression pattern of the most
# statistically significant genes across all samples.
#
# Row-wise Z-scores emphasize relative expression patterns.
#
############################################################

top50 <- res_df %>%
  
  filter(
    !is.na(padj)
  ) %>%
  
  arrange(
    padj
  ) %>%
  
  head(50) %>%
  
  pull(
    GeneID
  )


heatmap_genes <- intersect(
  
  top50,
  
  rownames(vsd_matrix)
)


heatmap_matrix <- vsd_matrix[
  
  heatmap_genes,
  
  ,
  
  drop = FALSE
]


############################################################
# ROW-WISE Z-SCORE
############################################################

heatmap_matrix <- t(
  
  scale(
    
    t(
      heatmap_matrix
    )
  )
)


############################################################
# SAMPLE ANNOTATION
############################################################

annotation_col <- data.frame(
  
  Condition =
    sample_info$condition
)


rownames(
  annotation_col
) <-
  rownames(
    sample_info
  )


annotation_colors <- list(
  
  Condition =
    condition_colors
)


############################################################
# HEATMAP
############################################################

png(
  
  "DESeq2_results/plots/Figure10_top50_DEG_heatmap.png",
  
  width = 3000,
  
  height = 3400,
  
  res = 300
)

pheatmap(
  
  heatmap_matrix,
  
  annotation_col =
    annotation_col,
  
  annotation_colors =
    annotation_colors,
  
  color =
    colorRampPalette(
      
      c(
        
        "#2166AC",
        
        "#F7F7F7",
        
        "#B2182B"
      )
      
    )(100),
  
  border_color =
    NA,
  
  fontsize_row =
    7,
  
  fontsize_col =
    11,
  
  cluster_rows =
    TRUE,
  
  cluster_cols =
    TRUE,
  
  main =
    "Top 50 differentially expressed genes"
)

dev.off()


############################################################
# ==========================================================
# 46. TOP GENES HIGHER IN CONTROL
# ==========================================================
#
# WHY?
#
# These genes have:
#
# log2FC >= 1
# padj < 0.05
#
# Therefore they have significantly higher expression
# in Control compared with Cancer.
#
############################################################

top_control <- control_genes %>%
  
  head(20)


p11 <- ggplot(
  
  top_control,
  
  aes(
    
    x =
      reorder(
        GeneID,
        log2FoldChange
      ),
    
    y =
      log2FoldChange
  )
) +
  
  geom_col(
    
    width = 0.7,
    
    fill =
      condition_colors[
        "Control"
      ]
  ) +
  
  coord_flip() +
  
  publication_theme +
  
  labs(
    
    title =
      "Top genes higher in Control",
    
    x = NULL,
    
    y =
      "log2 fold change"
  )


ggsave(
  
  "DESeq2_results/plots/Figure11_top_Control_genes.png",
  
  p11,
  
  width = 9,
  
  height = 7,
  
  dpi = 600
)


ggsave(
  
  "DESeq2_results/plots/Figure11_top_Control_genes.pdf",
  
  p11,
  
  width = 9,
  
  height = 7
)


############################################################
# ==========================================================
# 47. TOP GENES HIGHER IN CANCER
# ==========================================================
#
# WHY?
#
# These genes have:
#
# log2FC <= -1
# padj < 0.05
#
# Because the contrast is Control vs Cancer, negative
# log2FC indicates higher expression in Cancer.
#
############################################################

top_cancer <- cancer_genes %>%
  
  arrange(
    log2FoldChange
  ) %>%
  
  head(20)


p12 <- ggplot(
  
  top_cancer,
  
  aes(
    
    x =
      reorder(
        GeneID,
        log2FoldChange
      ),
    
    y =
      abs(log2FoldChange)
  )
) +
  
  geom_col(
    
    width = 0.7,
    
    fill =
      condition_colors[
        "Cancer"
      ]
  ) +
  
  coord_flip() +
  
  publication_theme +
  
  labs(
    
    title =
      "Top genes higher in Cancer",
    
    x = NULL,
    
    y =
      "|log2 fold change|"
  )


ggsave(
  
  "DESeq2_results/plots/Figure12_top_Cancer_genes.png",
  
  p12,
  
  width = 9,
  
  height = 7,
  
  dpi = 600
)


ggsave(
  
  "DESeq2_results/plots/Figure12_top_Cancer_genes.pdf",
  
  p12,
  
  width = 9,
  
  height = 7
)


############################################################
# ==========================================================
# 48. LOG2FC DISTRIBUTION
# ==========================================================
#
# WHY?
#
# Shows the overall distribution of gene expression
# changes between Control and Cancer.
#
############################################################

fc_df <- res_df %>%
  
  filter(
    
    !is.na(log2FoldChange)
  )


p13 <- ggplot(
  
  fc_df,
  
  aes(
    x = log2FoldChange
  )
) +
  
  geom_histogram(
    
    bins = 80,
    
    fill =
      "#6A51A3",
    
    color =
      "white",
    
    linewidth =
      0.25
  ) +
  
  geom_vline(
    
    xintercept = 0,
    
    linetype =
      "dashed",
    
    linewidth =
      0.7,
    
    color =
      "grey30"
  ) +
  
  publication_theme +
  
  labs(
    
    title =
      "Distribution of log2 fold changes",
    
    x =
      "log2 fold change (Control vs Cancer)",
    
    y =
      "Number of genes"
  )


ggsave(
  
  "DESeq2_results/plots/Figure13_log2FC_distribution.png",
  
  p13,
  
  width = 8,
  
  height = 6,
  
  dpi = 600
)


############################################################
# ==========================================================
# 49. ADJUSTED P-VALUE DISTRIBUTION
# ==========================================================
#
# WHY?
#
# Shows the distribution of multiple-testing-adjusted
# statistical significance across genes.
#
############################################################

padj_df <- res_df %>%
  
  filter(
    
    !is.na(padj),
    
    padj >= 0,
    
    padj <= 1
  )


p14 <- ggplot(
  
  padj_df,
  
  aes(
    x = padj
  )
) +
  
  geom_histogram(
    
    bins = 60,
    
    fill =
      "#756BB1",
    
    color =
      "white",
    
    linewidth =
      0.25
  ) +
  
  geom_vline(
    
    xintercept =
      0.05,
    
    linetype =
      "dashed",
    
    linewidth =
      0.7,
    
    color =
      "grey30"
  ) +
  
  publication_theme +
  
  labs(
    
    title =
      "Adjusted P-value distribution",
    
    x =
      "Adjusted P-value",
    
    y =
      "Number of genes"
  )


ggsave(
  
  "DESeq2_results/plots/Figure14_adjusted_Pvalue_distribution.png",
  
  p14,
  
  width = 8,
  
  height = 6,
  
  dpi = 600
)


############################################################
# ==========================================================
# 50. DEG SUMMARY
# ==========================================================
#
# WHY?
#
# Provides a simple overview of how many genes are:
#
# - Higher in Control
# - Higher in Cancer
# - Not significant
#
############################################################

summary_df <- data.frame(
  
  Category = c(
    
    "Higher in Control",
    
    "Higher in Cancer",
    
    "Not significant"
  ),
  
  Count = c(
    
    sum(
      res_df$significance ==
        "Higher in Control"
    ),
    
    sum(
      res_df$significance ==
        "Higher in Cancer"
    ),
    
    sum(
      res_df$significance ==
        "Not significant"
    )
  )
)


p15 <- ggplot(
  
  summary_df,
  
  aes(
    
    x = Category,
    
    y = Count,
    
    fill = Category
  )
) +
  
  geom_col(
    
    width = 0.65,
    
    color = "white",
    
    linewidth = 0.4
  ) +
  
  geom_text(
    
    aes(
      label = comma(Count)
    ),
    
    vjust = -0.4,
    
    size = 4.5,
    
    fontface =
      "bold"
  ) +
  
  scale_fill_manual(
    
    values = c(
      
      "Higher in Control" =
        condition_colors[
          "Control"
        ],
      
      "Higher in Cancer" =
        condition_colors[
          "Cancer"
        ],
      
      "Not significant" =
        "#BDBDBD"
    )
  ) +
  
  scale_y_continuous(
    
    labels = comma,
    
    expand =
      expansion(
        mult = c(
          0,
          0.1
        )
      )
  ) +
  
  publication_theme +
  
  labs(
    
    title =
      "Differential expression summary",
    
    x = NULL,
    
    y =
      "Number of genes",
    
    fill =
      NULL
  )


ggsave(
  
  "DESeq2_results/plots/Figure15_DEG_summary.png",
  
  p15,
  
  width = 8,
  
  height = 6,
  
  dpi = 600
)


ggsave(
  
  "DESeq2_results/plots/Figure15_DEG_summary.pdf",
  
  p15,
  
  width = 8,
  
  height = 6
)


############################################################
# ==========================================================
# 51. DISPERSION PLOT
# ==========================================================
#
# WHY?
#
# Dispersion describes the variability of counts relative
# to their mean expression.
#
# DESeq2 estimates dispersion before differential expression
# testing.
#
############################################################

png(
  
  "DESeq2_results/plots/Figure16_dispersion.png",
  
  width = 2600,
  
  height = 2200,
  
  res = 300
)

plotDispEsts(
  dds
)

dev.off()


############################################################
# ==========================================================
# 52. COOK'S DISTANCE
# ==========================================================
#
# WHY?
#
# Cook's distance identifies observations that may strongly
# influence individual gene-level statistical results.
#
# This is primarily a diagnostic rather than a biological
# result.
#
############################################################

cooks <- assays(
  dds
)[["cooks"]]


cooks_df <- data.frame(
  
  Sample = rep(
    
    colnames(cooks),
    
    each =
      nrow(cooks)
  ),
  
  CookDistance =
    as.vector(cooks)
)


p17 <- ggplot(
  
  cooks_df,
  
  aes(
    
    x = Sample,
    
    y = CookDistance
  )
) +
  
  geom_boxplot(
    
    linewidth =
      0.4,
    
    outlier.size =
      0.2
  ) +
  
  scale_y_log10() +
  
  publication_theme +
  
  labs(
    
    title =
      "Cook's distance",
    
    subtitle =
      "Influence diagnostics across samples",
    
    x = NULL,
    
    y =
      "Cook's distance"
  ) +
  
  theme(
    
    axis.text.x =
      element_text(
        
        angle = 45,
        
        hjust = 1
      )
  )


ggsave(
  
  "DESeq2_results/plots/Figure17_Cooks_distance.png",
  
  p17,
  
  width = 8,
  
  height = 6,
  
  dpi = 600
)


ggsave(
  
  "DESeq2_results/plots/Figure17_Cooks_distance.pdf",
  
  p17,
  
  width = 8,
  
  height = 6
)


############################################################
# ==========================================================
# 53. SAVE QC DATA
# ==========================================================
############################################################

dispersion_df <- data.frame(
  
  GeneID =
    rownames(dds),
  
  BaseMean =
    mcols(dds)$baseMean,
  
  Dispersion =
    mcols(dds)$dispersion
)


write.csv(
  
  dispersion_df,
  
  "DESeq2_results/qc/dispersion_values.csv",
  
  row.names = FALSE
)


############################################################
# ==========================================================
# 54. SAVE SESSION INFORMATION
# ==========================================================
#
# WHY?
#
# A journal submission should be reproducible.
#
# sessionInfo() records the R version and package versions
# used for the analysis.
#
############################################################

writeLines(
  
  capture.output(
    sessionInfo()
  ),
  
  "DESeq2_results/sessionInfo.txt"
)


############################################################
# ==========================================================
# 55. FINAL SUMMARY
# ==========================================================
############################################################

cat(
  
  "\n\n"
)

cat(
  
  "====================================================\n"
)

cat(
  
  "          RNA-seq ANALYSIS COMPLETED\n"
)

cat(
  
  "====================================================\n"
)

cat(
  
  "Comparison: CONTROL vs CANCER\n"
)

cat(
  
  "Adjusted P-value cutoff: < 0.05\n"
)

cat(
  
  "Fold-change cutoff: |log2FC| >= 1\n"
)

cat(
  
  "Total genes analyzed: ",
  
  nrow(res_df),
  
  "\n",
  
  sep = ""
)

cat(
  
  "Higher in Control: ",
  
  nrow(control_genes),
  
  "\n",
  
  sep = ""
)

cat(
  
  "Higher in Cancer: ",
  
  nrow(cancer_genes),
  
  "\n",
  
  sep = ""
)

cat(
  
  "Total significant DEGs: ",
  
  nrow(sig_genes),
  
  "\n",
  
  sep = ""
)

cat(
  
  "\nOutput directory:\n"
  
)

cat(
  
  "DESeq2_results/\n"
  
)

cat(
  
  "\n====================================================\n"
)

cat(
  
  "Publication-quality figures saved as:\n"
  
)

cat(
  
  "600 DPI PNG + PDF\n"
  
)

cat(
  
  "====================================================\n"
)