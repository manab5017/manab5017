#############################################
# VCF Mutation Analysis Script
# Author: You
#############################################

# ==============================
# 1. Install packages (run once)
# ==============================
install.packages(c("data.table", "ggplot2", "dplyr"))

# ==============================
# 2. Load libraries
# ==============================
library(data.table)
library(ggplot2)
library(dplyr)

# ==============================
# 3. Load VCF file
# ==============================
# Replace with your file
vcf_file <- "variants.vcf.gz"

# Read VCF (skip metadata lines)
vcf <- fread(vcf_file, skip = "#CHROM")

# Rename columns properly
colnames(vcf)[1:8] <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO")

cat("Total variants:", nrow(vcf), "\n")

# ==============================
# 4. SNP count per chromosome
# ==============================
snp_count <- vcf %>%
  group_by(CHROM) %>%
  summarise(count = n())

print(snp_count)

# Plot
p1 <- ggplot(snp_count, aes(x = CHROM, y = count)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "SNP Count per Chromosome",
       x = "Chromosome",
       y = "Number of SNPs")

print(p1)

# ==============================
# 5. SNP distribution plot
# ==============================
p2 <- ggplot(vcf, aes(x = POS)) +
  geom_histogram(bins = 100, fill = "steelblue") +
  facet_wrap(~CHROM, scales = "free_x") +
  theme_minimal() +
  labs(title = "SNP Distribution Across Genome",
       x = "Position",
       y = "Count")

print(p2)

# ==============================
# 6. SNP density plot
# ==============================
p3 <- ggplot(vcf, aes(x = POS)) +
  geom_density(fill = "red", alpha = 0.5) +
  facet_wrap(~CHROM, scales = "free_x") +
  theme_minimal() +
  labs(title = "SNP Density",
       x = "Position",
       y = "Density")

print(p3)

# ==============================
# 7. Filter high-quality SNPs
# ==============================
vcf_high <- vcf %>% filter(QUAL > 30)

cat("High-quality SNPs:", nrow(vcf_high), "\n")

p4 <- ggplot(vcf_high, aes(x = POS)) +
  geom_histogram(bins = 100, fill = "darkgreen") +
  facet_wrap(~CHROM, scales = "free_x") +
  theme_minimal() +
  labs(title = "High-quality SNP Distribution")

print(p4)

# ==============================
# 8. Mutation types
# ==============================
vcf$mutation <- paste(vcf$REF, vcf$ALT, sep = ">")

mut_table <- as.data.frame(table(vcf$mutation))

print(mut_table)

p5 <- ggplot(mut_table, aes(x = Var1, y = Freq)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Mutation Types",
       x = "Mutation",
       y = "Frequency") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p5)

# ==============================
# 9. Save plots
# ==============================
ggsave("snp_count.png", p1, width = 6, height = 4)
ggsave("snp_distribution.png", p2, width = 8, height = 6)
ggsave("snp_density.png", p3, width = 8, height = 6)
ggsave("high_quality_snp.png", p4, width = 8, height = 6)
ggsave("mutation_types.png", p5, width = 8, height = 6)

# ==============================
# 10. Save summary table
# ==============================
write.csv(snp_count, "snp_count_per_chromosome.csv", row.names = FALSE)

cat("Analysis complete! Files saved.\n")