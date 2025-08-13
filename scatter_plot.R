# Load required libraries
library(ggplot2)
library(dplyr)
library(gridExtra)
library(reshape2)

# Read the data
sample1 <- read.delim("~/Documents/18.05_studios/R24-5134_bwa_v1.tsv", header = TRUE, stringsAsFactors = FALSE)
sample2 <- read.delim("~/Documents/18.05_studios/r24-5135_bwa_v1.tsv", header = TRUE, stringsAsFactors = FALSE)

# Rename the columns for clarity
colnames(sample1) <- c("ref", "count1")
colnames(sample2) <- c("ref", "count2")

# Merge the datasets
merged_data <- merge(sample1, sample2, by = "ref")

# Get total mapped reads for normalization (excluding special categories)
real_refs <- !grepl("antisense|not_unique|\\*", merged_data$ref)
total_sample1 <- sum(merged_data$count1[real_refs])
total_sample2 <- sum(merged_data$count2[real_refs])

# Calculate normalization factors
norm_factor_1 <- 1000000 / total_sample1  # Normalize to CPM (counts per million)
norm_factor_2 <- 1000000 / total_sample2

# Add normalized counts
merged_data$norm_count1 <- merged_data$count1 * norm_factor_1
merged_data$norm_count2 <- merged_data$count2 * norm_factor_2

# Calculate correlation coefficient using normalized counts
correlation_raw <- cor(merged_data$count1, merged_data$count2)
correlation_norm <- cor(merged_data$norm_count1, merged_data$norm_count2)

# Filter for data points where at least one sample has a count > 0
filtered_data <- merged_data %>% filter(count1 > 0 | count2 > 0)

# Add a column to identify major chromosomes (CM patterns and mitochondrial)
filtered_data$is_major <- grepl("^CM|^J01415", filtered_data$ref)

# Create a scatterplot with log scale for raw counts
p1 <- ggplot(filtered_data, aes(x = count1, y = count2, color = is_major)) +
  geom_point(alpha = 0.7) +
  scale_x_log10(breaks = c(1, 10, 100, 1000, 10000), labels = scales::comma) +
  scale_y_log10(breaks = c(1, 10, 100, 1000, 10000), labels = scales::comma) +
  labs(
    title = "Raw Count Comparison Between Samples",
    subtitle = paste("Correlation coefficient:", round(correlation_raw, 3)),
    x = "R24-5134 Count (log scale)",
    y = "R24-5135 Count (log scale)",
    color = "Major Chromosome"
  ) +
  scale_color_manual(values = c("FALSE" = "grey", "TRUE" = "blue")) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  theme_light() +
  theme(legend.position = "bottom")

# Create a scatterplot with log scale for normalized counts
p2 <- ggplot(filtered_data, aes(x = norm_count1, y = norm_count2, color = is_major)) +
  geom_point(alpha = 0.7) +
  scale_x_log10(breaks = c(1, 10, 100, 1000, 10000), labels = scales::comma) +
  scale_y_log10(breaks = c(1, 10, 100, 1000, 10000), labels = scales::comma) +
  labs(
    title = "Normalized Count Comparison (CPM)",
    subtitle = paste("Correlation coefficient:", round(correlation_norm, 3)),
    x = "R24-5134 Normalized Count (log scale)",
    y = "R24-5135 Normalized Count (log scale)",
    color = "Major Chromosome"
  ) +
  scale_color_manual(values = c("FALSE" = "grey", "TRUE" = "blue")) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  theme_light() +
  theme(legend.position = "bottom")

# For non-log scale visualization, create a plot focusing on major chromosomes
major_chroms <- filtered_data %>% 
  filter(is_major == TRUE) %>%
  arrange(desc(count1 + count2)) %>%
  head(15)  # Top 15 chromosomes by total count

# Calculate fold change based on normalized counts
major_chroms$fold_change <- major_chroms$norm_count2 / major_chroms$norm_count1
major_chroms$log2_fold_change <- log2(major_chroms$fold_change)

# Create a bar plot for fold changes
p3 <- ggplot(major_chroms, aes(x = reorder(ref, log2_fold_change), y = log2_fold_change)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(
    title = "Log2 Fold Change (R24-5135/R24-5134)",
    subtitle = "After normalization by total mapped reads",
    x = "Reference",
    y = "Log2 Fold Change"
  ) +
  theme_light() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Create a bar plot for top chromosomes
bar_data <- melt(major_chroms[, c("ref", "count1", "count2")], id.vars = "ref")
bar_data$variable <- ifelse(bar_data$variable == "count1", "R24-5134", "R24-5135")

p4 <- ggplot(bar_data, aes(x = reorder(ref, -value), y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    title = "Top Chromosomes by Count",
    x = "Reference",
    y = "Raw Count",
    fill = "Sample"
  ) +
  theme_light() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save plots to PDF
pdf("gene_count_comparison.pdf", width = 10, height = 16)
grid.arrange(p1, p2, p3, p4, ncol = 1)
dev.off()

# Display plots in RStudio
plot_grid <- grid.arrange(p1, p2, p3, p4, ncol = 1)
print(plot_grid)

# Create a more detailed statistical summary
cat("\n--- Summary Statistics ---\n")

# Count of genes with non-zero reads in both samples
both_nonzero <- sum(merged_data$count1 > 0 & merged_data$count2 > 0)
cat(sprintf("References with reads in both samples: %d\n", both_nonzero))

# Count of genes unique to each sample
only_sample1 <- sum(merged_data$count1 > 0 & merged_data$count2 == 0)
only_sample2 <- sum(merged_data$count1 == 0 & merged_data$count2 > 0)
cat(sprintf("References with reads only in R24-5134: %d\n", only_sample1))
cat(sprintf("References with reads only in R24-5135: %d\n", only_sample2))

cat(sprintf("\nTotal reads mapped to references in R24-5134: %d\n", total_sample1))
cat(sprintf("Total reads mapped to references in R24-5135: %d\n", total_sample2))
cat(sprintf("Ratio of total reads (R24-5135/R24-5134): %.2f\n", total_sample2/total_sample1))

# Display fold changes with explanation
cat("\n--- What is Fold Change? ---\n")
cat("Fold change is the ratio of normalized counts between samples, showing how much more or less abundant\n")
cat("a reference is in one sample compared to the other. A fold change of 2 means the reference is twice as abundant.\n")
cat("We typically use log2 for fold change, where:\n")
cat(" - log2(fold change) = 0 means equal abundance\n")
cat(" - log2(fold change) = 1 means twice as abundant in R24-5135\n")
cat(" - log2(fold change) = -1 means twice as abundant in R24-5134\n\n")

# Show normalized fold changes for top chromosomes
cat("Normalized fold changes for top chromosomes (R24-5135/R24-5134):\n")
print(major_chroms[, c("ref", "norm_count1", "norm_count2", "fold_change", "log2_fold_change")])