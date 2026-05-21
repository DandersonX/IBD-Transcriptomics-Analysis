# ==============================================================================
# Pipeline: Transcriptomic Profiling of Active vs. Inactive IBD
# Dataset: GSE193677 (Adult IBD Cohort - 1,832 Samples)
# Description: Custom data engineering pipeline using limma for DGE modeling
# ==============================================================================

# Step 1: Load Core Dependencies
library(GEOquery)
library(limma)
library(ggplot2)

# Step 2: Fetch and Load GEO Core Objects
message("Loading dataset core matrix from memory/GEO...")
if(!exists("gse_new")) {
  gse_new <- getGEO("GSE193677", GSEMatrix = TRUE)
}
all_metadata <- pData(gse_new[[1]])

# Step 3: Extract and Filter Valid Clinical Cohorts
clean_conditions <- make.names(all_metadata$`ibd_clinicianmeasure_inactive_active:ch1`)
valid_indices    <- which(clean_conditions == "Active" | clean_conditions == "Inactive")

filtered_metadata <- all_metadata[valid_indices, ]
filtered_metadata$condition <- factor(clean_conditions[valid_indices])

# Step 4: Access and Parse Supplementary Expression Matrices
message("Accessing downloaded supplementary files...")
supp_files <- getGEOSuppFiles("GSE193677")

# Target the adjusted/normalized expression count matrix
file_path <- rownames(supp_files)[1]

# Use read.table with sep="" to safely handle any variable whitespace/tab setups
raw_counts <- read.table(file_path, header = TRUE, row.names = 1, check.names = FALSE, sep = "")

# Step 5: DATA ENGINEERING - Resolve Metadata vs Matrix Header Mismatch
# Strips the trailing strings (e.g. ", UC participants...") from metadata titles
short_titles <- sub(",.*$", "", filtered_metadata$title)
cleaned_metadata_titles <- make.names(short_titles)

# Step 6: Synchronize and Align Matrix Dimensions
common_in_counts   <- which(colnames(raw_counts) %in% cleaned_metadata_titles)
numeric_counts     <- as.matrix(raw_counts[, common_in_counts, drop = FALSE])
mode(numeric_counts) <- "numeric"

# Re-align sample metadata rows to match the matrix columns exactly
metadata_indices  <- match(colnames(numeric_counts), cleaned_metadata_titles)
filtered_metadata <- filtered_metadata[metadata_indices, ]

# Step 7: Low-Variance Gene Filtering (Remove non-informative features)
gene_variances <- apply(numeric_counts, 1, var)
numeric_counts <- numeric_counts[gene_variances > 0, , drop = FALSE]

# Print execution checks to console
print(paste("Verified Synchronized Genes:", nrow(numeric_counts)))
print(paste("Verified Synchronized Samples:", ncol(numeric_counts)))

# Step 8: Build Design Matrix and Fit Linear Models
design <- model.matrix(~ 0 + condition, data = filtered_metadata)
colnames(design) <- levels(filtered_metadata$condition)

fit <- lmFit(numeric_counts, design)

# Step 9: Compute Contrasts and Apply Empirical Bayes Moderation
contrast_matrix <- makeContrasts(Active_vs_Inactive = Active - Inactive, levels = design)
fit2  <- contrasts.fit(fit, contrast_matrix)
fit2  <- eBayes(fit2)

# Step 10: Extract Top 100 Differentially Expressed Genes
top_genes <- topTable(fit2, coef = "Active_vs_Inactive", number = 100, adjust.method = "BH")

# Save full results spreadsheet for portfolio tracking
write.csv(top_genes, "Top_IBD_Active_vs_Inactive_Genes.csv")
print("Top 10 Differentially Expressed Genes:")
print(head(top_genes, 10))

# Step 11: Render and Export Publication-Quality Volcano Plot
message("Generating and saving portfolio volcano plot...")
if(!dir.exists("figures")) dir.create("figures")

# Define thresholds for coloring
top_genes$Significance <- "Not Significant"
top_genes$Significance[top_genes$adj.P.Val < 0.05 & top_genes$logFC > 0.3] <- "Upregulated in Active"
top_genes$Significance[top_genes$adj.P.Val < 0.05 & top_genes$logFC < -0.3] <- "Downregulated in Active"

p <- ggplot(top_genes, aes(x = logFC, y = -log10(P.Value), color = Significance)) +
  geom_point(alpha = 0.8, size = 2) +
  scale_color_manual(values = c("Not Significant" = "grey60", 
                                "Upregulated in Active" = "firebrick3", 
                                "Downregulated in Active" = "dodgerblue3")) +
  theme_minimal(base_size = 14) +
  labs(title = "Volcano Plot: Active vs Inactive Clinical IBD States",
       subtitle = "Dataset: GSE193677 | 1,832 Intestinal Biopsies",
       x = "Log2 Fold Change",
       y = "-Log10 Raw P-Value") +
  theme(plot.title = element_text(face = "bold", size = 16),
        legend.position = "bottom")

ggsave("figures/02_ibd_volcano_plot.png", plot = p, width = 8, height = 6, dpi = 300)
message("Pipeline execution complete! Check your project folder for generated assets.")
