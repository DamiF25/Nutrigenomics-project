#install required packages
install.packages("DESeq2")
install.packages("tidyverse")
install.packages("ggplot2")
install.packages("dplyr")
install.packages("reshape2")
install.packages("pheatmap")
install.packages("vegan")
BiocManager::install()
install.packages("BiocManager")
BiocManager::install("ComplexHeatmap", version = "3.18")  # Bioconductor release for R 4.3
install.packages("circlize")
install.packages("grid")
install.packages("compositions")
install.packages("zCompositions")

#Load libraries
library(DESeq2)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(reshape2)
library(pheatmap)
library(vegan)
library(RColorBrewer)
library(BiocManager)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(compositions)
library(zCompositions)

#Set work directory

# THE GOAL HERE WAS TO EXPOLRE CORRELATIONS BETWEEN MICROBIOME OTU MATRIX AND DESeq2 MATRIX

# From the normalised "otu_table_tss.csv" file, we extracted the relatively  abundant (RA) taxa,
# the differentially abundant (DA) taxa, indicator taxa, and other taxa of interest into a 
# separate file. The data was renormalised using total sum scaling (TSS)

# Load the prepared OTU table 
sig_otu_table <- read.csv("sig_otu_taxa_tss.csv", row.names = 1, check.names = FALSE) 

# Make sure all columns are numeric
sig_otu_table[] <- lapply(sig_otu_table, as.numeric)

# Renormalise so each sample column sums to 1
sig_otu_renorm <- sweep(sig_otu_table, 2, colSums(sig_otu_table, na.rm = TRUE), FUN = "/")

# Check sums
colSums(sig_otu_renorm)

View(sig_otu_renorm)

# Save file
write.csv(sig_otu_renorm, "renorm_sig_otu_taxa_tss.csv")

# Transform the TSS-normalised OTU table using centered Log-Ratio (CLR)
# Load TSS-transformed OTU table
renorm_sig_otu_tss <- read.csv("renorm_sig_otu_taxa_tss.csv", row.names = 1, check.names = FALSE)

# Zero-replacement using multiplicative method
otu_nozeros <- cmultRepl(renorm_sig_otu_tss, label = 0, method ="CZM") # CZM = count zero multiplicative

# Check for zeros
sum(otu_nozeros == 0)  # Should be 0

# Check minimum value
min(otu_nozeros)       # Should be a small positive number

# Then warnings, if any, can be ignored for correlation analysis

# CLR transformation
clr(otu_nozeros)

# Convert to data frame and add sample
otu_clr <- clr(otu_nozeros)
otu_clr_df <- as.data.frame(otu_clr)
otu_clr_df$SampleID <- rownames(otu_clr_df)

# Save to file
write.csv(otu_clr_df, "clr_transformed_otu.csv", row.names = FALSE)

# Load CLR-transformed data
clr_df <- read.csv("clr_transformed_otu.csv", header = TRUE, stringsAsFactors = FALSE)

head(clr_df)

# Move SampleID (last column) to the first position 
clr_df <- clr_df[, c(ncol(clr_df), 1:(ncol(clr_df) - 1))]

# Rename first column to 'Taxon'
colnames(clr_df)[1] <- "Taxon"

# Preview structure and first few columns
str(clr_df)
head(clr_df[, 1:5])

# Save reordered data to overwrite the previous clr-transformed table
write.csv(clr_df, "clr_transformed_otu.csv", row.names = FALSE)

# OPTIONAL: view density plot of CLR-transformed otu values

# Reshape to long format
clr_long <- pivot_longer(otu_clr_df, -SampleID, names_to = "Taxon", values_to = "CLR")

# Density plot
ggplot(clr_long, aes(x = CLR)) +
  geom_density(fill = "steelblue", alpha = 0.6) +
  labs(title = "Density Plot of CLR Values",
       x = "CLR Value",
       y = "Density") +
  theme_minimal()

# Prepare variance-stablised transformed DEG matrix from DESeq2
#Load gene count data
count_data <- read.csv("Gene_count_matrix.csv",header=TRUE,row.names=1)

#Load gene sample information
sample_info <- read.csv("PHENO_DATA.csv",header=TRUE,row.name=1) 

#Confirm that sample names match in both count data and sample_info files. Should be TRUE
all(colnames(count_data) %in% rownames(sample_info))

#Confirm that they are also in the same order. Should be TRUE
all(colnames(count_data) == rownames(sample_info))

#Construct DESeq dataset object by linking both the count_data and the sample_info files to DESeq
dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = sample_info, design = ~ Diet)
nrow(dds)

#Data pre-filtering: remove lowly expressed genes
smallestGroupSize <-5
keep <- rowSums(counts(dds) >=10) >= smallestGroupSize

dds <- dds[keep,]

nrow(dds)

#Set factor levels
factors <- factor(sample_info$Diet)
groups <- unique(sample_info$Diet)

#To set FL as the reference (control)
dds$Diet <- relevel(dds$Diet, ref = "FL")
dds$Diet

#Run DESeq2 using Wald's tests
#Perform statistical tests
#Wald test for differential gene expression 
#poscount for estimates of size factor for normalisation
dds <- DESeq(dds,test="Wald",sfType='poscount')

#Fetch the DESeq result
res <- results(dds)

#Check DESeq result at p <.01
res <- results(dds, alpha = 0.01)
summary(res)

#Extract DE genes with padj < 0.01 and log2Foldchange <=-1 or >=1
deg <- subset(res, padj<0.01 & abs(log2FoldChange) >=1)
summary(deg)

#Remove NA from deg - OPTIONAL
deg <- na.omit(deg)
summary(deg)

#Change DESeq result into a dataframe
res <- as.data.frame(res)
class(res)

#Variance stabilising transformation
vst <- vst(dds,blind=FALSE)

# Get significant DEGs
deg_list <- rownames(res[which(res$padj < 0.01 & abs(res$log2FoldChange) >= 1), ])

# Extract VST counts
vst_counts <- assay(vst)

# Subset only DEGs
vst_deg_counts <- vst_counts[deg_list, ]

# Convert to data frame with gene names as a column
vst_deg_counts_df <- data.frame(Gene = rownames(vst_deg_counts), vst_deg_counts)

nrow(vst_deg_counts_df)

# Save to CSV
write.csv(vst_deg_counts_df, "vst_DEG_p0.01.csv", row.names = FALSE)

# Load CLR-transformed OTU matrix and VST DEG data
taxa <- read.csv("clr_transformed_otu.csv", row.names = 1, check.names = FALSE)
deg  <- read.csv("FL_vs_FF_Reproduction_vst_DEGs.csv", row.names = 1, check.names = FALSE)

# Align samples (columns) across both datasets
common_samples <- intersect(colnames(taxa), colnames(deg))
taxa_aligned <- taxa[, common_samples]
deg_aligned  <- deg[, common_samples]

# Transpose matrices so samples become rows
clr_mat <- t(taxa_aligned)
deg_mat <- t(deg_aligned)

# Compute Euclidean distances between samples
dist_taxa <- dist(clr_mat, method = "euclidean")
dist_deg  <- dist(deg_mat, method = "euclidean")

# Run Mantel test
mantel_result <- mantel(dist_taxa, dist_deg, method = "spearman", permutations = 999)
print(mantel_result)

# Visualise Mantel correlation in a scatter plot
mantel_df <- data.frame(
  microbiome = as.vector(dist_taxa),
  transcriptome = as.vector(dist_deg)
)

ggplot(mantel_df, aes(x = microbiome, y = transcriptome)) +
  geom_point(size = 3.0, alpha = 0.5) +
  geom_smooth(method = "lm", color = "blue", se = TRUE) +
  labs(
    title = paste0("Mantel Test (General): r = ", round(mantel_result$statistic, 3),
                   ", p = ", mantel_result$signif),
    x = "Taxa distance (CLR)",
    y = "DEG distance (VST)"
  ) +
  theme_minimal(base_size = 14) +  # Increase base font size
  theme(
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
    axis.title.x = element_text(face = "bold", size = 16),
    axis.title.y = element_text(face = "bold", size = 16),
    axis.text = element_text(size = 14)
  )

# Save plot
ggsave("Transcriptome_Microbiome_Mantel_correlation.png", width = 6.2, height = 6, dpi = 300)


# ==== SIMILAR APPROACHES WERE USED FOR FL_vs_FF ANALYSES ====

#Load gene count data
count_data <- read.csv("FL_vs_FF_Gene_count_matrix.csv",header=TRUE,row.names=1)

#Load gene sample information
sample_info <- read.csv("FL_vs_FF_PHENO_DATA.csv",header=TRUE,row.name=1) 

#Confirm that sample names match in both count data and sample_info files. Should be TRUE
all(colnames(count_data) %in% rownames(sample_info))

#Confirm that they are also in the same order. Should be TRUE
all(colnames(count_data) == rownames(sample_info))

#Construct DESeq dataset object by linking both the count_data and the sample_info files to DESeq
dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = sample_info, design = ~ Diet)
nrow(dds)

#Data pre-filtering: remove lowly expressed genes
smallestGroupSize <-5
keep <- rowSums(counts(dds) >=10) >= smallestGroupSize

dds <- dds[keep,]

nrow(dds)

#Set factor levels
factors <- factor(sample_info$Diet)
groups <- unique(sample_info$Diet)

#To set FL as the reference (control)
dds$Diet <- relevel(dds$Diet, ref = "FL")
dds$Diet

#Run DESeq2 using Wald's tests
#Perform statistical tests
#Wald test for differential gene expression 
#poscount for estimates of size factor for normalisation
dds <- DESeq(dds,test="Wald",sfType='poscount')

#Fetch the DESeq result
res <- results(dds)

#Check DESeq result at p <.01
res <- results(dds, alpha = 0.01)
summary(res)

#Extract DE genes with padj < 0.01 and log2Foldchange <=-1 or >=1
deg <- subset(res, padj<0.01 & abs(log2FoldChange) >=1)
summary(deg)

#Remove NA from deg - OPTIONAL
deg <- na.omit(deg)
summary(deg)

#Change DESeq result into a dataframe
res <- as.data.frame(res)
class(res)

#Variance stabilizing transformation
vst <- vst(dds,blind=FALSE)

# Get significant DEGs
deg_list <- rownames(res[which(res$padj < 0.01 & abs(res$log2FoldChange) >= 1), ])

# Extract VST counts
vst_counts <- assay(vst)

# Subset only DEGs
vst_deg_counts <- vst_counts[deg_list, ]

# Convert to data frame with gene names as a column
vst_deg_counts_df <- data.frame(Gene = rownames(vst_deg_counts), vst_deg_counts)

nrow(vst_deg_counts_df)

# Save to CSV
write.csv(vst_deg_counts_df, "FL_vs_FF_vst_p0.01.csv", row.names = FALSE)

# From the "FL_vs_FF_vst_p.01.csv" file, an "FL_vs_FF_Growth_vst_DEGs.csv" file, containing identified 
# growth-related genes, was prepared 

# Again, from the "sig_otu_taxa_tss.csv" OTU table, an "FL_vs_FF_sig_otu_taxa_tss.csv" file, containing OTU 
# for FL and FF samples, was prepared

# Read the prepared OTU table while keeping only numeric part
otu <- read.csv("FL_vs_FF_sig_otu_taxa_tss.csv", row.names = 1, check.names = FALSE)

# Make sure all columns are numeric
otu[] <- lapply(otu, as.numeric)

# Renormalise so each sample column sums to 1
otu_renorm <- sweep(otu, 2, colSums(otu, na.rm = TRUE), FUN = "/")

# Check sums
colSums(otu_renorm)

# Save file
write.csv(otu_renorm, "renorm_FL_vs_FF_sig_otu_taxa.tss.csv")

# Load TSS-transformed for FLvFF OTU table
FF_otu_table_tss <- read.csv("renorm_FL_vs_FF_sig_otu_taxa.tss.csv", row.names = 1, check.names = FALSE)

# Zero-replacement using multiplicative method
FF_otu_nozeros <- cmultRepl(FF_otu_table_tss, label = 0, method ="CZM") # CZM = count zero multiplicative

# Check for zeros
sum(FF_otu_nozeros == 0)  # Should be 0

# Check minimum value
min(FF_otu_nozeros)       # Should be a small positive number

# If the above are true, any warning can be ignored for correlation analysis

# CLR transformation
clr(FF_otu_nozeros)

# Convert to data frame and add sample
FF_otu_clr <- clr(FF_otu_nozeros)
FF_otu_clr_df <- as.data.frame(FF_otu_clr)
FF_otu_clr_df$SampleID <- rownames(FF_otu_clr_df)

# Save to file
write.csv(FF_otu_clr_df, "clr_transformed_otu_FFvFL.csv", row.names = FALSE)

# Load CLR-transformed data
FF_clr_df <- read.csv("clr_transformed_otu_FFvFL.csv", header = TRUE, stringsAsFactors = FALSE)

# Move SampleID (last column) to the first position 
FF_clr_df <- FF_clr_df[, c(ncol(FF_clr_df), 1:(ncol(FF_clr_df) - 1))]

# Rename first column to 'Taxon'
colnames(FF_clr_df)[1] <- "Taxon"

# Preview structure and first few columns
str(FF_clr_df)
head(FF_clr_df[, 1:5])

# Save reordered data to overwrite the previous clr-transformed table
write.csv(FF_clr_df, "clr_transformed_otu_FFvFL.csv", row.names = FALSE)

# Load transformed OTU and DEG datasets
deg <- read.csv("FL_vs_FF_vst_p0.01.csv", row.names = 1, check.names = FALSE) 
taxa <- read.csv("clr_transformed_otu_FFvFL.csv", row.names = 1, check.names = FALSE) 

# If "Group" column exists in deg, separate it
if ("Group" %in% colnames(deg)) {
  gene_groups <- deg$Group
  deg <- deg[ , !(colnames(deg) %in% "Group")]
} else {
  gene_groups <- rep("Unknown", nrow(deg))
}

# Ensure same sample order
samples <- intersect(colnames(deg), colnames(taxa))
deg <- deg[ , samples]
taxa <- taxa[ , samples]

# Mantel test
mantel_res <- mantel(
  dist(t(deg)), 
  dist(t(taxa)), 
  method = "spearman", permutations = 999
)

mantel_r <- round(mantel_res$statistic, 3)
mantel_p <- mantel_res$signif

cat("Mantel r =", mantel_r, "p =", mantel_p, "\n")

# Mantel test scatter plot
# Compute distance matrices
dist_deg <- dist(t(deg))
dist_taxa <- dist(t(taxa))

# Convert to vectors (upper triangular only, to avoid duplicates)
deg_vec <- as.vector(dist_deg)
taxa_vec <- as.vector(dist_taxa)

# Scatter plot of distances
mantel_df <- data.frame(
  DEG_dist = deg_vec,
  Taxa_dist = taxa_vec
)

# Plot
ggplot(mantel_df, aes(x = Taxa_dist, y = DEG_dist)) +
  geom_point(size = 3.0, alpha = 0.5) +
  geom_smooth(method = "lm", color = "blue", se = TRUE) +
  labs(
    title = paste0("Mantel Test (Growth): r = ", round(mantel_r, 3), 
                   ", p = ", signif(mantel_p, 3)),
    x = "Taxa distance (CLR)",
    y = "DEG distance (VST)"
  ) +
  theme_minimal(base_size = 14) + 
  theme(
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
    axis.title.x = element_text(face = "bold", size = 16),
    axis.title.y = element_text(face = "bold", size = 16),
    axis.text = element_text(size = 14)
  )

# Save plot
ggsave("Mantel_Gene-Taxon_correlation_Growth.png", width = 6.2, height = 6, dpi = 300)

# Perform pairwise correlations
cor_results <- expand.grid(Gene = rownames(deg), Taxon = rownames(taxa)) %>%
  mutate(
    corr = NA,
    pval = NA
  )

for (i in 1:nrow(cor_results)) {
  g <- cor_results$Gene[i]
  t <- cor_results$Taxon[i]
  test <- cor.test(as.numeric(deg[g, ]), as.numeric(taxa[t, ]),
                   method = "spearman")
  cor_results$corr[i] <- test$estimate
  cor_results$pval[i] <- test$p.value
}

# FDR correction
cor_results$padj <- p.adjust(cor_results$pval, method = "fdr")

# Create matrix of correlations
cor_matrix <- acast(cor_results, Gene ~ Taxon, value.var = "corr")

# Create matrix of significance labels
sig_labels <- acast(cor_results, Gene ~ Taxon, value.var = "padj")
sig_labels <- ifelse(sig_labels < 0.001, "***",
                  ifelse(sig_labels < 0.01, "**",
                     ifelse(sig_labels < 0.05, "*", "")))

# Assign colors for gene groups (if DEG table contains a Group column)
group_levels <- unique(gene_groups)
group_colors <- c(
  "PIP-associated" = "#000000", 
  "Myosin heavy chains" = "#0000FF",
  "Others" = "#FF0000",
  "Dynein chains" = "#B5B5B5",
  "Digestive enzymes" = "#ffb600"
  )

# Row annotation for gene groups
row_ha <- rowAnnotation(
  Group = gene_groups,
  col = list(Group = group_colors),
  annotation_legend_param = list(
    title = "Gene groups",
    title_gp = gpar(fontsize = 16, fontface = "bold"),  # legend title size/bold
    labels_gp = gpar(fontsize = 14),                     # legend labels size
    grid_height = unit(6, "mm")
  ),
  show_annotation_name = FALSE   # <---- hides the "Group" text
)

# Color function for correlation values
col_fun <- colorRamp2(c(-1, 1), c("#2ca8f4", "white"))  

# Mantel results
mantel_r <- 0.74 # add actual value
mantel_p <- 0.001 # add actual value

# Create the heatmap
ht <- Heatmap(
  cor_matrix,
  name = "Spearman's r",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  column_title = "FL v. FF growth gene-taxon correlation",
  column_title_gp = gpar(fontsize = 19, fontface = "bold"),  # column title size/bold
  column_names_gp = gpar(fontsize = 15),  # column label size
  column_names_rot = 55, # slant taxa name
  row_split = gene_groups,
  row_title = NULL,
  row_title_side = "left",
  left_annotation = row_ha,
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 16, fontface = "bold"),  # heatmap legend title size/bold
    labels_gp = gpar(fontsize = 14),               # heatmap legend labels size
    grid_height = unit(9, "mm")),  # increase box height/spacing
  # Add significance annotations inside heatmap cells
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (!is.na(sig_labels[i, j]) && sig_labels[i, j] != "") {
      grid.text(sig_labels[i, j], x, y, gp = gpar(fontsize = 14, col = "black"))
    }
  }
)

# Draw the heatmap
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")

# Add Mantel info to top-right under title
grid.text(
  paste0("Mantel r = ", round(mantel_r, 3), 
         ", p = ", signif(mantel_p, 3)),
  x = unit(0.98, "npc"),   # Right side
  y = unit(0.96, "npc"),   # Just below column title
  just = c("right", "top"),
  gp = gpar(fontsize = 14, fontface = "italic", col = "black")
)

# Export image: #1000W x 1100H (1000W x 900H for Reproduction heatmap)
