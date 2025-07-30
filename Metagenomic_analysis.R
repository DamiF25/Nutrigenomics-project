#install if required!

install.packages("dplyr")
install.packages("readr")
install.packages("purrr")
install.packages("stringr")
install.packages("tidyr")
install.packages("tidyverse")
install.packages("ggplot2")
install.packages("pheatmap")
install.packages("ape") # for PCoA
install.packages("rstatix") # for pairwise_wilcox_test
install.packages("ggpubr")  # for p-value annotation
install.packages("vegan") # for PERMANOVA
install.packages("indicspecies") # for multipatt
install.packages("forcats")

library(dplyr)
library(readr)
library(purrr)
library(stringr)
library(tidyr)
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(ape)
library(rstatix)
library(ggpubr)
library(vegan)
library(indicspecies)
library(forcats)

# Set working directory

# Bracken was performed using Galaxy.
# Galaxy Bracken result files are labelled as ".tabular" and often have no file extension even 
# though they contain TSV-formatted data. Renaming them manually to ".tsv" for readability 
# (e.g., in Excel) is fine locally, but R doesn't detect them unless they have the proper 
# .tsv extension in the file name.


# Recursively rename all .tabular files to .tsv in the Bracken reports folder
base_path <- "Bracken_reports"

tabular_files <- list.files(base_path, pattern = "\\.tabular$", recursive = TRUE, full.names = TRUE)

for (file in tabular_files) {
  new_name <- sub("\\.tabular$", ".tsv", file)
  file.rename(file, new_name)
}

# PREPARE OTU TABLE AND METADATA FOR ALPHA DIVERSITY AND BETA DIVERSITY

# Setup folder paths
base_path <- "Bracken_reports"  # Go to the parent folder
group_dirs <- list(
  FF_group = file.path(base_path, "FF_group"),
  FL_group = file.path(base_path, "FL_group"),
  SG_group = file.path(base_path, "SG_group")
)

# Read, filter, normalize Bracken files 
read_group_files <- function(path, group) {
  files <- list.files(path, pattern = "\\.tsv$", full.names = TRUE)
  
  dfs <- lapply(files, function(file) {
    df <- read.delim(file, sep = "\t", header = TRUE)
    
    # Remove "Mammalia" (case-insensitive)
    df <- df[trimws(tolower(df$name)) != "mammalia", ]
    
    # Normalize fraction_total_reads to sum to 1
    df$fraction_total_reads <- df$fraction_total_reads / sum(df$fraction_total_reads)
    
    # Keep only taxon name and normalized values
    sample_id <- tools::file_path_sans_ext(basename(file))
    df <- df[, c("name", "fraction_total_reads")]
    colnames(df)[2] <- sample_id
    return(df)
  })
  
  names(dfs) <- rep(group, length(dfs))
  return(dfs)
}

# Apply the function to all groups
bracken_list <- lapply(names(group_dirs), function(g) read_group_files(group_dirs[[g]], g))
all_dfs <- do.call(c, bracken_list)

# Merge all files by 'name' (taxon/class)
merged_df <- Reduce(function(x, y) full_join(x, y, by = "name"), all_dfs)
merged_df[is.na(merged_df)] <- 0  # Fill missing taxa with 0

# Set row names to taxa name and prepare OTU matrix
rownames(merged_df) <- merged_df$name
otu_table <- as.data.frame(merged_df[, -1])
otu_table <- as.data.frame(t(otu_table))  # Samples as rows

view(otu_table)
# Save OTU table (samples as rows, taxa as columns)
write.csv(otu_table, file = "otu_table_class_normalized.csv", quote = FALSE)

# Create metadata with group labels
sample_ids <- rownames(otu_table)
sample_groups <- names(all_dfs)  # Group name was used as list name
sample_metadata <- data.frame(SampleID = sample_ids, Group = sample_groups)
rownames(sample_metadata) <- sample_metadata$SampleID
view(sample_metadata)

# HEATMAP OF ALL TAXA

# Load OTU table with custom row names
otu_table <- read.csv("otu_table_class_normalized.csv", row.names = 1, check.names = FALSE) 
# There was a need to upload the OTU table ONLY because the original row names have been modified manually for cleaner labelling.

# Update row names of sample_metadata to match OTU table
sample_ids <- rownames(otu_table)
sample_groups <- names(all_dfs)  
sample_metadata <- data.frame(SampleID = sample_ids, Group = sample_groups)
rownames(sample_metadata) <- sample_metadata$SampleID
view(sample_metadata)

# Confirm row names
head(rownames(otu_table))         
head(rownames(sample_metadata))   # Must match above exactly

# Reorder metadata to match OTU table
sample_metadata_heatmap <- sample_metadata[rownames(otu_table), , drop = FALSE]

# Now annotate properly
annotation_row <- data.frame(Group = sample_metadata_heatmap$Group)
rownames(annotation_row) <- rownames(otu_table)

all(rownames(otu_table) == rownames(annotation_row))  # Should return TRUE

# Transpose the OTU table so taxa are in rows, samples in columns
otu_flipped <- t(scale(otu_table))

# Create annotation with renamed group labels
annotation_col <- data.frame(
  `Diet` = recode(sample_metadata_heatmap$Group,
                  "FF_group" = "FF",
                  "FL_group" = "FL",
                  "SG_group" = "SG"))
rownames(annotation_col) <- rownames(sample_metadata_heatmap) # Sample IDs

# Plot and save high-resolution heatmap
tiff("heatmap_all_otus.tif", 
     width = 3500, height = 9000, res = 400, compression = "lzw")

pheatmap(
  mat = otu_flipped,
  annotation_col = annotation_col,
  annotation_colors = list('Diet' = c(
    "FF" = "#990906",
    "FL" = "#000000",
    "SG" = "#0a5fb5")),
  show_rownames = TRUE,
  show_colnames = TRUE,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  fontsize = 14,
  colorRampPalette(c("black", "white", "orange"))(100))

dev.off()

# PERFORM ALPHA DIVERSITY

# Calculate diversity metrics using the vegan package
shannon <- diversity(otu_table, index = "shannon")
simpson <- diversity(otu_table, index = "simpson")
richness <- specnumber(otu_table)  # observed number of taxa (class level)

# Chao1 (requires estimateR)
# NOTE: Use raw or integer OTU instead of normalised OTU data
# If only normalised OTU data is availabele, multiply each value by a constant factor 
# (e.g., 10,000). This scales the relative abundances into pseudo-counts which can be used for Chao1
otu_counts <- round(otu_table * 10000)

# Confirm the new data is numeric and integer-like
str(otu_counts)
summary(rowSums(otu_counts))  # Will vary by sample

chao1_all <- estimateR(otu_counts)  # rows = metrics, cols = samples
chao1 <- chao1_all["S.chao1", ]     # Extract Chao1 only

# Combine all metrics
alpha_df <- data.frame(
  SampleID = rownames(otu_counts),
  Shannon = shannon,
  Simpson = simpson,
  Richness = richness,
  Chao1 = chao1)

# Merge with metadata
alpha_df <- left_join(alpha_df, sample_metadata, by = "SampleID")

# Format group labels. Ensure group is a factor with correct order
alpha_df$Group <- factor(alpha_df$Group, levels = c("FF_group", "FL_group", "SG_group"))

# Reshape to long format for faceting
alpha_long <- alpha_df %>%
  pivot_longer(cols = c(Shannon, Simpson, Richness, Chao1),
               names_to = "Metric",
               values_to = "Value") %>%
  mutate(
    Metric = factor(Metric,
                    levels = c("Shannon", "Simpson", "Richness", "Chao1"),
                    labels = c("Shannon Diversity", "Simpson Diversity", "Observed Richness", "Chao1 Richness")),
    Group = factor(Group, labels = c("FF", "FL", "SG"))
  )
# Pairwise Wilcoxon + clean p labels (e.g., *, **)
pairwise_alpha <- alpha_long %>%
  group_by(Metric) %>%
  pairwise_wilcox_test(Value ~ Group, p.adjust.method = "BH") %>%
  filter(p.adj < 0.05) %>%
  mutate(p.label = p_format(p.adj)) %>%
  ungroup()

print(pairwise_alpha)

# Add dynamic y-position only if results exist
if (nrow(pairwise_alpha) > 0) {
  y_max_lookup <- alpha_long %>%
    group_by(Metric) %>%
    summarise(max_val = max(Value, na.rm = TRUE), .groups = "drop")
  
  pairwise_alpha <- pairwise_alpha %>%
    left_join(y_max_lookup, by = "Metric") %>%
    group_by(Metric) %>%
    mutate(y.position = max_val * 1.05 + row_number() * 0.05) %>%
    ungroup() %>%
    select(-max_val)
}

# Facet box plot
ggplot(alpha_long, aes(x = Group, y = Value, fill = Group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.85, width = 0.4) +
  geom_jitter(width = 0.2, size = 1.6, alpha = 0.8) +
  facet_wrap(~ Metric, scales = "free_y") +
  {
    if (nrow(pairwise_alpha) > 0) {
      stat_pvalue_manual(
        pairwise_alpha,
        label = "p.label",
        y.position = "y.position",
        tip.length = 0.02,
        size = 5
      )
    } else {
      NULL
    }
  } +
  labs(x = "Diet Group", y = "Alpha Diversity Value") +
  scale_fill_manual(values = c("FF" = "#f00000", "FL" = "#000000", "SG" = "#0f00f7")) +
  theme_test() +
  theme(
    strip.text = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    axis.title.x = element_blank(),
    axis.text = element_text(size = 16),
    legend.position = "none"
  )

ggsave("alpha_diversity_facet_plot.tif", width = 7, height = 5, dpi = 300)
write.csv(alpha_df, "alpha_diversity_metrics.csv", row.names = FALSE)
write.csv(pairwise_alpha, "alpha_diversity_pairwise_wilcoxon.csv", row.names = FALSE) # If significance was observed

# PERFROM BETA-DIVERSITY ANALYSIS

# --------- STEP 6: Bray-Curtis Distance & Ordination ---------
bray_dist <- vegdist(otu_table, method = "bray")
pcoa_result <- cmdscale(bray_dist, eig = TRUE, k = 2)
pcoa_df <- as.data.frame(pcoa_result$points)
pcoa_df$SampleID <- rownames(pcoa_df)
pcoa_df <- left_join(pcoa_df, sample_metadata, by = "SampleID")

# Convert distance matrix to square matrix format for saving
bray_matrix <- as.matrix(bray_dist)
write.csv(bray_matrix, file = "bray_curtis_distance_matrix.csv", quote = FALSE)


# --------- STEP 7: PCoA Plot ---------
# Run PCoA with proper eigenvalue tracking
pcoa_result <- pcoa(bray_dist)

# Extract % variance explained
var_explained <- pcoa_result$values$Relative_eig * 100

view(var_explained)

# Create the plot data frame
pcoa_df <- as.data.frame(pcoa_result$vectors[, 1:2])
pcoa_df$SampleID <- rownames(pcoa_df)
pcoa_df <- left_join(pcoa_df, sample_metadata, by = "SampleID")

# Rename groups before plotting
pcoa_df$Group <- factor(pcoa_df$Group,
                        levels = c("FF_group", "FL_group", "SG_group"),
                        labels = c("FF", "FL", "SG"))

# Plot
ggplot(pcoa_df, aes(x = Axis.1, y = Axis.2, color = Group)) +
  geom_point(size = 4) +
  theme_bw() +
  labs(title = "Beta Diversity (Bray-Curtis index)", x = paste0("PCoA1 (", round(var_explained[1], 1), "%)"),
       y = paste0("PCoA2 (", round(var_explained[2], 1), "%)"),
       color = "Diet") +
  scale_color_manual(values = c(
    "FF" = "#f00000",
    "FL" = "#000000",
    "SG" = "#0f00f7"  
  )) +
  coord_fixed(ratio=1) + theme(panel.background = element_rect (),
                               panel.grid = element_line () ) +
  stat_ellipse(geom = "polygon", aes(fill = after_scale(alpha(colour, 0.3))),
               data = pcoa_df[pcoa_df$Group !="versicolor",]) +
  theme(plot.title = element_text(size = 20, face = "bold"),
        axis.title = element_text(size = 18, face = "bold"),
        axis.text = element_text(size = 16)) +
  theme(legend.text = element_text(size = 16), legend.title = element_text(size = 18, face = "bold"))

# Save the final figure
ggsave("PCoA (Bray-Curtis) on Class-Level Bracken Profiles.tif", width = 8, height = 6.43, dpi = 300)

# PERMANOVA
adonis_result <- adonis2(bray_dist ~ Group, data = sample_metadata)
print(adonis_result)

# Intergroup PERMANOVA
pairwise_permanova <- function(dist_matrix, metadata, group_col, n_perm = 999) {
  groups <- unique(metadata[[group_col]])
  results <- data.frame()
  
  for (i in 1:(length(groups) - 1)) {
    for (j in (i + 1):length(groups)) {
      group_i <- groups[i]
      group_j <- groups[j]
      
      # Subset metadata and distance matrix
      subset_metadata <- metadata[metadata[[group_col]] %in% c(group_i, group_j), ]
      subset_ids <- rownames(subset_metadata)
      subset_dist <- as.dist(as.matrix(dist_matrix)[subset_ids, subset_ids])
      
      # Create the formula dynamically
      formula <- as.formula(paste0("subset_dist ~ ", group_col))
      
      # Run PERMANOVA
      ad <- adonis2(formula, data = subset_metadata, permutations = n_perm)
      
      # Collect results
      results <- rbind(results, data.frame(
        Group1 = group_i,
        Group2 = group_j,
        R2 = ad$R2[1],
        F = ad$F[1],
        p.value = ad$`Pr(>F)`[1]
      ))
    }
  }
  
  results$p.adj <- p.adjust(results$p.value, method = "bonferroni")
  return(results)
}

pairwise_results <- pairwise_permanova(bray_dist, sample_metadata, "Group")
print(pairwise_results)

# Run a beta-dispersion test to validate PERMANOVA assumptions
dispersion <- betadisper(bray_dist, sample_metadata$Group)
anova(dispersion)

plot(dispersion)
boxplot(dispersion, main = "Group Dispersions", ylab = "Distance to Group Centroid")

# TO PERFORM STATISTICAL ANALYSIS FOR MOST SIGNIFICANTLY ABUNDANT TAXA ACROSS GROUPS

# Convert the bracken reports into long data format
base_dir <- "Bracken_reports"
groups <- c("FF_group", "FL_group", "SG_group")

# Function to read, clean, and normalize a Bracken file
read_bracken_long <- function(file, group_name) {
  sample_id <- tools::file_path_sans_ext(basename(file))
  
  df <- read_tsv(file, col_types = cols()) %>%
    select(name, taxonomy_id, taxonomy_lvl, new_est_reads) %>%
    filter(taxonomy_lvl == "C", name != "Mammalia") %>%
    mutate(class_name = name) %>%
    select(-name) %>%
    mutate(sample_id = sample_id, group = group_name) %>%
    
    # Normalize per sample_id
    group_by(sample_id) %>%
    mutate(relative_abundance = new_est_reads / sum(new_est_reads)) %>%
    ungroup() %>%
    select(sample_id, group, class_name, relative_abundance)
  
  return(df)
}

# Loop through all files
bracken_long <- map_dfr(groups, function(group_name) {
  files <- list.files(file.path(base_dir, group_name), pattern = "\\.tsv$", full.names = TRUE)
  map_dfr(files, read_bracken_long, group_name = group_name)
})

# Sanity check — y-values should sum to 1 per sample
check <- bracken_long %>%
  group_by(sample_id) %>%
  summarise(sum_abundance = sum(relative_abundance))

print(check)  # All should be ~1.0

# Run Kruskal–Wallis test
kw_results <- bracken_long %>%
  group_by(class_name) %>%
  summarise(
    p_value = tryCatch(
      kruskal.test(relative_abundance ~ group)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  filter(!is.na(p_value)) %>%
  mutate(adj_p = p.adjust(p_value, method = "BH")) %>%
  filter(adj_p < 0.05)

view(kw_results)

write.csv(kw_results, "Statistical data for significantly different classes.csv")

significant_taxa <- kw_results$class_name

filtered_data <- bracken_long %>%
  filter(class_name %in% significant_taxa)

# OPTIONAL, Re-normalize so that only the significant taxa per sample sum to 1.0. This gives
# stacked bars of equal heights. In this study, this option was skipped, and legend sorting was next.
# NOTE: If you use this option, you MUST follow the instruction in subsequent "IMPORTANT NOTE".
filtered_data <- filtered_data %>%
  group_by(sample_id) %>%
  mutate(relative_abundance = relative_abundance / sum(relative_abundance)) %>%
  ungroup()

# Sort legend by total abundance
class_order <- filtered_data %>%
  group_by(class_name) %>%
  summarise(total_abundance = sum(relative_abundance)) %>%
  arrange(desc(total_abundance)) %>%
  pull(class_name)

filtered_data$class_name <- factor(filtered_data$class_name, levels = class_order)

# Plot stacked barplot of significant taxa
# Customise bar colors
custom_colors <- c(
  "Bacilli" = "#2c1a57",
  "Flavobacteriia" = "#332b8e",
  "Cryptophyceae" = "#4e72da",
  "Clostridia" = "#8c564b",
  "Actinomycetes" = "#fffb00",
  "Epsilonproteobacteria" = "#ff7f0e",
  "Fusobacteriia" = "#0dd1f0",
  "Spirochaetia" = "#d10df0",
  "Thermotogae" = "#890df0",
  "Methanomicrobia" = "#d62728",
  "Tremellomycetes" = "#bbbbbe",
  "Caudoviricetes" = "#000000",
  "Desulfuromonadia" = "#e377c2",
  "Halobacteria" = "#bcbd22",
  "Holophagae" = "#318d3d",
  "Desulfobacteria" = "#aec7e8",
  "Deferribacteres" = "#9edae5",
  "Terriglobia" = "#98df8a",
  "Desulfobulbia" = "#ff9896",
  "Bacteriovoracia" = "#c5b0d5",
  "Fibrobacteria" = "#810000")

# Plot
ggplot(filtered_data, aes(x = group, y = relative_abundance, fill = class_name)) +
  stat_summary(fun = mean, geom = "bar", position = "stack") +
  scale_x_discrete(labels = c("FF_group" = "FF", "FL_group" = "FL", "SG_group" = "SG")) +
  labs(x = "Dietary group", y = "Mean relative abundance", fill = "Class") +
  scale_fill_manual(values = custom_colors) +
  theme_test() +
  theme(axis.title.y = element_text(size = 19, face = "bold"),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 17),
        legend.position = "right",
        legend.key.size = unit(0.65, "cm"),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 16, face = "bold"),
        panel.border = element_rect(linewidth = 1.0, fill = NA)) +
    guides(fill = guide_legend(ncol = 1))  # One column legend

# Save the final figure
ggsave("Relative microbial abundance.tif", width = 8, height = 6.43, dpi = 300)


# IMPORTANT NOTE: IF you ran the optional re-normalization step above, you MUST
# re-run the steps below before running the next intergroup significant difference
# in relative abundance step. 
filtered_data <- bracken_long %>%
  filter(class_name %in% significant_taxa)

class_order <- filtered_data %>%
  group_by(class_name) %>%
  summarise(total_abundance = sum(relative_abundance)) %>%
  arrange(desc(total_abundance)) %>%
  pull(class_name)

filtered_data$class_name <- factor(filtered_data$class_name, levels = class_order)

# INTERGROUP SIGNIFICANT DIFFERENCE IN RELATIVE ABUNDANCE

# Ensure required packages are loaded (tidyverse, ggpubr, rstatix)

# Create a named vector pattern for joining later
group_pattern <- c(rep("FF", 6), rep("FL", 5), rep("SG", 6)) # Integers rep number of samples per group

# Function to assign group within each class_name
assign_group <- function(n) {
  if (n == 17) {
    return(group_pattern)
  } else {
    # Truncate or recycle if n is different
    return(rep(group_pattern, length.out = n))
  }
}

# Apply the function per group using dplyr

filtered_data <- filtered_data %>%
  group_by(class_name) %>%
  mutate(group = assign_group(n())) %>%
  ungroup()

# Make sure group is factor
filtered_data$group <- factor(filtered_data$group, levels = c("FF", "FL", "SG"))

# Confirm grouping success
filtered_data %>%
  group_by(class_name, group) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = group, values_from = n, values_fill = 0)

# Ensure factor levels are consistent
filtered_data$group <- factor(filtered_data$group, levels = c("FF", "FL", "SG"))

# Compute per-class pairwise Wilcoxon tests
pairwise_results <- filtered_data %>%
  group_by(class_name) %>%
  pairwise_wilcox_test(
    relative_abundance ~ group,
    p.adjust.method = "BH"
  ) %>%
  ungroup() %>%
  filter(p.adj < 0.05) %>%
  mutate(p.label = paste0("p = ", signif(p.adj, 3)))

view(pairwise_results)

# Calculate per-class max y values for annotation placement
y_positions <- filtered_data %>%
  group_by(class_name) %>%
  summarise(y.max = max(relative_abundance, na.rm = TRUE), .groups = "drop")

# Join to p-values
pairwise_results <- pairwise_results %>%
  left_join(y_positions, by = "class_name") %>%
  group_by(class_name) %>%
  mutate(y.position = y.max + y.max * 0.1 * row_number()  # 10% spacing above max
  ) %>%
  ungroup()

# Plot
ggplot(filtered_data, aes(x = group, y = relative_abundance)) +
  geom_boxplot(aes(color = group), outlier.shape = NA, size = 1, alpha = 1) +
  facet_wrap(~ class_name, scales = "free_y") +
  stat_pvalue_manual(
    data = pairwise_results,
    mapping = aes(xmin = group1, xmax = group2, y.position = y.position, label = p.label),
    tip.length = 0.01,
    size = 4,
    inherit.aes = FALSE
  ) +
  scale_color_manual(values = c("FF" = "#00276f", "FL" = "#000000", "SG" = "#ffb600")) +
  labs(
    x = NULL,
    y = "Relative abundance" 
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 16, face = "bold"),  # X-axis title font size
    axis.title.y = element_text(size = 16, face = "bold"),  # Y-axis title font size
    legend.position = "none"                 # Remove legend
  )

# 4. Save the final figure
ggsave("Intergroup microbial abundance comparisons_3.tif", width = 12, height = 12, dpi = 300)

# INDICATOR TAXA ANALYSIS

# Re-scale (De-normalise) OTU table to Pseudo-counts
# First inspect your OTU table
head(otu_table)
summary(otu_table)

# Multiply each value by a constant factor (e.g., 10,000)
# This scales the relative abundances into pseudo-counts
otu_counts <- round(otu_table * 10000)

# Confirm the new data is numeric and integer-like
str(otu_counts)
summary(rowSums(otu_counts))  # Will vary by sample

# Check for zero-inflation
zero_proportion <- sum(otu_counts == 0) / length(otu_counts)
cat("Proportion of zeros:", round(zero_proportion, 3), "\n")

# Perform indicator species analysis
set.seed(123)
indval_result <- multipatt(otu_counts, sample_metadata$Group, func = "IndVal.g", control = how(nperm = 999))

# Capture printed output into a character vector
summary_text <- capture.output(summary(indval_result, indvalcomp = TRUE, alpha = 0.05))

# Write to file if needed
writeLines(summary_text, "Indicator_taxa_summary.txt")

# Convert the summary result into a dataframe manually
# `summary_result` is printed text; real values live in `indval_result$sign`
sign_df <- indval_result$sign

# Extract significant taxa info
sign_df$Taxon <- rownames(sign_df)

# Filter only significant taxa
sig_df <- sign_df %>%
  filter(!is.na(p.value), p.value < 0.05)

# Find group combinations per taxon
# We use max col index from `indval_result$str` to identify associated group/combo
group_names <- colnames(indval_result$str)
group_assignments <- apply(indval_result$str[sig_df$Taxon, , drop = FALSE], 1, function(row) {
  assigned <- which.max(row)
  group_names[assigned]
})

# Combine into a tidy data frame for plotting
plot_df <- sig_df %>%
  mutate(
    Group = group_assignments,
    Group = recode(Group, 
                   "FF_group" = "FF",
                   "FL_group" = "FL",
                   "SG_group" = "SG",
                   "FF_group+SG_group" = "FF+SG",
                   "FL_group+SG_group" = "FL+SG",
                   "FF_group+FL_group" = "FF+FL"),
    sig_label = case_when(
      p.value <= 0.001 ~ "***",
      p.value <= 0.01  ~ "**",
      p.value <= 0.05  ~ "*",
      TRUE ~ ""
    )
  )

# Order taxon names by effect size (stat)
plot_df$Taxon <- fct_reorder(plot_df$Taxon, plot_df$stat)

# Create bubble plot
ggplot(plot_df, aes(x = Group, y = Taxon)) +
  geom_point(aes(size = stat, fill = p.value), shape = 21, color = "NA") +
  scale_fill_gradient(low = "#e41a1c", high = "#fee0d2", name = "p-value") +
  scale_size_continuous(
    name = "IndVal stat",
    range = c(4, 10),
    limits = c(0.75, 1.0),
    breaks = seq(0.75, 1.0, length.out = 6),
    guide = guide_legend(title.position = "top", override.aes = list(fill = "black"))
  ) +
  labs(
    title = "Significant Indicator Taxa by Group",
    x = "Group",
    y = "Taxonomic class"
  ) +
  theme_bw  (base_size = 14) +
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_blank(), # Removes x-axis title
    axis.text.y = element_text(size = 16),
    axis.text.x = element_text(size = 16),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 16),
    plot.title = element_blank() # # Removes plot title 
  )
# Save the final figure
ggsave("Significant indicator taxa by group.tif", width = 8, height = 5, dpi = 300)
