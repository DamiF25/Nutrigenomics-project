# If required, install packages
install.packages("UpSetR")
install.packages("tidyverse")

# Load required libraries
library(UpSetR)
library(tidyverse)

# Set working directory

# Load the input CSV file
deg_data <- read.csv("Upset plot data.csv", stringsAsFactors = FALSE)

# Define the groups
groups <- c("FLvFF_up", "FLvFF_down", 
            "FLvSG_up", "FLvSG_down", 
            "SGvFF_up", "SGvFF_down")

# Clean and prepare gene list
gene_list <- lapply(groups, function(g) {
  genes <- deg_data[[g]]
  genes <- genes[!is.na(genes) & genes != "" & genes != "NA" & grepl("\\S", genes)]
  unique(genes)
})
names(gene_list) <- groups

# Create membership matrix
all_genes <- unique(unlist(gene_list))
membership_matrix <- sapply(groups, function(g) {
  all_genes %in% gene_list[[g]]
})

# View membership_matrix
head(membership_matrix)

# Convert TRUE/FALSE to integers (1 or 0)
membership_matrix <- apply(membership_matrix, 2, as.integer)

# View membership_matrix
head(membership_matrix)

# Create a dataframe and set row names
membership_df <- as.data.frame(membership_matrix)
rownames(membership_df) <- all_genes

# Now plot the upset plot
upset(
  membership_df,
  sets = groups,
  sets.bar.color = "#000000",
  order.by = "freq",
  keep.order = TRUE,
  main.bar.color = "#000000",
  point.size = 8,                   # Bigger points
  line.size = 1.5,                  # Thicker connection lines
  text.scale = c(5, 4, 5, 4, 4, 4)  # Text size tuning
)

# First text scale number = y-axis title (Intersection Size)
# second number = y-axis labels (Intersection size numbers)
# Third number = x-axis title (Set Size)
# Fourth number = Set size labels (on the left)
# Fifth number = y-axis labels (Groups)
# sixth number = For numeric values above vertical bars

#Export image size 2000W x 1608H 