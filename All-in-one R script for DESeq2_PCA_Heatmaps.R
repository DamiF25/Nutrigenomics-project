#install required packages

install.packages("DESeq2")
install.packages("tidyverse")
install.packages("ggplot2")
install.packages("ggrepel")
install.packages("dplyr")
install.packages("pheatmap")
install.packages("RColorBrewer")
BiocManager::install()


#Load libraries

library(DESeq2)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(BiocManager)


#Set work directory

#Load count data
count_data <- read.csv("FF_vs_FL_Gene_count_matrix.csv",header=TRUE,row.names=1)

#Check
colnames(count_data)
head(count_data)
dim(count_data)

#Load sample information
sample_info <- read.csv("FF_vs_FL_PHENO_DATA.csv",header=TRUE,row.name=1) 

#Check
colnames(sample_info)
head(sample_info)
dim(sample_info)

#Confirm that sample names match in both count data and sample_info files
all(colnames(count_data) %in% rownames(sample_info))

#Confirm that they are also in the same order
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

#Save DESeq result
write.csv(res, file = "FF_vs_FL_Gene_count_DE_final.csv")

#Check DESeq result at p <.01
res <- results(dds, alpha = 0.01)
summary(res)

#Check DESeq result at p <.05
res <- results(dds, alpha = 0.05)
summary(res)

#Extract DE genes with padj < 0.01 and log2Foldchange <=-1 or >=1
deg <- subset(res, padj<0.01 & abs(log2FoldChange) >=1)
summary(deg)

#Remove NA from deg - OPTIONAL
deg <- na.omit(deg)
summary(deg)

#Order deg result based on padj value
deg <- deg[order(deg$padj),]

#You can also order deg result based on log2FoldChange
deg <- deg[order(deg$log2FoldChange),]

head(deg)
dim(deg)

#Save
write.csv(deg, "FF_vs_FL_Gene_count_DE_pval_0.01_Final.csv")

#Explore result
summary(res)

#TO CHANGE THE ADJUSTED p-value IN THE RESULT SUMMARY
res0.05 <- results(dds, alpha = 0.05)
summary(res0.05)

#Change DESeq result into a dataframe
res <- as.data.frame(res)
class(res)

#To check the DESeq result components
resultsNames(dds)


#############DATA VISUALISATION#####################################

#Volcano plot for DEGs
res$diffexpressed <- 'NO'
res$diffexpressed[res$log2FoldChange > 1.0 & res$padj < 0.01] <- 'UP'
res$diffexpressed[res$log2FoldChange < -1.0 & res$padj < 0.01] <- 'DOWN'

ggplot(data = res, aes(x = log2FoldChange, y = -log10(padj), col = diffexpressed)) + 
  geom_vline(xintercept = c(-1.0, 1.0), col = 'gray', linetype = 'dashed') +
  geom_hline(yintercept = c(2), col = 'gray', linetype = 'dashed') +
  geom_point(size = 4.0, alpha = 0.7) +
  scale_colour_manual(values = c("#2ca8f4", "black", "orange"),
                      labels = c("Downregulated", "Not Significant", "Upregulated")) +
  theme_bw() +
  coord_cartesian(ylim = c(0, 100), xlim = c(-10, 20)) +
  scale_x_continuous(breaks = seq(-30, 30, 5)) +
  labs(color = 'FL vs FF',
       x = expression("log"[2]*"FC"), y = expression("-log"[10]*"padj.")) +
  theme(
    text = element_text(size = 15),  # Adjust the overall text size
    axis.text = element_text(size = 11),  # Adjust the axis text size
    axis.title = element_text(size = 13),  # Adjust the axis title size
    legend.text = element_text(size = 13),  # Adjust the legend text size
    legend.title = element_text(size = 15)  # Adjust the legend title size
  )

# Save the final figure
ggsave("Volcano plot_FFvsFL.png", width = 5.75, height = 4.62, dpi = 300)

#Plot dispersion estimates
plotDispEsts(dds, main="FF vs FL Dispersion Estimates") #800W x 664H

#MA plot can also be used to interpret result
plotMA(dds,ylim=c(-2,2))

#Remove the noise through LFC shrinkage
resLFC <- lfcShrink(dds,coef="Diet_FF._vs_FL",type="apeglm")

#Plot MA without noise
plotMA(resLFC,ylim=c(-5,5))


#Variance stabilizing transformation#Variancesample() stabilizing transformation
vsd <- vst(dds,blind=FALSE)

dim(vsd)

#Generate a PCA plot

#Use transformed values to generate a PCA plot
plotPCA(vsd,intgroup=c("Diet"))

pcaData <- plotPCA(vsd,intgroup=c("Diet"))


#Extract the data used in the PCA plot
pca_data <- pcaData$data

#Customise PCA plot using ggplot2
ggplot(pca_data, aes(PC1, PC2, colour = Diet)) +
  geom_point(size=3) + ggtitle("PCA plot with VST data") + # Delete ggtitle() to remove title
  xlab("PC1: insert % variance") + ylab("PC2: insert % variance") + 
  coord_fixed(ratio=1) + theme_light() +
  stat_ellipse(geom = "polygon", aes(fill = after_scale(alpha(colour, 0.3))),
               data = pca_data[pca_data$group !="versicolor",]) +
  theme(axis.title = element_text(size = 18, face = "bold"),
        axis.text = element_text(size = 16)) +
  theme(legend.text = element_text(size = 16), legend.title = element_text(size = 18))

# Save the final figure
ggsave("PCA plot_grey.tiff", width = 8, height = 6.43, dpi = 300)



#Heatmaps

#Heatmap of sample-to-sample distance matrix (with clustering) based on vsd

#Generate distance matrix (Euclidean Distance)
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
colnames(sampleDistMatrix)

#Set color scheme
colors <- colorRampPalette(rev(brewer.pal(9,"Blues")))(255)

#You can also use customised color scheme (if this preferred, change 'col=colors' to 'col=my_colors' below)
my_colors <- colorRampPalette(c("#3f447d","#ffffff"))(255)  


#Generate the heatmap
pheatmap(sampleDistMatrix, annotation_col = sample_info, clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists, col=colors)

#Or run
pheatmap(sampleDistMatrix, annotation_col = sample_info)



#To modify Diet color as preferred, using the annotation_colors parameter 
Diet_df <- data.frame(Diet= rep(c("FF", "FL"), c(6,5))) #Numeric values represent number of samples in each group

ann_colors = list(Diet = c("FF" = "#990906", "FL" = "#000000"))

row.names(Diet_df) <- colnames(sampleDistMatrix)

#Regenerate the heatmap
pheatmap(sampleDistMatrix, annotation_col = Diet_df, clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists, col=colors, annotation_colors = ann_colors, 
         fontsize = 10,           # overall font size
         fontsize_row = 10,        # font size for row labels
         fontsize_col = 10        # font size for column labels
)

#To generate heatmap for differentially expressed genes based on padj values

#Extract normalised counts from the DESeq result
normCounts <- counts(dds, normalized = T)
write.csv(normCounts, "FF vs FL_NormGene_count_matrix_final.csv")

#Log transformation (transform normalised counts) to avoid high expression values 
#to dominate the plot
transCounts <- log2(normCounts+1)

head(transCounts)

write.csv(transCounts, "FF vs FL_transGene_count_matrix_final.csv")

#HEATMAP FOR ALL DEGS

all_hits <- row.names(deg)

all_hits

#To extract respective transformed counts of all hits
all_hits <- transCounts[all_hits,]
head(all_hits)

#Generate heatmap (without clustering)
pheatmap(all_hits,cluster_rows=FALSE,cluster_cols=FALSE)

#To enable heatmap clustering
pheatmap(all_hits)

#To add annotation
pheatmap(all_hits, annotation_col = sample_info)

#OR add annotation without clustering
pheatmap(all_hits,cluster_rows=FALSE,cluster_cols=FALSE,annotation_col = sample_info)


#Generate heatmap of Z scores using the all hits
cal_z_score <- function(x) {(x-mean(x)) / sd(x)}

#First calculate z scores for all normalised counts
zscore_all <- t(apply(normCounts, 1, cal_z_score))

#generate heatmap
pheatmap(zscore_all,cluster_rows=TRUE,cluster_cols=TRUE,annotation_col = sample_info)

#To remove gene IDs in heatmap
pheatmap(zscore_all,cluster_rows=TRUE,cluster_cols=TRUE,show_rownames = FALSE, annotation_col = sample_info)

#To change heatmap colour
colors <- colorRampPalette(rev(brewer.pal(9,"RdGy")))(255)

#Regenerate heatmap
pheatmap(zscore_all,cluster_rows=TRUE,cluster_cols=TRUE,show_rownames = FALSE, annotation_col = sample_info, col=colors)


#To modify Diet color as preferred, using the annotation_colors parameter 

colors <- colorRampPalette(rev(brewer.pal(9,"RdYlBu")))(255)

Diet_df <- data.frame(Diet= rep(c("SG", "FF", "FL"), c(6,6,5)))

ann_colors = list(Diet = c("SG" = "#0a5fb5", "FF" = "#990906", "FL" = "#000000"))

row.names(Diet_df) <- colnames(zscore_all)


#Regenerate heatmap

pheatmap(zscore_all,cluster_rows=TRUE,cluster_cols=TRUE,show_rownames = FALSE,
         annotation_col = Diet_df,col=colors,annotation_colors = ann_colors,
         fontsize = 10,           # overall font size
         fontsize_row = 10,        # font size for row labels
         fontsize_col = 10        # font size for column labels
)


#HEATMAP FOR TOP 50 DEGs 

#The pheatmap package conventionally limits the number of ‘top hits’ subset to twenty when
#z score values are estimated and used. Therefore, I will use transformed normalised data (not Z score) to
#generate heatmaps for top DEGs. 

#Fetch top hits
top_hits <- row.names(deg[1:50, ])

top_hits

#To extract respective transformed counts of the top hits
top_hits <- transCounts[top_hits,]
head(top_hits)

#Generate heatmap (without clustering)
pheatmap(top_hits,cluster_rows=FALSE,cluster_cols=FALSE)

#To enable heatmap clustering
pheatmap(top_hits)

#To add annotation
pheatmap(top_hits, annotation_col = sample_info)


#To modify Diet color using the annotation_colors parameter 

colors <- colorRampPalette(rev(brewer.pal(9,"RdBu")))(255)

#To modify color based on personal preference
my_colors <- colorRampPalette(c("#2ca8f4", "black", "orange"))(255)

Diet_df <- data.frame(Diet= rep(c("FF", "FL"), c(6,5)))

ann_colors = list(Diet = c("FF" = "#990906", "FL" = "#000000"))

row.names(Diet_df) <- colnames(top_hits)


#Regenerate heatmap
pheatmap(top_hits, cluster_rows=TRUE,cluster_cols=TRUE,show_rownames = TRUE,
         annotation_col = Diet_df,col=my_colors,annotation_colors = ann_colors)

#SCALING THE HEATMAP
#Since we are interested in knowing the difference in each variable (DEG) across the snail
#samples, we will scale by row because the sample names (i.e. snail samples) are listed across 
#the columns. If the snail samples were listed across rows, we would scale by column.
pheatmap(top_hits, scale = "row", cluster_rows=TRUE,cluster_cols=TRUE,show_rownames = TRUE,
         annotation_col = Diet_df,col=my_colors,annotation_colors = ann_colors)

#Splitting the map into columns and rows
pheatmap(top_hits, scale = "row", cluster_rows=TRUE,cluster_cols=TRUE,show_rownames = TRUE,
         annotation_col = Diet_df,col=my_colors,annotation_colors = ann_colors,
         cutree_cols=2, cutree_rows=2)

#use 575W x 949H for image export

