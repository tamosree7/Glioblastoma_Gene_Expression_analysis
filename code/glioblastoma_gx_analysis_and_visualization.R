# Load necessary libraries
library(gplots)
library(RColorBrewer)
library(dplyr)
library(tidyr)
library(broom)
library(stats)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
# Load the dataset
df <- read.csv("glioblastoma.csv", row.names = 1)
# 1. Checking type of data ------------------------------------------------
mode(df) #list
genes = rownames(df)
df = apply(df, 2, as.numeric)
mode(df) #now it is numeric
rownames(df) = genes
# 2. Removing constant genes ----------------------------------------------
rowvar = apply(df, 1, var, na.rm = T)
constvar = (rowvar == 0 | is.na(rowvar))
sum(constvar)     #there is no constant genes
# 3. Handling NAs ---------------------------------------------------------
sum(is.na(df)) #there are no nulls
# 4. Log transformation ---------------------------------------------------
logged_data = log2(df +1)
#Visualizing data before and after transformation
par(mfrow = c(1, 2))
p1=plot(density(df), main = "Before Transformation")
p2=plot(density(logged_data), main = "After Transformation")
# 5. Scaling --------------------------------------------------------------
scaled_data = scale(t(logged_data), center = T, scale = T)
# 6. Visualization (Heatmaps) ---------------------------------------------
#############Gene only-Heatmap############
library(RColorBrewer)
# Set up color palette for heatmap
my_palette <- colorRampPalette(brewer.pal(11, "RdYlBu"))(n = 100)  # Diverging palette
# Heatmap 1: Cluster only genes (rows)
heatmap.2(as.matrix(t(scaled_data)),
          Rowv = TRUE,  # Cluster genes (rows)
          Colv = FALSE,  # Do not cluster samples (columns)
          col = my_palette,
          scale = "row",
          trace = "none",
          dendrogram = "row",  # Show dendrogram only for rows
          main = "Clustered Genes",
          density.info = "none",
          margins = c(10, 10),  # Increase margins for row/column names
          cexCol = 0.7,  # Adjust font size for column (sample) labels
          cexRow = 0.7,  # Adjust font size for row (gene) labels if needed)
          srtCol = 45,  # Rotate column labels by 45 degrees
          keysize = 1)
############# Sample only-Heatmap #####################
heatmap.2(as.matrix(t(scaled_data)),
          Rowv = FALSE,  # Cluster genes (rows)
          Colv = TRUE,  # Do not cluster samples (columns)
          col = my_palette,
          scale = "row",
          trace = "none",
          dendrogram = "column",  # Show dendrogram only for rows
          main = "Clustered Samples",
          density.info = "none",
          margins = c(10, 10),  # Increase margins for row/column names
          cexCol = 0.5,  # Adjust font size for column (sample) labels
          cexRow = 0.5,  # Adjust font size for row (gene) labels if needed)
          srtCol = 45,
          key = TRUE)  # Rotate column labels by 45 degrees
####################Both genes and samples################
# Set color palettes
diverging_palette <- brewer.pal(11, "RdBu")  # Diverging color palette
sequential_palette <- brewer.pal(9, "Greens")  # Sequential color palette
# Heatmap with clustering both genes and samples
heatmap.2(as.matrix(t(scaled_data)),
          Rowv = TRUE,
          Colv = TRUE,
          main = "Heatmap of Gene Expression",
          col = diverging_palette,
          scale = "row",  # Scale rows (genes)
          trace = "none",  # No trace lines
          dendrogram = "both",  # Cluster both rows and columns
          margins = c(10, 10),  # Margins for labels
          key = TRUE,  # Show color key
          keysize = 1.5,  # Size of the key
          density.info = "none",  # No density information
          cexRow = 0.5,  # Size of row labels
          cexCol = 0.5,  # Size of column labels
          srtCol = 45)
###########################divergent heatmap##################
heatmap.2(as.matrix(t(scaled_data)),
          Rowv = FALSE,  # Cluster genes (rows)
          Colv = FALSE,  # Do not cluster samples (columns)
          col = my_palette,
          scale = "row",
          trace = "none",
          dendrogram = "none",  # Show dendrogram only for rows
          main = "Diverging Heatmap",
          density.info = "none",
          margins = c(8, 10),  # Increase margins for row/column names
          cexCol = 0.7,  # Adjust font size for column (sample) labels
          cexRow = 0.7,  # Adjust font size for row (gene) labels if needed)
          srtCol = 45)  # Rotate column labels by 45 degrees
###############################Sequential palette Heatmap###############
# Sequential color palette
my_palette = colorRampPalette(brewer.pal(9, "YlOrRd"))(100)
heatmap.2(as.matrix(t(scaled_data)),
          Rowv = FALSE,  # Cluster genes (rows)
          Colv = FALSE,  # cluster samples (columns)
          col = my_palette,
          scale = "row",
          trace = "none",
          dendrogram = "both",  # Show dendrogram only for rows
          main = "Sequential Heatmap for clustered genes and samples",
          density.info = "none",
          margins = c(10, 10),  # Increase margins for row/column names
          cexCol = 0.7,  # Adjust font size for column (sample) labels
          cexRow = 0.7,
          srtCol = 45)  # Adjust font size for row (gene) labels if needed))
# 7. Clustering -----------------------------------------------------------
df= as.data.frame(df)
g1 <- df %>% dplyr::select(c("TCGA.06.2570.01A.01R.1849.01","TCGA.02.2483.01A.01R.1849.01","TCGA.06.5410.01A.01R.1849.01"))
g2 <- df %>% dplyr::select(c("TCGA.19.4065.02A.11R.2005.01","TCGA.19.0957.02A.11R.2005.01","TCGA.06.0152.02A.01R.2005.01","TCGA.14.1402.02A.01R.2005.01","TCGA.14.0736.02A.01R.2005.01","TCGA.19.5960.01A.11R.1850.01","TCGA.14.0781.01B.01R.1849.01"))
#Foldchange
fold_change <- rowMeans(g2) / rowMeans(g1)
fold_change <- log2(fold_change)
t_test_result <- t.test(g1, g2)
# Extract p-value
p_value <- t_test_result$p.value
# Combine results into a data frame
results <- data.frame(Gene = rownames(df), FoldChange = fold_change, PValue = p_value)
# Adjust p-values for multiple testing (optional)
results$AdjustedPValue <- p.adjust(results$PValue, method = "BH")
# View the results
head(results)
# thresholds
fc_threshold <- 2
pval_threshold <- 0.05
#===============================================================================
#subset upregulated genes
upregulated_genes <- df %>% filter(results$FoldChange > fc_threshold & results$PValue < pval_threshold)
# Subset downregulated genes
downregulated_genes <- df %>% filter(results$FoldChange < -fc_threshold+2 & results$PValue< pval_threshold)
# Upregulated pathways
Up_enrichment_results <- data.frame(
  Pathway = c("Loop of Henle development", " Negative reg. of myoblast differentiation","Adenylate cyclase-activating G protein-coupled receptor signaling pat", "Positive reg. of cytosolic calcium ion concentration", "Reg. of cytosolic calcium ion concentration"),
  GeneCount = c(3, 4, 9, 10, 11),
  PValue = c(3.4E-02, 3.2E-02	, 5.5E-03	, 3.8E-02	, 3.4E-02	)
)
# Calculate -log10(P-value) for better visualization of significance
Up_enrichment_results$logPValue <- -log10(Up_enrichment_results$PValue)
# Create a lollipop plot
# Up_enrichment pthway
ggplot(Up_enrichment_results, aes(x = reorder(Pathway, GeneCount), y = GeneCount)) +
  geom_segment(aes(xend = Pathway, yend = 0), color = "grey") +
  geom_point(aes(size = logPValue), color = "blue") +
  labs(title = "Top 5 Upregulated Enriched Pathways",
       x = "Pathway",
       y = "Number of Genes",
       size = "-log10(P-value)") +
  coord_flip() +  # Flips the axes for better readability
  theme_minimal()
#down_enrichment pthway
down_enrichment_results <- data.frame(
  Pathway = c("Neutrophil chemotaxis", "Neutrophil migration","Endocrine system development", "Granulocyte chemotaxis", "Granulocyte migration"),
  GeneCount = c(6, 7, 6, 6, 7),
  PValue = c(6.1E-05	, 2.7E-05	,1.0E-04, 1.0E-04, 4.5E-05	)
)
# Calculate -log10(P-value) for better visualization of significance
down_enrichment_results$logPValue <- -log10(down_enrichment_results$PValue)
# Create a lollipop plot
# Assuming down_enrichment_results is your data frame
ggplot(down_enrichment_results, aes(x = reorder(Pathway, GeneCount), y = GeneCount)) +
  geom_segment(aes(xend = Pathway, yend = 0), color = "grey") +
  geom_point(aes(size = logPValue), color = "blue") +
  labs(title = "Top 5 Downregulated Enriched Pathways",
       x = "Pathway",
       y = "Number of Genes",
       size = "-log10(P-value)") +
  coord_flip() +  # Flips the axes for better readability
  theme_minimal()