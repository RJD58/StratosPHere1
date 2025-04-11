# Hierarchical clustering of Stratosphere 1 qPCR and cell surface BMPR2 levels


library(stats)
library(pheatmap)
library(RColorBrewer)


data <- Biomarker_panel_cell_surface

# Columns 2-12 are relevant numeric. Column 1 = disease group; control, IPAH, BMPR2-PAH

set.seed(123)

# This is mixed protein and dCT data so needs adjusting so it is on the same scale

df<-data

df[, 2:12] <- lapply(df[, 2:12], log)
df[, 2:12] <- lapply(df[, 2:12], scale)
df <- as.data.frame(df)


# Separate Group column and numeric columns
group <- df$Group
numeric_data <- df[ , -1]  # Select only numeric data

# Set the order of group levels to the desired order
group <- factor(group, levels = c("PAH", "BMPR2", "Control"))

# Calculate Euclidean distance and hierarchical clustering for columns
dist_cols <- dist(t(numeric_data), method = "euclidean")
hclust_cols <- hclust(dist_cols, method = "complete")

# Order rows by grouping them together (to keep groups together on the Y-axis)
ordered_indices <- order(group)
numeric_data_ordered <- numeric_data[ordered_indices, ]
group_ordered <- group[ordered_indices]

# Create row annotations to indicate groups
annotation_row <- data.frame(Group = group_ordered)
rownames(annotation_row) <- rownames(numeric_data_ordered)

# Define a blue-white-red palette
color_palette <- colorRampPalette(rev(brewer.pal(n = 5, name = "RdBu")))(100)

# Assume `annotation_row` contains a "Group" column specifying group membership
group_boundaries <- which(diff(as.numeric(as.factor(annotation_row$Group))) != 0)

# Define breaks for data values (if needed)
breaks <- seq(min(numeric_data_ordered), max(numeric_data_ordered), length.out = 101)

# Generate the heatmap with gaps between groups
pheatmap(
  numeric_data_ordered, 
  cluster_cols = hclust_cols, 
  cluster_rows = FALSE, 
  annotation_row = annotation_row, 
  show_rownames = FALSE, 
  show_colnames = TRUE, 
  gaps_row = group_boundaries,       
  color = color_palette, 
  breaks = breaks
)



