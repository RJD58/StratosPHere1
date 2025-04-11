# Code used in 'BMPR-II Biomarkers for testing therapeutic efficacy in pulmonary arterial hypertension – novel findings from the StratosPHere 1 study.' Rowena J Jones

# CODE: DGE of RNAseq data by biomarker cluster using DESeq2

# Load the following

library(tidyverse)
library(dplyr)
library(DESeq2)
library(tibble)


# Prepare the data matrix and sample info ---------------------------------


#Load the file. 

counts <-cts

counts <-counts data

dim(counts)
# 508 patients with 56517 genes

#transpose counts and make a df - need to do this to merge with the cluster allocation file
counts <- as.data.frame(counts)
counts2 <- rbind(colnames(counts), counts)
colnames(counts2) <- NULL 
counts3 <-t(counts2)
counts3 <- as.data.frame(counts3)

#there are some duplicates in sampleID
counts3_unique <- counts3[!duplicated(counts3[, 1]), ]


#load in the list of patients with cluster allocation (KM3)

KM3<- biomarkerKM3_for_DGE
Biomark_counts <- merge (counts3_unique, KM3, by.x = '1', by.y = 'SampleID')



# Now need to make the Biomark_counts file ready for DESeq2

# Convert first column to row names 
Biomark_counts <- Biomark_counts %>% column_to_rownames(var = colnames(Biomark_counts)[1])


# remove superfluos rows

Biomark_counts <- Biomark_counts[,-56520]
Biomark_counts <- Biomark_counts[,-56519]
Biomark_counts <- Biomark_counts[,-56518]


# Convert all columns to numeric while preserving row names - otherwise lapply will delete row names
Biomark_counts2 <- Biomark_counts  
Biomark_counts2[] <- lapply(Biomark_counts2, as.numeric) 

# Transpose this so that the sample ID is along top and genes downside
Biomark_counts_matrix <- t(Biomark_counts2)
# Make a dataframe
Biomark_counts_matrix <- as.data.frame (Biomark_counts_matrix)
# 56517 variables (genes) for 356 patients - matches patients in paper

# Generate Sampleinfo meta data

SampleInfo <- KM3 %>% select(SampleID, KM3)
SampleInfo$KM3 <- as.factor(SampleInfo$KM3)


# It is important to be sure that the order of the samples in rows in the sample meta data table matches the order of the columns in the data matrix 

all(colnames(Biomark_counts_matrix) == SampleInfo$SampleID)

# if FALSE not re-order to match
df <- match(SampleInfo$SampleID, colnames(Biomark_counts_matrix))
Biomark_counts_matrix <- Biomark_counts_matrix[, df]

all(colnames(Biomark_counts_matrix) == SampleInfo$SampleID)
# TRUE

# Make sure the values are all integers

Biomark_counts_matrix <- round(Biomark_counts_matrix)

Biomark_counts_matrix <- as.matrix(Biomark_counts_matrix)



# Define the model --------------------------------------------------------

# Use the sample info meta data to define the model. We have one catergorical (KM3) with 3 levels in it.

# set KM=2 as the reference / intercept as this has the more 'normal' biomarker 
SampleInfo$KM3 <- relevel(SampleInfo$KM3, ref = "2")

model <- SampleInfo

model$KM3<- as.factor(model$KM3)

all(colnames(Biomark_counts_matrix) == model$SampleID)
# TRUE


# Make the SampleID to rownames for the next step
model <- model %>% column_to_rownames(var = "SampleID")



# Build a DESEq2 object ---------------------------------------------------



Biomark_dds_raw <- DESeqDataSetFromMatrix(countData = Biomark_counts_matrix,
                                          colData = model,
                                          design = ~ KM3)

# 56517 genes


# Perform filtering steps -------------------------------------------------


# Filter out low expressed genes

# Filter out a minimum number of counts in a minimum of 3 samples
minCountPerSample <- 10  # Minimum count in a sample to consider it
minSamples <- 3          # Minimum number of samples that should meet the threshold

# filter
Biomark_dds_filt <- Biomark_dds_raw[rowSums(counts(Biomark_dds_raw) >= minCountPerSample) >= minSamples,]

# 31775 variables



# Annotate and pull out protein coding only and removal of XY genes -------------------------------

# Pull out XY genes

XY <- mart_export.4.13.12.24

Biomark_dds_XYremoved <- Biomark_dds_filt[!rownames(Biomark_dds_filt) %in% XY$Gene.name, ]

# 31098 elements

# Pull out only Protein coding genes

biotype<- mart_export.biotype.13.12.24

unique(biotype$Gene.type)

# Subset rows where gene.type is "protein coding"
proteincoding <- biotype[biotype$Gene.type == "protein_coding", ]

Biomark_dds_proteincoding <- Biomark_dds_XYremoved[rownames(Biomark_dds_XYremoved) %in% proteincoding$Gene.name, ]

# 15073 protein coding genes

# can observe the genes we have kept

kept_genes <- rownames(Biomark_dds_proteincoding)




#  Run The DESeq command --------------------------------------------------


ddsObj <- DESeq(Biomark_dds_proteincoding)

# generate the normalised counts and save the file

normalized_counts <- counts(ddsObj, normalized=TRUE)

write.table(normalized_counts, "normalized_counts.txt", sep="\t", quote=F, col.names=NA)

norm.counts_KM3DGE <- as.data.frame(normalized_counts)

write.csv(norm.counts_KM3DGE, 'normalised counts of KM3 DGE.csv')


# view results


# Extracting specific contrasts -------------------------------------------


# Km3v2
results.simple <- results(ddsObj, alpha = 0.05)
results.simple

# Results. Simple is KM3v2
write.csv(results.simple, ' KM3v2 non shrunk.csv')

# This shows all the contrasts in the ddsObj
resultsNames(ddsObj)
# This only shows "KM3_1_vs_2" "KM3_3_vs_2" not KM1 vs 3 This is because KM=2 is set as the baseline. Have to do KM3v1 after all this is done as otherwise it messes with the contrasts


# Generate the KM1vs2 comparison
# Note the order of the names determines the direction of FC reported. The second element is the level that is used as the baseline 
contrast1v2 <- c("KM3", "1", "2")

results.KM1vs2 <- results(ddsObj, contrast=contrast1v2, alpha = 0.05)
results.KM1vs2
write.csv(results.KM1vs2, ' KM1v2 non shrunk.csv')

# The remaining contrast of KM3v1 needs to be performed seperately, otherwise there is a glitch.

# Perfrom LogFC shrinkage -------------------------------------------------

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("apeglm", force=TRUE)
library(apeglm)

# Use Apeglm for shrinkage, there are 3 options: normal (not advised due to bias), ashr and apeglm. I have selected apeglm for the purpose of ranking genes and reducing false positives
# Perform shrinkage 3v2

# TO AVOID CONFUSION - DO ONE SHRINKAGE AND VOLCANO ETC, DELETE AND MOVE ONTO THE NEXT

shrink_3v2_apeglm <- lfcShrink(ddsObj, 
                               coef="KM3_3_vs_2", 
                               type="apeglm")
shrink_3v2_apeglm

shrink_1v2_apeglm <- lfcShrink(ddsObj, 
                               coef="KM3_1_vs_2", 
                               type="apeglm")
shrink_1v2_apeglm

# Note: at this point, due to the design we only have these comparisons "Intercept"  "KM3_1_vs_2" "KM3_3_vs_2". We need to re-level to look at 3v1
resultsNames(ddsObj)


# Plot shrinkage of the above - we can see there is some shrinkage

plotMA(shrink_3v2_apeglm, ylim=c(-2,2))
# vs
plotMA(results.simple, ylim=c(-2,2))


plotMA(shrink_1v2_apeglm, ylim=c(-2,2))
#vs
plotMA(results.KM1vs2, ylim=c(-2,2))


# save these shrunken tables

write.csv(shrink_3v2_apeglm, ' KM3v2 shrunk_apeglm_15.2.25.csv')
write.csv(shrink_1v2_apeglm, ' KM1v2 shrunk_apeglm_15.2.25.csv')


# DE results table --------------------------------------------------------

shrink_3v2_apeglm<- as.data.frame(shrink_3v2_apeglm)

shrink_1v2_apeglm<- as.data.frame(shrink_1v2_apeglm)


library(tidyverse)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)


# Apply thresholds. Use p=0.01 and Log2FC of 0.59 = 1.5FC
shrink_3v2_vol <- shrink_3v2_apeglm %>% 
  mutate(threshold_OE = padj < 0.01 & abs(log2FoldChange) >= 0.59 )

shrink_1v2_vol <- shrink_1v2_apeglm %>% 
  mutate(threshold_OE = padj < 0.01 & abs(log2FoldChange) >= 0.59 )



# move the rownames to a new column 'GeneID'

shrink_3v2_vol$GeneID <- rownames(shrink_3v2_vol)
shrink_1v2_vol$GeneID <- rownames(shrink_1v2_vol)



# Volcano plot with annotations

# Filter for significant genes first (padj < 0.01), then select top up/downregulated genes

# DELETE TOP GENES AFTER EACH PASS

top_genes <- shrink_3v2_vol %>%
  filter(padj < 0.01) %>%  # Only keep genes with padj < 0.01
  arrange(desc(log2FoldChange)) %>% 
  slice_head(n = 5) %>%   # Top 5 upregulated
  bind_rows(shrink_3v2_vol %>%
              filter(padj < 0.01) %>%  # Apply padj filter again
              arrange(log2FoldChange) %>%
              slice_head(n = 5))   # Top 5 downregulated


### Set thresholds
padj.cutoff <- 0.01         
lfc.cutoff <- 0.59

# 0.59 = 1.5FC, 1 = 2FC

ggplot(shrink_3v2_vol) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = threshold_OE)) +
  ggtitle("KM3 v 2_1.5FC, p<0.01") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  scale_y_continuous(limits = c(0, 60)) +
  geom_text_repel(data = top_genes,  
                  aes(x = log2FoldChange, y = -log10(padj), label = GeneID),
                  size = 3,  # Adjust text size
                  max.overlaps = Inf,  # Ensure all labels are shown
                  box.padding = 0.5,   # Increase space around labels
                  point.padding = 0.3, # Increase space between points and labels
                  min.segment.length = 0) +  # Ensures label connection lines are always shown
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))


summary(shrink_1v2_vol$threshold_OE)

# Thresholds for KM3v2 1.5FC p<0.01 Mode   FALSE 13353    TRUE 1720

# save the lists that pass threholds

KM3v2_1.5Fc_p0.01 <- shrink_3v2_vol %>%
  filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff) 

write.csv(KM3v2_1.5Fc_p0.01, 'KM3v2_1.5Fc_p0.01_15.2.24.csv')


# Repeat for KM1v2 using same paramteters.
# Thresholds for KM1v2 1.5FC p<0.01 Mode   FALSE 14704    TRUE 369


# Gene ontology -----------------------------------------------------------

# Look at gene ontology ---------------------------------------------------


#switch on CLUSTERPROFILER, ANNOTATDBI and ORG.HS.eg.dg
library(org.Hs.eg.db)
library(clusterProfiler)
library(AnnotationDbi)

#load the significant lists from above

KM3v2 <- KM3v2_1_5Fc_p0_01


# 1. Gene Ontology
significant_genes <-KM3v2
significant_genes$entrez <- mapIds(org.Hs.eg.db, 
                                   keys = significant_genes$GeneID, 
                                   column = "ENTREZID", 
                                   keytype = "SYMBOL", 
                                   multiVals = "list")

# this raises a message to say thet 'many mapping between keys and columns
# run this table(sapply(significant_genes$entrez, length))
# shows 1 entrez ID maps twice - it is for gene HBD



#significant_genes <- significant_genes[!sapply(significant_genes$entrez, is.na), ]

gene_list <- significant_genes$entrez

print(gene_list)

go_results <- enrichGO(gene          = gene_list,
                       OrgDb         = org.Hs.eg.db,
                       keyType       = "ENTREZID",
                       ont           = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05)

barplot(go_results, showCategory = 10)
go_results_df <- as.data.frame(go_results)
dotplot(go_results, font.size = 8, showCategory = 20)

# save the dotplot


write.csv(go_results_df, 'GOntology of Biol.Proc.FC1.5_p<0.01 KM3v2.csv')



# Have a look and see which TGFb genes are changed ------------------------

# load list of TGFb genes 
Limited_TGF_genes_incl_BMP_16_2_25 <- TGFlist

# 128 genes - taken from KEGG pathway and annotated by hand

TGF <- Limited_TGF_genes_incl_BMP_16_2_25


# Now add a column to the df to say if they are on the TGF/ BMP genes

#TGF_genes <- significant_genes %>% filter(entrez %in% list_select)
TGF_select <- TGF$gene
shrink_3v2_vol$TGF <- ifelse(shrink_3v2_vol$GeneID %in% TGF_select, 1, 0)

# Reorder data so that TGF genes (1) come last in plotting order
shrink_3v2_vol3 <- shrink_3v2_vol %>%
  arrange(TGF)  # Ensures grey (0) is plotted first, then red (1)


# Select the top 5 upregulated and downregulated inflammatory genes
top_up2 <- shrink_3v2_vol3 %>%
  filter(TGF == 1) %>%
  arrange(desc(log2FoldChange)) %>%
  slice_head(n = 5)

top_down2 <- shrink_3v2_vol3 %>%
  filter(TGF == 1) %>%
  arrange(log2FoldChange) %>%
  slice_head(n = 5)

# Combine both for labeling
top_genes2 <- bind_rows(top_up2, top_down2)

# Volcano plot with labeled genes
ggplot(shrink_3v2_vol3, aes(x = log2FoldChange, y = -log10(padj), color = factor(TGF))) +
  geom_point(alpha = 0.7) +  # Adjust transparency for clarity
  scale_color_manual(values = c("0" = "grey", "1" = "blue")) +  # Custom colors
  geom_text_repel(  # Add labels
    data = top_genes2,
    aes(label = GeneID),  # Replace 'GeneID' with the column containing gene names
    size = 3,  # Adjust text size
    box.padding = 0.7,
    point.padding = 0.5,
    min.segment.length = 0.5,
    max.overlaps = Inf
  ) +
  labs(
    x = "log2 fold change",
    y = "-log10 adjusted p-value",
    color = "TGF/BMP"
  ) +
  theme_minimal() +  # Clean theme
  theme(legend.position = "right")  # Adjust legend position


# Additional amended code for including CIBERsort in model ----------------

# Load CIBERsort data and select cell types with detectable expression
CIBER <- CIBER %>% select(CohortID, ciber.B.cells.naive, ciber.B.cells.memory, ciber.Plasma.cells, ciber.T.cells.CD8, ciber.T.cells.CD4.naive, ciber.T.cells.CD4.memory.resting, ciber.T.cells.CD4.memory.activated, ciber.T.cells.follicular.helper, ciber.T.cells.regulatory..Tregs., ciber.T.cells.gamma.delta, ciber.NK.cells.resting, ciber.NK.cells.activated, ciber.Monocytes, ciber.Macrophages.M0, ciber.Macrophages.M1, ciber.Macrophages.M2, ciber.Dendritic.cells.resting, ciber.Dendritic.cells.activated, ciber.Mast.cells.resting, ciber.Mast.cells.activated, ciber.Eosinophils, ciber.Neutrophils)


# This becomes teh Design formula

design_formula <- ~ KM3 + ciber.B.cells.naive + ciber.B.cells.memory + ciber.Plasma.cells + ciber.T.cells.CD8 +  ciber.T.cells.CD4.naive + ciber.T.cells.CD4.memory.resting + ciber.T.cells.CD4.memory.activated  + ciber.NK.cells.resting + ciber.Monocytes   +  ciber.Dendritic.cells.activated + ciber.Mast.cells.resting +  ciber.Eosinophils + ciber.Neutrophils


Biomark_dds_raw <- DESeqDataSetFromMatrix(countData = Biomark_counts_matrix,
                                          colData = model,
                                          design =  design_formula)

# Continue DESEQ2 as normal
