# Correlation of Biomarker panel with original BMP gene list

library(rstatix)
library(ggplot2)
library(dplyr)
library(tidyr)
library(cocor)  
library(ggplot2) 
library(GGally)
library(tidyverse)
library(stats)
library(pheatmap)


# Import the original biomarkers file
# upload TPM data
TPM_rna_seq <- TPM


TPM <- TPM_rna_seq
TPM$GeneID <-rownames(TPM)
rownames(TPM) <- NULL
colnames(TPM)[duplicated(colnames(TPM))]
colnames(TPM) <- make.unique(colnames(TPM))
tpm_data <- TPM %>% filter(GeneID %in% c('ID2', 'ID3', 'NOTCH1', 'NOTCH2', 'SMAD1', 'SMAD5', 'PTGS2', 'ARL4C', 'APLN', 'CHSY3', 'DLX2', 'EML6', 'FGFR3', 'HEY1', 'HEY2', 'ID1', 'ID4', 'KIT',  'RAI2',  'SMAD6', 'SMAD7', 'SPRY1', 'IFNG'))
tmp <- as.data.frame(t(tpm_data))
colnames(tmp) <- tmp[509,]
tmp <- tmp[-509,]


Biomarkers.in.cohort_Feb24 <-file_rj


TPM_select <- tmp %>% mutate_all(as.numeric)
#make the rownames SampleID
TPM_select$SampleID <- rownames(TPM_select)  

keep <- Biomarkers.in.cohort_Feb24 %>% filter(group == 'PAH') %>% select(SampleID) %>% unique()
keep <- keep$SampleID

#Merge Biomarker Cohort with TPM Select by SampleID
TPM_BMP <- TPM_select %>% filter(SampleID %in% keep)

#356 are present

df_numeric <- TPM_BMP[, sapply(TPM_BMP, is.numeric)]  

# Run correlations in PAH
# Compute correlation matrix
cor_matrix <- cor(df_numeric, use = "pairwise.complete.obs", method = "spearman")

my_palette <- colorRampPalette(c("blue", "white", "red"))(100)
breaks <- seq(-1, 1, length.out = 101)


pheatmap(cor_matrix,
         clustering_distance_rows = "euclidean",
         color = my_palette,
         breaks = breaks,
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         display_numbers = TRUE,      
         main = "BMP genes correlation")


