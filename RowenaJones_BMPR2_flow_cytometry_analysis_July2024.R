# Code used in 'BMPR-II Biomarkers for testing therapeutic efficacy in pulmonary arterial hypertension – novel findings from the StratosPHere 1 study.' Rowena J Jones

# CODE: Analysis of cell surface BMPR-II in leukocyte populations

# Load the following

library(dplyr)
library(ggplot2)
library(stats)
library(ggResidpanel)
library(factoextra)
library(diceR)
library(rstatix)
library(Publish)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ConsensusClusterPlus")


#Import data - name 'flow'




# Make linear models and generate residuals of data -----------------------
# Linear models are generated to correct for antibody batch effect

lm_B <- lm(Bcell ~  Batch, data=flow)
lm_NCM <- lm(NC_Mono ~  Batch, data=flow)
lm_CM <- lm(CL_Mono ~  Batch, data=flow)
lm_DC <- lm(DC ~  Batch, data=flow)

anova(lm_B)
anova(lm_CM)
anova(lm_NCM)
anova(lm_DC)


# make a residuals file 

flow_resid <- flow


flow_resid$resid_Bcell <- residuals(lm_B)
resid_panel(lm_B,
            plots = c("resid", "qq", "ls", "cookd"),
            smoother = FALSE)


flow_resid$resid_clmono <- residuals(lm_CM)
resid_panel(lm_CM,
            plots = c("resid", "qq", "ls", "cookd"),
            smoother = FALSE)

flow_resid$resid_NCM <- residuals(lm_NCM)
resid_panel(lm_NCM,
            plots = c("resid", "qq", "ls", "cookd"),
            smoother = FALSE)

flow_resid$resid_DC <- residuals(lm_DC)
resid_panel(lm_DC,
            plots = c("resid", "qq", "ls", "cookd"),
            smoother = FALSE)

# compare histograms of raw data vs residuals
# example code:
#lapply(flow_resid[,4:7], hist)
#lapply(flow_resid[,9:12], hist)


# Continue with just residuals and relevant sample ID labelling
flow_resid2 <- flow_resid[, c(1,9:12)]

first_col_to_rowname =function(flow_resid2) {
  rownames(flow_resid2) <- flow_resid2[,1]
  flow_resid2 <- flow_resid2[,-1]
  return(flow_resid2)
}
flow_resid2<- first_col_to_rowname(as.data.frame(flow_resid2))



# Generate elbow and silhouette plots for PAM, Kmeans and HClust ----------


#PAM
fviz_nbclust(flow_resid2, cluster::pam, method="silhouette")+theme_classic() + ggtitle("PAM silhoutte plot factors")
fviz_nbclust(flow_resid2, cluster::pam, method="wss")+theme_classic() + ggtitle("PAM elbow plot factors")


#Kmeans
fviz_nbclust(flow_resid2, kmeans, method="wss")+theme_classic() + ggtitle("Km elbow plot factors")
fviz_nbclust(flow_resid2, kmeans, method="silhouette")+theme_classic() + ggtitle("Km  silhoutte plot factors")


#HCLUST
fviz_nbclust(flow_resid2, hcut, method="wss")+theme_classic() + ggtitle("Hcut elbow plot factors")
fviz_nbclust(flow_resid2, hcut, method="silhouette")+theme_classic() + ggtitle("Hcut  silhoutte plot factors")


# Need to save and reload the residuals file as can be glitchy

write.csv (flow_resid2, 'Flow_resid.csv')
resid_mat <- Flow_resid
resid_mat <- as.data.frame(resid_mat)
rownames(resid_mat) <- resid_mat[,1]
resid_mat <- resid_mat[,-1]
resid_mat <- as.matrix(resid_mat)


# Run consensus clustering using PAM, KMeans and HClust -------------------



CC_pam_resid <- ConsensusClusterPlus(t(resid_mat), maxK = 6, reps = 1000, pItem = 0.8, pFeature = 1, clusterAlg ="pam", distance = "euclidean", plot = "png", seed=123, title = "CC results factors pam euclidian Flow_resid")
CC_km_resid <- ConsensusClusterPlus(t(resid_mat), maxK = 6, reps = 1000, pItem = 0.8, pFeature = 1, clusterAlg ="km", distance = "euclidean", plot = "png", seed = 123, title = "CC results factors KM euclidian Flow_resid")
CC_hc_resid <- ConsensusClusterPlus(t(resid_mat), maxK = 6, reps = 1000, pItem = 0.8, pFeature = 1, clusterAlg ="hc", distance = "euclidean", plot = "png", seed = 123, title = "CC results factors Hclust euclidian Flow_resid")


# Check PAC scores

PAC(CC_km_resid[[2]][["consensusMatrix"]])
PAC(CC_km_resid[[3]][["consensusMatrix"]])
PAC(CC_km_resid[[4]][["consensusMatrix"]])
PAC(CC_pam_resid[[2]][["consensusMatrix"]])
PAC(CC_pam_resid[[3]][["consensusMatrix"]])
PAC(CC_pam_resid[[4]][["consensusMatrix"]])
PAC(CC_hc_resid[[2]][["consensusMatrix"]])
PAC(CC_hc_resid[[3]][["consensusMatrix"]])
PAC(CC_hc_resid[[4]][["consensusMatrix"]])


#DELETE THIS - NOTE to Eck - Km =3 has a slightly higher PAC score than PAM=3. However CDF plots etc look much better for Kmeans and clustering is more disctinct. I. have run PAM=3 and similar findings

#look at the clustering using optimal algorithm
KM_3_clust_resid <- eclust(resid_mat,FUNcluster="kmeans", k=3,hc_metric = "euclidean", seed=123)
fviz_cluster(KM_3_clust_resid, geom='point')


KM_3_clust_resid_df <- as.data.frame(KM_3_clust_resid$cluster)

# sort the row names back into a column for merging later on
KM_3_clust_resid_df <- KM_3_clust_resid_df
KM_3_clust_resid_df$SampleID <- row.names(KM_3_clust_resid_df) 
rownames(KM_3_clust_resid_df) <- NULL

KM_resid_merge <- merge(KM_3_clust_resid_df, flow, by = 'SampleID')

KM_resid_merge <- KM_resid_merge %>% 
  rename("KM3" = 'KM_3_clust_resid$cluster')

write.csv(KM_resid_merge, 'KM=3_residuals_merged.csv')




KM_resid_merge$KM3 <- as.factor(KM_resid_merge$KM3) 

# To be clear - this is Batch corrected but NOT log transformed scaled data clustered using BOTH controls and PAH. PAM=4 is semi-optimal and has been merged back with the basic data file

ggplot(KM_resid_merge,
       aes(x=KM3,y=Bcell)) +
  geom_boxplot() 
ggplot(KM_resid_merge,
       aes(x=KM3,y=CL_Mono)) +
  geom_boxplot() 
ggplot(KM_resid_merge,
       aes(x=KM3,y=NC_Mono)) +
  geom_boxplot() 
ggplot(KM_resid_merge,
       aes(x=KM3,y=DC)) +
  geom_boxplot() 


kruskal.test(Bcell ~ KM3, data = KM_resid_merge)
dunn_test(Bcell ~ KM3, data = KM_resid_merge)
kruskal.test(CL_Mono ~ KM3, data = KM_resid_merge)
dunn_test(CL_Mono ~ KM3, data = KM_resid_merge)
kruskal.test(NC_Mono ~ KM3, data = KM_resid_merge)
dunn_test(NC_Mono ~ KM3, data = KM_resid_merge)
kruskal.test(DC ~ KM3, data = KM_resid_merge)
dunn_test(DC ~ KM3, data = KM_resid_merge)


#Generate results tables
u1 <- univariateTable(KM3 ~ Group + Batch + Bcell + CL_Mono + NC_Mono + DC , data= KM_resid_merge, show.totals = TRUE, column.percent = TRUE, compare.groups = TRUE)
KM3_res<-  summary(u1)
print(u1)

#add Q(x) to pull out median and IQR to cell populations

u2 <- univariateTable(KM3 ~ Group + Batch + Q(Bcell) + Q(CL_Mono) + Q(NC_Mono) + Q(DC) , data= KM_resid_merge, , show.totals = TRUE, column.percent = TRUE, compare.groups = TRUE)
KM3_res_median<-  summary(u2)
print(u2)
write.csv(KM3_res_median, 'KM3_residuals_medians.csv')
write.csv(KM_resid_merge, 'KM3_residuals_clusters.csv')
