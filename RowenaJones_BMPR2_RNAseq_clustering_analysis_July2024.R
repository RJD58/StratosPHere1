# Code used in 'BMPR-II Biomarkers for testing therapeutic efficacy in pulmonary arterial hypertension – novel findings from the StratosPHere 1 study.' Rowena J Jones

# CODE: Clustering on RNAseq data

# Load the following

library(dplyr)
library(ggplot2)
library(stats)
library(ggResidpanel)
library(ConsensusClusterPlus)
library(factoextra)
library(data.table)
library(magrittr)
library(tidyverse) 
library(cluster)
library(Hmisc)
library(diceR)
library(FactoMineR)
library(tidyr)
library(Publish)


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ConsensusClusterPlus")


# load data as data frame - name 'biomark_clin'



# Log transform and then scale the RNAseq values

biom2 <- biomark_clin

lapply(biom2[,3:10], hist)

biom2[,3:10] <- lapply(biom2[,3:10], log)
biom2[,3:10] <- lapply(biom2[,3:10], scale)
biom2 <- biom2[,c(1,3:10)]
biom2 <- as.data.frame(biom2)

lapply(biom2[,2:9], hist)

# glitchy after lapply - rewrite csv
write.csv(biom2, 'biom2_transform.csv')

biom3 <- biom2_transform

biom3 <- biom3[,-1]

first_col_to_rowname =function(biom3) {
  rownames(biom3) <- biom3[,1]
  biom3 <- biom3[,-1]
  return(biom3)
}

biom3<- first_col_to_rowname(as.data.frame(biom3))


# Generate silhouette and elbow plots for PAM, Kmeans and HC --------------


#PAM
fviz_nbclust(biom3, cluster::pam, method="silhouette")+theme_classic() + ggtitle("PAM silhoutte plot factors")
fviz_nbclust(biom3, cluster::pam, method="wss")+theme_classic() + ggtitle("PAM elbow plot factors")

#Kmeans
fviz_nbclust(biom3, kmeans, method="wss")+theme_classic() + ggtitle("Km elbow plot factors")
fviz_nbclust(biom3, kmeans, method="silhouette")+theme_classic() + ggtitle("Km  silhoutte plot factors")

#HCLUST
fviz_nbclust(biom3, hcut, method="wss")+theme_classic() + ggtitle("Hcut elbow plot factors")
fviz_nbclust(biom3, hcut, method="silhouette")+theme_classic() + ggtitle("Hcut  silhoutte plot factors")



# Run consensus clustering ------------------------------------------------
#Glitchy after last step
write.csv (biom3, 'Biomarkers_log scaled.csv')
biom_mat <- Biomarkers_log_scaled
biom_mat <- as.data.frame(biom_mat)
rownames(biom_mat) <- biom_mat[,1]
biom_mat <- biom_mat[,-1]
biom_mat <- as.matrix(biom_mat)

CC_pam <- ConsensusClusterPlus(t(biom_mat), maxK = 6, reps = 1000, pItem = 0.8, pFeature = 1, clusterAlg ="pam", distance = "euclidean", plot = "png", seed=123, title = "CC results factors pam euclidian nonregressed Biomarkers")
CC_km <-CC_km <- ConsensusClusterPlus(t(biom_mat), maxK = 6, reps = 1000, pItem = 0.8, pFeature = 1, clusterAlg ="km", distance = "euclidean", plot = "png", seed = 123, title = "CC results factors KM euclidian nonregressed Biomarkers")
CC_hc <- ConsensusClusterPlus(t(biom_mat), maxK = 6, reps = 1000, pItem = 0.8, pFeature = 1, clusterAlg ="hc", distance = "euclidean", plot = "png", seed = 123, title = "CC results factors Hclust euclidian nonregressed Biomarkers")

PAC(CC_km[[2]][["consensusMatrix"]])
PAC(CC_km[[3]][["consensusMatrix"]])
PAC(CC_km[[4]][["consensusMatrix"]])
PAC(CC_pam[[2]][["consensusMatrix"]])
PAC(CC_pam[[3]][["consensusMatrix"]])
PAC(CC_pam[[4]][["consensusMatrix"]])

KM_3_clust <- eclust(biom_mat,FUNcluster="kmeans", k=3,hc_metric = "euclidean", seed=123)
fviz_cluster(KM_3_clust, geom='point')

KM_3_clust_df <- as.data.frame(KM_3_clust$cluster)

# sort the row names back into a column for merging later on
KM_3_clust_df2 <- KM_3_clust_df
KM_3_clust_df2$OpenClinicaID <- row.names(KM_3_clust_df2) 
rownames(KM_3_clust_df2) <- NULL

write.csv(KM_3_clust_df2, 'KM_3_clust_RNA.csv')




# Merge clusters with original data file ----------------------------------

#Make clin_clust file with cluster list and original data merged

clin_clust <- merge(KM_3_clust_df2, biomark_clin, by = 'OpenClinicaID')

clin_clust <- clin_clust %>% 
  rename("KM3" = 'KM_3_clust$cluster')


# Make sure cluster column is a factor and RNAseq data are numeric
# Example
clin_clust$KM3 <- as.factor(clin_clust$KM3)
clin_clust$ID3 <- as.numeric(clin_clust$ID3)
clin_clust$ID2 <- as.numeric(clin_clust$ID2)
clin_clust$NOTCH1 <- as.numeric(clin_clust$NOTCH1)
clin_clust$NOTCH2 <- as.numeric(clin_clust$NOTCH2)
clin_clust$SMAD1.x <- as.numeric(clin_clust$SMAD1.x)
clin_clust$SMAD5 <- as.numeric(clin_clust$SMAD5)
clin_clust$ARL4C <- as.numeric(clin_clust$ARL4C)
clin_clust$PTGS2 <- as.numeric(clin_clust$PTGS2)

clin_clust%>%
  count(KM3)

write.csv(clin_clust, 'KM=3 clin clust.csv')



# Generate data tables based on clusters ----------------------------------

u1 <- univariateTable(KM3 ~ ID2 + ID3 + NOTCH1 + NOTCH2 +  SMAD5 + SMAD1.x + PTGS2 + ARL4C , data= clin_clust, show.totals = TRUE, column.percent = TRUE, compare.groups = TRUE)
KM3_biomarkers_levels<-  summary(u1)
print(u1)

u2 <- univariateTable(KM3 ~ sex + age_diagnosis + BMPR2_gwas + hb_pawp_m + hb_pap_m + hb_rap_m + hb_cardiac_output_value_1 + hb_pvr_calc + hb_cardiac_index_value_1 + ep_1_distance_meters +pvr + cbt_card_ntprobnp_ngpl, data= clin_clust, show.totals = TRUE, column.percent = TRUE, compare.groups = TRUE)
KM3_biomarkers_clincharacters <-  summary(u2)
print(u2)

write.csv(KM3_biomarkers_clincharacters, 'KM3_biomarkers_clincharacters.csv')
write.csv(KM3_biomarkers_levels, 'KM3_biomarkers_nonregressed_level.csv')





# Perform Kaplan Meier on clustered data ----------------------------------

library(rstatix)
library(broom)
library(survminer)
library(data.table)
install.packages('survival')
library(survival)
library(stats)


## Cohort Study start date
start_of_study <- ymd("2013-06-01")
## Census date
date_of_census <- ymd("2022-07-01")

# something glitched in dates, added format to show programme what the original is

clin_clust$DOB <- as.Date(clin_clust$DOB, format="%d/%m/%Y")
clin_clust$sub_date <- as.Date(clin_clust$sub_date,  format="%d/%m/%Y")
#Generate diagnosis data = DOB + age of diagnosis
clin_clust$diagnosis_date <- clin_clust$DOB + (clin_clust$age_diagnosis*365.2422)
clin_clust$diagnosis_date <- as.Date(clin_clust$diagnosis_date)

clin_clust <- unique(clin_clust)

clin_clust$sub_date <- as.Date(clin_clust$sub_date)

clin_clust$sub_date2 <- fifelse(is.na(clin_clust$sub_date), ymd("2022-07-01"), clin_clust$sub_date)
#Determine survival time to census date
clin_clust$surv_time <- (clin_clust$sub_date2 - clin_clust$diagnosis_date)
clin_clust[c('surv_time', 'day')] <- str_split_fixed(clin_clust$surv_time, ' ', 2)
clin_clust$surv_time <- as.numeric(clin_clust$surv_time)
clin_clust$surv_time <- clin_clust$surv_time/365.2422
#Censor the information 0 = alive, lung transplant or dropout. NAs are assumed alive = 0
clin_clust$event <- ifelse(clin_clust$sub_cause == 'death', 1, 0)
clin_clust$event <- ifelse(is.na(clin_clust$event), 0, clin_clust$event)

# DELETE this has an NA and oftern causes a glitch
clin_clust <- clin_clust[-282,]

#Generate survival object. Type = right = right censoring
s_obj1 <- Surv(clin_clust$surv_time, clin_clust$event, type='right')

p <-list()

s_fit_KM3 <- survfit(s_obj1 ~ KM3, data=clin_clust)
p[[1]] <- ggsurvplot(s_fit_KM3, data=clin_clust, pval=TRUE, pval.method = TRUE, risk.table = TRUE, conf.int = TRUE, title='survival between KM=3 clusters non regressed in all cause PAH', xlim=c(0,15), legend.title='cluster', break.x.by=5, legend.labs=c('1', '2', '3'))
print(p[[1]])

#Generate a cox PH model

cox_biom <- coxph(s_obj1 ~ KM3 + sex + age_diagnosis, data=clin_clust)
ggforest(cox_biom, data=clin_clust)


# Run this below to change the reference cluster to cluster 2

clin_clust2 <- clin_clust

levels(clin_clust2$KM3) # [1]'1' '2' '3'

clin_clust2$KM3 <- relevel(clin_clust2$KM3, ref = '2')

cox_biom2 <- coxph(s_obj1 ~ KM3 + sex + age_diagnosis, data=clin_clust2)
ggforest(cox_biom2, data=clin_clust2)





# Bootstrap the survival finding ------------------------------------------

# Bootstrap the code ------------------------------------------------------


# https://bookdown.org/jgscott/DSGI/the-bootstrap.html
# https://stat.ethz.ch/R-manual/R-patched/library/boot/html/censboot.html
# https://journals.sagepub.com/doi/10.1177/0962280217751518

install.packages('boot')
library(boot)

# make dataframe of clin_clust and select relevant info OpenClinica, KM3, surv_time, event

df <- clin_clust %>%
  select(OpenClinicaID, KM3, surv_time, event)

df$KM3 <- as.factor(df$KM3)


rjd.fun <- function(df) {
  surv <- survfit(Surv(surv_time, event) ~ KM3, data = df) #makes survival model
  out <- NULL
  st <- 1
  for (s in 1:length(surv$strata)) {
    inds <- st:(st + surv$strata[s]-1)
    md <- min(surv$time[inds[1-surv$surv[inds] >= 0.5]])
    st <- st + surv$strata[s]
    out <- c(out, md)
  }
  out
}

# run the censboot resampling 1000 times (simple resampling with replacement from the observations)
rjd.case <- censboot(df, rjd.fun, R = 1000, strata = df$KM3)


# rjd.case is a list, the interesting data is 't'

median_surv <- as.data.frame(rjd.case$t)
# make numeric
median_surv[,1:3] <- lapply(median_surv[,1:3], as.numeric)


# This now needs making into a usable table

m2 <- median_surv
# first transpose
m2 <- t(m2)
m2 <-as.data.frame(m2)

#add a column for cluster ID
m2$cluster <- rownames(m2)

m2 <- m2 %>% pivot_longer(1:1000, values_to = 'median_surv', names_to = 'v')


# draw a box plot 
library(ggplot2)
library(stats)

ggplot(m2,
       aes(x=cluster,y=median_surv)) +
  geom_boxplot() 
+
geom_jitter()


kruskal.test(median_surv ~ cluster, data = m2)

dunn_test(median_surv ~ cluster, data = m2)



