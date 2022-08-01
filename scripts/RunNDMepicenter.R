source("scripts/functions/get_mean_rank_matrix.R")
source("scripts/functions/get_atrophy_t_vals_atlas.R")
source("scripts/functions/get_atrophy_t_vals_voxel.R")
source("scripts/functions/get_mean_hc_connectome.R")
source("scripts/functions/get_mean_hc_func_connectome.R")
source("scripts/functions/get_node_distance.R")
source("scripts/functions/get_nbs_connectome.R")
source("scripts/functions/apply_spin_nulls.R")
source("/scripts/functions/plot_on_brain_pysurf.R")
source("scripts/run_LDM.R")

#####################################
# compute spin nulls
######################################

type <- "illness"
timepoint <- "bl"
rank_index <- 1 ### 1=bl 2=3m 3=12m for illness effects or 1=3m 2=12m for medication effects
contrast <- ""
hc_connectome <- get_mean_hc_connectome(method = "traditionalConsistency", type = "001",sample = "gencog", age_match = T) #preserveDist #traditionalConsistency
hc_connectome_w<- get_mean_hc_connectome(method = "average", type = "0.001", sample = "gencog", age_match = T)
hc_connectome <- hc_connectome*hc_connectome_w
hc_funconnectome <- get_mean_hc_func_connectome(GSR = F, sample = 'gencog')

#observed mean t
mean_rank_matrix <- get_mean_rank_matrix(type = type, rank = FALSE, sample = 'gencog')
rankmean <- mean_rank_matrix[,rank_index] ### 1=bl 2=3m 3=12m
rankmean <- head(rankmean, -32)


#nulls
atrophy_nulls <- get_long_t_scores_dbm(smoothing = 3, fdr = F,
                                       timepoint = timepoint, 
                                       stat = "T",
                                       contrast = contrast, 
                                       p = 1, 
                                       remove.fov.rois = FALSE,
                                       zscore = T)

# remove subcortex  (no need to remove low fov rois as the apply_spin_nulls function will do that)
atrophy_nulls <- head(atrophy_nulls, -32)


nulls <- apply_spin_nulls(type = 'vasa', vector = as.matrix(atrophy_nulls), p =10000,remove.fov.rois = F)

# remove subcortex from sc and fc mat
hc_connectome <- hc_connectome[1:300,1:300]
hc_funconnectome <- hc_funconnectome[1:300,1:300]

# apply LDM to nulls to extract neighbor atrophy values
library(pbapply)
library(parallel)
cl <- makeCluster(8)
results_null <- pbapply(nulls, 2, run_LDM, 
                        sc = as.matrix(hc_connectome),
                        fc = hc_funconnectome, 
                        weight_by_sc = T,
                        weight_by_fc = F,
                        cor.type = "spearman",
                        null.type = NULL,
                        scale.nulls = F,
                        sample = "gencog", 
                        cl = cl)


# select neighbor atrophy values from the list of results 
neighbour_atrophy_nulls <- cbind(lapply(results_null,  "[[",2))

# compute mean node/neighbour atrophy for each spin
mean_nulls <- matrix(nrow = 300, ncol=10000)
for (s in 1:10000){
  mean_nulls[,s] <-  (nulls[,s] + neighbour_atrophy_nulls[[s]])/2
}

absmax <- function(x) { x[which.max( abs(x) )]}
pvec <- vector()
max_nulls <- vector()
for (s in 1:300){
  pvec[s] <- sum(abs(mean_nulls[s,]) >= abs(rankmean[s]))/length(mean_nulls[s,])
  max_nulls[s] <- absmax(mean_nulls[s,])
}

# get fwe corrected pvals
pvec_fwe <- vector()
for (s in 1:300){
  pvec_fwe[s] <- sum(abs(max_nulls) >= abs(rankmean[s]))/length(max_nulls)
}

# plot mean rank on brain 
# mask by pvals
pvec_fwe[pvec_fwe==1] <- NA
pvec_fwe[pvec_fwe<0.05] <- 1
pvec_fwe[pvec_fwe!=1] <- 0
pvec_fwe[is.na(pvec_fwe)]  <- 0

mean_rank_matrix <- get_mean_rank_matrix(type = type, rank = FALSE,sample = "gencog")
plot_vector <- mean_rank_matrix[,rank_index] 
plot_vector <- head(plot_vector, -32)
plot_vector <- c(plot_vector*pvec_fwe)
plot_vector <- round(plot_vector,3)
plot_vector <- c(plot_vector, rep(0,32)) # add 0s for subcortex
x <- plot_on_brain_pysurf(plot_vector = plot_vector, 
                          min = min(mean_rank_matrix[-301:-332,rank_index]), 
                          max = max(mean_rank_matrix[-301:-332,rank_index]),
                          colourscale = "coolwarm", 
                          jpeg = F, 
                          surf = "inflated", 
                          no_subcortex = T)


#Identify anatomical location names of epicenters 
source("scripts/functions/get_node_distance.R")
mean_rank_matrix <- get_mean_rank_matrix(type = "illness", rank = FALSE,sample = "gencog")
rankmean <- mean_rank_matrix[,1] ### 1=bl; 2=3m; 3=12m
x <- get_labs(remove_fov_rois = F)
row.names(mean_rank_matrix) <- x$ROI.Name

xyz <- get_COG(all_rois = T)
row.names(xyz) <- x$ROI.Name
xyz["7Networks_RH_Default_Temp_6",]
library("label4MRI")
label4MRI::mni_to_region_name()
