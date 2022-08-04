#Network Deformation Model 

apply_LDM_all <- function(weight_by_sc = NULL, weight_by_fc = NULL, cor.type="spearman", scale.nulls = F){
  source("scripts/functions/get_atrophy_t_vals_voxel.R")
  source("scripts/functions/get_mean_hc_connectome.R")
  source("scripts/functions/get_mean_hc_func_connectome.R")
  source("scripts/functions/get_node_distance.R")
  source("scripts/functions/apply_spin_nulls.R")
  source("scripts/run_LDM.R")
  hc_connectome <- get_mean_hc_connectome(method = "traditionalConsistency", type = "001") #type 001 = COMMIT2
  hc_connectome_w <- get_mean_hc_connectome(method = "average")
  hc_connectome <- as.matrix(hc_connectome*hc_connectome_w)
  hc_funconnectome <- get_mean_hc_func_connectome(GSR = F)
  
  if(scale.nulls ==T){
    temp <- hc_connectome != 0
    hc_connectome[temp] <- scale(hc_connectome[temp])
    
    temp <- hc_funconnectome != 0
    hc_funconnectome[temp] <- scale(hc_funconnectome[temp])
  }
  #load in parameter nulls
  #param_nulls <- list()
  #param_nulls[[5]] <- read.table("data/nulls/parameterized_nulls/parameterized_nulls_12m_medication.txt")
  #saveRDS(param_nulls, "data/nulls/parameterized_nulls/parameterized_nulls_list.RDS")
  param_nulls <- readRDS('data/nulls/parameterized_nulls/parameterized_nulls_list.RDS')
  if(scale.nulls == T) {
    param_nulls <- lapply(param_nulls, scale)
  }
  # Define list data structure to hold results and nulls
  All_results <- vector(mode = "list", length = 5)
  library(pbapply)
  library(parallel)
  cl <- makeCluster(7)
  index <- 1
  for (c in c("illness", "medication")){
    if (c =="illness") {time_vec <- c("bl", "bl_3m", "bl_12m")}
    if (c =="medication") {time_vec <- c("bl_3m", "bl_12m")}
    for (t in time_vec) {
      if(t == "bl") {contrast <- ""}
      if(all(c(t == "bl_3m" | t == "bl_12m"), c(c == "illness"))) {contrast <- "pipt_v_hc"}
      if(all(c(t == "bl_3m" | t == "bl_12m"), c(c == "medication"))) {contrast <- "mipt_v_pipt&hc"}
      atrophy <- get_long_t_scores_dbm(smoothing = 3, 
                                       fdr = F,
                                       timepoint = t, 
                                       stat = "T",
                                       contrast = contrast, 
                                       p = 1, 
                                       zscore = T, 
                                       remove.fov.rois = T)
      #run ldm + connectome nulls
      All_results[[index]] <- run_LDM(sc = as.matrix(hc_connectome), 
                                      fc = hc_funconnectome, 
                                      atrophy = atrophy,
                                      weight_by_sc = weight_by_sc,
                                      weight_by_fc = weight_by_fc,
                                      cor.type = cor.type,
                                      null.type = "bins_10_swap_50000", 
                                      scale.nulls = scale.nulls, 
                                      sample = "stages")
      
      
      #Param nulls
      pnulls  <- as.matrix(param_nulls[[index]][-c(50,51,52, 53, 86, 90, 197,199, 200, 201, 242, 243, 244),]) #remove low FOV rois
      results_null <- pbapply(pnulls, 2, run_LDM, 
                              sc = as.matrix(hc_connectome),
                              fc = hc_funconnectome, 
                              weight_by_sc = weight_by_sc,
                              weight_by_fc =weight_by_fc,
                              cor.type = cor.type,
                              null.type = NULL, 
                              scale.nulls = scale.nulls,
                              sample = "stages", 
                              cl = cl)
      
      All_results[[index]][[6]] <- as.vector(unlist((lapply(results_null,  "[[",3))))
      #two tailed p val
      All_results[[index]][[7]]  <- sum(abs(All_results[[index]][[6]])>abs(All_results[[index]][[3]]))/length(All_results[[index]][[6]])
      
      # Spin Nulls #
      
      #observed results for spin nulls (no subcortex)
      atrophy <- get_long_t_scores_dbm(smoothing = 3, 
                                       fdr = F,
                                       timepoint = t, 
                                       stat = "T",
                                       contrast = contrast,
                                       p = 1, 
                                       zscore = F)
      
      #remove subcortex
      atrophy <- head(atrophy, -32)
      hc_connectome_nosub <- hc_connectome[1:287, 1:287]
      hc_funconnectome_nosub <- hc_funconnectome[1:287, 1:287]
      #
      
      results_nosub <- run_LDM(sc = as.matrix(hc_connectome_nosub), 
                               fc = hc_funconnectome_nosub, 
                               atrophy = atrophy,
                               weight_by_sc = weight_by_sc,
                               weight_by_fc = weight_by_fc,
                               cor.type = cor.type)
      All_results[[index]][[8]] <- results_nosub[[3]] #r values when subcortex removed
      
      # compute spin nulls
      atrophy_nulls <- get_long_t_scores_dbm(smoothing = 3, 
                                             fdr = F,
                                             timepoint = t, 
                                             stat = "T",
                                             contrast = contrast, 
                                             p = 1, 
                                             remove.fov.rois = FALSE, 
                                             zscore = F)
      # remove subcortex
      atrophy_nulls <- head(atrophy_nulls, -32)
      nulls <- apply_spin_nulls(type = 'vasa', vector = as.matrix(atrophy_nulls), 1000)
      
      # apply LDM to nulls to extract neighbor atrophy values
      results_null <- pbapply(nulls, 2, run_LDM, 
                              sc = as.matrix(hc_connectome_nosub),
                              fc = hc_funconnectome_nosub, 
                              weight_by_sc = weight_by_sc,
                              weight_by_fc = weight_by_fc,
                              cor.type = cor.type,
                              null.type = NULL, cl = cl)
      
      
      All_results[[index]][[9]] <- as.vector(unlist((lapply(results_null,  "[[",3))))
      #two tailed
      All_results[[index]][[10]]  <- sum(abs(All_results[[index]][[9]])>abs(All_results[[index]][[8]]))/length(All_results[[index]][[9]])
      names(All_results[[index]]) <- c("ROI atrophy", "Neighbour atrophy", 
                                       "correlation", "connectome nulls", 
                                       "connectome p", "param nulls", "param p",
                                       "correlation nosub","spin nulls", "spin p")
      index <- index +1
    }
  }
  return(All_results)
}


#Results
No_weighting_spearman <- apply_LDM_all(weight_by_sc = F, weight_by_fc = F, cor.type = "spearman", scale.nulls  = T)
FC_weighting_spearman  <- apply_LDM_all(weight_by_sc = F, weight_by_fc = T, cor.type = "spearman", scale.nulls = T)
SC_weighting_spearman  <- apply_LDM_all(weight_by_sc = T, weight_by_fc = F, cor.type = "spearman", scale.nulls = T)


all_results <- list(No_weighting_spearman, FC_weighting_spearman ,SC_weighting_spearman)

#saveRDS(all_results, "output/results/all_results_LDM_scaled.RDS")



