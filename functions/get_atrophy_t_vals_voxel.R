get_long_t_scores_vbm <- function(smoothing =3, 
                                  p=1, 
                                  fdr=FALSE, 
                                  timepoint = NULL,
                                  contrast = NULL,
                                  stat = NULL,
                                  zscore = FALSE,
                                  remove.fov.rois = TRUE) {
  library(neurobase) 
  source("scripts/functions/nii_roi2ts.R")
  
  if(contrast != "medication") {
    stat_path <- list.files(path = paste0("data/vbm/SPM/s",smoothing ,"/",timepoint,"/",stat,"/",contrast,"/"), pattern = "stat_c", full.names = T)
    l10p_path <-  list.files(path = paste0("data/vbm/SPM/s",smoothing ,"/",timepoint,"/",stat,"/",contrast,"/"), pattern = "stat_lp", full.names = T)
    stat_nii <- fast_readnii(stat_path)
    if(zscore == T){stat_nii <- zscore_img(stat_nii)}
    l10p_nii <- fast_readnii(l10p_path)
    ifelse(p != 1,stat_nii <-   stat_nii*l10p_nii, stat_nii <-   stat_nii)
    
  }
  
  #extract mean zscore cortical
  cotrical_atlas <- fast_readnii(paste0("data/vbm/atlases/Schaefer2018_300Parcels_7Networks_order_FSLMNI152_1.5mm.nii.gz"))
  cotrical_atlas[is.nan(cotrical_atlas)] <- 0
  stat_nii[is.nan(stat_nii)] <- 0
  vector_of_mean_zscores_cortex <- nii_roi2ts(img = stat_nii, net = cotrical_atlas)
  
  #extract mean zscore subcortical
  subcotrical_atlas <- fast_readnii(paste0("data/vbm/atlases/Tian_Subcortex_S2_3T_1.5mm.nii.gz"))
  subcotrical_atlas[is.nan(subcotrical_atlas)] <- 0
  vector_of_mean_zscores_subcortex<- nii_roi2ts(img = stat_nii, net = subcotrical_atlas)
  
  all_zscores <- rbind(vector_of_mean_zscores_cortex, vector_of_mean_zscores_subcortex)
  if(remove.fov.rois==TRUE) {
    removed_roi <- c(50,51,52, 53, 86, 90, 197,199, 200, 201, 242, 243, 244)
    all_zscores  <- all_zscores[-c(removed_roi),]
  }
  all_zscores[is.infinite(all_zscores)] <- 0
  return(all_zscores)
  
}
#stat F is a chi square
get_long_t_scores_dbm <- function(smoothing =3, 
                                  p=1, fdr=FALSE, 
                                  timepoint = NULL,
                                  contrast = NULL,
                                  stat = NULL,
                                  zscore = FALSE,
                                  remove.fov.rois = TRUE,
                                  psychosis_only = FALSE) {
  library(neurobase)
  setwd("~/Dropbox/Sid/R_files/STAGES_difussion/")
  source("scripts/functions/nii_roi2ts.R")
  
  
  if(contrast != "medication") {
    if(psychosis_only==FALSE) {
      stat_path <- list.files(path = paste0("data/dbm/SPM/s",smoothing ,"/",timepoint,"/",stat,"/",contrast,"/"), pattern = "stat_c", full.names = T)
      l10p_path <-  list.files(path = paste0("data/dbm/SPM/s",smoothing ,"/",timepoint,"/",stat,"/",contrast,"/"), pattern = "stat_lp", full.names = T)
    }
    if(psychosis_only==TRUE) {
      stat_path <- list.files(path = paste0("data/dbm/SPM/s",smoothing ,"/psychosis_only/",timepoint,"/",stat,"/",contrast,"/"), pattern = "stat_c", full.names = T)
      l10p_path <-  list.files(path = paste0("data/dbm/SPM/s",smoothing ,"/psychosis_only/",timepoint,"/",stat,"/",contrast,"/"), pattern = "stat_lp", full.names = T)
    }
    stat_nii <- fast_readnii(stat_path)
    if(zscore == T){stat_nii <- zscore_img(stat_nii)}
    l10p_nii <- fast_readnii(l10p_path)
    ifelse(p != 1,stat_nii <-   stat_nii*l10p_nii, stat_nii <-   stat_nii)
    
    
  }
  
  
  #extract mean zscore cortical
  cotrical_atlas <- fast_readnii(paste0("data/vbm/atlases/Schaefer2018_300Parcels_7Networks_order_FSLMNI152_1.5mm.nii.gz"))
  cotrical_atlas[is.nan(cotrical_atlas)] <- 0
  stat_nii[is.nan(stat_nii)] <- 0
  vector_of_mean_zscores_cortex <- nii_roi2ts(img = stat_nii, net = cotrical_atlas)
  
  #extract mean zscore subcortical
  subcotrical_atlas <- fast_readnii(paste0("data/vbm/atlases/Tian_Subcortex_S2_3T_1.5mm.nii.gz"))
  subcotrical_atlas[is.nan(subcotrical_atlas)] <- 0
  vector_of_mean_zscores_subcortex<- nii_roi2ts(img = stat_nii, net = subcotrical_atlas)
  
  all_zscores <- rbind(vector_of_mean_zscores_cortex, vector_of_mean_zscores_subcortex)
  if(remove.fov.rois==TRUE) {
    removed_roi <- c(50,51,52, 53, 86, 90, 197,199, 200, 201, 242, 243, 244)
    all_zscores  <- all_zscores[-c(removed_roi),]
  }
  all_zscores[is.infinite(all_zscores)] <- 0
  return(all_zscores)
}
