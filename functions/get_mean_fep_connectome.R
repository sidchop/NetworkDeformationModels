get_mean_fep_connectome <- function(exclude_fov_rois = TRUE, type = "0.001", method = "preserveDist"){
  if(method == "average") {
    #load in data into 3d array 
    connectomes <- list.files(path = "data/connectomes/", pattern = paste0('*',type,'*'), full.names = TRUE, recursive = TRUE)  
    #exclude subs
    exclude.t1 <- c("011", "025")
    exclude.fov <- c("048","034c", "015", "060") 
    connectomes.qc <- grep(connectomes, 
                           pattern='*cont_*', inv=T, value=T)
    connectomes.qc <- grep(connectomes.qc, 
                           pattern='*034c_*|*011*|*025*|*048*|*034c*|*015|*060*', inv=F, value=T)
    
    #read into a list
    con_mat <- rep(NaN, 332*332*length(connectomes.qc));  
    con_mat <-array(con_mat, c(332, 332, length(connectomes.qc))); 
    for (i in 1:length(connectomes.qc)) {
      con_mat[, , i] <- as.matrix(read.table(connectomes.qc[[i]]))
    }
    
    ###Exclude ROIS due to poor FOV:
    #7Networks_LH_SomMot_26, #7Networks_LH_SomMot_27, #7Networks_LH_SomMot_28, #7Networks_LH_SomMot_29, #7Networks_LH_Limbic_OFC_1, #7Networks_LH_Limbic_TempPole_2, #7Networks_RH_SomMot_24, #7Networks_RH_SomMot_26, #7Networks_RH_SomMot_27, #7Networks_RH_SomMot_28, #7Networks_RH_Limbic_TempPole_1, #7Networks_RH_Limbic_TempPole_2, #7Networks_RH_Limbic_TempPole_3
    con_mat_fov <- rep(NaN, 319*319*length(connectomes.qc));  
    con_mat_fov <-array(con_mat, c(319, 319, length(connectomes.qc))); 
    removed_roi <- c(50,51,52, 53, 86, 90, 197,199, 200, 201, 242, 243, 244)
    for (i in 1:dim(con_mat)[3]) {
      temp <- as.data.frame(con_mat[, , i])
      temp <- temp[-c(removed_roi),-c(removed_roi)] 
      con_mat_fov[, , i] <- as.matrix(temp)
    }
    
    meanmat <- apply(con_mat_fov, c(1,2), mean)
    return(meanmat) 
  }
  if(method == "traditionalConsistency") {
    conmat <- read.csv("data/distanceDependent_betzel/HC_bin_traditional_consistency_based_consensus.txt", header = F)
    return(conmat)
    message("bin mat created using bentzel code `traditional_consistency_based_consensus`")
  }
}




