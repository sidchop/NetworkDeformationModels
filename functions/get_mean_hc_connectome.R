get_mean_hc_connectome <- function(type = "001", 
                                   method = "traditionalConsistency", 
                                   sample="gencog", 
                                   age_match = TRUE) {
  
  if (sample == "stages") {
    if(method == "average") {
      #load in data into 3d array 
      connectomes <- list.files(path = "data/connectomes/", pattern = paste0('*',type,'*'), full.names = T, recursive = TRUE)  
      #exclude subs
      exclude.t1 <- c("011", "025")
      exclude.fov <- c("048","034c", "015", "060") 
      connectomes.qc <- grep(connectomes, 
                             pattern='c/', inv=F, value=T)
      connectomes.qc <- grep(connectomes.qc, 
                             pattern='034c_', inv=T, value=T)
      
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
      # R.matlab::writeMat("scripts/distanceDependent_betzel/STAGES_HC_base_conmat_array.mat", A=con_mat_fov)
      meanmat <- apply(con_mat_fov, c(1,2), mean)
      return(meanmat) 
    }
    
    i
    if(method == "traditionalConsistency" &  type == "001") {
      conmat <- read.csv("scripts/distanceDependent_betzel/STAGES_HC_bin_traditional_consistency_based_consensus.txt", header = F)
      # pheatmap::pheatmap(conmat, cluster_cols = F, cluster_rows = F)
      return(conmat)
      message("bin mat created using bentzel code `traditional_consistency_based_consensus`")
    }
  
    if(method == "traditionalConsistency" &  type == "base") {
      conmat <- read.csv("scripts/distanceDependent_betzel/STAGES_HC_bin_base_traditional_consistency_based_consensus.txt", header = F)
      # pheatmap::pheatmap(conmat, cluster_cols = F, cluster_rows = F)
      return(conmat)
      message("bin mat created using bentzel code `traditional_consistency_based_consensus`")
    }
  }
  if (sample == "gencog") {
    #remove excluded subs
    connectomes <- list.files(path = "data/connectomes_GOC/", pattern = paste0('*',type,'*'), full.names = T, recursive = TRUE)  
    subjids <- stringr::str_split(connectomes, pattern = "/", simplify = T)[,4]
    excluded <- c("S103", "S11", "S13", "S147", "S156" ,"S206", "S24", 
                  "S45", "S271", "S275", "S343", "S362", "S377", "S400",
                  "S402", "S432", "S441", "S45", "S488", "S7", "S70", 
                  "S71") #scans exsluded due to QC;
    connectomes <- connectomes[-which(subjids %in% excluded)]
    subjids <- stringr::str_split(connectomes, pattern = "/", simplify = T)[,4] #update subj id
    
    if(method == "average") {
      if(age_match == T) {
        age_matched_list <- read.table("data/connectomes_GOC/age_matched_list.txt") #load in age matched list
        index <- which(subjids %in% age_matched_list$V1)
        connectomes_subset <- connectomes[index]
      } else { connectomes_subset <- connectomes }
      con_mat <- rep(NaN, 332*332*length(connectomes_subset))  
      con_mat <-array(con_mat, c(332, 332, length(connectomes_subset)))
      for (i in 1:length(connectomes_subset)) {
        con_mat[, , i] <- as.matrix(read.table(connectomes_subset[[i]]))
      }
      
      meanmat <- apply(con_mat, c(1,2), mean)
      return(meanmat)
    }
  
    if(method == "traditionalConsistency" &  type == "001" & age_match == T) {
      conmat <- read.csv("scripts/distanceDependent_betzel/GOC_HC_bin_traditional_consistency_based_consensus.txt", header = F)
      # pheatmap::pheatmap(conmat, cluster_cols = F, cluster_rows = F)
      return(conmat)
      message("bin mat created using bentzel code `traditional_consistency_based_consensus`")
    }

    if(method == "traditionalConsistency" &  type == "001" & age_match == F) {
      conmat <- read.csv("scripts/distanceDependent_betzel/GOC_allsubs_HC_bin_traditional_consistency_based_consensus.txt", header = F)
      # pheatmap::pheatmap(conmat, cluster_cols = F, cluster_rows = F)
      return(conmat)
      message("bin mat created using bentzel code `traditional_consistency_based_consensus`")
    }
  }
}

