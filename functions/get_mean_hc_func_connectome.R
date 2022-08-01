get_mean_hc_func_connectome <- function(GSR=FALSE, sample = NULL){
  if (sample == "stages") {
    if(GSR==FALSE) {filepath = "/data/ts_p300_n7/fmriprep_aroma_2p_hpfdt_noSmooth/"}
    if(GSR==TRUE) {filepath = "/data/ts_p300_n7/fmriprep_aroma_2p_hpfdt_gmr_NoSmooth/"}
    setwd(filepath)
    #load in data into 3d array 
    connectomes <- list.files(pattern = "c_1", full.names = TRUE, recursive = F)  
    qc_ts_a <- lapply(connectomes, function(i){
      data.table::fread(i)
    })
    rparc = 332
    nscan = 27
    all.ts.array.a <- array(numeric(),c(rparc,rparc,nscan))  
    for (i in 1:length(qc_ts_a)) {
      temp <- as.data.frame(qc_ts_a[i])
      all.ts.array.a[ , , i] <- cor(t(temp), use = "pairwise.complete.obs")
    }
    
    con_mat_fov <- rep(NaN, 319*319*length(connectomes));  
    con_mat_fov <-array(all.ts.array.a, c(319, 319, length(connectomes))); 
    removed_roi <- c(50,51,52, 53, 86, 90, 197,199, 200, 201, 242, 243, 244)
    for (i in 1:dim(all.ts.array.a)[3]) {
      temp <- as.data.frame(all.ts.array.a[, , i])
      temp <- temp[-c(removed_roi),-c(removed_roi)] 
      con_mat_fov[, , i] <- as.matrix(temp)
    }
    meanmat <- apply(con_mat_fov, c(1,2), mean)
  }
  if(sample== "gencog") {
    if(GSR==FALSE) {filepath = "/data/GOC_FC/conmat_fix_nogsr/"}
    if(GSR==TRUE) {filepath = "/data/GOC_FC/conmat_fix_gsr/"}
    setwd(filepath)
    #load in data into 3d array 
    connectomes <- list.files(pattern = ".txt", full.names = F, recursive = F) 
    subjid <- stringr::str_split(connectomes, pattern = ".txt", simplify = T)[,1]
    #select only age matched subjects
    age_matched_list <- read.table("~/Dropbox/Sid/R_files/STAGES_difussion/data/connectomes_GOC/age_matched_list.txt")
    index <- which(subjid %in% age_matched_list$V1)
    connectomes_subset <- connectomes[index]
    conmats <- lapply(connectomes_subset, function(i){
      data.table::fread(i)
    })
    rparc = 332
    nscan = length(conmats)
    conmat_array <- array(numeric(),c(rparc,rparc,nscan)) 
    for (i in 1:nscan) {
     conmat_array[ , ,i] <- as.matrix(conmats[[i]])
    }
    meanmat <- apply(conmat_array, c(1,2), mean)
  }
  return(meanmat)
}

