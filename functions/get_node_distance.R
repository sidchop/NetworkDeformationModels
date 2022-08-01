get_COG <- function(all_rois = F) {
  #get center of gravity for schafer300+tians2, minus excluded nodes
  cog <- read.table("data/COG.txt")
  if(all_rois==F) {
  removed_roi <- c(50,51,52, 53, 86, 90, 197,199, 200, 201, 242, 243, 244)
  cog <- cog[-c(removed_roi),]
  }
  
  message("get center of gravity for schafer300+tians2")
  return(cog)
}

get_node_dist  <- function(binned = F, percentile = F, all_rois = F) {
  #get distance for schaefer300+tians2, minus excluded nodes
  #Binned for QC checks
  library(bio3d)
  source("~/Dropbox/Sid/R_files/functions/vec_2_mat.R")
  cog <- read.table("data/COG.txt")
  
  if(all_rois==F) {
  removed_roi <- c(50,51,52, 53, 86, 90, 197,199, 200, 201, 242, 243, 244)
  cog <- cog[-c(removed_roi),]
  message("get distance matrix for schafer300+tians2, minus ecluded nodes: 50,51,52, 53, 86, 90, 197,199, 200, 201, 242, 243, 244")
  }
  xyz.dist <- dist.xyz(cog)  #calculate euclidean distance between each pair of centroids
  if(binned==T) {
    xyz.dist[xyz.dist<=50] <- 1
    xyz.dist[xyz.dist<=100 & xyz.dist>1] <- 2
    xyz.dist[xyz.dist>100] <- 3
    diag(xyz.dist) <- 0
    message("1=<50, 2=<100, 3=>100")
  }
  if(percentile==T) {
   dist_vec <- xyz.dist[upper.tri(xyz.dist)]
   dist_vec_pct <- fmsb::percentile(dist_vec)
   xyz.dist <- vec_2_mat(dist_vec_pct, 319, 0)
  }
  return(xyz.dist)
}

get_hemi_id <- function(all_rois = F) {
  #get id for hemisphere, 1 for l 2 for r
  if(all_rois==T) {
    hemi<-c(rep(1,150), rep(2,150), rep(2,16), rep(1,16))
  }
  if(all_rois==F) {
  labs <- readxl::read_excel("data/schaefer+aseg_labels.xlsx", sheet = "atlas")
  removed_roi <- c(50,51,52, 53, 86, 90, 197,199, 200, 201, 242, 243, 244)
  labs <- labs[-c(removed_roi),1]
  hemi<-c(rep(1,144), rep(2,143), rep(2,16), rep(1,16))
  }
  #View(cbind(labs, hemi))
  return(hemi)
  
}

get_labs <- function(remove_fov_rois = TRUE) {
  labs <- readxl::read_excel("~/Dropbox/Sid/R_files/STAGES_difussion/data/schaefer+aseg_labels.xlsx", sheet = "atlas")
  if(remove_fov_rois == TRUE) {
  removed_roi <- c(50,51,52, 53, 86, 90, 197,199, 200, 201, 242, 243, 244)
  labs <- as.data.frame(labs[-c(removed_roi),1])
  } 
  
  else if (remove_fov_rois == FALSE)  {
    labs <- as.data.frame(labs[,1])
  }
  return(labs)
}

get_network_ids <- function(remove_fov_rois = TRUE){
  labs <- readxl::read_excel("/data/schaefer+aseg_labels.xlsx", sheet = "atlas")
  if(remove_fov_rois == TRUE){
  removed_roi <- c(50,51,52, 53, 86, 90, 197,199, 200, 201, 242, 243, 244)
  labs <- as.data.frame(labs[-c(removed_roi),5])
  }
  return(labs)
}









