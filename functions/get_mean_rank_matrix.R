get_mean_rank_matrix <- function(type=NULL, rank = TRUE, sample = NULL, networks = F){
  
  if(sample == "stages") {all_results <- readRDS("/results/STAGES_all_results_LDM.RDS")}
  if(sample == "gencog") {all_results <- readRDS("/results/GOC_all_results_LDM.RDS")}
  
  message("higher rank = greater atrophy in patients")
  if (type == "illness") {
    if(sample == "stages") {mean_rank_matrix <- matrix(nrow = 319, ncol = 3)}
    if(sample == "gencog") {mean_rank_matrix <- matrix(nrow = 332, ncol = 3)}
    for (i in 1:3) {
      if(rank == TRUE) {
      atrophy_rank <- rank(all_results[[3]][[i]][["ROI atrophy"]])
      neigbour_rank <- rank(all_results[[3]][[i]][["Neighbour atrophy"]])
      mean_rank_matrix[,i] <- (atrophy_rank+neigbour_rank)/2
      } else {
        atrophy_rank <- all_results[[3]][[i]][["ROI atrophy"]]
        neigbour_rank <- all_results[[3]][[i]][["Neighbour atrophy"]]
        mean_rank_matrix[,i] <- (atrophy_rank+neigbour_rank)/2
      }
    }
  }
  
  if (type == "medication") {
    if(sample == "stages") {mean_rank_matrix <- matrix(nrow = 319, ncol = 2)}
    if(sample == "gencog") {mean_rank_matrix <- matrix(nrow = 332, ncol = 2)}
    for (i in 1:2) {
      if(rank == TRUE) {
      atrophy_rank <- rank(all_results[[3]][[3+i]][["ROI atrophy"]])
      neigbour_rank <- rank(all_results[[3]][[3+i]][["Neighbour atrophy"]])
      mean_rank_matrix[,i] <- (atrophy_rank+neigbour_rank)/2
      } else {
        atrophy_rank <- all_results[[3]][[3+i]][["ROI atrophy"]]
        neigbour_rank <- all_results[[3]][[3+i]][["Neighbour atrophy"]]
        mean_rank_matrix[,i] <- (atrophy_rank+neigbour_rank)/2
      }
    }
  }
  if(networks==T){
    source("scripts/functions/get_node_distance.R")
    netid <- get_network_ids(remove_fov_rois = F)[,5]
    netid <- as.factor(netid$netid)
    mean_rank_matrix <- aggregate(mean_rank_matrix, list(netid), FUN=mean) 
    mean_rank_matrix <- mean_rank_matrix[,-1]
  }
  return(mean_rank_matrix)
}

