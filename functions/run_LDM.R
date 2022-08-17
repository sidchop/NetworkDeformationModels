#Function to compute Linked deformation model AKA Network Deformation Models (NDM)
#See scripts/runNDM.R for use of this workhorse function

run_LDM <- function(sc = NULL, 
                    fc = NULL, 
                    atrophy=NULL, 
                    weight_by_fc = FALSE, 
                    weight_by_sc = FALSE, 
                    null.type = NULL, 
                    cor.type = "spearman",
                    scale.nulls = F,
                    sample = NULL) {
  
  neighbour_atrophy<-  r_obs <-  p <-  r_null <- NULL
  neighbour_atrophy <- vector()
  
  if(weight_by_fc == F & weight_by_sc == F) {
    for (i in 1:length(atrophy)) {
      neighbours <- which(sc[,i]!=0)
      if(length(neighbours) == 0 ){
        neighbour_atrophy[i]  <- NA
      } else {
        neighbour_atrophy[i] <- sum(atrophy[neighbours])*(1/length(atrophy[neighbours]))
      }
    }
  }
  
  if(weight_by_fc == T | weight_by_sc == T) {
    for (i in 1:length(atrophy)){
      neighbours <- which(sc[,i]!=0)
      if(weight_by_fc == T){weight_vec <- fc[i,c(neighbours)]}
      if(weight_by_sc == T){weight_vec <- sc[i,c(neighbours)]}
      if(length(neighbours) == 0 ){
        neighbour_atrophy[i]  <- NA
      } else {
        neighbour_atrophy[i] <- sum(atrophy[neighbours]*weight_vec)*(1/length(atrophy[neighbours]))
      }
    }
  }
  r_obs <-  c(cor(neighbour_atrophy, atrophy, method = cor.type))
  if(!is.null(null.type)){
    if(null.type != 'spin'){
      source("~/functions/get_connectome_nulls.R")
      nulls <- get_connectome_nulls(type = null.type, sample = sample)
      if (scale.nulls == T) {
        list.scale <- function(x) {
          x <- as.data.frame(x)
          x[x != 0] <- scale(x[x != 0])
          return(x)
        }
        nulls <-  lapply(nulls, list.scale)
      }
      neighbour_atrophy_null <- rnull <-  vector()
      
      for (n in 1:length(nulls)){
        null <- as.matrix(nulls[[n]])
        if(weight_by_fc == F & weight_by_sc == F) {
          for (i in 1:length(atrophy)) {
            neighbours <- which(null[,i]!=0)
            neighbour_atrophy_null[i] <- sum(atrophy[neighbours])*(1/length(atrophy[neighbours]))
          }
        }
        
        if(weight_by_fc == T | weight_by_sc == T) {
          for (i in 1:length(atrophy)){
            neighbours <- which(null[,i]!=0)
            if(weight_by_fc == T){weight_vec <- fc[i,c(neighbours)]}
            if(weight_by_sc == T){weight_vec <- null[i,c(neighbours)]}
            neighbour_atrophy_null[i] <- sum(atrophy[neighbours]*weight_vec)*(1/length(atrophy[neighbours]))
          }
        }
        r_null[n] <- cor(neighbour_atrophy_null, atrophy, method =  cor.type)
        message(n)
      }
      p <- sum(abs(r_null)>=abs(r_obs))/length(r_null) # two tailed
    }
  } 
  return(list(atrophy, neighbour_atrophy, r_obs, r_null, p))
}








#######################  #######################  #######################  #######################  #######################  
#######################  #######################  #######################  #######################  #######################  
#######################  #######################  #######################  #######################  #######################  
#######################  #######################  #######################  #######################  #######################  
#######################  #######################  #######################  #######################  #######################  
#######################  #######################  #######################  #######################  #######################  

#does neightbour atrophy predict atrophy at t2

run_LDM2 <- function(sc = NULL, fc = NULL, neighbour_atrophy_vec=NULL, actual_atrophy_vec=NULL, weight_by_fc = FALSE, weight_by_sc = FALSE, null.type = NULL, cor.type = "spearman") {
  
  r_obs <-  p <-  r_null <- NULL
  neighbour_atrophy <- vector()
  
  if(weight_by_fc == F & weight_by_sc == F) {
    for (i in 1:length(neighbour_atrophy_vec)) {
      neighbours <- which(sc[,i]!=0)
      neighbour_atrophy[i] <- sum(neighbour_atrophy_vec[neighbours])*(1/length(neighbour_atrophy_vec[neighbours]))
    }
  }
  
  if(weight_by_fc == T | weight_by_sc == T) {
    for (i in 1:length(neighbour_atrophy_vec)){
      neighbours <- which(sc[,i]!=0)
      if(weight_by_fc == T){weight_vec <- fc[i,c(neighbours)]}
      if(weight_by_sc == T){weight_vec <- sc[i,c(neighbours)]}
      neighbour_atrophy[i] <- sum(neighbour_atrophy_vec[neighbours]*weight_vec)*(1/length(neighbour_atrophy_vec[neighbours]))
    }
  }
  
  r_obs <-  cor(neighbour_atrophy, actual_atrophy_vec, method = cor.type )
  
  if(!is.null(null.type)){
    
    source("~/Dropbox/Sid/R_files/STAGES_difussion/scripts/functions/get_connectome_nulls.R")
    nulls <- get_connectome_nulls(type = null.type)
    neighbour_atrophy_null <- rnull <-  vector()
    for (n in 1:length(nulls)){
      null <- as.matrix(nulls[[n]])
      if(weight_by_fc == F & weight_by_sc == F) {
        for (i in 1:length(neighbour_atrophy_vec)) {
          neighbours <- which(null[,i]!=0)
          neighbour_atrophy_null[i] <- sum(neighbour_atrophy_vec[neighbours])*(1/length(neighbour_atrophy_vec[neighbours]))
        }
      }
      
      if(weight_by_fc == T | weight_by_sc == T) {
        for (i in 1:length(neighbour_atrophy_vec)){
          neighbours <- which(null[,i]!=0)
          if(weight_by_fc == T){weight_vec <- fc[i,c(neighbours)]}
          if(weight_by_sc == T){weight_vec <- null[i,c(neighbours)]}
          neighbour_atrophy_null[i] <- sum(neighbour_atrophy_vec[neighbours]*weight_vec)*(1/length(neighbour_atrophy_vec[neighbours]))
        }
      }
      r_null[n] <- cor(neighbour_atrophy_null, actual_atrophy_vec, method =  cor.type)
      message(n)
    }
   # p <- 1-ecdf(r_null)(r_obs)
    p <- NULL
  }
  return(list(neighbour_atrophy, r_obs, p, r_null))
}


#as a sensitivity test, see if only summing the edge weights (taking out the influence of neighbour atrophy is predictive)
run_LDM_edgeweightOnly <- function(sc = NULL, fc = NULL, atrophy=NULL, weight_by_fc = FALSE, weight_by_sc = FALSE, null.type = NULL, cor.type = "spearman") {
  
  neighbour_atrophy<-  r_obs <-  p <-  r_null <- NULL
  neighbour_atrophy <- vector()
  
  if(weight_by_fc == F & weight_by_sc == F) {
    for (i in 1:length(atrophy)) {
      neighbours <- which(sc[,i]!=0)
      neighbour_atrophy[i] <- length(atrophy[neighbours])
    }
  }
  
  if(weight_by_fc == T | weight_by_sc == T) {
    for (i in 1:length(atrophy)){
      neighbours <- which(sc[,i]!=0)
      if(weight_by_fc == T){weight_vec <- fc[i,c(neighbours)]}
      if(weight_by_sc == T){weight_vec <- sc[i,c(neighbours)]}
      neighbour_atrophy[i] <- sum(weight_vec)*(1/length(atrophy[neighbours]))
    }
  }
  r_obs <-  cor(neighbour_atrophy, atrophy, method = cor.type )
  if(!is.null(null.type)){
    source("~/Dropbox/Sid/R_files/STAGES_difussion/scripts/functions/get_connectome_nulls.R")
    nulls <- get_connectome_nulls(type = null.type)
    neighbour_atrophy_null <- rnull <-  vector()
    
    for (n in 1:length(nulls)){
      null <- as.matrix(nulls[[n]])
      if(weight_by_fc == F & weight_by_sc == F) {
        for (i in 1:length(atrophy)) {
          neighbours <- which(null[,i]!=0)
          neighbour_atrophy_null[i] <- length(atrophy[neighbours])
        }
      }
      
      if(weight_by_fc == T | weight_by_sc == T) {
        for (i in 1:length(atrophy)){
          neighbours <- which(null[,i]!=0)
          if(weight_by_fc == T){weight_vec <- fc[i,c(neighbours)]}
          if(weight_by_sc == T){weight_vec <- null[i,c(neighbours)]}
          neighbour_atrophy_null[i] <- sum(weight_vec)*(1/length(atrophy[neighbours]))
        }
      }
      r_null[n] <- cor(neighbour_atrophy_null, atrophy, method =  cor.type)
      message(n)
    }
    p <- 1-ecdf(r_null)(r_obs)
  }
  
  
  return(list(neighbour_atrophy, r_obs, p, r_null))
  
  
}



#how does the correlation change when you only include short, med, and long range connections 
run_LDM_dist_binned <- function(sc = NULL, 
                                dist_binned= NULL,
                                bin=NULL,
                                fc = NULL, 
                                atrophy=NULL, 
                                weight_by_fc = FALSE, 
                                weight_by_sc = FALSE, 
                                cor.type = "spearman") {
  
  neighbour_atrophy<-  r_obs <-  p <-  r_null <- NULL
  neighbour_atrophy <- vector()
  
  if(weight_by_fc == F & weight_by_sc == F) {
    for (i in 1:length(atrophy)) {
      neighbours <- which(sc[,i]>0 & dist_binned[,i]==bin)
      if(length(neighbours) == 0 ){
        neighbour_atrophy[i]  <- NA
      } else {
        neighbour_atrophy[i] <- sum(atrophy[neighbours])*(1/length(atrophy[neighbours]))
      }
    }
  }
  
  
  if(weight_by_fc == T | weight_by_sc == T) {
    for (i in 1:length(atrophy)){
      neighbours <- which(sc[,i]>0 & dist_binned[,i]==bin)
      if(weight_by_fc == T){weight_vec <- fc[i,c(neighbours)]}
      if(weight_by_sc == T){weight_vec <- sc[i,c(neighbours)]}
      if(length(neighbours) == 0 ){
        neighbour_atrophy[i]  <- NA
      } else {
        neighbour_atrophy[i] <- sum(atrophy[neighbours]*weight_vec)*(1/length(atrophy[neighbours]))
      }
    }
  }
  r_obs <-  cor(neighbour_atrophy, atrophy, method = cor.type, use = "pairwise.complete.obs")
  return(list(neighbour_atrophy, r_obs))
}



