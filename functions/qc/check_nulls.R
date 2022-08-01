setwd("~/Dropbox/Sid/R_files/STAGES_difussion/")
source("scripts/functions/get_atrophy_t_vals_atlas.R")
source("scripts/functions/get_atrophy_t_vals_voxel.R")
source("scripts/functions/get_mean_hc_connectome.R")
source("scripts/functions/get_mean_hc_func_connectome.R")
source("scripts/functions/get_node_distance.R")
source("scripts/functions/get_nbs_connectome.R")
source("scripts/run_LDM.R")
source("scripts/functions/get_connectome_nulls.R")


bins=30
swaps=25000
nulls <- get_connectome_nulls(type = paste0("bins_", bins, "_swap_", swaps))
nullmat <- as.matrix(nulls[[99]])
convec <- c(nullmat[upper.tri(nullmat)])
convec <- convec[convec != 0] 

#is degree dist the same
hist(log10(convec), breaks = 20)

#is edge distance the same
dist <- get_node_dist()
nullmatbin <-  nullmat
nullmatbin[nullmatbin!=0] <- 1
nullmatbin <- nullmatbin*dist
nullmatbin <- c(nullmatbin[upper.tri(nullmatbin)])
nullmatbin <- nullmatbin[nullmatbin != 0] 
hist(nullmatbin, breaks = 20)

#degree distance relationship 
lim <- quantile(convec, .85)
nullmat[nullmat<lim] <- 0
x_igraph <- igraph::graph_from_adjacency_matrix(as.matrix(nullmat)) #convert connectivity matrix into an graph object.
d <- igraph::degree(x_igraph)

atlas <- as.data.frame(xlsx::read.xlsx("data/schaefer+aseg_labels.xlsx", sheetName = "atlas"))
removed_roi <- c(50,51,52, 53, 86, 90, 197,199, 200, 201, 242, 243, 244)
atlas <- atlas[-c(removed_roi),]
atlas[,2:4] <- sapply(atlas[,2:4], as.integer)


brainconn(atlas, conmat = nullmat,
          node.color = "purple", 
          view = "top", 
          node.size = log(d)/1.3, 
          #  edge.width = 0.5,
          edge.alpha = 0.9, 
          background.alpha = 0.3,
          all.nodes = T,
          scale.edge.width = c(0,10)
) 


brainconn(atlas, conmat = hc_connectome,
          node.color = "purple", 
          view = "ortho", 
          node.size = log(d)/1.3, 
          #  edge.width = 0.5,
          edge.alpha = 0.9, 
          background.alpha = 0.3,
          all.nodes = T,
          scale.edge.width = c(0,10)
) 


#plot weight-length distribution 
# Observerd 
hc_connectome <- as.matrix(get_mean_hc_connectome(method = "traditionalConsistency", type = "0.001")) #preserveDist #traditionalConsistency
degree <- rowSums(hc_connectome)
hc_connectome_w <- as.matrix(get_mean_hc_connectome(method = "average", type = "001"))
hc_connectome <- hc_connectome*hc_connectome_w
weighteddegree <- rowSums(hc_connectome)
weight_vec <- hc_connectome[upper.tri(hc_connectome)]
weight_vec <- weight_vec[weight_vec != 0]

hc_connectome_dist <-  as.matrix(get_mean_hc_connectome(method = "traditionalConsistency", type = "0.001")) #preserveDist #traditionalConsistency
dist  <- as.matrix(get_node_dist())
hc_connectome_dist <- hc_connectome_dist*dist

dist_vec <- hc_connectome_dist[upper.tri(hc_connectome_dist)]
dist_vec <- dist_vec[dist_vec!=0]

plot(dist_vec,weight_vec)
temp <- as.data.frame(cbind(log10(weight_vec), dist_vec))
temp$bins <- cut(temp$dist_vec,breaks = 15)
ggplot(temp, aes(x = bins, y =V1)) +   
  geom_boxplot(alpha=0.7) + 
  ylab("log10(weight)") + 
  xlab("Distance bins") +
  ggtitle("Observerd log10(weight) and distance relationship") +  
  scale_x_discrete(labels=1:15)



#Nulls
bins=30
swaps=25000
nulls <- get_connectome_nulls(type = paste0("bins_", bins, "_swap_", swaps))
nullmat <- as.matrix(nulls[[69]])
weighteddegreeNull <- rowSums(nullmat)
weight_vec <- nullmat[upper.tri(nullmat)]
weight_vec <- weight_vec[weight_vec != 0]

nullmat[nullmat !=0] <- 1
degreeNull <- rowSums(nullmat)
dist  <- as.matrix(get_node_dist())
hc_connectome_dist <- nullmat*dist

dist_vec <- hc_connectome_dist[upper.tri(hc_connectome_dist)]
dist_vec <- dist_vec[dist_vec!=0]

plot(dist_vec,weight_vec)
temp <- as.data.frame(cbind(log10(weight_vec), dist_vec))
temp$bins <- cut(temp$dist_vec,breaks = 15)
ggplot(temp, aes(x = bins, y =V1)) +   
  geom_boxplot(alpha=0.7) + 
  ylab("log10(weight)") + 
  xlab("Distance bins") +
  ggtitle("Random Null log10(weight) and distance relationship") +  
  scale_x_discrete(labels=1:15)

hist(weighteddegree, breaks = 30, xlim = c(0,170), ylim = c(0,80))
hist(weighteddegreeNull, breaks = 30, xlim = c(0,170), ylim = c(0,80))

hist(degree, breaks = 30, xlim = c(0,250), ylim = c(0,30))
hist(degreeNull, breaks = 30, xlim = c(0,250), ylim = c(0,30))
