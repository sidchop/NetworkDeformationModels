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

#No_weighting_pearson <- apply_LDM_all(weight_by_sc = F, weight_by_fc = F, cor.type = "pearson")
#FC_weighting_pearson  <- apply_LDM_all(weight_by_sc = F, weight_by_fc = T, cor.type = "pearson")
#SC_weighting_pearson  <- apply_LDM_all(weight_by_sc = T, weight_by_fc = F, cor.type = "pearson")

all_results <- list(No_weighting_spearman, FC_weighting_spearman ,SC_weighting_spearman)
#    No_weighting_pearson,FC_weighting_pearson,SC_weighting_pearson)

#saveRDS(all_results, "output/results/all_results_LDM_scaled.RDS")




##Plot all results
setwd("~/Dropbox/Sid/R_files/STAGES_difussion/")
all_results <- readRDS("output/results/STAGES_all_results_LDM_noScaing.RDS")



#print neat results for table
library(weights)
for (m in 1:3) {
  for (c in 1:5) {
    r = unlist(all_results[[m]][[c]]["correlation"])
    p = unlist(all_results[[m]][[c]]["connectome p"])
    print(paste0(rd(r,3)," (", rd(p,3),")"))
  }
}

for (m in 1:3) {
  for (c in 1:5) {
    r = unlist(all_results[[m]][[c]]["correlation nosub"])
    p = unlist(all_results[[m]][[c]]["spin p"])
    print(paste0(rd(r,3)," (", rd(p,3),")"))
  }
}

for (m in 1:3) {
  for (c in 1:5) {
    r = unlist(all_results[[m]][[c]]["correlation"])
    p = unlist(all_results[[m]][[c]]["param p"])
    print(paste0(rd(r,3)," (", rd(p,3),")"))
  }
}

#plot nulls and results
# Merge all nulls for illness 
library(reshape2)
library(ggplot2)

null_list <- list()
null_list_all <- list()
index2 <- 1
for (t in c("Baseline","Δ3 Months",  "Δ12 Months")) {
  index <- 1
  for (w in c("None", "FC", "SC")){
    nulls <- qpcR:::cbind.na(all_results[[index]][[index2]][["connectome nulls"]], 
                             all_results[[index]][[index2]][["param nulls"]])
    nulls <- qpcR:::cbind.na(nulls, all_results[[index]][[index2]][["spin nulls"]])
    colnames(nulls) <- c("Cn", "Pn", "Sn") 
    nulls <- melt(nulls, na.rm = F)
    nulls$Var3 <- rep(paste(w,t), length(nulls$Var1))
    null_list[[index]] <- nulls
    index <- index + 1
  }
  nulls <- do.call("rbind", null_list)
  nulls$Var4 <- rep(t, length(nulls$Var1))
  null_list_all[[index2]] <- nulls
  index2 <- index2 + 1
}
nulls <- do.call("rbind", null_list_all)

nulls$Var3 <- as.factor(nulls$Var3)
nulls$Var3  <- factor(nulls$Var3,
                      levels = c("None Baseline","SC Baseline",  "FC Baseline",
                                 "None Δ3 Months","SC Δ3 Months",  "FC Δ3 Months",
                                 "None Δ12 Months","SC Δ12 Months",  "FC Δ12 Months"),
                      ordered = TRUE)


          
library(gghalves)

theme = theme_set(theme_minimal())

theme = theme_update(legend.position="right", 
                     legend.title = element_text(color = "black", size = 14), 
                     panel.grid.major.x=element_blank(),
                     panel.grid.minor.x = element_blank(),
                     axis.title.x = element_blank(),
                     axis.text.x = element_blank(),
                     axis.text.y =element_text(size=12),
                     axis.title.y =element_text(size=16),
                     legend.text=element_text(size=12))


p <- ggplot(nulls, aes(x=Var3, y=value, fill = Var2)) +
  geom_boxplot(aes(fill=Var2), outlier.shape = NA, width=0.2,
               alpha= 0) +
  geom_point(aes(colour=Var2, group=Var2), 
             position = position_jitterdodge(jitter.width = 0.2,dodge.width = 0.2), 
             size = .1, alpha = 0.01) + 
  geom_half_violin(position = position_nudge(x = -.2, y =0), 
                   alpha = .5, trim = T, adjust = 2)   + 
  ylab("Null distribution") + ylim(c(-0.3,0.65)) +
  annotate(geom = "text", x = 1:9, y = -0.3, 
           label = rep(c("SC (No Weighting)", "SC Weighted","FC Weighted"),3), size = 4) + 
  annotate(geom = "text", x = 2 + 3 * (0:2), y = -0.32, 
           label = c("Baseline", "Δ3 Months", "Δ12 Months"), size = 5) +
  scale_fill_brewer(palette = "Dark2",
                     labels= c("Connectome nulls", "Parameterised nulls", "Spin nulls"),
                     name= "Null type") +
  scale_colour_brewer(palette = "Dark2",
                       labels= c("Connectome nulls", "Parameterised nulls", "Spin nulls"),
                       name= "Null type") 


p
#annotate r and p
source("scripts/functions/extract_r_and_p.R")
pv <-  extract_r_and_p(value = "p", all_results)
r <- extract_r_and_p(value = "r", all_results)
pv[pv==0] <- 2  # 2 = <0.001, 1 = < 0.01
pv[pv<0.001] <- 1
ind <- 0
for (i in 1:9) {
  ind <- ind + 1
  if(pv[ind]==2){
    p <- p +  annotate("point", x = i-0.06, y = r[ind],  size = 2,
                       alpha = 1, 
                       colour = "red") 
  }
  ind <- ind + 1
  if(pv[ind]==2){
    p <- p +  annotate("point", x = i, y = r[ind],  size = 2,
                       alpha = 1, 
                       colour = "red") 
  }
  ind <- ind + 1
  if(pv[ind]==2){
    p <- p +  annotate("point", x = i+0.06, y = r[ind],  size = 2,
                       alpha = 1, 
                       colour = "red") 
  }
}

p


ind <- 0
for (i in 1:9) {
  ind <- ind + 1
  p <- p +  annotate("point", x = i-0.06, y = r[ind],  size = 1,
                     alpha = 0.8, 
                     colour = "black") 
  ind <- ind + 1
  p <- p +  annotate("point", x = i, y = r[ind],  size = 1,
                     alpha = 1, 
                     colour = "black") 
  ind <- ind + 1
  p <- p +  annotate("point", x = i+0.06, y = r[ind],  size = 1,
                     alpha = 0.8, 
                     colour = "black") 
}
p


# scatter plots
library(ggplot2)
plots <- list()
for (i in 1:3){
  data <- as.data.frame(cbind(all_results[[3]][[i]][["ROI atrophy"]],
                              all_results[[3]][[i]][["Neighbour atrophy"]]))
  ggplot(data, aes(x = V2, y = V1)) + 
    geom_point(colour="#13697E", size = 5) + 
    geom_point(shape = 1,size = 5,colour = "black", stroke = 0.3, alpha = 0.5) +
    geom_smooth(method = "lm", colour = "black")  + 
    xlab("SC-weighted mean neighbour deformatiom") +
    ylab("ROI deformaion") + 
    theme_classic() +
    theme(axis.text=element_text(size=16),axis.title=element_text(size=16,face="bold"))
  
}


# add line for oberved values instead of dots
#p + annotate("segment", x = 1-0.2, xend =1+0.2, y = r[1], yend = r[1], 
#             size = 0.4,
#             alpha = 0.8, 
#             colour = "red", linetype=2)


#p <- ggplot(data = nulls, aes(x=Var3, y=value, fill = Var2)) + 
#  geom_boxplot(aes(fill=Var2), 
#               outlier.shape = NA, width=0.5) + 
#  ylab("Null distribution") + ylim(c(-0.3,0.7)) +
#  annotate(geom = "text", x = 1:9, y = -0.27, 
#           label = rep(c("-", "FC","SC"),3), size = 4) + 
#  annotate(geom = "text", x = 2 + 3 * (0:2), y = -0.30, 
#           label = c("Baseline", "Δ3 Months", "Δ12 Months"), size = 5) +
#  scale_fill_viridis(discrete = T, option = "plasma") 
#p
#
#
#
##obs_r and p_vals
#source("scripts/functions/extract_r_and_p.R")
#r <- extract_r_and_p(value = "r")
## add observerd values
#ind <- 0
#for (i in 1:9) {
#  ind <- ind + 1
#  p <- p +  annotate("point", x = i-0.16, y = r[ind],  size = 2,
#                     alpha = 1, 
#                     colour = "black") 
#  ind <- ind + 1
#  p <- p +  annotate("point", x = i, y = r[ind],  size = 2,
#                     alpha = 1, 
#                     colour = "black") 
#  ind <- ind + 1
#  p <- p +  annotate("point", x = i+0.16, y = r[ind],  size = 2,
#                     alpha = 1, 
#                     colour = "black") 
#}



#### Plot medication results
null_list <- list()
null_list_all <- list()
index2 <- 4
for (t in c("Δ3 Months",  "Δ12 Months")) {
  index <- 1
  for (w in c("None", "FC", "SC")){
    nulls <- qpcR:::cbind.na(all_results[[index]][[index2]][["connectome nulls"]], 
                             all_results[[index]][[index2]][["param nulls"]])
    nulls <- qpcR:::cbind.na(nulls, all_results[[index]][[index2]][["spin nulls"]])
    colnames(nulls) <- c("Cn", "Pn", "Sn") 
    nulls <- melt(nulls, na.rm = F)
    nulls$Var3 <- rep(paste(w,t), length(nulls$Var1))
    null_list[[index]] <- nulls
    index <- index + 1
  }
  nulls <- do.call("rbind", null_list)
  nulls$Var4 <- rep(t, length(nulls$Var1))
  null_list_all[[index2]] <- nulls
  index2 <- index2 + 1
}

nulls <- do.call("rbind", null_list_all)

nulls$Var3 <- as.factor(nulls$Var3)
nulls$Var3  <- factor(nulls$Var3,
                      levels = c("None Δ3 Months","SC Δ3 Months",  "FC Δ3 Months",
                                 "None Δ12 Months","SC Δ12 Months",  "FC Δ12 Months"),
                      ordered = TRUE)



library(gghalves)

theme = theme_set(theme_minimal())

theme = theme_update(legend.position="right", 
                     legend.title = element_text(color = "black", size = 14), 
                     panel.grid.major.x=element_blank(),
                     panel.grid.minor.x = element_blank(),
                     axis.title.x = element_blank(),
                     axis.text.x = element_blank(),
                     axis.text.y =element_text(size=12),
                     axis.title.y =element_text(size=16),
                     legend.text=element_text(size=12))


p <- ggplot(nulls, aes(x=Var3, y=value, fill = Var2)) +
  geom_boxplot(aes(fill=Var2), outlier.shape = NA, width=0.2,
               alpha= 0) +
  geom_point(aes(colour=Var2, group=Var2), 
             position = position_jitterdodge(jitter.width = 0.2,dodge.width = 0.2), 
             size = .1, alpha = 0.01) + 
  geom_half_violin(position = position_nudge(x = -.2, y =0), 
                   alpha = .5, trim = T, adjust = 2)   + 
  ylab("Null distribution") + ylim(c(-0.3,0.65)) +
  annotate(geom = "text", x = 1:6, y = -0.3, 
           label = rep(c("SC (No Weighting)", "SC Weighted","FC Weighted"),2), size = 4) + 
 # annotate(geom = "text", x = 2 + 3* (0:2), y = -0.32, 
  #         label = c("Δ3 Months", "Δ12 Months"), size = 5) +
  scale_fill_brewer(palette = "Dark2",
                    labels= c("Connectome nulls", "Parameterised nulls", "Spin nulls"),
                    name= "Null type") +
  scale_colour_brewer(palette = "Dark2",
                      labels= c("Connectome nulls", "Parameterised nulls", "Spin nulls"),
                      name= "Null type") 


p
#annotate r and p
source("scripts/functions/extract_r_and_p.R")
pv <-  extract_r_and_p(value = "p", all_results, medication = T)
r <- extract_r_and_p(value = "r", all_results, medication = T)
pv[pv==0] <- 2  # 2 = <0.001, 1 = < 0.01
pv[pv<0.001] <- 1
ind <- 0
for (i in 1:6) {
  ind <- ind + 1
  if(pv[ind]==2){
    p <- p +  annotate("point", x = i-0.06, y = r[ind],  size = 2,
                       alpha = 1, 
                       colour = "red") 
  }
  ind <- ind + 1
  if(pv[ind]==2){
    p <- p +  annotate("point", x = i, y = r[ind],  size = 2,
                       alpha = 1, 
                       colour = "red") 
  }
  ind <- ind + 1
  if(pv[ind]==2){
    p <- p +  annotate("point", x = i+0.06, y = r[ind],  size = 2,
                       alpha = 1, 
                       colour = "red") 
  }
}

p


ind <- 0
for (i in 1:6) {
  ind <- ind + 1
  p <- p +  annotate("point", x = i-0.06, y = r[ind],  size = 1,
                     alpha = 0.8, 
                     colour = "black") 
  ind <- ind + 1
  p <- p +  annotate("point", x = i, y = r[ind],  size = 1,
                     alpha = 1, 
                     colour = "black") 
  ind <- ind + 1
  p <- p +  annotate("point", x = i+0.06, y = r[ind],  size = 1,
                     alpha = 0.8, 
                     colour = "black") 
}
p


#plot medication scatter plots
library(ggplot2)
plots <- list()
for (i in 4:5){
  data <- as.data.frame(cbind(all_results[[3]][[i]][["ROI atrophy"]],
                              all_results[[3]][[i]][["Neighbour atrophy"]]))
  ggplot(data, aes(x = V2, y = V1)) + 
    geom_point(colour="#13697E", size = 5) + 
    geom_point(shape = 1,size = 5,colour = "black", stroke = 0.3, alpha = 0.5) +
    geom_smooth(method = "lm", colour = "black")  + 
    xlab("SC-weighted mean neighbour deformatiom") +
    ylab("ROI deformaion") + 
    theme_classic() +
    theme(axis.text=element_text(size=16),axis.title=element_text(size=16,face="bold"))
  
}













##does neighbour atrophy at T1, predict atrophy at T2>

atrophy <- get_long_t_scores_dbm(smoothing = 3, fdr = F,
                                 timepoint = "bl", stat = "T",
                                 contrast = "", p = 1)


atrophy2 <- get_long_t_scores_dbm(smoothing = 3, fdr = F,
                                  timepoint = "bl_3m", stat = "T",
                                  contrast = "pipt_v_hc", p = 1)

rho = cor(atrophy, atrophy2, method = "spearman")
plot(atrophy, atrophy2)
results <- run_LDM2(sc = as.matrix(hc_connectome), 
                    fc = hc_funconnectome, 
                    neighbour_atrophy_vec=atrophy,
                    actual_atrophy_vec=atrophy2,
                    weight_by_sc = T,
                    weight_by_fc = F,
                    cor.type = "spearman",
                    null.type = NULL)

rho = cor(results[[1]], atrophy2, method = "spearman")
plot(x = results[[1]], y = atrophy2, xlab = "baseline neighbour atrophy", 
     ylab = "longitudinal atrophy", main = paste0("rho = ", round(rho,3)))

abline(lm(atrophy2~results[[1]]))


#### What if we put SC and atrophy on same scale

atrophy <- get_long_t_scores_dbm(smoothing = 3, fdr = F,
                                 timepoint = "bl_3m", stat = "T",
                                 contrast = "pipt_v_hc", p = 1, zscore = T)

#scale/zscore mat
x <- as.matrix(hc_connectome)
x1 <- x!=0
x[x1] <- scale(x[x1])
results <- run_LDM(sc = x, 
                   fc = hc_funconnectome, 
                   atrophy = atrophy,
                   weight_by_sc = T,
                   weight_by_fc = F,
                   cor.type = "spearman",
                   null.type = NULL)

plot(results[[2]], atrophy)


###########################################
# LDM with paramiterised nulls
###########################################
