
##Plot all results
all_results <- readRDS("output/results/STAGES_all_results_LDM.RDS")

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



#Plot Main Figure for illness-related results (Figure 2)
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

#annotate r and p values
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




#### Plot main medication related results (Figure 3)
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






