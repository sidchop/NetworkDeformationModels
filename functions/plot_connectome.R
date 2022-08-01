#Takes in csv containing conectivity weights, and outputs .png heatmap.
#Requires `pheatmap`, `Matrix`, `RColorBrewer`
plot_connectome <- function(csv, 
                            binary = FALSE, 
                            thr = NULL, 
                            ulim = NULL, 
                            llim=NULL,
                            breaks=NULL,
                            labs=NULL,
                            title=csv,
                            out="connectome.png"){
  x <- read.table(csv, header = F)
  x <- as.matrix(Matrix::forceSymmetric(as.matrix(x)))
  if(!is.null(labs)){labs <- read.csv(labs, header = F)
  colnames(x) <- labs$V1
  rownames(x) <- labs$V1}
  if(binary==TRUE){x <- as.matrix((x > thr) + 0)}
  library(RColorBrewer)
  if(!is.null(ulim) & !is.null(llim)){breaksList = seq(llim,ulim, by = breaks)
  color = colorRampPalette(brewer.pal(n = 7, name = "Blues"))(length(breaksList)) # Defines the vector of colors for the legend (it has to be of the same length of breaksList)
  png(out, width = 4, height = 4, units = 'in', res = 200)
  pheatmap::pheatmap(as.matrix(x), treeheight_row = 0, treeheight_col = 0,
                     cluster_rows=F, cluster_cols=F, fontsize=5, breaks = breaksList,
                     color = color,  main = title)
  dev.off()
  }
  
  if(is.null(ulim) & is.null(llim)){
    png(out, width = 4, height = 4, units = 'in', res = 200)
    pheatmap::pheatmap(as.matrix(x), treeheight_row = 0, treeheight_col = 0,
                       cluster_rows=F, cluster_cols=F, fontsize=5, main = title,
                       color =c("#000000", "#FFFFFF"), border_color = NA)
    dev.off()
  }
  
}
