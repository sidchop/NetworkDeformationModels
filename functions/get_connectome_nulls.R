get_connectome_nulls <- function(type = "bins_10_swap_50000", sample = NULL) {
if(sample == "stages") {filepath = paste0("~/Dropbox/Sid/R_files/STAGES_difussion/data/nulls/STAGES_",type,"/")}
if(sample == "gencog") {filepath = paste0("~/Dropbox/Sid/R_files/STAGES_difussion/data/nulls/GOC_",type,"/")}
  
  
connectomes <- list.files(path = filepath, full.names = TRUE, recursive = F)  
  nulls <- lapply(connectomes, function(i){
    data.table::fread(i)
  })
  return(nulls)
}

#x <- get_connectome_nulls("preserve_weight_degree_length_weighted")
#pheatmap::pheatmap(x[[1]], cluster_rows = F, cluster_cols = F)

