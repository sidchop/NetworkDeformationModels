get_connectome_nulls <- function(type = "bins_10_swap_50000", sample = NULL) {
if(sample == "stages") {filepath = paste0("~/data/nulls/STAGES_",type,"/")}
if(sample == "gencog") {filepath = paste0("~/data/nulls/GOC_",type,"/")}
  
connectomes <- list.files(path = filepath, full.names = TRUE, recursive = F)  
  nulls <- lapply(connectomes, function(i){
    data.table::fread(i)
  })
  return(nulls)
}


