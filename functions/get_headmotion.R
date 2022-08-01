#This displacement is calculated for each volume with respect to both the previous volume, 
#and to the first volume; from these, we take the mean across all volumes to obtain scalar
#estimates of relative (RELall) and absolute (ABSall) head motion respectively. For each 
#pipeline, the corresponding motion parameter estimates from EDDY1 or EDDY2 were utilized.
#In this article we focus on results obtained using ABSall; the results for RELall, along 
#with five other measures that have also been used to characterise motion (Baum et al., 2018; 
#                                                                          Roalf et al., 2016), 
#are presented in the supplementary material (Table S1). We focus on ABSall as a measure of 
#motion as it is an established way that quality control can be performed for motion in DWI 
#data (Bastiani et al., 2019), and this measure can capture both gradual and sudden movements.


 get_headmotion <- function(type = 'absolute') {
  #read subj list
  subjlist <- read.table("/data/squad_headmotion/subjlist.txt")
  squadjson <- jsonlite::read_json("data/squad_headmotion/group_db.json")
  
  if (type =='absolute') {
    headmotion <- unlist(lapply(squadjson$qc_motion, `[[`, 1))
  }
  if (type =='relative') {
    headmotion <- unlist(lapply(squadjson$qc_motion, `[[`, 2))
  }
data <- cbind(subjlist, headmotion)
data <- data[order(data$V1),]
return(data)
}


