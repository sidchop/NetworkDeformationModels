# This script visualizes different stages of the processing of dwi images, by outputting .png and .gif 
# files for each scan at each stage of processing. Visualizing each scan at each processing step using 
# gifs or representative slices (.png) can make quality checks for large datasets fast and efficient. 

# Currently the code is designed to be executed locally, downloading each scan from massive, before 
# creating the .png/.gif output. This is slow and inefficient, and the reason for doing this is that
# I had difficulty getting the `magick` package working smoothly on massive. However, for a dataset as 
# large as the biobank, running a script like this on massive may be your only option. Try the code
# out on massive, maybe the latest version of `magick` will work fine, or the massive team will be able 
# to assist. Also, you can probably just remove the dependence on 'magick' and just use the display 
# functions offered by `neurobase`, (e.g. slice, image, ortho2, overlay). See https://github.com/sidchop/neuRo-vis-BrainHack-2021/blob/main/nifti.Rmd 
# for a intro to visualizing voxel-level data in R. 

# As it stands, this code visualizes for following voxel-level images (however, the code and framework
# can be applied to any processing step):

# Raw T1
# Mean B0
# Raw diffusion 
# Denoised diffusion 
# Distortion corrected diffusion (i.e. post top-up)
# Eddy corrected difussion 
# T1 - tractoflow/mrtix processing (BET + n4 + crop + normalise + resample)
# T1 - T1 -> DWI ANTS registration
# DTI metrics  (FA)
# fODF metrics, first 10 directions

# Not included in this script but would be useful for QC:

# - gm/wm/csf/sc Tissue maps used for anatomically constrained tractography (see check_civet_segmentation.R)
# - Visualizing the shape of the diffusion ellipsoids, as opposed to just the voxel level fODF metrics
# - Visualize the tractogram. I have not found a good way to view .trk/.tck files in R, so this is best 
# done in python. See DIPY (https://dipy.org/documentation/1.4.1./examples_index/#visualization) OR 
# https://github.com/GuillaumeTh/dmriqcpy (I would recommend the scripts in dmriqcpy, as they are
# specifically designed for QC)
# - Registration of any atlas into dwi space (see make_atlas_overlay.R)

# The output from qc scripts like this is best complied and visualized using R-markdown. See qc_report.rmd
# for and example .rmd file and qc_report.html for example output. Note that these examples probably overuse
# gifs, when you preferred method for visualizing may be multiple representation slices, rather than gifs.
# Also, in the example report, you will see differences in the intensity ranges of some images. This can be fixed by
# using the `robust_window()` function, as implement in the script below. 

# Please contact me at sid.chopra@monash.edu if you have any questions. 

library(oro.nifti)
library(neurobase)
library(magick)
library(magrittr)
library(ssh)
library(magick)
library(RNifti)

#list of subject IDs (using BIDS format)
sublist <- read.table("~/Dropbox/Sid/R_files/STAGES_difussion/scripts/list_fmri_ses-1.txt", col.names = F)
setwd("~/Dropbox/Sid/R_files/STAGES_difussion/temp")

# Where the .png/.gif outputs will be stored
outdir="~/Dropbox/Sid/R_files/STAGES_difussion/output/qc_report/"

for (s in 1:dim(sublist)[1]){
  sub=sublist[s,]
  #======================================
  # 0. T1
  #=====================================
  session <- ssh_connect("scho0011@m3.massive.org.au")
  scp_download(session, 
               paste0("/home/scho0011/kg98/Sid/STAGES/STAGES_dti/data/",sub,"/ses-1/anat/",sub,"_ses-1_T1w.nii"), 
               to = "~/Dropbox/Sid/R_files/STAGES_difussion/temp/")
  
  
  t1 <- readNIfTI2(paste0(sub,"_ses-1_T1w.nii"))
  
  png(file="%02d.png", width=200, height=200)
  for (i in c(100,130,140,150,160,175,190)){
    oro.nifti::image(robust_window(t1) ,z = i, w = 1,
                     plot.type = "single",
                     zlim = c(0,800))
    
  }
  dev.off()

  
  list.files(pattern = '*.png', full.names = TRUE ,) %>% 
    gtools::mixedsort(decreasing = F) %>%
    image_read() %>% # reads each path file
    image_append() %>% # animates, can opt for number of loops
    image_annotate(paste0(sub), size = 15, gravity = "south", color = "white") %>%
    image_write(paste0(outdir,"0_t1/",sub,"_t1.gif"))# write to current dir
  
  file.remove(list.files(pattern=".png"))
  unlink("*") #clear temp folder
  
  #======================================
  # 1. Mean b0 
  #=====================================
  session <- ssh_connect("scho0011@m3.massive.org.au")
  
  scp_download(session, 
               paste0("/home/scho0011/kg98/Sid/STAGES/STAGES_dti/data/",sub,"/ses-1/dwi/syn/b0.nii.gz"), 
               to = "~/Dropbox/Sid/R_files/STAGES_difussion/temp/")
  meanB0 <- readNIfTI2("b0.nii.gz")
  
  png(file="Meanb0%02d.png", width=200, height=200)
  for (i in c(10, 20, 30, 40,50)){
    oro.nifti::image(robust_window(meanB0) ,z = i, w = 1,
                     plot.type = "single",
                     zlim = c(0,1300))
    
  }
  
  dev.off()

  
  list.files(pattern = '*.png', full.names = TRUE) %>% 
    gtools::mixedsort(decreasing = F) %>%
    image_read() %>% # reads each path file
    image_append() %>% # animates, can opt for number of loops
    image_annotate(paste0(sub), size = 15, gravity = "south", color = "white") %>%
    image_write(paste0(outdir,"1_B0/",sub,"_B0.png")) # write to current dir
  
  file.remove(list.files(pattern=".png"))
  unlink("*") #clear temp folder
  
  #======================================
  # 2. raw dwi
  #=====================================
  session <- ssh_connect("scho0011@m3.massive.org.au")
  
  scp_download(session, 
               paste0("/home/scho0011/kg98/Sid/STAGES/STAGES_dti/data/",sub,"/ses-1/dwi/",sub,"_ses-1_dwi.nii"), 
               to = "~/Dropbox/Sid/R_files/STAGES_difussion/temp/")
  
  rawDWI <- readNIfTI2(paste0(sub,"_ses-1_dwi.nii"))
  
  png(file="rawDWI%02d.png", width=200, height=200)
   for (i in 11:rawDWI@dim_[5]){
    orthographic(rawDWI, crosshairs = F, w = i, useRaster = TRUE,  zlim = c(0,180)
                 , text = i-10)
  }
  
  dev.off()
  list.files(pattern = '*.png', full.names = TRUE) %>% 
    gtools::mixedsort(decreasing = F) %>%
    image_read() %>% # reads each path file
    image_append(stack = TRUE) %>% #
    image_annotate(paste0(sub), size = 15, gravity = "west", color = "white") %>%
    image_write(paste0(outdir,"2_raw_dwi/",sub,"_dwi.png")) # write to current dir
  
  file.remove(list.files(pattern=".png"))
  unlink("*") #clear temp folder
  
  
  #======================================
  # 3. denoised dwi 
  #=====================================
  session <- ssh_connect("scho0011@m3.massive.org.au")
  
  scp_download(session, 
               paste0("/home/scho0011/kg98/Sid/STAGES/STAGES_dti/data/",sub,"/ses-1/dwi/",sub,"_ses-1_dwi_denoised.nii"), 
               to = "~/Dropbox/Sid/R_files/STAGES_difussion/temp/")
  
  DWI_denoised <- readNIfTI2(paste0(sub,"_ses-1_dwi_denoised.nii"))
  png(file="DWI_denoised%02d.png", width=200, height=200)
  for (i in 11:DWI_denoised@dim_[5]){
      orthographic(DWI_denoised, crosshairs = F, w = i, useRaster = TRUE,  zlim = c(0,180)
                   , text = i-10)
    }
  
  
  dev.off()
  list.files(pattern = '*.png', full.names = TRUE) %>% 
    gtools::mixedsort(decreasing = F) %>%
    image_read() %>% # reads each path file
    image_append(stack = TRUE) %>% #
    image_annotate(paste0(sub), size = 15, gravity = "north", color = "white") %>%
    image_write(paste0(outdir,"3_denoised_dwi/",sub,"_dwi.png")) # write to current dir
  
  file.remove(list.files(pattern=".png"))
  unlink("*")
  #======================================
  # 3. Synb0-DISCO (SDC)
  #=====================================
  session <- ssh_connect("scho0011@m3.massive.org.au")
  
  scp_download(session, 
               paste0("/home/scho0011/kg98/Sid/STAGES/STAGES_dti/data/",sub,"/ses-1/dwi/b0_all.nii.gz"), 
               to = "~/Dropbox/Sid/R_files/STAGES_difussion/temp/")
  
  #mean B0 to synb0
  synb0 <- readNIfTI2("b0_all.nii.gz")
  #t1 <- readNIfTI2("ses-1/dwi/syn/T1.nii.gz")
  png(file="3_synb0%02d.png", width=300, height=300)
  for (i in 1:synb0@dim_[5]){
    oro.nifti::orthographic(synb0, w = i,
                            zlim = c(0,1300),
                            crosshairs = FALSE)
  }
  dev.off()
  
  
  
  b0 <- image_read("3_synb001.png")
  b0 <- image_annotate(b0, paste0(sub, " B0 raw"), size = 15, gravity = "south", color = "white")
  synb0 <- image_read("3_synb002.png")
  synb0 <- image_annotate(synb0, paste0(sub, " B0 corrected"), size = 15, gravity = "south", color = "white")
  new <- image_resize(c(b0, synb0), '200x200')
  image_morph(new) %>%
    image_animate(optimize = TRUE,delay = 20 ) %>%
    image_write(paste0(outdir,"4_synB0/",sub,"_synb0.gif"))
  
  file.remove(list.files(pattern=".png"))
  unlink("*")
  
  #======================================
  # 4. eddied + top-up dwi image
  #=====================================
  session <- ssh_connect("scho0011@m3.massive.org.au")
  
  scp_download(session, 
               paste0("/home/scho0011/kg98/Sid/STAGES/STAGES_dti/data/",
                      sub,"/ses-1/dwi/",sub,"_ses-1_dwi_denoised_eddy.nii.gz"), 
               to = "~/Dropbox/Sid/R_files/STAGES_difussion/temp/")
  
  
  
  DWI_eddy <- readNIfTI2(paste0(sub,"_ses-1_dwi_denoised_eddy.nii.gz"))
  png(file="4_DWI_eddy0%02d.png", width=200, height=200)
  for (i in 11:DWI_eddy@dim_[5]){
    orthographic(DWI_eddy, crosshairs = F, w = i, useRaster = TRUE,  zlim = c(0,180)
                 , text = i-10)
    }
    
  
  
  dev.off()
  list.files(pattern = '*.png', full.names = TRUE) %>% 
    gtools::mixedsort(decreasing = F) %>%
    image_read() %>% # reads each path file
    image_append(stack = TRUE) %>% #
    image_annotate(paste0(sub), size = 15, gravity = "north", color = "white") %>%
    image_write(paste0(outdir,"5_eddy/",sub,"_eddy.gif")) # write to current dir
  
  file.remove(list.files(pattern=".png"))
  unlink("*")
  
  
  #======================================
  # 5. Tractoflow  - BET + n4 + crop + normalise DWI + resample
  #=====================================
  session <- ssh_connect("scho0011@m3.massive.org.au")
  
  scp_download(session, 
               paste0("/home/scho0011/kg98/Sid/STAGES/STAGES_dti/scripts/TractoFlow/results/",
                      sub,"/Resample_DWI/",sub,"__dwi_resampled.nii.gz"), 
               to = "~/Dropbox/Sid/R_files/STAGES_difussion/temp/")
  
  
  DWI_resample <- readNIfTI2(paste0(sub,"__dwi_resampled.nii.gz"))
  
  png(file="5_DWI_resample%02d.png", width=200, height=200)

  for (i in 11:DWI_resample@dim_[5]) {
    orthographic(DWI_resample, crosshairs = F, w = i, useRaster = TRUE,  zlim = c(0,600)
                 , text = i-10)
  }
  
dev.off()
  
  list.files(pattern = '*.png', full.names = TRUE) %>% 
    gtools::mixedsort(decreasing = F) %>%
    image_read() %>% # reads each path file
    image_append(stack = TRUE) %>% #
    image_annotate(paste0(sub), size = 15, gravity = "north", color = "white") %>%
    image_write(paste0(outdir,"6_tracto_dwi/",sub,"_dwi.png")) # write to current dir
  
  file.remove(list.files(pattern=".png"))

  
  #======================================
  # 8. Tractoflow  - T1-denoise+N4+Resample+BET+crop+register
  #=====================================
  session <- ssh_connect("scho0011@m3.massive.org.au")
  scp_download(session, 
               paste0("/home/scho0011/kg98/Sid/STAGES/STAGES_dti/scripts/TractoFlow/results/",
                      sub,"/Register_T1/",sub,"__t1_warped.nii.gz"), 
               to = "~/Dropbox/Sid/R_files/STAGES_difussion/temp/")
  
  t1 <- readNIfTI2(paste0(sub,"__t1_warped.nii.gz"))
  
  png(file="warpedT1%02d.png", width=200, height=200)
  for (i in c(20,40,50,60,75,90,100)){
    oro.nifti::image(robust_window(t1) ,z = i, w = 1,
                     plot.type = "single",
                     zlim = c(0,700))
    
  }
  dev.off()
  
  list.files(pattern = '*.png', full.names = TRUE ,) %>% 
    gtools::mixedsort(decreasing = F) %>%
    image_read() %>% # reads each path file
    image_append() %>% # animates, can opt for number of loops
    image_annotate(paste0(sub), size = 15, gravity = "south", color = "white") %>%
    image_write(paste0(outdir,"7_tracto_t1/",sub,"_t1.png")) # write to current dir
  
  file.remove(list.files(pattern=".png"))
  

  
  #======================================
  # 9. Tractoflow  - T1 to DWI registration
  #=====================================
  session <- ssh_connect("scho0011@m3.massive.org.au")
  
  #mean B0 to synb0
  
  png(file="7_tf_t1.png", width=300, height=300)
    oro.nifti::slice(robust_window(t1), z=60,
                            zlim = c(0,700),
                            crosshairs = FALSE)

  
  dev.off()
  
  png(file="7_tf_dwi.png", width=300, height=300)

  oro.nifti::slice(robust_window(DWI_resample), w=1,z=60,
                   zlim = c(0,4200),
                   crosshairs = FALSE)
  
  dev.off()

  dwi <- image_read("7_tf_dwi.png")
  dwi <- image_annotate(dwi, paste0(sub, " b0"), size = 12, gravity = "south", color = "white")
  t1 <- image_read("7_tf_t1.png")
  t1 <- image_annotate(t1, paste0(sub, " t1"), size = 12, gravity = "south", color = "white")
  new <- image_resize(c(t1, dwi), '200x200')
  image_morph(new) %>%
    image_animate(optimize = TRUE,delay = 25 ) %>%
    image_write(paste0(outdir,"8_tracto_b0_2_t1/",sub,"_t12b0.gif"))
  
  file.remove(list.files(pattern=".png"))
  unlink("*")
  
  
  #======================================
  # 9. Tractoflow  - DTI metrics  (only FA)
  #=====================================
  session <- ssh_connect("scho0011@m3.massive.org.au")
  scp_download(session, 
               paste0("/home/scho0011/kg98/Sid/STAGES/STAGES_dti/scripts/TractoFlow/results/",
                      sub,"/DTI_Metrics/",sub,"__fa.nii.gz"), 
               to = "~/Dropbox/Sid/R_files/STAGES_difussion/temp/")
  
  fa <- readNIfTI2(paste0(sub,"__fa.nii.gz"))
  
  png(file="fa%02d.png", width=200, height=200)
  
  for (i in c(20,40,50,60,75,90,100)){
    oro.nifti::slice(robust_window(fa) ,z = i, 
                     zlim = c(0,1))
    
  }
  
  
  
  dev.off()
  
  list.files(pattern = '*.png', full.names = TRUE) %>% 
    gtools::mixedsort(decreasing = F) %>%
    image_read() %>% # reads each path file
    image_append() %>% # animates, can opt for number of loops
    image_annotate(paste0(sub), size = 11, gravity = "southwest", color = "white") %>%
    image_write(paste0(outdir,"9_tracto_fa/",sub,"_fa.png")) # write to current dir
 
  
  file.remove(list.files(pattern=".png"))
  
  
  #======================================
  # 8. Tractoflow  - fODF metrics, first 10 directions
  #=====================================
  session <- ssh_connect("scho0011@m3.massive.org.au")
  scp_download(session, 
               paste0("/home/scho0011/kg98/Sid/STAGES/STAGES_dti/scripts/TractoFlow/results/",
                      sub,"/DTI_Metrics/",sub,"__tensor.nii.gz"), 
               to = "~/Dropbox/Sid/R_files/STAGES_difussion/temp/")
  
  system("fsleyes render sub-076c__fa.nii.gz sub-076c__tensor.nii.gz -ot tensor")
  system("pwd")
  fODF <- readNIfTI2(paste0(sub,"__tensor.nii.gz"))

}


