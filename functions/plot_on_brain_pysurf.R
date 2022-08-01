plot_on_brain_pysurf <- function(plot_vector = NULL,
                                 temppath='/output/figures/temp/',
                                 min=NULL,
                                 max=NULL,
                                 colourscale='viridis', 
                                 jpeg = FALSE,
                                 clean_files = TRUE,
                                 surf = "inflated",
                                 no_subcortex = F) {
  
  library(cowplot)
  library(ggplot2)
  library(reticulate)
  library(magick)
  library(RColorBrewer)
  library(rstudioapi)
  
  if(is.null(min)) {
    min <- min(plot_vector)
  }
  
  if(is.null(max)) {
    max <- max(plot_vector)
  }
  
  #write out degree vector for cortex
  plot_vector <- c(plot_vector)
  if(length(plot_vector)!=332) {
    stop("Plot vector should be 332 long")
  }
  
  
  #split plot vector into lh_cortex,rh_cortex and subcortex
  lh_plotvec <- plot_vector[1:150]
  rh_plotvec <- plot_vector[151:300]
  sub_plotvec <- plot_vector[301:332]
  
  #add  -1 to cortical files with for medial wall 
  lh_plotvec <- c(0,lh_plotvec) 
  rh_plotvec <- c(0,rh_plotvec) 
  write.table(lh_plotvec, paste0(temppath, "temp_degree_lh.txt"), quote = F, row.names = F, col.names = F)
  write.table(rh_plotvec, paste0(temppath, "temp_degree_rh.txt"),quote = F, row.names = F, col.names = F)
  
  pyactivate <- "source activate /opt/anaconda3/envs/pysurfer \n"
  pyexecute <- paste0("python /functions/make_pysurf.py -r ",
                      paste0(temppath, "temp_degree_rh.txt"),
                      " -l ", paste0(temppath, "temp_degree_lh.txt"),
                      " -m ", min,
                      " -a ", max,
                      " -c ", colourscale, 
                      " -s ", surf, " \n")
  
  myTerm <- terminalCreate()
  Sys.sleep(2)
  terminalSend(myTerm, pyactivate)
  Sys.sleep(2)
  terminalSend(myTerm, pyexecute)
  Sys.sleep(2)
  repeat{
    Sys.sleep(1)
    if(rstudioapi::terminalBusy(myTerm) == FALSE){
      print("Surfaces made. ")
      rstudioapi::terminalKill(myTerm)
      break
    }
  }
  
  #trim and remove background
  
  list_c <- list.files(path = temppath,
                       pattern = 'surf_*', full.names = T) #get list of file names you want
  
  for (x in list_c) {
    pic  <- image_read(x)
    tpic <- image_transparent(pic, 'white')
    tpic_c <- image_trim(tpic)
    image_write(tpic_c, path = x, format = "png") # new file is t_$file_name
  }
  
  if (no_subcortex == F) {
    ### Addin subcortex
    #use pyvista to make a png with subcortical surface mesh
    #source_python("/Users/sidchopra/Dropbox/Sid/R_files/STAGES_difussion/scripts/functions/get_subcortex_mesh_tian2.py")
    #get_subcortex_mesh_tian2(numvec = sub_plotvec, min=min,max=max,colourmap=colourscale,outputfolder=temppath)
    
    write.table(sub_plotvec, paste0(temppath, "sub_degree.txt"), quote = F, row.names = F, col.names = F)
    
    pyexecute <- paste0("python /functions/get_subcortex_mesh_tian2.py -n ",
                        paste0(temppath, "sub_degree.txt"),
                        " -m ", min,
                        " -a ", max,
                        " -c ", colourscale, 
                        " -o ", temppath, " \n")
    
    myTerm <- terminalCreate()
    Sys.sleep(2)
    terminalSend(myTerm, pyactivate)
    Sys.sleep(2)
    terminalSend(myTerm, pyexecute)
    Sys.sleep(2)
    repeat{
      Sys.sleep(1)
      if(rstudioapi::terminalBusy(myTerm) == FALSE){
        print("Subcortex made.")
        rstudioapi::terminalKill(myTerm)
        break
      }
    }
    
    
    
    #Trim , flip and remove background from subcortex pngs
    list_s <- list.files(path = temppath,
                         pattern = 'sub_temp*', full.names = T) #get list of file names you want
    
    
    #make bg transpatent
    for (x in list_s) {
      pic  <- image_read(x)
      tpic <- image_transparent(pic, 'white')
      tpic_c <- image_crop(tpic, "810x700+250")
      tpic_c <- image_trim(tpic)
      image_write(tpic_c, path = x, format = "png") # new file is t_$file_name
    }
    
    #flip subcortex images
    for (x in list_s) {
      pic  <- image_read(x)
      tpic <- image_flop(pic)
      tpic_c <- image_crop( tpic, "810x500")
      image_write(tpic_c, path = x, format = "png")
    }
    
    
    #add subcortex to medial images
    lm <- image_read(list_c[2])
    ls <- image_read(list_s[2]) #backwards on purpose, because the images have been fliped
    lm_s <- image_composite(lm, image_scale(ls,'490'), offset = "+190+222")
    
    
    rm <- image_read(list_c[4])
    rs <- image_read(list_s[1]) #backwards on purpose, because the images have been fliped
    rm_s <- image_composite(rm, image_scale(rs,'490'), offset = "+110+222")
    
    
    ll <-  image_read(list_c[1])
    rl <- image_read(list_c[3])
    
    
    rl_ll <- image_append(c(rl, ll))
    rm_lm <- image_append(c(rm_s, lm_s))
    all <- image_append(c(rl_ll, rm_lm), stack = T)
    
  } else {
    list_c <- list.files(path = temppath,
                         pattern = 'surf_*', full.names = T)
    lm <- image_read(list_c[2])
    rm <- image_read(list_c[4])
    ll <-  image_read(list_c[1])
    rl <- image_read(list_c[3])
    
    
    rl_ll <- image_append(c(rl, ll))
    rm_lm <- image_append(c(rm, lm))
    all <- image_append(c(rl_ll, rm_lm), stack = T)
  }
  
  #addcolour bar
  mat = matrix(rnorm(400,1,1), ncol=20)
  
  #change colour bar here
  #grad = viridis::viridis(n=100)
  #grad = RColorBrewer::brewer.pal(9,"coolwarm")
  #grad = pals::coolwarm(n=100)
  grad =  pals::brewer.reds(100)
  
  
  tiff(paste0(temppath, "temp_cb.tiff"), width=5, height=5, res=300, units = "in")
  print({lattice::levelplot(mat, col.regions=grad, colorkey = list(at=seq(min,max,length.out=100)))})
  dev.off()
  cb <- image_read(paste0(temppath,"temp_cb.tiff"))
  cb <- image_chop(cb, "1300x150")
  cb <- image_trim(cb)
  cb <- image_transparent(cb, 'white')
  
  final <- image_composite(image_scale(all, '1500'),image_scale(cb,'60'), offset = "+1410+0")
  if(jpeg==T){ image_write(final, paste0(temppath, "final_output.jpeg"), quality = 100) }
  
  if(clean_files==T) {
    remove <- list.files(temppath, pattern = "*temp*", full.names = T)
    file.remove(remove)
    remove <- list.files(temppath, pattern = "*surf*", full.names = T)
    file.remove(remove)
    remove <- list.files(temppath, pattern = "*degree*", full.names = T)
    file.remove(remove)
  }
  #
  return(final)
}
