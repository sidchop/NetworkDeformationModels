# This function is basically a copy of roi2ts from the XCP engine:
# https://xcpengine.readthedocs.io/utils/roiquants.html
# All credit to the authors
# https://github.com/PennBBL/xcpEngine/blob/master/utils/roi2ts.R
# It is designed to mimic the functionality of fsl_meants, but in R. 
# Note that fsl_meants includes all voxels, this function only includes non-zero-voxels (akin to fslstats -K)

# Edited to work only with 3D images (sid - 09/04/2021)

nii_roi2ts <- function(img, net) {
library("RNifti")
library("pracma")


labs    <- NA 
###################################################################
# 3. Compute the network timeseries. This functionality is based
#    on the matrix2timeseries function from ANTsR, written by
#    Shrinidhi KL.
#
# First, obtain all unique nonzero values in the mask.
###################################################################
labels                  <- sort(unique(net[net > 0]))
if (is.na(labs[1])) {
  labs                 <- labels
}
###################################################################
# Create a logical over all voxels, indicating whether each
# voxel has a nonzero mask value.
###################################################################
logmask                 <- (net > 0)
###################################################################
# Use the logical mask to subset the 4D data. Determine the
# dimensions of the timeseries matrix: they should equal
# the number of voxels in the mask by the number of time points.
###################################################################
mat                     <- img[logmask]
#dim(mat)                <- c(sum(logmask), dim(img)[length(dim(img))]) #for3d data need to edit
dim(mat)                <- c(sum(logmask), 1) 
mat                     <- t(mat) 
###################################################################
# Determine how many unique values are present in the RoI map.
#  * If only one unique value is present, then the desired
#    output is a voxelwise timeseries matrix, which has already
#    been computed.
#  * If multiple unique values are present in the map, then
#    the map represents a network, and the desired output is
#    a set of mean node timeseries.
###################################################################
if (length(labs) == 1) {
  mmat                 <- mat
} else {
  mmat                 <- zeros(dim(mat)[1],length(labs))
  ################################################################
  # If the script enters this statement, then there are multiple
  # unique values in the map, indicating multiple mask RoIs: a
  # network timeseries analysis should be prepared.
  #
  # Prime the modified matrix. Extract the timeseries of all
  # voxels in the first RoI submask. If only one voxel is in the
  # RoI, then the extracted timeseries will lack dimension
  # according to R; it must be made into a column vector so that
  # it can be appended to the modified matrix. The user is warned,
  # as singleton voxels are more susceptible to artefactual
  # influence. If multiple voxels are in the RoI, then the mean
  # RoI timeseries is computed and added to the model.
  ################################################################
  nodevec              <- net[logmask]
  if (labs[1] %in% labels) {
    voxelwise         <- mat[, nodevec == labs[1]]
    if (is.null(dim(voxelwise)) && !is.null(length(voxelwise))) {
      warning("Warning: node 1 contains one voxel\n")
      dim(voxelwise) <-c(length(voxelwise),1)
    }
    mmat[,1]          <- matrix(mean(voxelwise[voxelwise!=0]), ncol = 1)
  }
  ################################################################
  # Repeat for all remaining RoIs.
  ################################################################
  for (i in 2:length(labs)) {
    if (! labs[i] %in% labels) { next }
    voxelwise         <-mat[, nodevec == labs[i]]
    if (is.null(dim(voxelwise)) && !is.null(length(voxelwise))) {
      warning(paste("Warning: node ", labs[i], " contains one voxel\n"))
      dim(voxelwise) <-c(length(voxelwise),1)
    }
    mmat[,i]          <- matrix(mean(voxelwise[voxelwise!=0]), ncol = 1)
  }
  colnames(mmat)       <- paste0("L_", labs)
}

return(t(mmat))
}
