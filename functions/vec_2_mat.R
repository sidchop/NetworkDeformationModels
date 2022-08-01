# input = the vector you want to transform into a symetric square matrix
# length = length of rows/column of matrix
# diag = what value you want in the diagonal 
# lower.tri = should the lower triange be filled first [Always check that the matrix has been filled correctly]
# to visualise the matrix: pheatmap::pheatmap(martix, cluster_rows = F, cluster_cols = F)

vec_2_mat <- function(vec, length, diag, lower.tri = FALSE, supress = F) {
  if (supress == F) {message("always check upper or lower triangle.")}
  temp <- matrix(nrow=length, ncol=length)
  if(lower.tri == TRUE) {temp[lower.tri(temp)] <- as.matrix(vec) 
  temp <- Matrix::forceSymmetric(temp, uplo = "L")} else 
  { temp[upper.tri(temp)] <- as.matrix(vec)
  temp <- Matrix::forceSymmetric(temp, uplo = "U")}
   diag(temp) <- diag
return(as.matrix(temp))
}
