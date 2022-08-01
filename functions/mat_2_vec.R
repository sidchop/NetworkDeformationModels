mat_2_vec <- function(mat) {
 vec <-  t(mat[upper.tri(mat, diag = FALSE)])
  return(vec)
}
