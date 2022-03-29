
Check_Sigma <- function(Sigma, q){

  if(q == 0){
    stop(paste0("The argument q = 0, it should be larger than 0. Notably, Sigma (and Drift) is a q x q matrix."))
  }

  # Checks on Sigma
  if(length(Sigma) != 1){
    if(length(dim(Sigma)) < 2){
      stop(paste0("The residual covariance matrix Sigma should be a square matrix. It should be of size q x q, with q = ", q))
    }else if(length(dim(Sigma)) > 2){
      stop(paste0("The residual covariance matrix Sigma should, like Drift (or Phi), be a q x q matrix, with q = ", q, "."))
    }else if(dim(Sigma)[1] != dim(Sigma)[2]){
      stop(paste0("The residual covariance matrix Sigma should, like Drift (or Phi), be a square matrix of size q x q, with q = ", q, ". Currently, it is of size ", dim(Sigma)[1], " x ", dim(Sigma)[2]))
    }else if(dim(Sigma)[1] != q){
      stop(paste0("The residual covariance matrix Sigma should, like Drift (or Phi), be a matrix of size q x q, with q = ", q, "."))
    }
  }else if(q != 1){
    stop(paste0("The residual covariance matrix Sigma is a scalar, but it should be of size q x q, with q = ", q))
  }
  if(!is.numeric(Sigma)){
    stop("Sigma should be a numerical matrix.")
  }

}
