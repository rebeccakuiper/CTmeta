
Check_Sigma <- function(Sigma, q){

  if(q == 0){
    print(paste0("The argument q = 0, it should be larger than 0. Notably, Sigma (and Drift) is a q times q matrix."))
    stop()
  }

  # Checks on Sigma
  if(length(Sigma) != 1){
    if(length(dim(Sigma)) < 2){
      print(paste0("The residual covariance matrix Sigma should be a square matrix. It should be of size q times q, with q = ", q))
      stop()
    }else if(length(dim(Sigma)) > 2){
      print(paste0("The residual covariance matrix Sigma should, like Drift (or Phi), be an q times q matrix, with q = ", q, ". Currently, it is of size ", dim(Sigma)))
      stop()
    }else if(dim(Sigma)[1] != dim(Sigma)[2]){
      print(paste0("The residual covariance matrix Sigma should, like Drift (or Phi), be a square matrix of size q times q, with q = ", q, ". Currently, it is of size ", dim(Sigma)[1], " times ", dim(Sigma)[2]))
      stop()
    }else if(dim(Sigma)[1] != q){
      print(paste0("The residual covariance matrix Sigma should, like Drift (or Phi), be a matrix of size q times q, with q = ", q, ". Currently, it is of size ", dim(Sigma)))
      stop()
    }
  }else if(q != 1){
    print(paste0("The residual covariance matrix Sigma is a scalar, but it should be of size q times q, with q = ", q))
    stop()
  }

}
