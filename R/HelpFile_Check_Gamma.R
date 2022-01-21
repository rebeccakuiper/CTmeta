
Check_Gamma <- function(Gamma, q){

  if(q == 0){
    print(paste0("The argument q = 0, but q should be larger than 0."))
    stop()
  }

  # Checks on Gamma
  if(length(Gamma) != 1){
    if(length(dim(Gamma)) < 2){
      print(paste0("The stationary covariance matrix Gamma should be a square matrix. It should be of size q x q, with q = ", q))
      stop()
    }else if(length(dim(Gamma)) > 2){
      print(paste0("The stationary covariance matrix Gamma should be a q times q matrix, with q = ", q, "."))
      stop()
    }else if(dim(Gamma)[1] != dim(Gamma)[2]){
      print(paste0("The stationary covariance matrix Gamma should, like Phi, be a square matrix of size q x q, with q = ", q, ". Currently, it is of size ", dim(Gamma)[1], " x ", dim(Gamma)[2]))
      stop()
    }else if(dim(Gamma)[1] != q){
      print(paste0("The stationary covariance matrix Gamma should, like Phi, be a matrix of size q x q, with q = ", q, "."))
      stop()
    }
  }else if(q != 1){
    print(paste0("The stationary covariance matrix Gamma is a scalar, but it should be of size q x q, with q = ", q))
    stop()
  }

}
