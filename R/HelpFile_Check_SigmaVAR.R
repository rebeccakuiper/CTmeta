
Check_SigmaVAR <- function(SigmaVAR, q){

  # Checks on SigmaVAR
  if(length(SigmaVAR) != 1){
    if(length(dim(SigmaVAR)) < 2){
      print(paste0("The residual covariance matrix SigmaVAR should be a square matrix. It should be of size q times q, with q = ", q))
      stop()
    }else if(length(dim(SigmaVAR)) > 2){
      print(paste0("The residual covariance matrix SigmaVAR should, like Phi, be an q times q matrix, with q = ", q, ". Currently, it is of size ", dim(SigmaVAR)))
      stop()
    }else if(dim(SigmaVAR)[1] != dim(SigmaVAR)[2]){
      print(paste0("The residual covariance matrix SigmaVAR should, like Phi, be a square matrix of size q times q, with q = ", q, ". Currently, it is of size ", dim(SigmaVAR)[1], " times ", dim(SigmaVAR)[2]))
      stop()
    }else if(dim(SigmaVAR)[1] != q){
      print(paste0("The residual covariance matrix SigmaVAR should, like Phi, be a matrix of size q times q, with q = ", q, ". Currently, it is of size ", dim(SigmaVAR)))
      stop()
    }
  }else if(q != 1){
    print(paste0("The residual covariance matrix SigmaVAR is a scalar, but it should be of size q times q, with q = ", q))
    stop()
  }

}
