
Check_SigmaVAR <- function(SigmaVAR, q){

  if(q == 0){
    print(paste0("q = 0, but we should have q > 0. SigmaVAR and Phi are q x q matrices."))
    stop()
  }

  # Checks on SigmaVAR
  if(length(SigmaVAR) != 1){
    if(length(dim(SigmaVAR)) < 2){
      print(paste0("The residual covariance matrix SigmaVAR should be a square matrix. It should be with dimensions q x q, with q = ", q))
      stop()
    }else if(length(dim(SigmaVAR)) > 2){
      print(paste0("The residual covariance matrix SigmaVAR should, like Phi, be a q x q matrix, with q = ", q, "."))
      stop()
    }else if(dim(SigmaVAR)[1] != dim(SigmaVAR)[2]){
      print(paste0("The residual covariance matrix SigmaVAR should, like Phi, be a square matrix with dimensions q x q, with q = ", q, ". In the given input, it has dimensions ", dim(SigmaVAR)[1], " x ", dim(SigmaVAR)[2]))
      stop()
    }else if(dim(SigmaVAR)[1] != q){
      print(paste0("The residual covariance matrix SigmaVAR should, like Phi, be a matrix with dimensions q x q, with q = ", q, "."))
      stop()
    }
  }else if(q != 1){
    print(paste0("The residual covariance matrix SigmaVAR is a scalar, but it should be a matrix with dimensions q x q, with q = ", q))
    stop()
  }

}
