
Check_SigmaVAR <- function(SigmaVAR, q){

  if(q == 0){
    stop(paste0("q = 0, but we should have q > 0. SigmaVAR and Phi are q x q matrices."))
  }

  # Checks on SigmaVAR
  if(length(SigmaVAR) != 1){
    if(length(dim(SigmaVAR)) < 2){
      stop(paste0("The residual covariance matrix SigmaVAR should be a square matrix. It should have dimensions q x q, with q = ", q))
    }else if(length(dim(SigmaVAR)) > 2){
      stop(paste0("The residual covariance matrix SigmaVAR should, like Phi, be a q x q matrix, with q = ", q, "."))
    }else if(dim(SigmaVAR)[1] != dim(SigmaVAR)[2]){
      stop(paste0("The residual covariance matrix SigmaVAR should, like Phi, be a square matrix with dimensions q x q, with q = ", q, ". In the given input, it has dimensions ", dim(SigmaVAR)[1], " x ", dim(SigmaVAR)[2]))
    }else if(dim(SigmaVAR)[1] != q){
      stop(paste0("The residual covariance matrix SigmaVAR should, like Phi, be a matrix with dimensions q x q, with q = ", q, "."))
    }
  }else if(q != 1){
    stop(paste0("The residual covariance matrix SigmaVAR is a scalar, but it should be a matrix with dimensions q x q, with q = ", q))
  }

}
