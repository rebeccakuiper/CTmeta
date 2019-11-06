
#' Calculate the stationary covariance matrix Gamma
#'
#' Calculates the stationary covariance matrix (Gamma) from the discrete-time (un)standardized lagged effects matrix (Phi) and residual covariance matrix (SigmaVAR). There is also an interactive web application on my website: Standardizing and/or transforming lagged regression estimates (\url{https://www.uu.nl/staff/RMKuiper/Websites\%20\%2F\%20Shiny\%20apps}).
#'
#' @param Phi (Un)standardized lagged effects matrix.
#' @param SigmaVAR Residual covariance matrix of the first-order discrete-time vector autoregressive (DT-VAR(1)) model.
#'
#' @return This function returns the stationary covariance matrix, that is, the contemporaneous covariance matrix of the data.
#' @export
#' @examples
#'
#' Phi <- myPhi[1:2,1:2]
#' #Phi <- matrix(c(0.50, 0.15, 0.25, 0.40), byrow=T, ncol = 2)
#' SigmaVAR <- diag(q) # for ease
#' calc.Gamma.fromVAR(Phi, SigmaVAR)


calc.Gamma.fromVAR <- function(Phi, SigmaVAR) {

  # Check on Phi
  if(length(Phi) == 1){
    q <- 1
  }else{
    #
    if(is.null(dim(Phi))){
      if(!is.null(length(Phi))){
        print(paste("The argument Phi is not a matrix of size q times q."))
        stop()
      }else{
        print(paste("The argument Phi is not found: The continuous-time lagged effects matrix Phi is unknown, but should be part of the input."))
        stop()
      }
    }else{
      if(dim(Phi)[1] != dim(Phi)[2] | length(dim(Phi)) != 2){
        print(paste("The argument Phi is not a matrix of size q times q."))
        stop()
      }
      q <- dim(Phi)[1]
    }
  }

  # Checks on SigmaVAR
  if(length(SigmaVAR) != 1){
    if(dim(SigmaVAR)[1] != dim(SigmaVAR)[2]){ # Should be square
      print(paste("The residual covariance matrix SigmaVAR should be a square matrix of size q times q, with q = ", q, "."))
      stop()
    }
    if(dim(SigmaVAR)[1] != q){ # Should have same dimension as Phi
      print(paste("The residual covariance matrix SigmaVAR should, like Phi, be a matrix of size q times q, with q = ", q, "."))
      stop()
    }
    if(length(dim(SigmaVAR)) > 2){ # Should be matrix, not array
      print(paste("The residual covariance matrix SigmaVAR should be an q times q matrix, with q = ", q, "."))
      stop()
    }
  }else if(q != 1){
    print(paste("The residual covariance matrix SigmaVAR should, like Phi, be a scalar."))
    stop()
  }



if(q > 1){
  vecS <- as.vector(SigmaVAR)
  PhiKronPhi <- kronecker(Phi, Phi)
  vecGamma <- solve(diag(q*q) - PhiKronPhi) %*% vecS

  Gamma_fromVAR <- matrix(vecGamma, nrow=q)
} else{
  Gamma_fromVAR <- SigmaVAR / (1 - Phi^2)
}

return(Gamma_fromVAR)
}
