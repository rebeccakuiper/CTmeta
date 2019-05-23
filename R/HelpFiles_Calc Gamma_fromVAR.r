
#' Calculates the stationary covariance matrix from the discrete-time (un)standardized lagged effects matrix and residual covariance matrix 
#'
#' @param Phi (Un)standardized lagged effects matrix. 
#' @param SigmaVAR Residual covariance matrix of the first-order discrete-time vector autoregressive (DT-VAR(1)) model.
#' 
#' @return Stationary covariance matrix, that is, the contemporaneous covariance matrix of the data.
#' @export
#' @examples
#' 
#' q = 2
#' Phi <- matrix(c(0.50, 0.15, 0.25, 0.40), byrow=T, ncol = q)
#' SigmaVAR <- diag(q) # for ease
#' calc.Gamma.fromVAR(Phi, SigmaVAR)


calc.Gamma.fromVAR <- function(Phi, SigmaVAR) {
  
q <- dim(Phi)[1]

vecS <- as.vector(SigmaVAR)
PhiKronPhi <- kronecker(Phi, Phi)
vecGamma <- solve(diag(q*q) - PhiKronPhi) %*% vecS

Gamma_fromVAR <- matrix(vecGamma, nrow=q)


return(Gamma_fromVAR)
}