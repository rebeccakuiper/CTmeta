
#' Calculate the stationary covariance matrix Gamma
#'
#' Calculates the stationary covariance matrix (Gamma) from the discrete-time (un)standardized lagged effects matrix (Phi) and residual covariance matrix (SigmaVAR). There is also an interactive web application on my website: Standardizing and/or transforming lagged regression estimates (\url{https://www.uu.nl/staff/RMKuiper/Websites\%20\%2F\%20Shiny\%20apps}).
#'
#' @param Phi (Un)standardized lagged effects matrix.
#' It also takes a fitted object from the classes "varest" (from the VAR() function in vars package) and "ctsemFit" (from the ctFit() function in the ctsem package); see example below. From such an object, the Phi and SigmaVAR matrices are calculated/extracted.
#' @param SigmaVAR Residual covariance matrix of the first-order discrete-time vector autoregressive (DT-VAR(1)) model.
#'
#' @return This function returns the stationary covariance matrix, that is, the contemporaneous covariance matrix of the data.
#' @importFrom expm expm
#' @export
#' @examples
#'
#' # library(CTmeta)
#'
#' ## Example 1 ##
#' #
#' Phi <- myPhi[1:2,1:2]
#' #Phi <- matrix(c(0.50, 0.15, 0.25, 0.40), byrow=T, ncol = 2)
#' SigmaVAR <- diag(q) # for ease
#' Gamma.fromVAR(Phi, SigmaVAR)
#'
#'
#' ## Example 2: input from fitted object of class "varest" ##
#' #
#' data <- myData
#' if (!require("vars")) install.packages("vars")
#' library(vars)
#' out_VAR <- VAR(data, p = 1)
#' Gamma.fromVAR(out_VAR)
#'


Gamma.fromVAR <- function(Phi, SigmaVAR) {

  # Check on Phi
  if(any(class(Phi) == "varest")){
    SigmaVAR <- cov(resid(Phi))
    Phi <- Acoef(Phi)[[1]]
    #
    if(length(Phi) == 1){
      q <- 1
    }else{
      q <- dim(Phi)[1]
    }
  } else if(any(class(Phi) == "ctsemFit")){
    B <- -1 * summary(Phi)$DRIFT
    Sigma <- summary(Phi)$DIFFUSION
    #
    VarEst <- VARparam(DeltaT, -B, Sigma)
    Phi <- VarEst$Phi
    SigmaVAR <- VarEst$SigmaVAR
    #
    if(length(Phi) == 1){
      q <- 1
    }else{
      q <- dim(Phi)[1]
    }
  } else{
    #
    if(length(Phi) == 1){
      q <- 1
    }else{
      Check_Phi(Phi)
      q <- dim(Phi)[1]
    }

    # Checks on SigmaVAR
    Check_SigmaVAR(SigmaVAR, q)

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
