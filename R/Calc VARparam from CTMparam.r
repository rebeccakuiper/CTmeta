#' Discrete-time estimates from continuous-time estimates
#'
#' The discrete-time lagged-effects model matrices corresponding to the continuous-time ones. The interactive web application 'Phi-and-Psi-Plots and Find DeltaT' also contains this functionality, you can find it on my website: \url{https://www.uu.nl/staff/RMKuiper/Websites\%20\%2F\%20Shiny\%20apps}.
#'
#' @param DeltaT Optional. The time interval used. By default, DeltaT = 1.
#' @param Drift Matrix of size q times q of (un)standardized continuous-time lagged effects, called drift matrix. Note that Phi(DeltaT) = expm(Drift*DeltaT). By default, input for Phi is used; only when Phi = NULL, Drift will be used.
#' It also takes a fitted object from the class "ctsemFit" (from the ctFit() function in the ctsem package); see example below. From such an object, the (standardized) Drift and Sigma matrices are extracted.
#' @param Sigma Residual covariance matrix of the first-order continuous-time (CT-VAR(1)) model, that is, the diffusion matrix. By default, input for Sigma is used; only when Sigma = NULL, Gamma will be used.
#' @param Gamma Optional (either Sigma or Gamma). Stationary covariance matrix, that is, the contemporaneous covariance matrix of the data. By default, input for Sigma is used; only when Sigma = NULL, Gamma will be used.
#' Note that if Drift and Sigma are known, Gamma can be calculated; hence, either Sigma or Gamma is needed as input.
#'
#' @return The output renders the discrete-times equivalent matrices of the continuous-time ones.
#' @importFrom expm expm
#' @export
#' @examples
#'
#' # library(CTmeta)
#'
#' ## Example 1 ##
#'
#' ##################################################################################################
#' # Input needed in examples below with q=2 variables.
#' # I will use the example matrix stored in the package:
#' DeltaT <- 1
#' Drift <- myDrift
#' q <- dim(Drift)[1]
#' Sigma <- diag(q) # for ease
#' #
#' Gamma <- Gamma.fromCTM(Drift, Sigma)
#' ##################################################################################################
#'
#' DeltaT <- 1
#' VARparam(DeltaT, Drift, Sigma)
#' # or
#' VARparam(DeltaT, Drift, Gamma = Gamma)
#'
#'
#' ## Example 2: input from fitted object of class "ctsemFit" ##
#' #
#' #data <- myData
#' #if (!require("ctsemFit")) install.packages("ctsemFit")
#' #library(ctsemFit)
#' #out_CTM <- ctFit(...)
#' #
#' #DeltaT <- 1
#' #VARparam(DeltaT, out_CTM)
#'


VARparam <- function(DeltaT = 1, Drift, Sigma = NULL, Gamma = NULL) {
# DeltaT = time interval
# q = Number of dimensions, that is, number of dependent variables

  #DeltaT <- 1
  #kronsum <- kronecker(diag(q),B) + kronecker(B,diag(q))
  #SigmaVAR <- matrix((solve(kronsum) %*% (diag(q*q) - expm(-kronsum * DeltaT)) %*% as.vector(Sigma)), ncol=q, nrow=q)



  #  #######################################################################################################################
  #
  #  #if (!require("expm")) install.packages("expm")
  #  library(expm)
  #
  #  #######################################################################################################################

  # Checks:
  if(length(DeltaT) != 1){
    ErrorMessage <- (paste0("The argument DeltaT should be a scalar (i.e., one number or a vector with one element)."))
    stop(ErrorMessage)
  }
  #
  #
  if(any(class(Drift) == "ctsemFit")){
    B <- -1 * summary(Drift)$DRIFT
    Sigma <- summary(Drift)$DIFFUSION
  }else{
    # Drift = A = -B
    # B is drift matrix that is pos def, so Phi(DeltaT) = expm(-B*DeltaT)
    B <- -Drift
    #
    # Check on B
    if(length(Drift) > 1){
      Check_B(B)
      if(all(Re(eigen(Drift)$val) > 0)){
        cat("All (the real parts of) the eigenvalues of the drift matrix Drift are positive. Therefore. I assume the input for Drift was B = -A instead of A. I will use Drift = -B = A.")
        cat("Note that Phi(DeltaT) = expm(-B*DeltaT) = expm(A*DeltaT) = expm(Drift*DeltaT).")
        Drift = -Drift
      }
      if(any(Re(eigen(Drift)$val) > 0)){
        #ErrorMessage <- ("The function stopped, since some of (the real parts of) the eigenvalues of the drift matrix Drift are positive.")
        #return(ErrorMessage)
        #stop(ErrorMessage)
        cat("If the function stopped, this is because some of (the real parts of) the eigenvalues of the drift matrix Drift are positive.")
      }
    }
  }
  #
  if(length(B) == 1){
    q <- 1
  }else{
    q <- dim(B)[1]
  }
  #
  # Check on Sigma, and Gamma
  if(is.null(Gamma) & is.null(Sigma)){ # Both unknown
    ErrorMessage <- (paste0("The arguments Sigma and Gamma are NULL. One of them must be part of the input. If specifying Gamma, do not forget the name the argument (i.e., specify 'Gamma = yourGamma'."))
    stop(ErrorMessage)
  }else if(is.null(Gamma)){ # Gamma unknown, calculate Gamma from Drift & Sigma

    if(!is.null(Sigma)){ # Sigma known, calculate Gamma from Drift=-B & Sigma

      # Checks on Sigma
      Check_Sigma(Sigma, q)

      # Calculate Gamma
      Gamma <- Gamma.fromCTM(-B, Sigma)

    }

  }else if(!is.null(Gamma)){ # Gamma known, check on Gamma needed and calculate Sigma

    # Checks on Gamma
    Check_Gamma(Gamma, q)

    # Calculate Sigma
    if(q == 1){
      Sigma <- B * Gamma + (B * Gamma)
    }else{
      Sigma <- B %*% Gamma + t(B %*% Gamma)
    }

  }else{

    # Checks on Sigma
    Check_Sigma(Sigma, q)

    # Checks on Gamma
    Check_Gamma(Gamma, q)

  }



# Phi = exp{-B * DeltaT} = V exp{-Eigenvalues * DeltaT} V^-1 = V D_Phi V^-1
if(q == 1){
  ParamVAR <- exp(-B*DeltaT)

  kronsum <- kronecker(diag(q),B) + kronecker(B,diag(q))
  Sigma_VAR <- matrix((solve(kronsum) %*% (diag(q*q) - exp(-kronsum * DeltaT)) %*% as.vector(Sigma)), ncol=q, nrow=q)
  #if(all(abs(Im(Sigma_VAR)) < 0.0001) == TRUE){Sigma_VAR <- Re(Sigma_VAR)}

}else{
  ParamVAR <- expm(-B*DeltaT)

  kronsum <- kronecker(diag(q),B) + kronecker(B,diag(q))
  Sigma_VAR <- matrix((solve(kronsum) %*% (diag(q*q) - expm(-kronsum * DeltaT)) %*% as.vector(Sigma)), ncol=q, nrow=q)
  #if(all(abs(Im(Sigma_VAR)) < 0.0001) == TRUE){Sigma_VAR <- Re(Sigma_VAR)}

}
#if(all(abs(Im(ParamVAR)) < 0.0001) == TRUE){ParamVAR <- Re(ParamVAR)}



# Stable process info
Decomp_ParamVAR <- eigen(Phi, symmetric=F)
Eigen_ParamVAR <- Decomp_ParamVAR$val
Eigen_ParamCTM <- diag(eigen(B)$val) # = -log(Eigen_ParamVAR) / DeltaT
StableProcess = FALSE
if(all(Re(Eigen_ParamCTM) > 0)){
  if(all(abs(Eigen_ParamVAR) < 1)){
    StableProcess = TRUE
  }
}

# Standardized matrices
Sxy <- sqrt(diag(diag(Gamma)))
Gamma_s <- solve(Sxy) %*% Gamma %*% solve(Sxy)
ParamVAR_s <- solve(Sxy) %*% ParamVAR %*% Sxy
Sigma_VAR_s <- solve(Sxy) %*% Sigma_VAR %*% solve(Sxy)



############################################################################################################

#final <- list(Phi = ParamVAR, SigmaVAR = Sigma_VAR, Gamma = Gamma)
final <- list(eigenvalueDrift = -Eigen_ParamCTM,
              eigenvaluePhi_DeltaT = Eigen_ParamVAR,
              StableProcess = StableProcess,
              Phi = ParamVAR,
              SigmaVAR = Sigma_VAR,
              Gamma = Gamma,
              standPhi = ParamVAR_s,
              standSigmaVAR = Sigma_VAR_s,
              standGamma = Gamma_s)
return(final)

}

