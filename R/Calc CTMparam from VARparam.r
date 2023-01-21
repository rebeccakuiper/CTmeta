#' Continuous-time estimates from discrete-time estimates
#'
#'
#' Transforms discrete-time parameters into their continuous-time (CT) equivalent ones and checks CT matrices (whether covariance matrices are positive definite, whether process is stable, etc.). The interactive web application 'Phi-and-Psi-Plots and Find DeltaT' can also do this. It is available here: \url{https://www.uu.nl/staff/RMKuiper/Websites\%20\%2F\%20Shiny\%20apps}.
#'
#' @param DeltaT Optional. The time interval used in the discrete-time matrices. By default, DeltaT = 1.
#' @param Phi (Un)standardized lagged effects matrix. If necessary, it can be standardized. The covariance matrix is determined for the standardized and vectorized Phi.
#' It also takes a fitted object from the class "varest" (from the VAR() function in vars package); see example below. The Phi and SigmaVAR matrices are extracted from this object.
#' @param SigmaVAR Optional (either SigmaVAR or Gamma). Residual covariance matrix of the first-order discrete-time vector autoregressive (DT-VAR(1)) model.
#' @param Gamma Optional (either SigmaVAR or Gamma). Stationary covariance matrix, i.e., the contemporaneous covariance matrix of the data.
#' Note that if Phi and SigmaVAR are known, Gamma can be calculated. Therefore, only one of SigmaVAR or Gamma is needed as input.
#'
#' @return The output renders the continuous-time matrices equivalent to the discrete-times ones. Some checks are also performed (namely, the stability of the process and uniqueness of the solution; see the function 'ChecksCTM' for checks on all matrices).
#' @importFrom expm expm
#' @importFrom expm logm
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
#' Phi <- myPhi[1:2, 1:2]
#' #
#' q <- dim(Phi)[1]
#' SigmaVAR <- diag(q) # for ease
#' #
#' Gamma <- Gamma.fromVAR(Phi, SigmaVAR)
#' ##################################################################################################
#'
#' CTMparam(DeltaT, Phi, SigmaVAR)
#' # or
#' CTMparam(DeltaT, Phi, Gamma = Gamma)
#'
#'
#' ## Example 2: input from fitted object of class "varest" ##
#' #
#' DeltaT <- 1
#' data <- myData
#' if (!require("vars")) install.packages("vars")
#' library(vars)
#' out_VAR <- VAR(data, p = 1)
#' CTMparam(DeltaT, out_VAR)
#'


CTMparam <- function(DeltaT = 1, Phi, SigmaVAR = NULL, Gamma = NULL) {
# DeltaT = time interval (although in VAR itself 1 is used!)
# q = Number of dimensions, that is, number of dependent variables
# Phi = matrix of size qxq with autoregressive and cross-lagged effects, matrix with lagged coefficients
# SigmaVAR = covariance matrix of size qxq, covariance matrix of innovations

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
  if(is.na(DeltaT)) {
    stop("Please specify DeltaT. This should be a scalar (i.e., one number or a vector with one element).")
  }
  if(DeltaT == 0){
    stop("Please specify a non-zero value for DeltaT.")
  }
  if(anyNA(SigmaVAR)){
    stop("There are NA values in SigmaVAR.")
  }
  if(!is.null(SigmaVAR) && !is.numeric(SigmaVAR)){
    stop("There are non-numerical values in SigmaVAR.")
  }
  if(!is.null(Gamma) && !is.numeric(Gamma)){
    stop("There are non-numerical values in Gamma.")
  }
  if(!is.numeric(Phi)){
    stop("There are non-numerical values in Phi.")
  }
  #
  B <- NULL
  Sigma <- NULL
  # Check on Phi
  if(any(class(Phi) == "varest")){
    SigmaVAR <- cov(resid(Phi))
    Phi <- Acoef(Phi)[[1]]
    # Calculate Gamma
    Gamma <- Gamma.fromVAR(Phi, SigmaVAR)
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
    Gamma <- Gamma.fromCTM(-B, Sigma)
    #
    #VarEst <- VARparam(DeltaT, -B, Sigma)
    #Phi <- VarEst$Phi
    #SigmaVAR <- VarEst$SigmaVAR
    #
    if(length(B) == 1){
      q <- 1
    }else{
      q <- dim(B)[1]
    }
  } else{
    #
    if(length(Phi) == 1){
      q <- 1
    }else{
      Check_Phi(Phi)
      q <- dim(Phi)[1]
    }
    #
    #Check on Phi
    Check_Phi(Phi)
    #
    # Check on SigmaVAR and Gamma
    if(!is.null(SigmaVAR) & !is.null(Gamma)){
      Gamma.from.VAR <- Gamma.fromVAR(Phi, SigmaVAR)
      if(any(Gamma.from.VAR != Gamma)) {
        warning("One of Phi, SigmaVAR, or Gamma, is likely mismatched from the other two.")
      }
    }
    if(is.null(SigmaVAR) & is.null(Gamma)){ # Both SigmaVAR and Gamma unknown
      #ErrorMessage <- (paste0("The arguments SigmaVAR and Gamma are not found: Both SigmaVAR and Gamma are unknown; either one (or both) should be part of the input. In case of first matrix, specify 'SigmaVAR = SigmaVAR'."))
      #return(ErrorMessage)
      #stop(ErrorMessage)
      #
      #print(paste0("Note: Both SigmaVAR and Gamma are unknown; Hence, the continuous-time residual covariance matrix Sigma and the the standardized parameter matrices cannot be calculated."))
    }else if(is.null(Gamma)){ # Gamma unknown, calculate Gamma from SigmaVAR and Phi

      # Check on SigmaVAR
      if (!is.null(try(Check_SigmaVAR(SigmaVAR, q), silent = TRUE)) &&
        grepl("The residual covariance matrix SigmaVAR should, like Phi, be a matrix with dimensions q x q, with q = ",
               as.character(try(Check_SigmaVAR(SigmaVAR, q), silent = TRUE)),
               fixed = TRUE)) {
        stop("SigmaVAR and Phi have different dimensions, but should both be square matrices with dimensions q x q.")
      } else {
        Check_SigmaVAR(SigmaVAR, q)
      }

      # Calculate Gamma
      Gamma <- Gamma.fromVAR(Phi, SigmaVAR)

    }else if(is.null(SigmaVAR)){ # SigmaVAR unknown, calculate SigmaVAR from Gamma and Phi

      # Checks on Gamma
      if (!is.null(try(Check_Gamma(Gamma, q), silent = TRUE)) &&
          grepl("The stationary covariance matrix Gamma should, like Phi, be a matrix of size q x q, with q = ",
                as.character(try(Check_Gamma(Gamma, q), silent = TRUE)),
                fixed = TRUE)) {
        stop("Gamma and Phi have different dimensions, but should both be square matrices with dimensions q x q.")
      } else {
        Check_Gamma(Gamma, q)
      }

      # Calculate SigmaVAR
      if(q == 1){
        SigmaVAR <- Gamma - Phi * Gamma * Phi
      }else{
        SigmaVAR <- Gamma - Phi %*% Gamma %*% t(Phi)
      }

    }else{ # Both SigmaVAR and Gamma known

      # Check on SigmaVAR
      Check_SigmaVAR(SigmaVAR, q)

      # Checks on Gamma
      Check_Gamma(Gamma, q)

    }
  }


if(is.null(B)){
  Decomp_ParamVAR <- eigen(Phi, symmetric=F)
  Eigen_ParamVAR <- Decomp_ParamVAR$val
  #
  #if (any(is.na(logm(Phi)) == TRUE)){ # In that case, there does not exist a solution A=-B for Phi # Note: logm(Phi) can (sometimes) exist when an EV < 0... so, I use the next check:
  if (any(is.na(log(Eigen_ParamVAR)))){ # In that case, there does not exist a solution A=-B for Phi
    ErrorMessage <- "'Phi' does not have a CTM-equivalent drift matrix. That is, there is no positive autocorrelation (as in the first-order continuous-time models), since one or more eigenvalues (have real parts which) are negative."
    stop(ErrorMessage)
  }
  # I at first wanted to filter out those where the real part of the eigenvalue of Phi is negative, but I don't think that is correct in all complex cases, so I decided to do the check above.
  #if (any(Re(Eigen_ParamVAR) <= 0)){ #
  #  ErrorMessage <- "At least one of the eigenvalues of 'Phi' is (has a real part)  smaller than or equal to 0; so, no positive autocorrelation (like CTM models)."
  #  final <- list(ErrorMessage = ErrorMessage)
  #  return(final)
  #  stop(ErrorMessage)
  #}


  # B = - log(Phi) / DeltaT = V log(-Eigenvalues / DeltaT) V^-1
  if(q == 1){
    B <- -log(Phi)/DeltaT
  }else{
    B <- -logm(Phi)/DeltaT
  }
}
#if(all(abs(Im(B)) < 0.0001) == TRUE){B <- Re(B)}

if(is.null(SigmaVAR) & is.null(Gamma)){
  Sigma = "Input for SigmaVAR(DeltaT) or Gamma is required to calculate the corresponding Sigma."
  Gamma <- "Input for SigmaVAR(DeltaT) or Gamma is required to calculate the corresponding Gamma."
  Gamma_s <- "Input for SigmaVAR(DeltaT) or Gamma is required to calculate the corresponding standGamma."
  Drift_s <- "Input for SigmaVAR(DeltaT) or Gamma is required to calculate the corresponding standDrift."
  Sigma_s <- "Input for SigmaVAR(DeltaT) or Gamma is required to calculate the corresponding standSigma."
}else{
  if(is.null(Sigma)){
    kronprod <- kronecker(Phi,Phi)
    if(q == 1){
      Sigma <- (-1/DeltaT) * matrix((solve(diag(q*q) - kronprod) %*% log(kronprod) %*% as.vector(SigmaVAR)), ncol=q, nrow=q)
    }else{
      Sigma <- (-1/DeltaT) * matrix((solve(diag(q*q) - kronprod) %*% logm(kronprod) %*% as.vector(SigmaVAR)), ncol=q, nrow=q)
      #if(all(abs(Im(Sigma)) < 0.0001) == TRUE){Sigma <- Re(Sigma)}
    }
  }

  Sxy <- sqrt(diag(diag(Gamma)))
  Gamma_s <- solve(Sxy) %*% Gamma %*% solve(Sxy)
  Drift_s <- solve(Sxy) %*% (-B) %*% Sxy
  Sigma_s <- solve(Sxy) %*% Sigma %*% solve(Sxy)

} # end 'else' belonging to Svar not NULL



Eigen_ParamCTM <- eigen(B)$val # = -log(Eigen_ParamVAR) / DeltaT

StableProcess <- FALSE
StableProcess_message <- "The process is NOT stable: it is exploding (since not all (real parts of) the eigenvalues of Drift are positive or, equivalently, not all absolute values of the eigenvalues of Phi are smaller than one)."
if(all(Re(Eigen_ParamCTM) > 0)){
  if(all(abs(Eigen_ParamVAR) < 1)){
    StableProcess <- TRUE
    StableProcess_message <- "The process is stable, it is returning to its equilibrium (since all (real parts of) the eigenvalues of Drift are positive or, equivalently, all absolute values of the eigenvalues of Phi are smaller than one)."
  }
}

UniqueSolution <- FALSE
UniqueSolution_message <- "The resulting drift matrix Drift is NOT unique, multiple solutions exist (since the eigenvalues of Drift (and thus also Phi) are complex, that is, have a non-zero imaginary part)."
if(all(Im(Eigen_ParamCTM) == 0)){
  UniqueSolution <- TRUE
  UniqueSolution_message <- "The resulting drift matrix Drift is unique (since the eigenvalues of Drift (and thus also Phi) are real, that is, have a zero imaginary part)."
}else if(all(abs(Im(Eigen_ParamCTM)) < base::pi/DeltaT)){
  UniqueSolution <- TRUE
  UniqueSolution_message <- "The resulting drift matrix Drift is unique for your sampling frequency, that is, for the DeltaT you used (since the imaginary part of the complex eigenvalues of Drift lie in (-pi/DeltaT, pi/DeltaT)); i.e., you measured frequently enough to know that this drift matrix is the solution, and not another (which do exist, as can be seen from the plots)."
}
#
##Multiple solutions - only if q=2!
#N = 1
#for(N in 1:3){
#  im <- complex(real = 0, imaginary = 1)
#  diagN <- diag(2)
#  diagN[1,1] <- -1
#  diagN <- N * diagN
#  Drift_N = (-B) + (2 * base::pi * im / DeltaT) * V_ParamVAR_1 %*% diagN %*% solve(V_ParamVAR_1)
#  Drift_N
#  print(Drift_N)
#}
#all(abs((2 * base::pi * im / DeltaT) * V_ParamVAR_1 %*% diagN %*% solve(V_ParamVAR_1)) == 0)




############################################################################################################

final <- list(Drift = -B,
              Sigma = Sigma,
              Gamma = Gamma,
              standDrift = Drift_s,
              standSigma = Sigma_s,
              standGamma = Gamma_s,
              UniqueSolutionDrift = UniqueSolution,
              UniqueSolutionDrift_message = UniqueSolution_message,
              StableProcess = StableProcess,
              StableProcess_message = StableProcess_message,
              eigenvalueDrift = -Eigen_ParamCTM,
              eigenvaluePhi_DeltaT = Eigen_ParamVAR)

return(final)



}
