#' Checks on CT estimates
#'
#' Checks on the continuous-time lagged-effects model matrices (like check whether matrices are positive-definite). The interactive web application 'Phi-and-Psi-Plots and Find DeltaT' also contains this functionality, you can find it on my website: \url{https://www.uu.nl/staff/RMKuiper/Websites\%20\%2F\%20Shiny\%20apps}.
#'
#' @param Drift Matrix of size q times q of (un)standardized continuous-time lagged effects, called drift matrix. Note that Phi(DeltaT) = expm(Drift*DeltaT). By default, input for Phi is used; only when Phi = NULL, Drift will be used.
#' It also takes a fitted object from the class "ctsemFit" (from the ctFit() function in the ctsem package); see example below. From such an object, the (standardized) Drift matrix is extracted.
#' @param Sigma Residual covariance matrix of the first-order continuous-time (CT-VAR(1)) model, that is, the diffusion matrix. By default, input for SigmaVAR is used; only when SigmaVAR = NULL, Sigma will be used.
#' @param Gamma Optional (either Sigma or Gamma). Stationary covariance matrix, that is, the contemporaneous covariance matrix of the data. By default, input for SigmaVAR is used; only when SigmaVAR = NULL, Gamma will be used.
#' Note that if Drift and Sigma are known, Gamma can be calculated; hence, only one out of SigmaVAR, Sigma, and Gamma is needed as input.
#'
#' @return The output renders the conclusions from the checks on the continuous-time lagged-effects model matrices.
#' @importFrom expm expm
#' @export
#' @examples
#'## Example 1 ##
#'
#' ##################################################################################################
#' # Input needed in examples below with q=2 variables.
#' # I will use the example matrices stored in the package:
#' Phi <- myPhi[1:2, 1:2]
#' if (!require("expm")) install.packages("expm") # Use expm package for function logm()
#' library(expm)
#' Drift <- logm(Phi)/DeltaT
#' #
#' Sigma <- diag(2) # for ease
#' #
#' Gamma <- Gamma.fromCTM(Drift, Sigma)
#' ##################################################################################################
#'
#' ChecksCTM(Drift, Sigma)
#' # or
#' ChecksCTM(Drift, Gamma = Gamma)
#'
#'
#' ## Example 2: input from fitted object of class "ctsemFit" ##
#' #
#' #data <- myData
#' #if (!require("ctsemFit")) install.packages("ctsemFit")
#' #library(ctsemFit)
#' #out_CTM <- ctFit(...)
#' #ChecksCTM(out_CTM)
#'


ChecksCTM <- function(B, Sigma = NULL, Gamma = NULL) {

  #  #######################################################################################################################
  #
  #  #if (!require("expm")) install.packages("expm")
  #  library(expm)
  #
  #  #######################################################################################################################

  # Checks:

  if(any(class(Drift) == "ctsemFit")){
    B <- -1 * summary(Drift)$DRIFT
    Sigma <- summary(Drift)$DIFFUSION
  }else{
    # Drift = A = -B
    # B is drift matrix that is pos def, so Phi(DeltaT) = expm(-B*DeltaT)
    B <- -Drift
    if(all(eigen(B)$val < 0)){
      #("All the eigenvalues of the drift matrix B are negative; therefore. I assume the input was A=-B instead of B. I will use -A=B in the calculation.")
      #("Note that Phi(DeltaT) = expm(-B*DeltaT).")
      ("All the eigenvalues of the drift matrix Drift are positive. Therefore. I assume the input was B=-A instead of A. I will use -B=A in the calculation.")
      B = -B
    }
    #
    # Check on B
    if(any(eigen(B)$val <= 0)){
      #("The function stopped, since some of the eigenvalues of the drift matrix B are negative or zero.")
      ("The function stopped, since some of the eigenvalues of the drift matrix Drift are positive or zero.")
      stop()
    }
    if(dim(B)[1] != dim(B)[2]){
      print(paste("The matrix (Drift or Phi) should be a square (q times q) matrix."))
      stop()
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
    print(paste("The arguments Sigma or Gamma are not found: one should be part of the input. Notably, in case of the last matrix, specify 'Gamma = Gamma'."))
    stop()
  }else if(is.null(Gamma)){ # Gamma unknown, calculate Gamma from Drift & Sigma

    if(!is.null(Sigma)){ # Sigma known, calculate Gamma from Drift & Sigma

      # Checks on Sigma
      if(length(Sigma) != 1){
        if(dim(Sigma)[1] != dim(Sigma)[2]){
          print(paste("The residual covariance matrix Sigma should be a square matrix of size q times q, with q = ", q, "."))
          stop()
        }
        if(dim(Sigma)[1] != q){
          print(paste("The residual covariance matrix Sigma should, like Phi (or Drift), be a matrix of size q times q, with q = ", q, "."))
          stop()
        }
        if(length(dim(Sigma)) > 2){
          print(paste("The residual covariance matrix Sigma should be an q times q matrix, with q = ", q, "."))
          stop()
        }
      }else if(q != 1){
        print(paste("The residual covariance matrix Sigma should, like Phi (or Drift), be a scalar."))
        stop()
      }

      # Calculate Gamma
      Gamma <- Gamma.fromCTM(Drift, Sigma)

    }

  }else if(!is.null(Gamma)){ # Gamma known, only check on Gamma needed and calculate Sigma

    # Checks on Gamma
    if(length(Gamma) != 1){
      if(dim(Gamma)[1] != dim(Gamma)[2]){
        print(paste("The stationary covariance matrix Gamma should be a square matrix of size q times q, with q = ", q, "."))
        stop()
      }
      if(dim(Gamma)[1] != q){
        print(paste("The stationary covariance matrix Gamma should, like Phi (or Drift), be a matrix of size q times q, with q = ", q, "."))
        stop()
      }
      if(length(dim(Gamma)) > 2){
        print(paste("The stationary covariance matrix Gamma should be an q times q matrix, with q = ", q, "."))
        stop()
      }
    }else if(q != 1){
      print(paste("The stationary covariance matrix Gamma should, like Phi (or Drift), be a scalar."))
      stop()
    }

    # Calculate Sigma
    Sigma <- B %*% Gamma + t(B %*% Gamma)

  }


ChecksAreFine <- TRUE
error <- list()
append(error, "Checks on the matrices Drift, Sigma, and Gamma are inspected, the following problems exist:")


#BplusTransB = B + t(B) # Namely, B niet per se symmetric, hence positive eigenvalues is not a suficient condition for positive (semi-)definiteness
##Decomp_BB <- eigen(BplusTransB, sym=TRUE)
##diagD_BB <- Decomp_BB$values
#diagD_BB <- eigen(BplusTransB, symmetric=FALSE)$val
#if (any(diagD_BB <= 0)){
#  ("'B' is not positive definite")
#  ChecksAreFine = FALSE
#}
# Eigenvalues positive (or at least real part)
diagD_B <- eigen(B, symmetric=FALSE)$val
if (any(Re(diagD_B) <= 0)){
  cat("The (real parts of the) eigenvalues of the drift matrix 'Drift' are positive. \n Hence, the process is not stable. \n")
  ChecksAreFine <- FALSE
  error <- append(error, "Not all (real parts of the) eigenvalues of the drift matrix 'Drift' are negative. \n Hence, the process is not stable.")
}
#
# CHECK on complex
Eigen_ParamCTM <- eigen(B)$val
if (is.complex(Eigen_ParamCTM) == TRUE){
  if(all(Im(Eigen_ParamCTM) > -base::pi/1)){
    if(all(Im(Eigen_ParamCTM) < base::pi/1)){
      cat("The drift matrix 'Drift' has complex eigenvalues and -vectors, but is a unique solution for Phi(1). \n (nl imaginary part of the complex eigenvalues of drift matrix lie in (-pi/1, pi/1)). \n")
      ChecksAreFine <- FALSE
      error <- append(error, "Note: The drift matrix 'Drift' has complex eigenvalues and -vectors, but is a unique solution for Phi(1). \n (nl imaginary part of the complex eigenvalues of drift matrix lie in (-pi/1, pi/1)).")
    }
  } else{
    cat("The drift matrix 'Drift' has complex eigenvalues and -vectors and is not a unique solution for this Phi(1). \n (nl imaginary part of the complex eigenvalues of drift matrix do NOT lie in (-pi/1, pi/1)). \n")
    ChecksAreFine <- FALSE
    error <- append(error, "Note: The drift matrix 'Drift'  has complex eigenvalues and -vectors and is not a unique solution for this Phi(1). \n (nl imaginary part of the complex eigenvalues of drift matrix do NOT lie in (-pi/1, pi/1)).")
  }
}


# CHECK pos def
#Decomp_Sigma <- eigen(Sigma, sym=TRUE)
#diagDSigma <- Decomp_Sigma$values
diagDSigma <- eigen(Sigma, symmetric=T)$val
if (any(diagDSigma <= 0)){
  cat("'Sigma' is not positive definite. \n")
  ChecksAreFine <- FALSE
  error <- append(error, "'Sigma' is not positive definite.")
}


# Determine corresponding Gamma (if q < 4)
# Sigma <- B %*% Gamma + Gamma %*% t(B)
# Determine via other function!
if(q == 1){
  # Sigma = B %*% Gamma + Gamma %*% t(B) = 2 * beta * gamma, dus gamma = sigma / (2 * beta)
  Gamma <- Sigma / (2 * B)
}else if(q > 1){
Gamma <- Gamma.fromCTM(B, Sigma)
}
#
#if(q > 3){
#  cat("'Gamma', corresponding to 'B' and 'SigmaCTM', cannot be calculated for q > 3 (in this sofware). \n")
#  cat("Thus, positive definiteness of 'Gamma' is not checked (in generation of CTM data). \n")
#}else{
  # CHECK pos def
  #Decomp_Gamma <- eigen(Gamma, sym=TRUE)
  #diagDGamma <- Decomp_Gamma$values
  diagDGamma <- eigen(Gamma, symmetric=T)$val
  if (any(diagDGamma <= 0)){
    cat("The stationary covariance matrix ('Gamma'), corresponding to 'Drift' and 'Sigma', is not positive definite. \n")
    ChecksAreFine <- FALSE
    error <- append(error, "'Gamma' (stationary covariance matrix) is not positive definite.")
  }
#}


  if(ChecksAreFine){
    error <- "All matrices (Drift, Sigma, and Gamma) are fine."
  }


############################################################################################################

final <- list(ChecksAreFine = ChecksAreFine, Gamma = Gamma, error = error)
return(final)

}



