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
#'
#' # library(CTmeta)
#'
#' ## Example 1 ##
#'
#' ##################################################################################################
#' # Input needed in examples below with q=2 variables.
#' # I will use the example matrix stored in the package:
#' Drift <- myDrift
#' #
#' q <- dim(Drift)[1]
#' Sigma <- diag(q) # for ease
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


ChecksCTM <- function(Drift, Sigma = NULL, Gamma = NULL) {

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
    #
    if(length(B) > 1){
      eigenvalDrift <- eigen(B)$val
    }else{
      eigenvalDrift <- B
    }
    if(length(B) == 1){
      q <- 1
    }else{
      q <- dim(B)[1]
    }
  }else{
    # Drift = A = -B
    # B is drift matrix that is pos def, so Phi(DeltaT) = expm(-B*DeltaT)
    if(is.null(Drift)){
      ErrorMessage <- ("The drift matrix Drift should be input to the function.")
      #("Note that Phi(DeltaT) = expm(-B*DeltaT).")
      stop(ErrorMessage)
    }else{ # !is.null(Drift)
      B <- -Drift
      #
      if(length(B) == 1){
        q <- 1
      }else{
        q <- dim(B)[1]
      }
    }

    # Check on B
    if(length(B) > 1){
      Check_B(B)
      eigenvalDrift <- eigen(B)$val
    }else{
      eigenvalDrift <- B
    }
    #
    if(all(Re(eigen(B)$val) < 0)){
      cat("All (the real parts of) the eigenvalues of the drift matrix Drift are positive. Therefore. I assume the input for Drift was B = -A instead of A. I will use Drift = -B = A.")
      cat("Note that Phi(DeltaT) = expm(-B*DeltaT) = expm(A*DeltaT) = expm(Drift*DeltaT).")
      B = -B
    }
    #if(any(Re(eigen(B)$val) < 0)){
    #  #ErrorMessage <- ("The function stopped, since some of (the real parts of) the eigenvalues of the drift matrix Drift are positive.")
    #  #stop(ErrorMessage)
    #  cat("If the function stopped, this is because some of (the real parts of) the eigenvalues of the drift matrix Drift are positive.")
    #}

    # Check on Sigma and Gamma - need Sigma
    if(is.null(Sigma) & is.null(Gamma)){ # Both Sigma and Gamma unknown
      ErrorMessage <- (paste0("The arguments Sigma and Gamma are not found: Both Sigma and Gamma are unknown; one should be part of the input."))
      stop(ErrorMessage)
    }
    if(!is.null(Sigma)){ # Sigma known
      # Check on Sigma
      Check_Sigma(Sigma, q)
    }else{ # Sigma unknown and Gamma known, calculate Sigma
      # Calculate Sigma
      if(q != 1){
        Sigma <- B %*% Gamma + t(B %*% Gamma)
      }else{
        Sigma <- B * Gamma + (B * Gamma)
      }
    }
  }
  #
  # Check on Gamma
  if(!is.null(Gamma)){ # Gamma known
    # Checks on Gamma
    Check_Gamma(Gamma, q)
  }else{ # Gamma unknown, calculate Gamma from Sigma and -B=Drift
    # Calculate Gamma
    # Sigma <- B %*% Gamma + Gamma %*% t(B)
    Gamma <- Gamma.fromCTM(-B, Sigma)
  }


ChecksAreFine <- TRUE
error <- list()
append(error, "Checks on the matrices Drift, Sigma, and Gamma are inspected, the following problems exist:")


# Eigenvalues positive (or at least real part)
if (any(Re(eigenvalDrift) <= 0)){
  cat("The (real parts of the) eigenvalues of the drift matrix 'Drift' are positive. \n Hence, the process is not stable. \n")
  ChecksAreFine <- FALSE
  error <- append(error, "Not all (real parts of the) eigenvalues of the drift matrix 'Drift' are negative. \n Hence, the process is not stable.")
}
#
# CHECK on complex
if (is.complex(eigenvalDrift) == TRUE){
  if(all(Im(eigenvalDrift) > -base::pi/1)){
    if(all(Im(eigenvalDrift) < base::pi/1)){
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
#eigenvalSigma <- Decomp_Sigma$values
eigenvalSigma <- eigen(Sigma, symmetric=T)$val
if (any(eigenvalSigma <= 0)){
  cat("'Sigma' is not positive definite. \n")
  ChecksAreFine <- FALSE
  error <- append(error, "'Sigma' is not positive definite.")
}

# CHECK pos def
#Decomp_Gamma <- eigen(Gamma, sym=TRUE)
#eigenvalGamma <- Decomp_Gamma$values
eigenvalGamma <- eigen(Gamma, symmetric=T)$val
if (any(eigenvalGamma <= 0)){
  cat("The stationary covariance matrix ('Gamma'), corresponding to 'Drift' and 'Sigma', is not positive definite. \n")
  ChecksAreFine <- FALSE
  error <- append(error, "'Gamma' (stationary covariance matrix) is not positive definite.")
}


  if(ChecksAreFine){
    error <- "No error: All matrices (Drift, Sigma, and Gamma) are fine."
  }


############################################################################################################

final <- list(ChecksAreFine = ChecksAreFine, error = error,
              Drift = -B, Sigma = Sigma, Gamma = Gamma,
              EigenVal_Drift = eigenvalDrift, EigenVal_Sigma = eigenvalSigma, EigenVal_Gamma = eigenvalGamma)
return(final)

}



