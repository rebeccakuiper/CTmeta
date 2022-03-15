#' Checks on CT estimates
#'
#' Various checks on the continuous-time lagged-effects model matrices (e.g. check whether matrices are positive-definite). The interactive web application 'Phi-and-Psi-Plots and Find DeltaT': \url{https://www.uu.nl/staff/RMKuiper/Websites\%20\%2F\%20Shiny\%20apps}.
#'
#' @param Drift Matrix of size q x q of (un)standardized continuous-time lagged effects, called the drift matrix. Can also be a fitted object from the class "ctsemFit" (from the ctFit() function in the ctsem package); see example below. The (standardized) Drift matrix is extracted from such an object.
#' @param Sigma Optional (either Sigma or Gamma). Residual covariance matrix (of size q x q) of the first-order continuous-time (CT-VAR(1)) model, i.e. the diffusion matrix.
#' Note that if Drift and Sigma are known, Gamma can be calculated: hence, only one out of Sigma and Gamma is needed as input.
#' @param Gamma Optional (either Sigma or Gamma). Stationary covariance matrix (of size q x q), that is, the contemporaneous covariance matrix of the data.
#' Note that if Drift and Sigma are known, Gamma can be calculated: hence, only one out of Sigma and Gamma is needed as input. If all three are provided, their compatibility is checked.
#'
#' @return The output provides conclusions from the checks on the continuous-time lagged-effects model matrices: whether any problems were found, the exact errors that were found, the Drift, Sigma and Gamma matrices and their eigenvalues.
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
#' # We use the example matrix stored in the package:
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
#' if (!require("ctsem")) install.packages("ctsem")
#' library(ctsem)
#' library(ctsemOMX)
#' #
#' ############ adapted from https://rdrr.io/cran/ctsemOMX/man/ctFit.html ############
#' data(ctExample1)
#' model <- ctModel(n.manifest=2, n.latent=2, Tpoints=6, LAMBDA=diag(2),
#'                  manifestNames=c('LeisureTime', 'Happiness'),
#'                  latentNames=c('LeisureTime', 'Happiness'), TRAITVAR="auto")
#' out_Drift <- ctFit(dat=ctExample1, ctmodelobj=model)
#' ##################################################################################
#' #
#' ChecksCTM(out_Drift, Sigma)
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
      ErrorMessage <- ("The drift matrix Drift should be part of the input.")
      #("Note that Phi(DeltaT) = expm(-B*DeltaT).")
      stop(ErrorMessage)
    }else{ # !is.null(Drift)
      if(anyNA(Drift)) {
        stop("There are missing values in Drift.")
      }
      if(!is.numeric(Drift)) {
        stop("There are non-numerical values in Drift.")
      }
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
      cat("All (the real parts of) the eigenvalues of the drift matrix Drift are positive. Therefore, I assume the input for Drift was B = -A instead of A. I will use Drift = -B = A.")
      cat("Note that Phi(DeltaT) = expm(-B*DeltaT) = expm(A*DeltaT) = expm(Drift*DeltaT).")
      B = -B
    }
    #if(any(Re(eigen(B)$val) < 0)){
    #  #ErrorMessage <- ("The function stopped, since some of (the real parts of) the eigenvalues of the drift matrix Drift are positive.")
    #  #return(ErrorMessage)
    #  #stop(ErrorMessage)
    #  cat("If the function stopped, this is because some of (the real parts of) the eigenvalues of the drift matrix Drift are positive.")
    #}

    # Check on Sigma and Gamma - need Sigma
    if(is.null(Sigma) & is.null(Gamma)){ # Both Sigma and Gamma unknown
      ErrorMessage <- (paste0("Both Sigma and Gamma are NULL. One of them must be part of the input."))
      stop(ErrorMessage)
    }
    if(!is.null(Sigma)){ # Sigma known
      if(anyNA(Sigma)) {
        stop("There are NA values in Sigma.")
      }
      if(!is.numeric(Sigma)) {
        stop("There are non-numerical values in Sigma.")
      }
      # Check on Sigma
      Check_Sigma(Sigma, q)
    }else{ # Sigma unknown and Gamma known, calculate Sigma
      Check_Gamma(Gamma, q)
      if(anyNA(Gamma)) {
        stop("There are missing values in Gamma.")
      }
      if(!is.numeric(Gamma)) {
        stop("There are non-numerical values in Gamma.")
      }
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
error <- append(error, "Checks on the matrices Drift, Sigma, and Gamma are inspected, the following problems exist:")

if(all(Re(eigen(B)$val) < 0)){
  error <- append(error, "All (the real parts of) the eigenvalues of the drift matrix Drift are positive. Therefore, B = -A was used for Drift instead of A. I will use Drift = -B = A.")
  error <- append(error, "Note that Phi(DeltaT) = expm(-B*DeltaT) = expm(A*DeltaT) = expm(Drift*DeltaT).")
  B = -B
}

# Check that the three matrices are compatible if all three are provided
if(!is.null(Sigma) & !is.null(Gamma)) {
  if(q != 1){
    Sigma.from.Gamma <- B %*% Gamma + t(B %*% Gamma)
  }else{
    Sigma.from.Gamma <- B * Gamma + (B * Gamma)
  }
  Gamma.from.Sigma <- Gamma.fromCTM(-B, Sigma)
  if(any(Gamma != Gamma.from.Sigma) | any(Sigma != Sigma.from.Gamma)) {
    cat("The Gamma and/or Sigma matrices obtained from the Drift and Sigma and/or Gamma matrices do not match the provided Sigma and Gamma.\n One of the given matrices is likely wrong.\n")
    ChecksAreFine <- FALSE
    error <- append(error, "The provided Sigma and/or Gamma do not match the computed Sigma and/or Gamma.")
  }
}


# Eigenvalues positive (or at least real part)
if (any(Re(eigenvalDrift) <= 0)){
  cat("The (real parts of the) eigenvalues of the drift matrix 'Drift' are positive. \n Hence, the process is not stable. \n")
  ChecksAreFine <- FALSE
  error <- append(error, "Not all (real parts of the) eigenvalues of the drift matrix 'Drift' are negative. Hence, the process is not stable.")
}
#
# CHECK on complex
if (is.complex(eigenvalDrift) == TRUE){
  if(all(Im(eigenvalDrift) > -base::pi/1)){
    if(all(Im(eigenvalDrift) < base::pi/1)){
      cat("The drift matrix 'Drift' has complex eigenvalues and eigenvectors, but is a unique solution for Phi(1). \n (namely imaginary part of the complex eigenvalues of drift matrix lie in (-pi/1, pi/1)). \n")
      ChecksAreFine <- FALSE
      error <- append(error, "Note: The drift matrix 'Drift' has complex eigenvalues and eigenvectors, but is a unique solution for Phi(1). (namely imaginary part of the complex eigenvalues of drift matrix lie in (-pi/1, pi/1)).")
    }
  } else{
    cat("The drift matrix 'Drift' has complex eigenvalues and -vectors and is not a unique solution for this Phi(1). \n (namely imaginary part of the complex eigenvalues of drift matrix do NOT lie in (-pi/1, pi/1)). \n")
    ChecksAreFine <- FALSE
    error <- append(error, "Note: The drift matrix 'Drift'  has complex eigenvalues and -vectors and is not a unique solution for this Phi(1). \n (namely imaginary part of the complex eigenvalues of drift matrix do NOT lie in (-pi/1, pi/1)).")
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



