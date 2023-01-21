
#' Calculate the stationary covariance matrix Gamma
#'
#' Calculates the stationary covariance matrix (Gamma) from the continuous-time (un)standardized lagged effects matrix (Drift) and residual covariance matrix (Sigma). There is also an interactive web application on my website: Standardizing and/or transforming lagged regression estimates (\url{https://www.uu.nl/staff/RMKuiper/Websites\%20\%2F\%20Shiny\%20apps}).
#'
#' @param Drift (Un)standardized lagged effects matrix of the first-order continuous-time (CT-VAR(1)) model.
#' Can also take a fitted object from the class "ctsemFit" (from the ctFit() function in the ctsem package); see example below. If such an object is used, the Drift and Sigma matrices are extracted from its.
#' @param Sigma Residual covariance matrix of the first-order continuous-time (CT-VAR(1)) model, that is, the diffusion matrix.
#'
#' @return This function returns the stationary covariance matrix, that is, the contemporaneous covariance matrix of the data.
#' @export
#' @examples
#'
#' # library(CTmeta)
#'
#' ## Example 1 ##
#' #
#' Drift <- myDrift
#' #
#' q <- dim(Drift)[1]
#' Sigma <- diag(q) # for ease
#' #
#' Gamma.fromCTM(Drift, Sigma)
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
#'                                   latentNames=c('LeisureTime', 'Happiness'), TRAITVAR="auto")
#' out_CTM <- ctFit(dat=ctExample1, ctmodelobj=model)
#' ##################################################################################
#' #
#' Gamma.fromCTM(out_CTM)
#'

Gamma.fromCTM <- function(Drift, Sigma) {

  # Checks:
  if(anyNA(Drift)){
    stop("There are NA values in the Drift matrix.")
  }
  if(any(class(Drift) == "ctsemFit")){
    B <- -1 * summary(Drift)$DRIFT
    Sigma <- summary(Drift)$DIFFUSION
    #
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
      B <- -Drift
      #
      if(length(B) == 1){
        q <- 1
      }else{
        q <- dim(B)[1]
      }
    }
    
    # Check on compatible dimensions
    if (!is.null(try(Check_Sigma(Sigma, q), silent = TRUE)) &&
        grepl("The residual covariance matrix Sigma should, like Drift (or Phi), be a matrix of size q x q, with q = ",
              as.character(try(Check_Sigma(Sigma, q), silent = TRUE)),
              fixed = TRUE)) {
      stop("Sigma and Drift have different dimensions, but should both be square matrices with dimensions q x q.")
    } else {
      Check_Sigma(Sigma, q)
    }

    # Check on B
    if(length(B) > 1){
      Check_B(B)
      if(all(Re(eigen(B)$val) < 0)){
        cat("All (the real parts of) the eigenvalues of the drift matrix Drift are positive. Therefore, I assume the input was B=-A instead of A. I will use -B=A in the calculation.")
        B = -B
      }
      if(any(Re(eigen(B)$val) < 0)){
        #ErrorMessage <- ("The function stopped, since some of (the real parts of) the eigenvalues of the drift matrix Drift are positive.")
        #return(ErrorMessage)
        #stop(ErrorMessage)
        cat("Some of (the real parts of) the eigenvalues of the drift matrix Drift are positive.")
      }
    }

    # Check on Sigma
    if(is.null(Sigma)){ # Sigma unknown
      ErrorMessage <- (paste0("The argument Sigma is NULL, but should be part of the input."))
      stop(ErrorMessage)
    }else if(!is.null(Sigma)){ # Sigma known
      # Check on Sigma
      Check_Sigma(Sigma, q)
    }
  }


    # Check on B
    if(length(Drift) > 1){
      Check_B(B)
      if(all(Re(eigen(B)$val) > 0)){
        cat("All (the real parts of) the eigenvalues of the drift matrix Drift are positive. Therefore, I assume the input for Drift was B = -A instead of A. I will use Drift = -B = A.")
        cat("Note that Phi(DeltaT) = expm(-B*DeltaT) = expm(A*DeltaT) = expm(Drift*DeltaT).")
        Drift = -B
      }
      if(any(Re(eigen(B)$val) > 0)){
        #ErrorMessage <- ("The function stopped, since some of (the real parts of) the eigenvalues of the drift matrix Drift are positive.")
        #return(ErrorMessage)
        #stop(ErrorMessage)
        cat("Some of (the real parts of) the eigenvalues of the drift matrix Drift are positive.")
      }
    }


  #Gamma <- matrix((solve(kronecker(diag(q),B) + kronecker(B,diag(q))) %*% as.vector(Sigma)), ncol=q, nrow=q)
  #
  R <- kronecker(diag(q),B) + kronecker(B,diag(q))
  #
  if(q == 1){
    # Sigma = B %*% Gamma + Gamma %*% t(B) = 2 * beta * gamma, dus gamma = sigma / (2 * beta)
    Gamma_q <- Sigma / (2 * B)
  }else if(abs(det(R)) > 0.0001){
    Gamma_q <- matrix((solve(R) %*% as.vector(Sigma)), ncol=q, nrow=q)
  } else{
    Gamma_q = NULL
  }


  return(Gamma_q)
}
