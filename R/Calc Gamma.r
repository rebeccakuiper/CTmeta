
#' Calculate the stationary covariance matrix Gamma
#'
#' Calculates the stationary covariance matrix (Gamma) from the continuous-time (un)standardized lagged effects matrix (Drift) and residual covariance matrix (Sigma). There is also an interactive web application on my website: Standardizing and/or transforming lagged regression estimates (\url{https://www.uu.nl/staff/RMKuiper/Websites\%20\%2F\%20Shiny\%20apps}).
#'
#' @param Drift (Un)standardized lagged effects matrix of the first-order continuous-time (CT-VAR(1)) model.
#' It also takes a fitted object from the class "ctsemFit" (from the ctFit() function in the ctsem package); see example below. From such an object, the Drift and Sigma matrices are extracted.
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
#' #data <- myData
#' #if (!require("ctsemFit")) install.packages("ctsemFit")
#' #library(ctsemFit)
#' #out_CTM <- ctFit(...)
#' #Gamma.fromCTM(out_CTM)
#'

Gamma.fromCTM <- function(Drift, Sigma) {

  # Checks:
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
      ErrorMessage <- ("The drift matrix Drift should be input to the function.")
      #("Note that Phi(DeltaT) = expm(-B*DeltaT).")
      return(ErrorMessage)
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
      if(all(Re(eigen(B)$val) < 0)){
        cat("All (the real parts of) the eigenvalues of the drift matrix Drift are positive. Therefore. I assume the input was B=-A instead of A. I will use -B=A in the calculation.")
        B = -B
      }
      if(any(Re(eigen(B)$val) < 0)){
        #ErrorMessage <- ("The function stopped, since some of (the real parts of) the eigenvalues of the drift matrix Drift are positive.")
        #return(ErrorMessage)
        #stop(ErrorMessage)
        cat("If the function stopped, this is because some of (the real parts of) the eigenvalues of the drift matrix Drift are positive.")
      }
    }

    # Check on Sigma
    if(is.null(Sigma)){ # Sigma unknown
      ErrorMessage <- (paste0("The argument Sigma is NULL, but should be part of the input."))
      return(ErrorMessage)
      stop(ErrorMessage)
    }else if(!is.null(Sigma)){ # Sigma known
      # Check on Sigma
      Check_Sigma(Sigma, q)
    }
  }


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
