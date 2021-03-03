
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
#' ## Example 1 ##
#' #
#' Phi <- myPhi[1:2,1:2]
#' #Phi <- matrix(c(0.50, 0.15, 0.25, 0.40), byrow=T, ncol = 2)
#' if (!require("expm")) install.packages("expm") # Use expm package for function logm()
#' library(expm)
#' DeltaT <- 1
#' Drift <- logm(Phi)/DeltaT
#' #
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

  if(any(class(Drift) == "ctsemFit")){
    B <- -1 * summary(Drift)$DRIFT
    Sigma <- summary(Drift)$DIFFUSION

  } else{

    # Check on Drift
    if(length(Drift) == 1){
      q <- 1
    }else{
      #
      if(is.null(dim(Drift))){
        if(!is.null(length(Drift))){
          print(paste("The argument Drift is not a matrix of size q times q."))
          stop()
        }else{
          print(paste("The argument Drift is not found: The continuous-time lagged effects matrix Drift is unknown, but should be part of the input."))
          stop()
        }
      }else{
        if(dim(Drift)[1] != dim(Drift)[2] | length(dim(Drift)) != 2){
          print(paste("The argument Drift is not a matrix of size q times q."))
          stop()
        }
        q <- dim(Drift)[1]
      }
    }

    # Checks on Sigma
    if(length(Sigma) != 1){
      if(dim(Sigma)[1] != dim(Sigma)[2]){ # Should be square
        print(paste("The residual covariance matrix Sigma should be a square matrix of size q times q, with q = ", q, "."))
        stop()
      }
      if(dim(Sigma)[1] != q){ # Should have same dimension as Drift
        print(paste("The residual covariance matrix Sigma should, like Drift, be a matrix of size q times q, with q = ", q, "."))
        stop()
      }
      if(length(dim(Sigma)) > 2){ # Should be matrix, not array
        print(paste("The residual covariance matrix Sigma should be an q times q matrix, with q = ", q, "."))
        stop()
      }
    }else if(q != 1){
      print(paste("The residual covariance matrix Sigma should, like Drift, be a scalar."))
      stop()
    }

  }

  # Drift = A = -Drift
  # B is drift matrix that is pos def, so Phi(DeltaT) = expm(-B*DeltaT)
  B <- -Drift

  #Gamma <- matrix((solve(kronecker(diag(q),B) + kronecker(B,diag(q))) %*% as.vector(Sigma)), ncol=q, nrow=q)
  #
  R <- kronecker(diag(q),B) + kronecker(B,diag(q))
  #
  if(abs(det(R)) > 0.0001){
    Gamma_q <- matrix((solve(R) %*% as.vector(Sigma)), ncol=q, nrow=q)
  } else{
    Gamma_q = NULL
  }

  return(Gamma_q)
}
