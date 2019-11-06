
#' calc.CovMxStandPhi
#'
#' This function calculates the (vectorized) standardized lagged effects matrix, their covariance matrix, and corresponding elliptical 95\% confidence interval (CI). There is also an interactive web application on my website: Standardizing and/or transforming lagged regression estimates (\url{https://www.uu.nl/staff/RMKuiper/Websites\%20\%2F\%20Shiny\%20apps}).
#'
#' @param N Number of persons (panel data) or number of measurement occasions - 1 (time series data). This is used in determining the covariance matrix of the vectorized standardized lagged effects.
#' @param Phi (Un)standardized lagged effects matrix. If necessary, it is standardized and for the standardized and vectorized Phi the covariance matrix is determined.
#' @param Gamma Stationary covariance matrix, that is, the contemporaneous covariance matrix of the data.
#' Note that if Phi and SigmaVAR are known, Gamma can be calculated; hence, only SigmaVAR or Gamma is needed as input (if only Gamma, then use 'Gamma = Gamma' or set SigmaVAR to NULL, see examples below).
#' @param SigmaVAR Residual covariance matrix of the first-order discrete-time vector autoregressive (DT-VAR(1)) model.
#' Note that if Phi and SigmaVAR are known, Gamma can be calculated; hence, only SigmaVAR or Gamma is needed as input (if only Gamma, then use 'Gamma = Gamma' or set SigmaVAR to NULL, see examples below).
#' @param alpha The alpha level in determining the (1-alpha)*100\% CI. By default, alpha is set to 0.05, resulting in a 95\% CI.
#'
#' @return This function returns the vectorized standardized lagged effects, their covariance matrix, and the corresponding elliptical/multivariate 95\% CI.
#' @export
#' @examples
#'
#' # Input for examples below
#' N <- 643
#' Phi <- myPhi[1:2,1:2]
#' #Phi <- matrix(c(0.25, 0.10,
#' #                0.20, 0.36), byrow=T, ncol = 2)
#' SigmaVAR <- diag(2) # for ease
#' # Calculate the Gamma corresponding to Phi and SigmaVAR - used in the second example
#' Gamma <- calc.Gamma.fromVAR(Phi, SigmaVAR) # ?calc.Gamma.fromVAR
#'
#' #Example where only SigmaVAR is known and not Gamma
#' calc.CovMxStandPhi(N, Phi, NULL, SigmaVAR)
#' # or
#' calc.CovMxStandPhi(N, Phi, SigmaVAR = SigmaVAR)
#'
#' #Example where only Gamma is known and not SigmaVAR
#' calc.CovMxStandPhi(N, Phi, Gamma)
#' # or
#' calc.CovMxStandPhi(N, Phi, Gamma, NULL)
#'


calc.CovMxStandPhi <- function(N, Phi, Gamma = NULL, SigmaVAR = NULL, alpha=0.05) {

  # Checks:
  if(length(N) != 1){
    print(paste("The argument N should be a scalar, that is, one number, that is, a vector with one element."))
    stop()
  }
  #
  # Check on Phi
  if(length(Phi) != 1){
    if(is.null(dim(Phi))){
      if(!is.null(length(Phi))){
        print(paste("The argument Phi is not a matrix of size q times q."))
        stop()
      }else{
        print(paste("The argument Phi is not found: The lagged effects matrix Phi is unknown, but should be part of the input."))
        stop()
      }
    }
    q <- dim(Phi)[1]
  } else{
    q <- 1
  }

  # Check on SigmaVAR and Gamma
  if(is.null(SigmaVAR) & is.null(Gamma)){ # Both SigmaVAR and Gamma unknown
    print(paste("The arguments SigmaVAR and Gamma are not found: Both SigmaVAR and Gamma are unknown; either one (or both) should be part of the input. In case of first matrix, specify 'SigmaVAR = SigmaVAR'."))
    stop()
  }else if(is.null(Gamma)){ # Gamma unknown, calculate Gamma from SigmaVAR and Phi

    # Checks on SigmaVAR
    if(length(SigmaVAR) != 1){
      if(dim(SigmaVAR)[1] != dim(SigmaVAR)[2]){
        print(paste("The residual covariance matrix SigmaVAR should be a square matrix of size q times q, with q = ", q, "."))
        stop()
      }
      if(dim(SigmaVAR)[1] != q){
        print(paste("The residual covariance matrix SigmaVAR should, like Phi, be a matrix of size q times q, with q = ", q, "."))
        stop()
      }
      if(length(dim(SigmaVAR)) > 2){
        print(paste("The residual covariance matrix SigmaVAR should be an q times q matrix, with q = ", q, "."))
        stop()
      }
    }else if(q != 1){
      print(paste("The residual covariance matrix SigmaVAR should, like Phi, be a scalar."))
      stop()
    }

    # Calculate Gamma
    Gamma <- calc.Gamma.fromVAR(Phi, SigmaVAR)

  }else if(is.null(SigmaVAR)){ # SigmaVAR unknown, calculate SigmaVAR from Gamma and Phi


    # Checks on Gamma
    if(length(Gamma) != 1){
      if(dim(Gamma)[1] != dim(Gamma)[2]){
        print(paste("The stationary covariance matrix Gamma should be a square matrix of size q times q, with q = ", q, "."))
        stop()
      }
      if(dim(Gamma)[1] != q){
        print(paste("The stationary covariance matrix Gamma should, like Phi, be a matrix of size q times q, with q = ", q, "."))
        stop()
      }
      if(length(dim(Gamma)) > 2){
        print(paste("The stationary covariance matrix Gamma should be an q times q matrix, with q = ", q, "."))
        stop()
      }
    }else if(q != 1){
      print(paste("The stationary covariance matrix Gamma should, like Phi, be a scalar."))
      stop()
    }

    # Calculate SigmaVAR
    SigmaVAR <- Gamma - Phi %*% Gamma %*% t(Phi)

  }else{ # Both SigmaVAR and Gamma known

    # Checks on SigmaVAR
    if(length(SigmaVAR) != 1){
      if(dim(SigmaVAR)[1] != dim(SigmaVAR)[2]){
        print(paste("The residual covariance matrix SigmaVAR should be a square matrix of size q times q, with q = ", q, "."))
        stop()
      }
      if(dim(SigmaVAR)[1] != q){
        print(paste("The residual covariance matrix SigmaVAR should, like Phi, be a matrix of size q times q, with q = ", q, "."))
        stop()
      }
      if(length(dim(SigmaVAR)) > 2){
        print(paste("The residual covariance matrix SigmaVAR should be an q times q matrix, with q = ", q, "."))
        stop()
      }
    }else if(q != 1){
      print(paste("The residual covariance matrix SigmaVAR should, like Phi, be a scalar."))
      stop()
    }

    # Checks on Gamma
    if(length(Gamma) != 1){
      if(dim(Gamma)[1] != dim(Gamma)[2]){
        print(paste("The stationary covariance matrix Gamma should be a square matrix of size q times q, with q = ", q, "."))
        stop()
      }
      if(dim(Gamma)[1] != q){
        print(paste("The stationary covariance matrix Gamma should, like Phi, be a matrix of size q times q, with q = ", q, "."))
        stop()
      }
      if(length(dim(Gamma)) > 2){
        print(paste("The stationary covariance matrix Gamma should be an q times q matrix, with q = ", q, "."))
        stop()
      }
    }else if(q != 1){
      print(paste("The stationary covariance matrix Gamma should, like Phi, be a scalar."))
      stop()
    }

  }


Sxy <- sqrt(diag(diag(Gamma)))
Gamma_s <- solve(Sxy) %*% Gamma %*% solve(Sxy)
Phi_s <- solve(Sxy) %*% Phi %*% Sxy
SigmaVAR_s <- solve(Sxy) %*% SigmaVAR %*% solve(Sxy)

vecPhi <- as.vector(t(Phi_s))

CovMx <- kronecker(SigmaVAR_s, solve(Gamma_s)) / (N-q)



# Determine points on 95% LL contour
mu_Phi <- vecPhi
CovMx_Phi <- CovMx
eigenCovMx <- eigen(CovMx_Phi)
lambda <- eigenCovMx$val
E <- eigenCovMx$vec
#df1F <- q*q*qf(p=alpha, df1=q*q, df2=(N-q*q), lower.tail=FALSE) # geeft bij mij 25 getallen?! TO DO
Chi2 <- qchisq(p=alpha, df=(q*q), lower.tail=FALSE) # for large N df1F goes to Chi2
LB_vecPhi <- matrix(NA, nrow=q*q, ncol =q*q)
UB_vecPhi <- matrix(NA, nrow=q*q, ncol =q*q)
LL <- matrix(NA, nrow=q*q, ncol=2)
teller = 0
for(row in 1:q){
  for(column in 1:q){
    teller = teller + 1
    #LB_vecPhi[teller,] <- matrix(mu_Phi - sqrt(df1F * lambda[teller]) * E[,teller], nrow = 1)
    #UB_vecPhi[teller,] <- matrix(mu_Phi + sqrt(df1F * lambda[teller]) * E[,teller], nrow = 1)
    LB_vecPhi[teller,] <- matrix(mu_Phi - sqrt(Chi2 * lambda[teller]) * E[,teller], nrow = 1)
    UB_vecPhi[teller,] <- matrix(mu_Phi + sqrt(Chi2 * lambda[teller]) * E[,teller], nrow = 1)
    LL[teller,1] <- t(LB_vecPhi[teller,]-mu_Phi) %*% solve(CovMx_Phi) %*% (LB_vecPhi[teller,]-mu_Phi)
    LL[teller,2] <- t(UB_vecPhi[teller,]-mu_Phi) %*% solve(CovMx_Phi) %*% (UB_vecPhi[teller,]-mu_Phi)
  }
}
minPhi <- apply(rbind(LB_vecPhi, UB_vecPhi), 2, min)
maxPhi <- apply(rbind(LB_vecPhi, UB_vecPhi), 2, max)
multiCI <- rbind(minPhi, maxPhi)
rownames(multiCI) <- c("LB", "UB")
sub = NULL
for(i in 1:q){
  sub = c(sub, paste0("Phi", i, 1:q, sep=""))
}
colnames(multiCI) <- sub


############################################################################################################

final <- list(vecStandPhi_DeltaT = vecPhi, CovMx_vecStandPhi_DeltaT = CovMx, multiCI_vecStandPhi_DeltaT = multiCI,
              standSigmaVAR_DeltaT = SigmaVAR_s, standGamma = Gamma_s)

return(final)
}
