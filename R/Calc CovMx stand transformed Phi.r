
#' StandTransPhi
#'
#' This function calculates the (vectorized) transformed standardized Phi, their covariance matrix and elliptical 95\% confidence interval (CI). There is also an interactive web application on my website: Standardizing and/or transforming lagged regression estimates (\url{https://www.uu.nl/staff/RMKuiper/Websites\%20\%2F\%20Shiny\%20apps}).
#'
#' @param DeltaTStar The time interval to which the (un)standardized lagged effects matrix (Phi) should be transformed to.
#' @param DeltaT Optional. The time interval used. Hence, Phi(DeltaT) will be transformed to Phi(DeltaTStar) and standardized. By default, DeltaT = 1.
#' @param N Optional. Number of persons (panel data) or number measurement occasions - 1 (time series data). This is used in determining the covariance matrix of the vectorized standardized lagged effects. By default, N = NULL.
#' @param Phi (Un)standardized lagged effects matrix. If necessary, it is standardized, then it is transformed and for this vectorized Phi the covariance matrix is determined.
#' It also takes a fitted object from the classes "varest" (from the VAR() function in vars package) and "ctsemFit" (from the ctFit() function in the ctsem package); see example below. From such an object, the Phi, SigmaVAR, and Gamma matrices are calculated/extracted.
#' @param SigmaVAR Residual covariance matrix of the first-order discrete-time vector autoregressive (DT-VAR(1)) model.
#' @param Gamma Optional (either SigmaVAR or Gamma). Stationary covariance matrix, that is, the contemporaneous covariance matrix of the data.
#' Note that if Phi and SigmaVAR are known, Gamma can be calculated; hence, only SigmaVAR or Gamma is needed as input (if only Gamma, then use 'Gamma = Gamma' or set SigmaVAR to NULL, see examples below).
#' @param alpha Optional. The alpha level in determining the (1-alpha)*100\% CI. By default, alpha = 0.05; resulting in a 95\% CI
#'
#' @return This function returns the vectorized transformed (standardized) lagged effects and - if N is part of input - their covariance matrix and the corresponding elliptical/multivariate 95\% CI. If both Sigma and Gamma are not part of the input, only the transformed Phi is calculated.
#' @importFrom expm expm
#' @export
#' @examples
#' ## Example 1 ##
#'
#' # Input for examples below
#' DeltaTStar <- 1
#' DeltaT <- 2
#' N <- 643
#' # Phi(DeltaT)
#' Phi <- myPhi[1:2,1:2]
#' #Phi <- matrix(c(0.25, 0.10,
#' #               0.20, 0.36), byrow=T, ncol = 2)
#' # SigmaVAR(DeltaT)
#' SigmaVAR <- diag(q) # for ease
#' # Calculate the Gamma corresponding to Phi and SigmaVAR - used in the second example
#' Gamma <- Gamma.fromVAR(Phi, SigmaVAR) # ?Gamma.fromVAR
#'
#' #Example where only SigmaVAR is known and not Gamma
#' StandTransPhi(DeltaTStar, DeltaT, N, Phi, SigmaVAR)
#'
#' #Example where only Gamma is known and not SigmaVAR
#' StandTransPhi(DeltaTStar, DeltaT, N, Phi, NULL, Gamma)
#' # or
#' StandTransPhi(DeltaTStar, DeltaT, N, Phi, Gamma = Gamma)
#'
#'
#' ## Example 2: input from fitted object of class "varest" ##
#' #
#' DeltaTStar <- 1
#' DeltaT <- 2
#' N <- 643
#' data <- myData
#' if (!require("vars")) install.packages("vars")
#' library(vars)
#' out_VAR <- VAR(data, p = 1)
#' StandTransPhi(DeltaTStar, DeltaT, N, out_VAR)
#'
#'
#' ## Example 3: obtain only (un)standardized transformed lagged effects ##
#' # Note: Use Phi and SigmaVAR from Example 1
#' DeltaTStar <- 1
#' DeltaT <- 2
#' StandTransPhi(DeltaTStar, DeltaT, N = NULL, Phi, SigmaVAR)
#' # or
#' StandTransPhi(DeltaTStar, DeltaT, Phi = Phi, SigmaVAR = SigmaVAR)
#'
#'
#' ## Example 4: obtain only (unstandardized) transformed lagged effects ##
#' # Note: Use Phi from Example 1
#' DeltaTStar <- 1
#' DeltaT <- 2
#' StandTransPhi(DeltaTStar, DeltaT, N = NULL, Phi)
#' # or
#' StandTransPhi(DeltaTStar, DeltaT, Phi = Phi)
#'


StandTransPhi <- function(DeltaTStar, DeltaT = 1, N = NULL, Phi, SigmaVAR = NULL, Gamma = NULL, alpha=0.05) {

  # Checks:
  if(length(DeltaTStar) != 1){
    print(paste("The argument DeltaTStar should be a scalar, that is, one number, that is, a vector with one element."))
    stop()
  }
  if(length(DeltaT) != 1){
    print(paste("The argument DeltaT should be a scalar, that is, one number, that is, a vector with one element."))
    stop()
  }
  if(!is.null(N) & length(N) != 1){
    print(paste("The argument N should be a scalar, that is, one number, that is, a vector with one element."))
    stop()
  }
  #
  # Check on Phi
  if(any(class(Phi) == "varest")){
    SigmaVAR <- cov(resid(Phi))
    Phi <- Acoef(Phi)[[1]]
    if(length(Phi) == 1){
      q <- 1
    }else{
      q <- dim(Phi)[1]
    }
  } else if(any(class(Phi) == "ctsemFit")){
    B <- -1 * summary(Phi)$DRIFT
    Sigma <- summary(Phi)$DIFFUSION
    #Phi <- summary(Phi)$discreteDRIFT # Is no longer output in ctsem...
    Phi <- expm(-B*DeltaT)
    #source("Calc VARparam from CTMparam.r") # werkt zo niet, moet er dan ws ook net als andere files package fie van maken
    #VarEst <- VARparam(DeltaT, B, Sigma)
    #Phi <- VarEst$Phi
    if(length(Phi) == 1){
      q <- 1
    }else{
      q <- dim(Phi)[1]
    }
    #SigmaVAR <- VarEst$SigmaVAR
    kronsum <- kronecker(diag(q),B) + kronecker(B,diag(q))
    SigmaVAR <- matrix((solve(kronsum) %*% (diag(q*q) - expm(-kronsum * DeltaT)) %*% as.vector(Sigma)), ncol=q, nrow=q)
  } else{
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
      #
      if(dim(Phi)[1] != dim(Phi)[2]){
        print(paste("The lagged effects matrix Phi should be a square matrix of size q times q, with q = ", q, "."))
        stop()
      }
      if(length(dim(Phi)) > 2){
        print(paste("The lagged effects matrix Phi should be an q times q matrix, with q = ", q, "."))
        stop()
      }
    } else{
      q <- 1
    }
  }

  # Check on SigmaVAR and Gamma
  if(is.null(SigmaVAR) & is.null(Gamma)){ # Both SigmaVAR and Gamma unknown
    print("Just to be sure:")
    print(paste("The arguments SigmaVAR and Gamma are not found: Both SigmaVAR and Gamma are NULL."))
    print(paste("Now, only the (unstandardized) transformed Phi is calculated."))
    print(paste("In case you want a standarized matrix or more: either SigmaVAR and/or Gamma should be part of the input. In case of first matrix, specify 'SigmaVAR = <insert your matrix name>'."))
    print("")
    #stop()
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
    Gamma <- Gamma.fromVAR(Phi, SigmaVAR)

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



ratioDeltaT <- DeltaTStar / DeltaT

if(q > 1){
  eigenPhi <- eigen(Phi)
  V <- eigenPhi$vectors
  D <- diag(eigenPhi$values)
  Phi_DeltaT <- V %*% D^ratioDeltaT %*% solve(V)
} else{
  Phi_DeltaT <- Phi^ratioDeltaT
}

warning <- "No warnings (since there are no complex eigenvalues)"
Warning <- 0
if(any(is.complex(eigenPhi$values))){
  if(DeltaTStar%%DeltaT == 0){
    warning <- "There is at least one pair of complex eigenvalues and the ratio of DeltaTs (i.e., DeltaT*/DeltaT) is not an integer, so the solution for Phi(DeltaT*) is NOT unique."
    Warning <- 1
  }else{
    warning <- "There is at least one pair of complex eigenvalues, but the ratio of DeltaTs (i.e., DeltaT*/DeltaT) is an integer, so the solution for Phi(DeltaT*) is unique."
  }
}


if(!(is.null(SigmaVAR) & is.null(Gamma))){

  SigmaVAR_DeltaT <- Gamma - Phi_DeltaT %*% Gamma %*% t(Phi_DeltaT)

  #CTMparam <- CTMparam(ratioDeltaT, Phi, SigmaVAR)
  #CTMparam$WarningPhi
  #if (!is.null(CTMparam$WarningPhi)){
  #  WarningPhi = "'Phi' does not have a CTM-equivalent drift matrix; so, no positive autocorrelation (like CTM models). For example, one or more real eigenvalues is negative."
  #  final <- list(WarningPhi = WarningPhi)
  #  return(final)
  #  stop(WarningPhi)
  #}
  #VARparam <- VARparam(1, CTMparam$A , CTMparam$SigmaCTM)
  #Phi_DeltaT <- VARparam$Phi_DeltaT
  #SigmaVAR_DeltaT <- VARparam$SigmaVAR_DeltaT


  Sxy <- sqrt(diag(diag(Gamma)))
  Gamma_s <- solve(Sxy) %*% Gamma %*% solve(Sxy)
  Phi_DeltaT_s <- solve(Sxy) %*% Phi_DeltaT %*% Sxy
  SigmaVAR_DeltaT_s <- solve(Sxy) %*% SigmaVAR_DeltaT %*% solve(Sxy)


  vecPhi <- as.vector(t(Phi_DeltaT_s))


  if(!is.null(N)){

    CovMx = kronecker(SigmaVAR_DeltaT_s, solve(Gamma_s)) / (N-q)

    # Determine points on 95% LL contour
    mu_Phi <- vecPhi
    CovMx_Phi <- CovMx
    eigenCovMx <- eigen(CovMx_Phi)
    lambda <- eigenCovMx$val
    E <- eigenCovMx$vec
    #df1F <- q*q*qf(p=alpha, df1=q*q, df2=(N-q*q), lower.tail=FALSE)
    Chi2 <- qchisq(p=alpha, df=(q*q), lower.tail=FALSE) # for large N, df1F goes to Chi2
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
  }

}


############################################################################################################

if(is.null(SigmaVAR) & is.null(Gamma)){
  final <- list(Phi_DeltaTStar = Phi_DeltaT,
                Warning = Warning, warning = warning)
}else if(!is.null(N)){
  final <- list(Phi_DeltaTStar = Phi_DeltaT,
                SigmaVAR_DeltaTStar = SigmaVAR_DeltaT, Gamma = Gamma,
                standPhi_DeltaTStar = Phi_DeltaT_s,
                vecStandPhi_DeltaTStar = vecPhi, CovMx_vecStandPhi_DeltaTStar = CovMx, multiCI_vecStandPhi_DeltaTStar = multiCI,
                standSigmaVAR_DeltaTStar = SigmaVAR_DeltaT_s, standGamma = Gamma_s,
                Warning = Warning, warning = warning)
} else{
  final <- list(Phi_DeltaTStar = Phi_DeltaT,
                SigmaVAR_DeltaTStar = SigmaVAR_DeltaT, Gamma = Gamma,
                standPhi_DeltaTStar = Phi_DeltaT_s,
                standSigmaVAR_DeltaTStar = SigmaVAR_DeltaT_s, standGamma = Gamma_s,
                Warning = Warning, warning = warning)
}

return(final)
}
