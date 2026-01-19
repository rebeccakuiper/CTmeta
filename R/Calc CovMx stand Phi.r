
#' Calculate standardized lagged effects matrix and accompanying matrices
#'
#' This function calculates the (vectorized) standardized lagged effects matrix, their covariance matrix, and corresponding elliptical 95\% confidence interval (CI). There is also an interactive web application on my website: Standardizing and/or transforming lagged regression estimates (\url{https://www.uu.nl/staff/RMKuiper/Websites\%20\%2F\%20Shiny\%20apps}).
#'
#' @param N Optional. Number of persons (panel data) or number of measurement occasions - 1 (time series data). This is used in determining the covariance matrix of the vectorized standardized lagged effects. By default, N = NULL.
#' @param Phi (Un)standardized lagged effects matrix. If necessary, it is standardized and for the standardized and vectorized Phi the covariance matrix is determined.
#' It also takes a fitted object from the classes "varest" (from the VAR() function in vars package) and "ctsemFit" (from the ctFit() function in the ctsem package); see example below. From such an object, the Phi, SigmaVAR, and Gamma matrices are calculated/extracted.
#' @param SigmaVAR Residual covariance matrix of the first-order discrete-time vector autoregressive (DT-VAR(1)) model.
#' @param Gamma Optional (either SigmaVAR or Gamma). Stationary covariance matrix, that is, the contemporaneous covariance matrix of the data.
#' Note that if Phi and SigmaVAR are known, Gamma can be calculated; hence, only SigmaVAR or Gamma is needed as input (if only Gamma, then use 'Gamma = Gamma' or set SigmaVAR to NULL, see examples below).
#' @param alpha Optional. The alpha level in determining the (1-alpha)*100\% CI. By default, alpha is set to 0.05, resulting in a 95\% CI.
#'
#' @return This function returns the vectorized standardized lagged effects and - if N is part of input - their covariance matrix and the corresponding elliptical/multivariate 95\% CI.
#' @importFrom expm expm
#' @export
#' @examples
#'
#' # library(CTmeta)
#'
#' ## Example 1 ##
#'
#' # Input for examples below
#' N <- 643
#' Phi <- myPhi[1:2,1:2]
#' #Phi <- matrix(c(0.25, 0.10,
#' #                0.20, 0.36), byrow=T, ncol = 2)
#' q <- dim(Phi)[1]
#' SigmaVAR <- diag(q) # for ease
#' # Calculate the Gamma corresponding to Phi and SigmaVAR - used in the second example
#' Gamma <- Gamma.fromVAR(Phi, SigmaVAR) # ?Gamma.fromVAR
#'
#' #Example where SigmaVAR is known and Gamma unknown
#' StandPhi(N, Phi, SigmaVAR)
#'
#' #Example where Gamma is known and SigmaVAR unknown
#' StandPhi(N, Phi, NULL, Gamma)
#' # or
#' StandPhi(N, Phi, Gamma = Gamma)
#'
#'
#' ## Example 2: input from fitted object of class "varest" ##
#' #
#' N <- 643
#' data <- myData
#' if (!require("vars")) install.packages("vars")
#' library(vars)
#' out_VAR <- VAR(data, p = 1)
#' StandPhi(N, out_VAR)
#'
#'
#' ## Example 3: obtain only standardized lagged effects ##
#' Phi <- myPhi[1:2,1:2]
#' q <- dim(Phi)[1]
#' SigmaVAR <- diag(q) # for ease
#' #
#' StandPhi(N = NULL, Phi, SigmaVAR)
#' # or
#' StandPhi(Phi = Phi, SigmaVAR = SigmaVAR)
#'


StandPhi <- function(N = NULL, Phi, SigmaVAR = NULL, Gamma = NULL, alpha = 0.05) {

  # Checks:
  if(!is.null(N) & length(N) != 1){
    ErrorMessage <- (paste0("The argument N should be a scalar, that is, one number, that is, a vector with one element. Currently, N = ", N))
    return(ErrorMessage)
    stop(ErrorMessage)
  }
  #
  if(length(alpha) != 1){
    ErrorMessage <- (paste0("The argument alpha should be a scalar, that is, one number, that is, a vector with one element. Currently, alpha = ", alpha))
    return(ErrorMessage)
    stop(ErrorMessage)
  }
  #
  # Check on Phi
  if(any(class(Phi) == "varest")){
    SigmaVAR <- cov(resid(Phi))
    Phi <- Acoef(Phi)[[1]]
  } else if(any(class(Phi) == "ctsemFit")){
    B <- -1 * summary(Phi)$DRIFT
    Sigma <- summary(Phi)$DIFFUSION
    #
    VarEst <- VARparam(DeltaT, -B, Sigma)
    Phi <- VarEst$Phi
    SigmaVAR <- VarEst$SigmaVAR
  } else{
    #
    if(length(Phi) != 1){
      Check_Phi(Phi)
      q <- dim(Phi)[1]
    } else{
      q <- 1
    }
    #
    # Check on SigmaVAR and Gamma
    if(is.null(SigmaVAR) & is.null(Gamma)){ # Both SigmaVAR and Gamma unknown
      ErrorMessage <- (paste0("Both SigmaVAR and Gamma are NULL; either one (or both) should be part of the input. In case of first matrix, specify 'SigmaVAR = SigmaVAR'."))
      return(ErrorMessage)
      stop(ErrorMessage)
    }else if(!is.null(SigmaVAR)){ # SigmaVAR known
      # Check on SigmaVAR
      Check_SigmaVAR(SigmaVAR, q)
    }else if(!is.null(Gamma)){ # Gamma known
      # Checks on Gamma
      Check_Gamma(Gamma, q)
    }
  }
  #
  if(length(Phi) == 1){
    q <- 1
  }else{
    q <- dim(Phi)[1]
  }

  # Calculate Gamma and Sigma - if necessary
  if(!is.null(Gamma)){ # Gamma known
    if(is.null(SigmaVAR)){ # SigmaVAR unknown, calculate it
      # Calculate SigmaVAR
      if(q != 1){
        SigmaVAR <- Gamma - Phi %*% Gamma %*% t(Phi)
      }else{
        SigmaVAR <- Gamma - Phi * Gamma * Phi
      }
    }
  }else if(!is.null(SigmaVAR)){  # Gamma unknown and SigmaVAR known, calculate Gamma from SigmaVAR and Phi
    # Calculate Gamma
    Gamma <- Gamma.fromVAR(Phi, SigmaVAR)
  }


if(q > 1){
  Sxy <- sqrt(diag(diag(Gamma)))
  Gamma_s <- solve(Sxy) %*% Gamma %*% solve(Sxy)
  Phi_s <- solve(Sxy) %*% Phi %*% Sxy
  SigmaVAR_s <- solve(Sxy) %*% SigmaVAR %*% solve(Sxy)
  #
  S <- sqrt(diag(diag(SigmaVAR)))
  ResidCorrMx <- solve(S) %*% SigmaVAR %*% solve(S)
  #
  vecPhi <- as.vector(t(Phi_s))
}else{
  Sxy <- sqrt(diag(diag(Gamma)))
  Gamma_s <- solve(Sxy) * Gamma * solve(Sxy)
  Phi_s <- solve(Sxy) * Phi * Sxy
  SigmaVAR_s <- solve(Sxy) * SigmaVAR * solve(Sxy)
  #
  S <- sqrt(diag(diag(SigmaVAR)))
  ResidCorrMx <- solve(S) * SigmaVAR * solve(S)
  #
  vechi <- Phi_s
}


if(!is.null(N)){

  CovMx <- kronecker(SigmaVAR_s, solve(Gamma_s)) / (N-q)

  # Determine points on 95% LL contour
  mu_Phi <- vecPhi
  CovMx_Phi <- CovMx
  eigenCovMx <- eigen(CovMx_Phi)
  lambda <- eigenCovMx$val
  if(any(lambda < 0)){
    message("Some of the eigenvalues of the covariance matrix of Phi are negative. \n",
            "The corresponding confidence interval(s), in $multiCI_vecStandPhi_DeltaT, are NA."
    )
    which_lambda <- which(eigenCovMx$val < 0)
    lambda[which_lambda] <- 0
  }
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
  if(any(lambda <= 0)){
    minPhi[which_lambda] <- NA
    maxPhi[which_lambda] <- NA
  }
  multiCI <- rbind(minPhi, maxPhi)
  rownames(multiCI) <- c("LB", "UB")
  sub = NULL
  for(i in 1:q){
    sub = c(sub, paste0("Phi", i, 1:q, sep=""))
  }
  colnames(multiCI) <- sub
}


############################################################################################################

if(!is.null(N)){
  final <- list(Phi_DeltaT = Phi, StandPhi_DeltaT = Phi_s,
                vecStandPhi_DeltaT = vecPhi, CovMx_vecStandPhi_DeltaT = CovMx,
                multiCI_vecStandPhi_DeltaT = multiCI,
                SigmaVAR_DeltaT = SigmaVAR, standSigmaVAR_DeltaT = SigmaVAR_s,
                ResidCorrMx = ResidCorrMx,
                Gamma = Gamma, standGamma = Gamma_s)
  #N = N, alpha = alpha
}else{
  final <- list(Phi_DeltaT = Phi, StandPhi_DeltaT = Phi_s,
                SigmaVAR_DeltaT = SigmaVAR, standSigmaVAR_DeltaT = SigmaVAR_s,
                ResidCorrMx = ResidCorrMx,
                Gamma = Gamma, standGamma = Gamma_s)
  #N = N
}

return(final)
}
