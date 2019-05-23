
#' Calculates the (vectorized) standardized lagged effects matrix, their covariance matrix, and corresponding elliptical 95 percent confidence interval
#'
#' @param N Number of persons (panel data) or measurement occasions - 1 (time series data). Used in determining the covariance matrix of the vectorized standardized lagged effects.
#' @param Phi (Un)standardized lagged effects matrix. If necessary, it is standardized and for the standardized and vectorized Phi the covariance matrix is determined.
#' @param Gamma Stationary covariance matrix, that is, the contemporaneous covariance matrix of the data.
#' @param SigmaVAR Residual covariance matrix of the first-order discrete-time vector autoregressive (DT-VAR(1)) model
#' @param alpha The alpha level in determining the (1-alpha)*100 percent confidence interval (CI). By default, alpha is set to 0.05, resulting in a 95& CI
#'
#' @return vectorized standardized lagged effects, their covariance matrix, and the corresponding elliptical/multivariate 95 percent CI.
#' @export
#' @examples
#'
#' # Example where Phi and SigmaVAR are known
#' N <- 643
#' q <- 2
#' Phi <- matrix(c(0.25, 0.10,
#'                 0.20, 0.36), byrow=T, ncol = q)
#' SigmaVAR <- diag(q) # for ease
#' Gamma <- calc.Gamma.fromVAR(Phi, SigmaVAR)
#' calc.CovMxStandPhi(N, Phi, Gamma, SigmaVAR)
#'
#' # Example where Phi and Gamma are known
#' N <- 643
#' q <- 2
#' Phi <- matrix(c(0.25, 0.10,
#'                 0.20, 0.36), byrow=T, ncol = q)
#' # Use Gamma from above
#' SigmaVAR <- Gamma - Phi %*% Gamma %*% t(Phi)
#' calc.CovMxStandPhi(N, Phi, Gamma, SigmaVAR)



calc.CovMxStandPhi <- function(N, Phi, Gamma, SigmaVAR, alpha=0.05) {

q <- dim(Phi)[1]

Sxy <- sqrt(diag(diag(Gamma)))
Gamma_s <- solve(Sxy) %*% Gamma %*% solve(Sxy)
Phi_s <- solve(Sxy) %*% Phi %*% Sxy
SigmaVAR_s <- solve(Sxy) %*% SigmaVAR %*% solve(Sxy)

vecPhi <- as.vector(t(Phi_s))

CovMx = kronecker(SigmaVAR_s, solve(Gamma_s)) / (N-q)



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

final <- list(vecStandPhi_DeltaT = vecPhi, CovMx_vecStandPhi_DeltaT = CovMx, multiCI_vecStandPhi_DeltaT = multiCI)

return(final)
}
