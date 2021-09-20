
#' TransPhi_Corr
#'
#' Calculates the (vectorized) transformed standardized Phi, their covariance matrix, the corresponding elliptical 95\% confidence interval (CI) from a correlation matrix with contemporaneous and lagged correlations. There is also an interactive web application on my website: Standardizing and/or transforming lagged regression estimates (\url{https://www.uu.nl/staff/RMKuiper/Websites\%20\%2F\%20Shiny\%20apps}).
#'
#' @param DeltaTStar The time interval to which the standardized lagged effects matrix (Phi) should be transformed to.
#' @param DeltaT Optional. The time interval used. Hence, Phi(DeltaT) will be transformed to Phi(DeltaTStar) and standardized. By default, DeltaT = 1.
#' @param N Optional. Number of persons (panel data) or number of measurement occasions - 1 (time series data). This is used in determining the covariance matrix of the vectorized standardized lagged effects. By default, N = NULL.
#' @param corr_YXYX The correlation matrix of the variables and the lagged variables (of size 2q x 2q). The upper q x q matrix is the correlation matrix between the (q) variables and the lower q x q matrix is the correlation matrix between the (q) lagged variables.
#' @param alpha Optional. The alpha level in determining the (1-alpha)*100\% CI. By default, alpha = 0.05; resulting in a 95\% CI.
#'
#' @return This function returns the vectorized transformed standardized lagged effects (i.e., for DeltaTStar) and - if N is part of input - their covariance matrix and the corresponding elliptical/multivariate 95\% CI; SigmaVAR: residual covariance matrix for DeltaTStar; and Gamma: stationary covariance matrix.
#' @export
#' @examples
#'
#' # library(CTmeta)
#'
#' # In the examples below, the following values are used:
#' DeltaTStar <- 12
#' DeltaT <- 24
#' N <- 2235
#'
#' # Example with full correlation matrix
#' corr_YXYX <- matrix(c(1.00, 0.40, 0.63, 0.34,
#' 0.40, 1.00, 0.31, 0.63,
#' 0.63, 0.31, 1.00, 0.41,
#' 0.34, 0.63, 0.41, 1.00), byrow = T, ncol = 2*2)
#' # Run function
#' TransPhi_Corr(DeltaTStar, DeltaT, N, corr_YXYX)
#'
#' # Example with vector of lower triangular correlation matrix
#' LT <- c(0.40, 0.63, 0.34, 0.31, 0.63, 0.41) # corr_YXYX[lower.tri(corr_YXYX,diag = F)]
#' # Make full correlation matrix of size 2*q times 2*q, with q=2 and thus 2*q=4
#' corr_YXYX <- diag(4) # As check: length(LT) = 4*(4-1)/2
#' corr_YXYX[lower.tri(corr_YXYX,diag = F)] <- LT
#' corr_YXYX[upper.tri(corr_YXYX,diag = F)] <- t(corr_YXYX)[upper.tri(t(corr_YXYX),diag = F)]
#' # Run function
#' TransPhi_Corr(DeltaTStar, DeltaT, N, corr_YXYX)
#'
#' # Example with vector of lower triangular correlation matrix including diagonals
#' LTD <- c(1.00, 0.40, 0.63, 0.34, 1.00, 0.31, 0.63, 1.00, 0.41, 1.00) # corr_YXYX[lower.tri(corr_YXYX,diag = T)]
#' # Make full correlation matrix of size 2*q times 2*q, with q=2 and thus 2*q=4
#' corr_YXYX <- matrix(NA, nrow=(4), ncol=(4))  # As check: length(LTD) = 4*(4+1)/2
#' corr_YXYX[lower.tri(corr_YXYX,diag = T)] <- LTD
#' corr_YXYX[upper.tri(corr_YXYX,diag = F)] <- t(corr_YXYX)[upper.tri(t(corr_YXYX),diag = F)]
#' # Run function
#' TransPhi_Corr(DeltaTStar, DeltaT, N, corr_YXYX)
#'
#' # The output (standPhi_DeltaTStar, standSigmaVAR_DeltaTStar, and standGamma) can be used to make stacked matrices or arrays which can serve as input for continuous-time meta-analysis CTmeta (using the function CTmeta; see ?CTmeta).
#'


TransPhi_Corr <- function(DeltaTStar, DeltaT = 1, N = NULL, corr_YXYX, alpha = 0.05) {

  # Checks:
  if(length(DeltaTStar) != 1){
    ErrorMessage <- (paste0("The argument DeltaTStar should be a scalar, that is, one number, that is, a vector with one element. Currently, DeltaTStar = ", DeltaTStar))
    return(ErrorMessage)
    stop(ErrorMessage)
  }
  if(length(DeltaT) != 1){
    ErrorMessage <- (paste0("The argument DeltaT should be a scalar, that is, one number, that is, a vector with one element. Currently, DeltaT = ", DeltaT))
    return(ErrorMessage)
    stop(ErrorMessage)
  }
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
  # Check on corr_YXYX
  if(is.null(dim(corr_YXYX))){
    if(!is.null(length(corr_YXYX))){ # Should be matrix
      ErrorMessage <- (paste0("The argument corr_YXYX is not a matrix. It should be a matrix of size 2q times 2q."))
      return(ErrorMessage)
      stop(ErrorMessage)
    }else{
      ErrorMessage <- (paste0("The argument corr_YXYX is not found: The (lagged) correlation matrix corr_YXYX is unknown, but should be part of the input."))
      return(ErrorMessage)
      stop(ErrorMessage)
    }
  }
  # Should be square matrix
    #if(dim(corr_YXYX)[1] != dim(corr_YXYX)[2] | length(dim(corr_YXYX)) != 2){
    #  ErrorMessage <- (paste0("The argument corr_YXYX is not a matrix of size 2q times 2q. Currently, it is of size ", dim(corr_YXYX)))
    #  stop(ErrorMessage)
    #}
    if(length(dim(corr_YXYX)) < 2){
      ErrorMessage <- (paste0("The argument corr_YXYX is not a matrix. It should be a matrix of size 2q times 2q."))
      return(ErrorMessage)
      stop(ErrorMessage)
    }
    if(length(dim(corr_YXYX)) > 2){
      ErrorMessage <- (paste0("The argument corr_YXYX is not a matrix of size 2q times 2q. Currently, it is of size ", dim(corr_YXYX)))
      return(ErrorMessage)
      stop(ErrorMessage)
    }
    if(dim(corr_YXYX)[1] != dim(corr_YXYX)[2]){
      ErrorMessage <- (paste0("The argument corr_YXYX should be a square matrix. It should be a matrix of size 2q times 2q. Currently, it is of size ", dim(corr_YXYX)[1], " times ", dim(corr_YXYX)[2]))
      return(ErrorMessage)
      stop(ErrorMessage)
    }
    #
    if((dim(corr_YXYX)[1] %% 2) != 0) { # 2q should be an even number
      ErrorMessage <- (paste0("The argument corr_YXYX is not a matrix of size 2q times 2q, since the number of rows/columns is not even but an odd number, namely ", dim(corr_YXYX)[1]))
      return(ErrorMessage)
      stop(ErrorMessage)
    }


    q <- dim(corr_YXYX)[1]/2

    ratioDeltaT <- DeltaTStar / DeltaT

    RYY <- corr_YXYX[1:q,1:q]
    RYX <- corr_YXYX[1:q,(q+1):(q+q)]
    RXY <- corr_YXYX[(q+1):(q+q),1:q]
    RXX <- corr_YXYX[(q+1):(q+q),(q+1):(q+q)]

    Phi <- t(solve(RXX) %*% RXY)

    eigenPhi <- eigen(Phi)
    V <- eigenPhi$vectors
    D <- diag(eigenPhi$values)
    TransPhi <- V %*% D^ratioDeltaT %*% solve(V)

    warning = "No warnings (since there are no complex eigenvalues)"
    if(any(is.complex(eigenPhi$values))){
      if(DeltaTStar%%DeltaT == 0){
        warning <- "There is at least one pair of complex eigenvalues and the ratio of DeltaTs (i.e., DeltaT*/DeltaT) is not an integer, so the solution for Phi(DeltaT*) is NOT unique."
      }else{
        warning <- "There is at least one pair of complex eigenvalues, but the ratio of DeltaTs (i.e., DeltaT*/DeltaT) is an integer, so the solution for Phi(DeltaT*) is unique."
      }
      if(all(Im(TransPhi) == 0)){
        TransPhi <- Re(TransPhi)
      }
    }


    RXY <- RXX %*% t(TransPhi)
    RYX <- t(RXY)
    # RYY and RXX remain the same!


    vecPhi <- as.vector(t(TransPhi))

    SigmaVAR <- RYY - RYX  %*%  solve(RXX) %*% RXY
    Gamma <- RXX


    if(!is.null(N)){

      CovMx <- kronecker(SigmaVAR, solve(Gamma)) / (N-q)

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


    ############################################################################################################

    if(!is.null(N)){
      final <- list(standPhi_DeltaTStar = TransPhi,
                    vecStandPhi_DeltaTStar = vecPhi, CovMx_vecStandPhi_DeltaTStar = CovMx, multiCI_vecStandPhi_DeltaT = multiCI,
                    standSigmaVAR_DeltaTStar = SigmaVAR, standGamma = Gamma,
                    warning = warning)
    }else{
      final <- list(standPhi_DeltaTStar = TransPhi,
                    standSigmaVAR_DeltaTStar = SigmaVAR, standGamma = Gamma,
                    warning = warning)
    }

    return(final)
  }

