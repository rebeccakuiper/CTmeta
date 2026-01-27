
#' TransPhi_Corr
#'
#' Calculates the (vectorized) transformed standardized Phi, their covariance matrix, the corresponding elliptical 95\% confidence interval (CI) from a correlation matrix with contemporaneous and lagged correlations. There is also an interactive web application on my website: Standardizing and/or transforming lagged regression estimates (\url{https://www.uu.nl/staff/RMKuiper/Websites\%20\%2F\%20Shiny\%20apps}).
#'
#' @param DeltaTStar The time interval to which the standardized lagged effects matrix (Phi) should be transformed to. This should be a vector with one number.
#' @param DeltaT Optional. The time interval used. Hence, Phi(DeltaT) will be transformed to Phi(DeltaTStar) and standardized. This should be a vector with one or ('waves'-1) numbers. By default, DeltaT = 1.
#' @param N Optional. Number of persons (panel data) or number of measurement occasions - 1 (time series data). This is used in determining the covariance matrix of the vectorized standardized lagged effects. By default, N = NULL.
#' @param waves Optional. Number of waves, that is, number of measurement occasions; which should be an integer of at least 2. If waves > 3, then both a) wave-independent and b) wave-specific transformed standardized order-1 lagged effects are calculated. By default, waves = 2.
#' @param corr_lagged The (lagged) correlation matrix of the q variables of size (q x waves) times (q x waves). The upper left q x q matrix is the correlation matrix between the q variables on the first time point; the upper right q x q matrix is the lagged correlation matrix between the q variables for the first and last time point (btw this block matrix and the other ones for non-subsequent time points will not be used to determine the order-1 relationships). An example can be found below. By default, corr_lagged = NULL.
#' @param corr_YXYX Used if not NULL (default) & corr_lagged = NULL. Also the (lagged) correlation matrix of the q variables of size (q x waves) times (q x waves) but, in this matrix, the time points are ordered from high to low. Therefore, the upper left q x q matrix is the correlation matrix between the q variables on the last time point; the upper right q x q matrix is the lagged correlation matrix between the q variables for the last and first time point. An example can be found below. By default, corr_YXYX = NULL.
#' @param alpha Optional. The alpha level in determining the (1-alpha)*100\% CI. By default, alpha = 0.05; resulting in a 95\% CI.
#'
#' @return This function returns i) the vectorized transformed standardized order-1 lagged effects (i.e., for DeltaTStar), both a) wave-independent and b) wave-specific ones; and ii) if N is part of input, their covariance matrix and the corresponding elliptical/multivariate 95\% CI. Additionally, it renders iii) SigmaVAR, the residual covariance matrix for DeltaTStar, and iv) Gamma, the stationary covariance matrix, also both the a) wave-independent and b) wave-specific ones.
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
#' ## 2 waves ##
#'
#' # Example with full (lagged) correlation matrix corr_lagged
#' q <- 2 # number of variables, so here we assume/use a bivariate process; say, Y1 and Y2.
#' waves <- 2 # number of waves, that is, number of time points; say, T0 and T1.
#' corr_lagged <- matrix(c(1.00, 0.41, 0.63, 0.34,
#' 0.41, 1.00, 0.31, 0.63,
#' 0.63, 0.31, 1.00, 0.40,
#' 0.34, 0.63, 0.40, 1.00), byrow = T, ncol = q*waves)
#' # corr_lagged contains the following (lagged) correlations:
#' # r(Y1_T0, Y1_T0), r(Y1_T0, Y2_T0), r(Y1_T0, Y1_T1), r(Y1_T0, Y2_T1),
#' # r(Y2_T0, Y1_T0), r(Y2_T0, Y2_T0), r(Y2_T0, Y1_T1), r(Y2_T0, Y2_T1),
#' # r(Y1_T1, Y1_T0), r(Y1_T1, Y2_T0), r(Y1_T1, Y1_T1), r(Y1_T1, Y2_T1),
#' # r(Y2_T1, Y1_T0), r(Y2_T1, Y2_T0), r(Y2_T1, Y1_T1), r(Y2_T1, Y2_T1),
#' #
#' # Run function
#' TransPhi_Corr(DeltaTStar, DeltaT, N, corr_lagged = corr_lagged)
#'
#'
#' # Example with full (lagged) correlation matrix corr_YXYX
#' corr_YXYX <- matrix(c(1.00, 0.40, 0.63, 0.34,
#' 0.40, 1.00, 0.31, 0.63,
#' 0.63, 0.31, 1.00, 0.41,
#' 0.34, 0.63, 0.41, 1.00), byrow = T, ncol = 2*2)
#' # Run function
#' TransPhi_Corr(DeltaTStar, DeltaT, N, corr_YXYX = corr_YXYX)
#'
#' # Example with vector of lower triangular correlation matrix
#' LT <- c(0.40, 0.63, 0.34, 0.31, 0.63, 0.41) # corr_YXYX[lower.tri(corr_YXYX,diag = F)]
#' # Make full correlation matrix of size 2xq times 2xq, with q=2 and thus 2xq=4 (note: waves = 2).
#' corr_YXYX <- diag(4) # As check: length(LT) = 4*(4-1)/2
#' corr_YXYX[lower.tri(corr_YXYX,diag = F)] <- LT
#' corr_YXYX[upper.tri(corr_YXYX,diag = F)] <- t(corr_YXYX)[upper.tri(t(corr_YXYX),diag = F)]
#' # Run function
#' TransPhi_Corr(DeltaTStar, DeltaT, N, corr_YXYX = corr_YXYX)
#'
#' # Example with vector of lower triangular correlation matrix including diagonals
#' LTD <- c(1.00, 0.40, 0.63, 0.34, 1.00, 0.31, 0.63, 1.00, 0.41, 1.00) # corr_YXYX[lower.tri(corr_YXYX,diag = T)]
#' # Make full correlation matrix of size 2*q times 2*q, with q=2 and thus 2*q=4
#' corr_YXYX <- matrix(NA, nrow=(4), ncol=(4))  # As check: length(LTD) = 4*(4+1)/2
#' corr_YXYX[lower.tri(corr_YXYX,diag = T)] <- LTD
#' corr_YXYX[upper.tri(corr_YXYX,diag = F)] <- t(corr_YXYX)[upper.tri(t(corr_YXYX),diag = F)]
#' # Run function
#' TransPhi_Corr(DeltaTStar, DeltaT, N, corr_YXYX = corr_YXYX)
#'
#'
#' # NOTE: The output (standPhi_DeltaTStar, standSigmaVAR_DeltaTStar, and standGamma) can be used to make stacked matrices or arrays which can serve as input for continuous-time meta-analysis CTmeta (using the function CTmeta; see ?CTmeta).
#'
#'
#' ## 3 waves ##
#'
#' # In case of waves > 2, then both the wave-independent and wave-specific parameter estimates are determined.
#'
#' # Example with full (lagged) correlation matrix corr_lagged
#' q <- 2 # number of variables, so here we assume/use a bivariate process; say, Y1 and Y2.
#' waves <- 3 # number of waves, that is, number of time points; say, T0, T1 and T2.
#' corr_lagged <- matrix(c(
#' 1.00, 0.41, 0.63, 0.34, 0.55, 0.28,
#' 0.41, 1.00, 0.31, 0.63, 0.28, 0.60,
#' 0.63, 0.31, 1.00, 0.40, 0.38, 0.45,
#' 0.34, 0.63, 0.40, 1.00, 0.40, 0.66,
#' 0.55, 0.28, 0.38, 0.40, 1.00, 0.39,
#' 0.28, 0.60, 0.45, 0.66, 0.39, 1.00
#' ), byrow = T, ncol = q*waves)
#' # corr_lagged contains the following (lagged) correlations:
#' # r(Y1_T0, Y1_T0), r(Y1_T0, Y2_T0), r(Y1_T0, Y1_T1), r(Y1_T0, Y2_T1), r(Y1_T0, Y1_T2), r(Y1_T0, Y2_T2),
#' # r(Y2_T0, Y1_T0), r(Y2_T0, Y2_T0), r(Y2_T0, Y1_T1), r(Y2_T0, Y2_T1), r(Y2_T0, Y1_T2), r(Y2_T0, Y2_T2),
#' # r(Y1_T1, Y1_T0), r(Y1_T1, Y2_T0), r(Y1_T1, Y1_T1), r(Y1_T1, Y2_T1), r(Y1_T1, Y1_T2), r(Y1_T1, Y2_T2),
#' # r(Y2_T1, Y1_T0), r(Y2_T1, Y2_T0), r(Y2_T1, Y1_T1), r(Y2_T1, Y2_T1), r(Y2_T1, Y1_T2), r(Y2_T1, Y2_T2),
#' # r(Y1_T2, Y1_T0), ...
#  # r(Y2_T2, Y1_T0), ...
#' #
#' # Run function
#' TransPhi_Corr(DeltaTStar, DeltaT, N, waves = waves, corr_lagged = corr_lagged)
#'
#'
#' # Example with full (lagged) correlation matrix corr_lagged and multiple DeltaT
#' q <- 2 # number of variables, so here we assume/use a bivariate process.
#' waves <- 3 # number of waves, that is, number of time points.
#' corr_lagged <- matrix(c(
#' 1.00, 0.41, 0.63, 0.34, 0.55, 0.28,
#' 0.41, 1.00, 0.31, 0.63, 0.28, 0.60,
#' 0.63, 0.31, 1.00, 0.40, 0.38, 0.45,
#' 0.34, 0.63, 0.40, 1.00, 0.40, 0.66,
#' 0.55, 0.28, 0.38, 0.40, 1.00, 0.39,
#' 0.28, 0.60, 0.45, 0.66, 0.39, 1.00
#' ), byrow = T, ncol = q*waves)
#' # Run function
#' DeltaT_2 <- c(24, 12)
#' TransPhi_Corr(DeltaTStar, DeltaT = DeltaT_2, N, waves = waves, corr_lagged = corr_lagged)
#'


TransPhi_Corr <- function(DeltaTStar, DeltaT = 1, N = NULL, waves = 2, corr_lagged = NULL, corr_YXYX = NULL, alpha = 0.05) {

  # Checks:
  #
  if(!is.null(waves) & length(waves) != 1){
    ErrorMessage <- (paste0("The argument 'waves' should be an integer, that is, one number, that is, a vector with one element. Currently, waves = ", waves))
    return(ErrorMessage)
    stop(ErrorMessage)
  }
  if(waves < 2){
    ErrorMessage <- (paste0("The argument 'waves' should have an integer value of at least two (to have repeated measures). Currently, waves = ", waves))
    return(ErrorMessage)
    stop(ErrorMessage)
  }
  #
  # The length of the DeltaTStar should be 1
  if(length(DeltaTStar) != 1){
    ErrorMessage <- (paste0("The argument 'DeltaTStar' should be a vector of length 1. Currently, length(DeltaTStar) = ", length(DeltaTStar)))
    return(ErrorMessage)
    stop(ErrorMessage)
  }
  # The length of DeltaT should be either 1 or (waves-1)
  if(length(DeltaT) != 1 & length(DeltaT) != (waves-1)){
    ErrorMessage <- (paste0("The argument 'DeltaT' should be a vector of length 1 or ('waves'-1). Currently, length(DeltaT) = ", length(DeltaT)))
    return(ErrorMessage)
    stop(ErrorMessage)
  }
  #
  if(!is.null(N) & length(N) != 1){
    ErrorMessage <- (paste0("The argument 'N' should be an integer, that is, one number, that is, a vector with one element. Currently, N = ", N))
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
  if(is.null(corr_lagged) & is.null(corr_YXYX)){
    ErrorMessage <- (paste0("Either the argument corr_lagged or corr_YXYX should be entered (they should not both be NULL). ",
                            "It should be a matrix of size (q x waves) times (q x waves), with q the number of variables and waves the number of time points."))
    return(ErrorMessage)
    stop(ErrorMessage)
  }


  if(is.null(corr_lagged)){

  # Check on corr_YXYX
  if(is.null(dim(corr_YXYX))){
    if(!is.null(length(corr_YXYX))){ # Should be matrix
      ErrorMessage <- (paste0("The argument 'corr_YXYX' is not a matrix. ",
                              "It should be a matrix of size (q x waves) times (q x waves), with q the number of variables and waves the number of time points."))
      return(ErrorMessage)
      stop(ErrorMessage)
    }else{
      ErrorMessage <- (paste0("The argument 'corr_YXYX' is not found: The (lagged) correlation matrix corr_YXYX is unknown, but should be part of the input."))
      return(ErrorMessage)
      stop(ErrorMessage)
    }
  }
  # Should be square matrix
    #if(dim(corr_YXYX)[1] != dim(corr_YXYX)[2] | length(dim(corr_YXYX)) != 2){
    #  ErrorMessage <- (paste0("The argument 'corr_YXYX' is not a matrix of size (q x waves) times (q x waves), with q the number of variables and waves the number of time point. Currently, it is of size ", dim(corr_YXYX)))
    #  stop(ErrorMessage)
    #}
    if(length(dim(corr_YXYX)) < 2){
      ErrorMessage <- (paste0("The argument 'corr_YXYX' is not a matrix. ",
                              "It should be a matrix of size (q x waves) times (q x waves), with q the number of variables and waves the number of time points."))
      return(ErrorMessage)
      stop(ErrorMessage)
    }
    if(length(dim(corr_YXYX)) > 2){
      ErrorMessage <- (paste0("The argument 'corr_YXYX' is not a matrix of (q x waves) times (q x waves), with q the number of variables and waves the number of time points. Currently, it is of size ", dim(corr_YXYX)))
      return(ErrorMessage)
      stop(ErrorMessage)
    }
    if(dim(corr_YXYX)[1] != dim(corr_YXYX)[2]){
      ErrorMessage <- (paste0("The argument 'corr_YXYX' should be a square matrix. It should be a matrix of size (q x waves) times (q x waves), with q the number of variables and waves the number of time points. Currently, it is of size ", dim(corr_YXYX)[1], " times ", dim(corr_YXYX)[2]))
      return(ErrorMessage)
      stop(ErrorMessage)
    }
    #
    if((dim(corr_YXYX)[1] %% 2) != 0) { # 2q should be an even number
      ErrorMessage <- (paste0("The argument 'corr_YXYX' is not a matrix of (q x waves) times (q x waves), with q the number of variables and waves the number of time points. Namely, the number of rows/columns is not even but an odd number, namely ", dim(corr_YXYX)[1]))
      return(ErrorMessage)
      stop(ErrorMessage)
    }
    if((dim(corr_YXYX)[1] %% waves) != 0) { # 2q should be an even number
      ErrorMessage <- (paste0("The argument 'corr_YXYX' is not a matrix of (q x waves) times (q x waves), with q the number of variables and waves the number of time points. Namely, the number of rows/columns is not a multiple of 'waves', namely ", dim(corr_YXYX)[1]))
      return(ErrorMessage)
      stop(ErrorMessage)
    }


    q <- dim(corr_YXYX)[1]/waves

    RYY <- NULL
    RYX <- NULL
    RXY <- NULL
    RXX <- NULL
    for(w in 2:waves){
      #w <- 2
      RYY[[w]] <- corr_lagged[(1+(w-2)*q):((w-1)*q), (1+(w-2)*q):((w-1)*q)]
      RYX[[w]] <- corr_lagged[(1+(w-2)*q):((w-1)*q), (1+(w-1)*q):(w*q)]
      RXY[[w]] <- corr_lagged[(1+(w-1)*q):(w*q), (1+(w-2)*q):((w-1)*q)]
      RXX[[w]] <- corr_lagged[(1+(w-1)*q):(w*q), (1+(w-1)*q):(w*q)]
    }


  } else{ # thus, if corr_lagged # So, this is used when both exist

    # Check on corr_lagged
    if(is.null(dim(corr_lagged))){
      if(!is.null(length(corr_lagged))){ # Should be matrix
        ErrorMessage <- (paste0("The argument 'corr_lagged' is not a matrix. ",
                                "It should be a matrix of size (q x waves) times (q x waves), with q the number of variables and waves the number of time points."))
        return(ErrorMessage)
        stop(ErrorMessage)
      }else{
        ErrorMessage <- (paste0("The argument 'corr_lagged' is not found: The (lagged) correlation matrix corr_lagged is unknown, but should be part of the input."))
        return(ErrorMessage)
        stop(ErrorMessage)
      }
    }
    # Should be square matrix
    #if(dim(corr_lagged)[1] != dim(corr_lagged)[2] | length(dim(corr_lagged)) != 2){
    #  ErrorMessage <- (paste0("The argument 'corr_lagged' is not a matrix of size (q x waves) times (q x waves), with q the number of variables and waves the number of time point. Currently, it is of size ", dim(corr_lagged)))
    #  stop(ErrorMessage)
    #}
    if(length(dim(corr_lagged)) < 2){
      ErrorMessage <- (paste0("The argument 'corr_lagged' is not a matrix. It should be a matrix of size (q x waves) times (q x waves), with q the number of variables and waves the number of time point."))
      return(ErrorMessage)
      stop(ErrorMessage)
    }
    if(length(dim(corr_lagged)) > 2){
      ErrorMessage <- (paste0("The argument 'corr_lagged' is not a matrix of size (q x waves) times (q x waves), with q the number of variables and waves the number of time point. Currently, it is of size ", dim(corr_lagged)))
      return(ErrorMessage)
      stop(ErrorMessage)
    }
    if(dim(corr_lagged)[1] != dim(corr_lagged)[2]){
      ErrorMessage <- (paste0("The argument 'corr_lagged' should be a square matrix. It should be a matrix of size (q x waves) times (q x waves), with q the number of variables and waves the number of time point. Currently, it is of size ", dim(corr_lagged)[1], " times ", dim(corr_lagged)[2]))
      return(ErrorMessage)
      stop(ErrorMessage)
    }
    #
    if((dim(corr_lagged)[1] %% 2) != 0) { # 2q should be an even number
      ErrorMessage <- (paste0("The argument 'corr_lagged' is not a matrix of size (q x waves) times (q x waves), with q the number of variables and waves the number of time point, since the number of rows/columns is not even but an odd number, namely ", dim(corr_lagged)[1]))
      return(ErrorMessage)
      stop(ErrorMessage)
    }
    if((dim(corr_lagged)[1] %% waves) != 0) { # 2q should be an even number
      ErrorMessage <- (paste0("The argument 'corr_lagged' is not a matrix of (q x waves) times (q x waves), with q the number of variables and waves the number of time points. Namely, the number of rows/columns is not a multiple of 'waves', namely ", dim(corr_lagged)[1]))
      return(ErrorMessage)
      stop(ErrorMessage)
    }


    q <- dim(corr_lagged)[1]/waves


    RYY_w <- NULL
    RYX_w <- NULL
    RXY_w <- NULL
    RXX_w <- NULL
    for(w in 2:waves){
      #w <- 2
      RXX_w[[w]] <- corr_lagged[(1+(w-2)*q):((w-1)*q), (1+(w-2)*q):((w-1)*q)]
      RXY_w[[w]] <- corr_lagged[(1+(w-2)*q):((w-1)*q), (1+(w-1)*q):(w*q)]
      RYX_w[[w]] <- corr_lagged[(1+(w-1)*q):(w*q), (1+(w-2)*q):((w-1)*q)]
      RYY_w[[w]] <- corr_lagged[(1+(w-1)*q):(w*q), (1+(w-1)*q):(w*q)]
    }
    #RXX <- corr_lagged[1:q,1:q]
    #RXY <- corr_lagged[1:q,(q+1):(q+q)]
    #RYX <- corr_lagged[(q+1):(q+q),1:q]
    #RYY <- corr_lagged[(q+1):(q+q),(q+1):(q+q)]

  }


  # Calculate wave-specific ratioDeltaT
  ratioDeltaT <- DeltaTStar / DeltaT
  if(length(ratioDeltaT) == 1 & waves > 2){
    ratioDeltaT <- rep(ratioDeltaT, (waves-1))
  }else if(length(ratioDeltaT) != (waves-1)){
    ErrorMessage <- (paste0("The arguments 'DeltaTStar' and 'DeltaT' should either have length 1 or ('waves'-1) such that ratioDeltaT = DeltaTStar / DeltaT has length 1 or ('waves'-1), currently ratioDeltaT has length", length(ratioDeltaT)))
    return(ErrorMessage)
    stop(ErrorMessage)
  }


  warning <- "No warnings (since there are no complex eigenvalues)"
  if(waves > 2){
    # Transform wave-specific Phi(DeltaT) to wave-specific Phi(DeltaTStar)
    Phi_w <- NULL
    TransPhi_w <- NULL
    vecTransPhi_w <- NULL
    SigmaVAR_w <- NULL
    Gamma_w <- NULL
    CovMx_w <- NULL
    multiCI_w <- NULL
    for(w in 2:waves){
      #w <- 2
      #
      # Calculate wave-specific Phi(DeltaT)
      Phi_w[[w]] <- t(solve(RXX_w[[w]]) %*% RXY_w[[w]])
      #
      if(ratioDeltaT[(w-1)] != 1){
        #
        # Transform wave-specific Phi(DeltaT) to wave-specific Phi(DeltaTStar)
        eigenPhi <- eigen(Phi_w[[w]])
        V <- eigenPhi$vectors
        D <- diag(eigenPhi$values)
        TransPhi_w[[w]] <- V %*% D^ratioDeltaT[(w-1)] %*% solve(V)
        #
        if(any(is.complex(eigenPhi$values))){
          if(DeltaTStar%%DeltaT == 0){
            warning <- "There is at least one pair of complex eigenvalues and the ratio of DeltaTs (i.e., DeltaT*/DeltaT) is not an integer, so the solution for Phi(DeltaT*) is NOT unique."
          }else{
            warning <- "There is at least one pair of complex eigenvalues, but the ratio of DeltaTs (i.e., DeltaT*/DeltaT) is an integer, so the solution for Phi(DeltaT*) is unique."
          }
          if(all(Im(TransPhi) == 0)){
            TransPhi_w[[w]] <- Re(TransPhi_w[[w]])
          }
        }
        #
        RXY_w[[w]] <- RXX_w[[w]] %*% t(TransPhi_w[[w]])
        RYX_w[[w]] <- t(RXY_w[[w]])
        # RYY and RXX remain the same!
      } else{
        TransPhi_w[[w]] <- Phi_w[[w]]
        #
        warning <- "No warnings (since there are no complex eigenvalues)"
      }
      vecTransPhi_w[[w]] <- as.vector(t(TransPhi_w[[w]]))

      SigmaVAR_w[[w]] <- RYY_w[[w]] - RYX_w[[w]]  %*%  solve(RXX_w[[w]]) %*% RXY_w[[w]]
      Gamma_w[[w]] <- RXX_w[[w]]


      if(!is.null(N)){

        CovMx_w[[w]] <- kronecker(SigmaVAR_w[[w]], solve(Gamma_w[[w]])) / (N-q)

        # Determine points on 95% LL contour
        mu_Phi <- vecTransPhi_w[[w]]
        CovMx_Phi <- CovMx_w[[w]]
        eigenCovMx <- eigen(CovMx_Phi)
        lambda <- eigenCovMx$val # elt.wise[(eigenvalSigmaVAR / eigenvalGamma)] / (N-q)
        if(any(is.complex(lambda))){
          if(warning == "No warnings (since there are no complex eigenvalues)"){
            warning <- NULL
          }
          warning <- c(warning,
                       paste("Some of the eigenvalues of the covariance matrix of Phi are complex (i.e., not real).",
                             "The corresponding confidence interval(s), in $multiCI_vecStandPhi_DeltaT, are NA.")
          )
          which_lambda <- which(Im(eigenCovMx$val) != 0)
          lambda[which_lambda] <- 0
          lambda <- Re(lambda)
        }
        if(any(lambda < 0)){
          if(warning == "No warnings (since there are no complex eigenvalues)"){
            warning <- NULL
          }
          warning <- c(warning,
                       paste("Some of the eigenvalues of the covariance matrix of Phi are negative.",
                             "The corresponding confidence interval(s), in $multiCI_vecStandPhi_DeltaT, are NA.")
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
        teller <- 0
        for(row in 1:q){
          for(column in 1:q){
            teller <- teller + 1
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
        multiCI_w[[w]] <- rbind(minPhi, maxPhi)
        rownames(multiCI) <- c("LB", "UB")
        sub = NULL
        for(i in 1:q){
          sub = c(sub, paste0("Phi", i, 1:q, sep=""))
        }
        colnames(multiCI) <- sub
      }
    }


    # Determine wave-independent Phi based on averaging correlation matrices
    RYY <- apply(simplify2array(RYY_w[-1]), c(1,2), mean)
    RYX <- apply(simplify2array(RYX_w[-1]), c(1,2), mean)
    RXX <- apply(simplify2array(RXX_w[-1]), c(1,2), mean)
    RXY <- apply(simplify2array(RXY_w[-1]), c(1,2), mean)
  }else if (waves == 2){
    RYY <- RYY_w[[2]]
    RYX <- RYX_w[[2]]
    RXX <- RXX_w[[2]]
    RXY <- RXY_w[[2]]
  }

  # (vectorized) transformed Phi, so (vectorized) Phi(DeltaTStar)
  TransPhi <- t(solve(RXX) %*% RXY)
  vecTransPhi <- as.vector(t(TransPhi))

  SigmaVAR <- RYY - RYX  %*%  solve(RXX) %*% RXY
  Gamma <- RXX


  if(!is.null(N)){

    CovMx <- kronecker(SigmaVAR, solve(Gamma)) / (N-q)

    # Determine points on 95% LL contour
    mu_Phi <- vecTransPhi
    CovMx_Phi <- CovMx
    eigenCovMx <- eigen(CovMx_Phi)
    lambda <- eigenCovMx$val # elt.wise[(eigenvalSigmaVAR / eigenvalGamma)] / (N-q)
    if(any(is.complex(lambda))){
      if(warning == "No warnings (since there are no complex eigenvalues)"){
        warning <- NULL
      }
      warning <- c(warning,
                   paste("Some of the eigenvalues of the covariance matrix of Phi are complex (i.e., not real).",
                         "The corresponding confidence interval(s), in $multiCI_vecStandPhi_DeltaT, are NA.")
      )
      which_lambda <- which(Im(eigenCovMx$val) != 0)
      lambda[which_lambda] <- 0
      lambda <- Re(lambda)
    }
    if(any(lambda < 0)){
      if(warning == "No warnings (since there are no complex eigenvalues)"){
        warning <- NULL
      }
      warning <- c(warning,
                   paste("Some of the eigenvalues of the covariance matrix of Phi are negative.",
                         "The corresponding confidence interval(s), in $multiCI_vecStandPhi_DeltaT, are NA.")
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
    teller <- 0
    for(row in 1:q){
      for(column in 1:q){
        teller <- teller + 1
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

    if(waves == 2){
      if(!is.null(N)){
        final <- list(DeltaTStar = DeltaTStar, DeltaT = DeltaT,
                      nr_variables = q, nr_waves = waves, N = N,
                      alpha = alpha,
                      #corr_lagged = corr_lagged, corr_YXYX = corr_YXYX,
                      #
                      standPhi_DeltaTStar = TransPhi, vecStandPhi_DeltaTStar = vecTransPhi,
                      CovMx_vecStandPhi_DeltaTStar = CovMx,
                      multiCI_vecStandPhi_DeltaTStar = multiCI,
                      standSigmaVAR_DeltaTStar = SigmaVAR, standGamma = Gamma,
                      #
                      warning = warning)
      } else{
        final <- list(DeltaTStar = DeltaTStar, DeltaT = DeltaT,
                      nr_variables = q, nr_waves = waves,
                      #corr_lagged = corr_lagged, corr_YXYX = corr_YXYX,
                      #
                      standPhi_DeltaTStar = TransPhi, vecStandPhi_DeltaTStar = vecTransPhi,
                      standSigmaVAR_DeltaTStar = SigmaVAR, standGamma = Gamma,
                      #
                      warning = warning)
      }
    } else if(waves > 2){
      if(!is.null(N)){
        final <- list(DeltaTStar = DeltaTStar, DeltaT = DeltaT,
                      nr_variables = q, nr_waves = waves, N = N,
                      alpha = alpha,
                      #corr_lagged = corr_lagged, corr_YXYX = corr_YXYX,
                      #
                      standPhi_DeltaTStar = TransPhi, vecStandPhi_DeltaTStar = vecTransPhi,
                      CovMx_vecStandPhi_DeltaTStar = CovMx,
                      multiCI_vecStandPhi_DeltaTStar = multiCI,
                      standSigmaVAR_DeltaTStar = SigmaVAR, standGamma = Gamma,
                      #
                      standPhi_DeltaTStar_w = TransPhi_w[-1], vecStandPhi_DeltaTStar_w = vecTransPhi_w[-1],
                      CovMx_vecStandPhi_DeltaTStar_w = CovMx_w[-1],
                      multiCI_vecStandPhi_DeltaTStar_w = multiCI_w[-1],
                      standSigmaVAR_DeltaTStar_w = SigmaVAR_w[-1], standGamma_w = Gamma_w[-1],
                      #
                      warning = warning)
      } else{
        final <- list(DeltaTStar = DeltaTStar, DeltaT = DeltaT,
                      nr_variables = q, nr_waves = waves,
                      #corr_lagged = corr_lagged, corr_YXYX = corr_YXYX,
                      #
                      standPhi_DeltaTStar = TransPhi, vecStandPhi_DeltaTStar = vecTransPhi,
                      standSigmaVAR_DeltaTStar = SigmaVAR, standGamma = Gamma,
                      #
                      standPhi_DeltaTStar_w = TransPhi_w[-1], vecStandPhi_DeltaTStar_w = vecTransPhi_w[-1],
                      standSigmaVAR_DeltaTStar_w = SigmaVAR_w[-1], standGamma_w = Gamma_w[-1],
                      #
                      warning = warning)
      }
    }

    return(final)
  }

