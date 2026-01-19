#' @exportS3Method print CTmeta
#' @export print.CTmeta
#' @export

print.CTmeta <- function(x, digits = NULL)
{

  x <- as.list(x)

  if(is.null(digits)){
    NrDigits <- 3
  }else{
    NrDigits <- digits
  }

  DF <- data.frame(Overall_Phi = x$Overall_vecStandPhi,
                   LB = x$LB,
                   UB = x$UB)

  cat("\n")
  if(!is.null(x$LB_elliptical)){
    cat(paste0("The overall estimates obtained with CTmeta and their ", (1-x$alpha)*100, "% elliptical/multivariate confidence interval: \n"))
  }else{
    cat(paste0("The overall estimates obtained with CTmeta and their ", (1-x$alpha)*100, "% (univariate) confidence interval: \n"))
  }
  cat("\n")
  print(DF, digits = NrDigits, right = F)
  cat("\n")

  if(!is.null(x$tau2)){
    cat(paste0("Note: A random-effects model is used. The tau^2 values can be obtained via '$tau2'; see '$summaryMetaAnalysis' for more information."))
    cat("\n")
    cat("\n")
  }

  cat("CTmeta Messages: \n")
  cat("- ")
  cat(x$messageTrans)
  cat("\n")
  cat("- ")
  cat(x$messageMultivar)
  cat("\n")

  if(!is.null(x$StudiesComplexEV)){
    cat("\n")
    cat("The studies with complex eigenvalues: \n")
    cat(x$StudiesComplexEV)
    cat("\n")
  }
  if(!is.null(x$StudiesNegEV)) {
    cat("\n")
    cat("The studies with negative (real parts of the) eigenvalues: \n")
    cat(x$StudiesNegEV)
    cat("\n")
  }
  if(!is.null(x$StudiesCovMxNotPosDef)) {
    cat("\n")
    cat("The studies with not-positive definite (i.e., negative semidefinite) covariance matrices: \n")
    cat(x$StudiesCovMxNotPosDef)
    cat("\n")
  }




  return(invisible(x))

}

