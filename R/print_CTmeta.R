#' @S3method print CTmeta
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
  if(!is.null(x$tau2)){
    cat(paste0("An random-effects model is used and tau2 = ", x$tau2, ". See '$summaryMetaAnalysis' for more information."))
  }

  cat("\n")
  cat("CTmeta Messages: \n")
  cat("- ")
  cat(x$messageTrans)
  cat("\n")
  cat("- ")
  cat(x$messageMultivar)

  if(!is.null(x$StudiesComplexEV)){
    cat("\n")
    cat(x$StudiesComplexEV)
  }
  if(!is.null(x$StudiesNegEV)) {
    cat("\n")
    cat(x$StudiesNegEV)
  }
  if(!is.null(x$StudiesCovMxNotPosDef)) {
    cat("\n")
    cat(x$StudiesCovMxNotPosDef)
  }


  return(invisible(x))

}

