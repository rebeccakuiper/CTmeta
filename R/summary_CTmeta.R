#' @importFrom jtools md_table
#' @exportS3Method summary CTmeta
#' @export summary.CTmeta
#' @export

summary.CTmeta <- function(x, digits = NULL)
{

  x <- as.list(x)

  if(is.null(digits)){
    NrDigits <- 3
    sig.digits <- TRUE
    align <- NULL
  }else{
    NrDigits <- digits
    sig.digits <- FALSE
    align <- 'c'
  }

  DF <- data.frame(Overall_Phi = x$Overall_vecStandPhi,
                   LB = x$LB,
                   UB = x$UB)

  md_table(DF, digits = NrDigits, sig.digits = sig.digits, align = align)

}
