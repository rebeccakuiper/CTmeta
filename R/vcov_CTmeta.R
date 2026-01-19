#' @exportS3Method vcov CTmeta
#' @export vcov.CTmeta
#' @export

vcov.CTmeta <- function(x)
{

  x <- as.list(x)

  if(!is.null(x$CovMx)){
    VCOV <- x$CovMx
    message <- NULL
  }else{
    VCOV <- diag(x$se)
    message <- "VCOV is a diagonal matrix, since only the standard errors are known"
  }

  final <- list(VCOV = VCOV, message = message)
  #return(final)
  #VCOV
  #return(invisible(final))
  return(VCOV)

}
