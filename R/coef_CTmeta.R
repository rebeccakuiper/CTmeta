#' @S3method coef CTmeta
#' @export coef.CTmeta
#' @export

coef.CTmeta <- function(x)
{

  x <- as.list(x)

  est <- x$Overall_vec

  return(est)

}
