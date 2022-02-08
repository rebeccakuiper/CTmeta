#' Area under the curves in the Phi(DeltaT)-plot
#'
#' Area under the curves in the Phi(DeltaT)-plot. Area is taken for time intervals between 0 and infinity and optionally for a user-specified range. The interactive web application 'Phi-and-Psi-Plots and Find DeltaT' also contains this functionality: \url{https://www.uu.nl/staff/RMKuiper/Websites\%20\%2F\%20Shiny\%20apps}.
#'
#' @param DeltaT Optional. The time interval used (only of interest if Phi is also part of the input). By default, DeltaT = 1.
#' @param Phi Matrix of size q x q of (un)standardized lagged effects.
#' Can also be a fitted object of the class "varest" (from the VAR() function in vars package) or "ctsemFit" (from the ctFit() function in the ctsem package); see example below. The (standardized) Drift matrix is calculated/extracted from these objects.
#' @param Drift Optional (either Phi or Drift). Matrix of size q x q of (un)standardized continuous-time lagged effects (called the drift matrix). Note that Phi(DeltaT) = expm(Drift*DeltaT). If both Phi and Drift are specified, Phi is used and Drift is ignored.
#' @param t_min Optional. The time from which the area is taken into account. By default, t_min = 0.
#' @param t_max Optional. The time until which the area is taken into account. By default, t_max = Inf.
#'
#' @return For each (i,j), the output renders the area under the curve for Phi_ij. Area is the total area under the curve; and Area_range is the area under the specifed range of time.
#' @importFrom expm expm
#' @export
#' @examples
#'
#' # library(CTmeta)
#'
#' ## Example 1 ##
#'
#' ##################################################################################################
#' # Input needed in the examples below; q=2 variables.
#' # We use the example matrices stored in the package:
#' # specifying Phi and not Drift
#' DeltaT <- 1
#' Phi <- myPhi[1:2, 1:2]
#' # specifying Drift and not Phi
#' DeltaT <- 1
#' Drift <- myDrift
#' # note that in this example, Phi and Drift are taken from the same dataset and thus specifying one or the other gives the same results
#' ##################################################################################################
#'
#' # Input: Phi
#' Area(DeltaT = DeltaT, Phi = Phi)
#' # or
#' Area(DeltaT, Phi)
#' # or, since DeltaT = 1, as per default:
#' Area(Phi = Phi)
#'
#' # Input: Drift
#' Area(DeltaT, Drift = Drift)
#'
#' # If inspecting the time interval from e.g. 1 to 2 (instead of 0 to infinity):
#' Area(DeltaT, Phi, t_min = 1, t_max = 2)
#'
#'
#' # The function 'PhiPlot' can help for visualization of the curves in the Phi(DeltaT)-plot.
#' PhiPlot(DeltaT, Phi)
#'
#'
#' ## Example 2: input from fitted object of class "varest" ##
#' #
#' data <- myData
#' if (!require("vars")) install.packages("vars")
#' library(vars)
#' 
#' out_VAR <- VAR(data, p = 1)
#' DeltaT <- 1
#' 
#' Area(DeltaT = DeltaT, Phi = out_VAR)
#'
#' ## Example 3: input from fitted object of class "CTsemFit" ##
#' #
#' library(ctsem)
#' library(ctsemOMX)
#' #
#' ############ adapted from https://rdrr.io/cran/ctsemOMX/man/ctFit.html ############
#' data(ctExample1)
#' model <- ctModel(n.manifest=2, n.latent=2, Tpoints=6, LAMBDA=diag(2),
#'                  manifestNames=c('LeisureTime', 'Happiness'),
#'                  latentNames=c('LeisureTime', 'Happiness'), TRAITVAR="auto")
#' out_Phi <- ctFit(dat=ctExample1, ctmodelobj=model)
#' ##################################################################################
#' #
#' Area(Phi = out_Phi)
#' Area(Phi = out_Phi, t_min = 1, t_max = 2)
#' PhiPlot(DeltaT = 1, out_Phi)
#' #
#' # Four plots are produced since the eigenvalues of the Drift matrix are complex, so there are multiple solutions for Phi.
#' 

Area <- function(DeltaT = 1, Phi = NULL, Drift = NULL, t_min = 0, t_max = Inf) {

  # Checks:
  if(!is.null(Drift) & DeltaT != 1) {
    warning("The input for DeltaT was ignored as Drift was part of the input.")
    DeltaT <- 1 # this way the function can run even if DeltaT is not valid input
  }
  #
  if(length(DeltaT) != 1 | !is.numeric(DeltaT)){
    ErrorMessage <- (paste0("The argument DeltaT should be a scalar (i.e., one number or a vector with one numerical element)."))
    stop(ErrorMessage)
  }
  if(!(DeltaT > 0)) {
    stop("DeltaT must be a positive number.")
  }
  if(DeltaT == Inf) {
    stop("DeltaT must be a finite number.")
  }
  #
  if(!is.null(Phi) & !is.null(Drift)) {
    warning("Both Phi and Drift are specified. Drift is ignored.")
  }
  #
  if(!is.numeric(Drift) & !is.null(Drift)) {
    stop("Drift contains non-numerical values.")
  }
  #
  if(length(t_min) > 1 | !is.numeric(t_min)) {
    stop("t_min should be a scalar (i.e., a real number or a vector with one numerical element).")
  }
  if(!(t_min >= 0)) {
    stop("t_min must be >= 0.")
  }
  #
  if(length(t_max) > 1 | !is.numeric(t_max)) {
    stop("t_max should be a number.")
  }
  if (!(t_max > 0)) {
    stop("t_max must be strictly positive.")
  }
  if(t_min > t_max) {
    warning("t_min is larger than t_max. The values in Area_range will be negative.")
  }
  # Check on Phi
  if(anyNA(Phi) | anyNA(Drift)) {
    stop("There are missing values in Phi or Drift.")
  }
  if(!is.numeric(Phi) & !is.null(Phi) & !(any(class(Phi) %in% c("varest", "ctsemFit")))) {
    stop("Phi contains non-numerical values.")
  }
  if(any(class(Phi) == "varest")){
    Phi <- Acoef(Phi)[[1]]
    CTMp <- CTMparam(DeltaT, Phi)
    if(is.null(CTMp$ErrorMessage)){
      B <- -CTMp$Drift  # Drift <- logm(Phi)/DeltaT  # Phi <- expm(Drift * DeltaT)
    }else{
      ErrorMessage <- CTMp$ErrorMessage
      stop(ErrorMessage)
    }
  } else if(any(class(Phi) == "ctsemFit")){
    B <- -1 * summary(Phi)$DRIFT
  } else{

    # Drift = A = -B
    # B is drift matrix that is pos def, so Phi(DeltaT) = expm(-B*DeltaT)
    if(is.null(Drift)){
      if(!is.null(Phi)){
        CTMp <- CTMparam(DeltaT, Phi)
        if(is.null(CTMp$ErrorMessage)){
          B <- -CTMp$Drift  # Drift <- logm(Phi)/DeltaT  # Phi <- expm(Drift * DeltaT)
        }else{
          ErrorMessage <- CTMp$ErrorMessage
          stop(ErrorMessage)
        }
      }else{ # is.null(Phi)
        ErrorMessage <- ("Either the drift matrix Drift or the autoregressive matrix Phi should be part of the input.")
        #("Note that Phi(DeltaT) = expm(-B*DeltaT).")
        stop(ErrorMessage)
      }
    }else{ # !is.null(Drift)
      B <- -Drift
    }
    # Check on B
    if(length(B) > 1){
      Check_B_or_Phi(B)
      if(all(Re(eigen(B)$val) < 0)){
        cat("All (the real parts of) the eigenvalues of the drift matrix Drift are positive. Therefore, it is assumed the input for Drift was B = -A instead of A (or -Phi instead of Phi). Drift = -B = A is used.")
        cat(" Note that Phi(DeltaT) = expm(-B*DeltaT) = expm(A*DeltaT) = expm(Drift*DeltaT).")
        B = -B
      }
      if(any(Re(eigen(B)$val) < 0)){
        #ErrorMessage <- ("The function stopped, since some of (the real parts of) the eigenvalues of the drift matrix Drift are positive.")
        #return(ErrorMessage)
        #return(ErrorMessage)
        #stop(ErrorMessage)
        cat("If the function stopped, this is because some of (the real parts of) the eigenvalues of the drift matrix Drift are positive.")
      }
    }
  }
  #
  if(length(B) == 1){
    q <- 1
  }else{
    q <- dim(B)[1]
  }


  Area_full <- solve(B)
  if(t_max == Inf){
    expB_max <- matrix(0, nrow = q, ncol = q)
  }else{
    expB_max <- expm(-B*t_max)
  }
  Area_range <- - solve(B) %*% ( expB_max - expm(-B*t_min) )


  ############################################################################################################

  final <- list(Area = Area_full,
                Area_range = Area_range)

  return(final)

}
