#' Area under the curves in the Phi(DeltaT)-plot
#'
#' Area under the curves in the Phi(DeltaT)-plot, for time-intervals from 0 to infinity and optionally for a user-specified range. The interactive web application 'Phi-and-Psi-Plots and Find DeltaT' also contains this functionality, you can find it on my website: \url{https://www.uu.nl/staff/RMKuiper/Websites\%20\%2F\%20Shiny\%20apps}.
#'
#' @param DeltaT Optional. The time interval used (only of interest if Phi is part of input). By default, DeltaT = 1.
#' @param Phi Matrix of size q times q of (un)standardized lagged effects.
#' It also takes a fitted object from the classes "varest" (from the VAR() function in vars package) and "ctsemFit" (from the ctFit() function in the ctsem package); see example below. From such an object, the (standardized) Drift matrix is calculated/extracted.
#' @param Drift Optional (either Phi or Drift). Matrix of size q times q of (un)standardized continuous-time lagged effects, called drift matrix. Note that Phi(DeltaT) = expm(Drift*DeltaT). By default, input for Phi is used; only when Phi = NULL, Drift will be used.
#'
#' @return The output renders, per element (i,j), the area under the curve for Phi_ij.
#' @importFrom expm expm
#' @export
#' @examples
#'
#' # library(CTmeta)
#'
#' ## Example 1 ##
#'
#' ##################################################################################################
#' # Input needed in examples below with q=2 variables.
#' # I will use the example matrices stored in the package:
#' DeltaT <- 1
#' Phi <- myPhi[1:2, 1:2]
#' # or: Drift
#' DeltaT <- 1
#' Drift <- myDrift
#' ##################################################################################################
#'
#' Area(DeltaT = DeltaT, Phi = Phi)
#' # or
#' Area(DeltaT, Phi)
#' # or, since DeltaT = 1
#' Area(Phi = Phi)
#'
#' # If you would use the drift matrix Drift as input, then use:
#' Area(DeltaT, Drift = Drift)
#'
#' # If, for instance, the time-interval range from 1 to 2 should be inspected (and not 0 to infinity), then use:
#' Area(DeltaT, Phi, t_min = 1, t_max = 2)
#'
#'
#' # Note that the function 'PhiPlot' can help for visualization of the curves in the Phi(DeltaT)-plot.
#' PhiPlot(DeltaT, Phi)
#'
#'
#' ## Example 2: input from fitted object of class "varest" ##
#' #
#' data <- myData
#' if (!require("vars")) install.packages("vars")
#' library(vars)
#' out_VAR <- VAR(data, p = 1)
#' DeltaT <- 1
#' Area(DeltaT, out_VAR)
#'

Area <- function(DeltaT = 1, Phi = NULL, Drift = NULL, t_min = 0, t_max = "inf") {

  # Checks:
  if(length(DeltaT) != 1){
    ErrorMessage <- (paste0("The argument DeltaT should be a scalar, that is, one number, that is, a vector with one element. Currently, DeltaT = ", DeltaT))
    return(ErrorMessage)
    stop(ErrorMessage)
  }
  #
  # Check on Phi
  if(any(class(Phi) == "varest")){
    Phi <- Acoef(Phi)[[1]]
    CTMp <- CTMparam(DeltaT, Phi)
    if(is.null(CTMp$ErrorMessage)){
      B <- -CTMp$Drift  # Drift <- logm(Phi)/DeltaT  # Phi <- expm(Drift * DeltaT)
    }else{
      return(ErrorMessage)
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
          return(ErrorMessage)
          stop(ErrorMessage)
        }
      }else{ # is.null(Phi)
        ErrorMessage <- ("Either the drift matrix Drift or the autoregressive matrix Phi should be input to the function.")
        #("Note that Phi(DeltaT) = expm(-B*DeltaT).")
        return(ErrorMessage)
        stop(ErrorMessage)
      }
    }else{ # !is.null(Drift)
      B <- -Drift
    }
    # Check on B
    if(length(B) > 1){
      Check_B_or_Phi(B)
      if(all(Re(eigen(B)$val) < 0)){
        cat("All (the real parts of) the eigenvalues of the drift matrix Drift are positive. Therefore. I assume the input for Drift was B = -A instead of A (or -Phi instead of Phi). I will use Drift = -B = A.")
        cat("Note that Phi(DeltaT) = expm(-B*DeltaT) = expm(A*DeltaT) = expm(Drift*DeltaT).")
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
  if(t_max == "inf"){
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
