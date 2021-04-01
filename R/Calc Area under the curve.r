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
#' @importFrom expm logm
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
#' Phi <- myPhi[1:2, 1:2]
#'
#' # If you would use the drift matrix Drift as input, then use:
#' if (!require("expm")) install.packages("expm") # Use expm package for function logm()
#' library(expm)
#' Drift <- logm(Phi)/DeltaT
#' ##################################################################################################
#'
#' DeltaT <- 1
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
    print(paste0("The argument DeltaT should be a scalar, that is, one number, that is, a vector with one element. Currently, DeltaT = ", DeltaT))
    stop()
  }
  #
  # Check on Phi
  if(any(class(Phi) == "varest")){
    Phi <- Acoef(Phi)[[1]]
    B <- -logm(Phi)/DeltaT # Phi = expm(Drift * DeltaT)
  } else if(any(class(Phi) == "ctsemFit")){
    B <- -1 * summary(Phi)$DRIFT
  } else{

    # Drift = A = -B
    # B is drift matrix that is pos def, so Phi(DeltaT) = expm(-B*DeltaT)
    if(is.null(Drift)){
      if(!is.null(Phi)){
        B <- -logm(Phi)/DeltaT
      }else{ # is.null(Phi)
        ("Either the drift matrix Drift or the autoregressive matrix Phi should be input to the function.")
        #("Note that Phi(DeltaT) = expm(-B*DeltaT).")
        stop()
      }
    }else{ # !is.null(Drift)
      B <- -Drift
    }
    # Check on B
    if(length(B) > 1){
      Check_B_or_Phi(B)
      if(all(eigen(B)$val < 0)){
        #("All the eigenvalues of the drift matrix B are negative; therefore. I assume the input was A=-B instead of B. I will use -A=B in the calculation.")
        #("Note that Phi(DeltaT) = expm(-B*DeltaT).")
        ("All the eigenvalues of the drift matrix Drift are positive. Therefore. I assume the input for Drift was B = -A instead of A. I will use Drift = -B = A.")
        ("Note that Phi(DeltaT) = expm(-B*DeltaT) = expm(A*DeltaT) = expm(Drift*DeltaT).")
        B = -B
      }
      if(any(eigen(B)$val <= 0)){
        #("The function stopped, since some of the eigenvalues of the drift matrix B are negative or zero.")
        ("The function stopped, since some of the eigenvalues of the drift matrix Drift are positive or zero.")
        stop()
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
