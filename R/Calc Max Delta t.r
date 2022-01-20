#' Time-interval (DeltaT) for which Phi_ij(DeltaT) reaches its minimum or maximum (together with that minimum or maximum)
#'
#' Time-interval (DeltaT) for which Phi_ij(DeltaT) reaches its minimum or maximum (together with that minimum or maximum). The interactive web application 'Phi-and-Psi-Plots and Find DeltaT' also contains this functionality, you can find it on my website: \url{https://www.uu.nl/staff/RMKuiper/Websites\%20\%2F\%20Shiny\%20apps}.
#'
#' @param DeltaT Optional. The time interval used. By default, DeltaT = 1.
#' @param Phi Matrix of size q times q of (un)standardized lagged effects.
#' It also takes a fitted object from the classes "varest" (from the VAR() function in vars package) and "ctsemFit" (from the ctFit() function in the ctsem package); see example below. From such an object, the (standardized) Drift matrix is calculated/extracted.
#' @param Drift Optional (either Phi or Drift). Matrix of size q times q of (un)standardized continuous-time lagged effects, called drift matrix. Note that Phi(DeltaT) = expm(Drift*DeltaT). By default, input for Phi is used; only when Phi = NULL, Drift will be used.
#'
#' @return The output renders, per element (i,j), the time-interval for which Phi_ij reaches its minimum/maximum together with this minimum/maximum Phi_ij. Note that even though a matrix is presented, the elements in it refer to different time-intervals when DeltaT differs per element (see in examples below).
#' @importFrom expm expm
#' @importFrom nleqslv nleqslv
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
#' Drift <- myDrift
#' ##################################################################################################
#'
#' MaxDeltaT(DeltaT = DeltaT, Phi = Phi)
#' # or
#' MaxDeltaT(DeltaT, Phi)
#'
#' # Note that the DeltaT for which Phi_ij reaches its maximum or minimum ('DeltaT_MinOrMaxPhi') differs per Phi_ij.
#' # Therefore, the matrix 'MinOrMaxPhi' is not a Phi-matrix, but each element should be inspected separately.
#' # To obtain the full Phi-matrix for a specific DeltaT one can use:
#' DeltaT_MinOrMaxPhi <- MaxDeltaT(DeltaT, Phi)$DeltaT_MinOrMaxPhi
#' StandTransPhi(DeltaTStar = DeltaT_MinOrMaxPhi[1,2], DeltaT, N = NULL, Phi)
#'
#' # If you would use the drift matrix Drift as input, then use:
#' MaxDeltaT(DeltaT, Drift = Drift)
#'
#'
#' # Note that the function 'PhiPlot' can help to see (per element) whether a minimum or maximum is reached.
#' PhiPlot(DeltaT, Phi)
#' # or:
#' ggPhiPlot(DeltaT, Phi)
#'
#'
#' ## Example 2: input from fitted object of class "varest" ##
#' #
#' data <- myData
#' if (!require("vars")) install.packages("vars")
#' library(vars)
#' out_VAR <- VAR(data, p = 1)
#' DeltaT <- 1
#' MaxDeltaT(DeltaT, out_VAR)
#' #
#' ggPhiPlot(DeltaT, out_VAR)
#'

MaxDeltaT <- function(DeltaT = 1, Phi = NULL, Drift = NULL) {

  # Checks:
  if(length(DeltaT) != 1){
    ErrorMessage <- (paste0("The argument DeltaT should be a scalar (i.e., one number or a vector with one element."))
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
        cat("All (the real parts of) the eigenvalues of the drift matrix Drift are positive. Therefore. I assume the input for Drift was B = -A instead of A (or -Phi instead of Phi). I will use Drift = -B = A.")
        cat("Note that Phi(DeltaT) = expm(-B*DeltaT) = expm(A*DeltaT) = expm(Drift*DeltaT).")
        B = -B
      }
      if(any(Re(eigen(B)$val) < 0)){
        #ErrorMessage <- ("The function stopped, since some of (the real parts of) the eigenvalues of the drift matrix Drift are positive.")
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


  message <- "There is a DeltaT such that the Phi(DeltaT) functions reach a minimum or maximum."


  SolveForMaxDelta_ij_fie <- function(q, B, i, j) {
    function(x) {
      y <- numeric(1)
      y <- -((-B) %*% expm(-B*x))[i,j]
    }
  }

  MaxDeltaT_mx <- matrix(NA, nrow = q, ncol = q)
  # If single i and j
  #i <- 1
  #j <- 2
  # If loop
  for(i in 1:q){
    for(j in 1:q){
      #
      SolveForMaxDelta_ij <- SolveForMaxDelta_ij_fie(q, B, i, j)
      #
      #
      xstart_ij <- 1
      fstart_ij <- SolveForMaxDelta_ij(xstart_ij)
      # Check
      if(all(is.nan(fstart_ij) == FALSE) == FALSE){
        ("There is no DeltaT such that the Phi(DeltaT) functions reach a minimum or maximum!")
        message <- "There is no DeltaT such that the Phi(DeltaT) functions reach a minimum or maximum."
      }
      #
      #
      #
      sol_ij <- nleqslv(xstart_ij, SolveForMaxDelta_ij, control=list(btol=.0000001, allowSingular = TRUE)) #, method="Newton")
      MaxDeltaT_ij <- sol_ij$x
      #MaxDeltaT_ij
      # Function terminated if sol$termcd does not equal 1 (and 1 is "Function criterion is near zero. Convergence of function values has been achieved.")
      sol_ij$message
      sol_ij$fvec
      if(sol_ij$termcd != 1){
        cat("The nleqslv-function terminated.")
        cat("Hence, there is no DeltaT such that the Phi(DeltaT) functions reach a minimum or maximum.")
        message <- "There is no DeltaT such that the Phi(DeltaT) functions reach a minimum or maximum."
        #stop(message)
      }
      #
      #
      MaxDeltaT_mx[i,j] <- MaxDeltaT_ij
    } #end loop j
  } #end loop i

  ## For a single i and j
  #MaxDeltaT_ij
  #cat("Phi(Max DeltaT)_ij = ")
  #cat((expm(-B*MaxDeltaT_mx[i,j]))[i,j])
  #cat("Check on zero for 1st order deriv = ")
  #cat((B %*% expm(-B*MaxDeltaT_mx[i,j]))[i,j])
  #cat("")

  # If loop is in place
  MaxDeltaT_mx
  Phi_MaxDeltaT_mx <- matrix(NA, nrow = q, ncol = q)
  for(i in 1:q){
    for(j in 1:q){
      Phi_MaxDeltaT_mx[i,j] <- (expm(-B*MaxDeltaT_mx[i,j]))[i,j]
      #cat("Phi(Max DeltaT)_ij = ")
      #cat((expm(-B*MaxDeltaT_mx[i,j]))[i,j])
      #cat("Check on zero for 1st order deriv = ")
      #cat((B %*% expm(-B*MaxDeltaT_mx[i,j]))[i,j])
      #cat("")
    }
  }

  ############################################################################################################

  final <- list(DeltaT_MinOrMaxPhi = MaxDeltaT_mx,
                MinOrMaxPhi = Phi_MaxDeltaT_mx)

  return(final)

}
