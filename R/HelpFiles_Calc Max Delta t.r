#' Time-interval (DeltaT) for which Phi_ij(DeltaT) reaches its minimum or maximum (together with that minimum or maximum)
#'
#' Time-interval (DeltaT) for which Phi_ij(DeltaT) reaches its minimum or maximum (together with that minimum or maximum). The interactive web application 'Phi-and-Psi-Plots and Find DeltaT' also contains this functionality, you can find it on my website: \url{https://www.uu.nl/staff/RMKuiper/Websites\%20\%2F\%20Shiny\%20apps}.
#'
#' @param Drift Matrix of size q times q of (un)standardized continuous-time lagged effects, called drift matrix. Note that Phi(DeltaT) = expm(Drift*DeltaT).
#' @param Phi Optional. Matrix of size q times q of (un)standardized lagged effects. By default, input for Drift is used; only when Drift = NULL, Phi will be used (to determine the corresponding Drift).
#'
#' @return The output renders, per element (i,j), the time-interval for which Phi_ij reaches its minimum/maximum together with this minimum/maximum.
#' @importFrom expm expm
#' @importFrom expm logm
#' @importFrom nleqslv nleqslv
#' @export
#' @examples
#'
#' ##################################################################################################
#' # Input needed in examples below with q=2 variables.
#' # I will use the example matrices stored in the package:
#' Phi <- myPhi[1:2, 1:2]
#' ##################################################################################################
#'
#' calc.MaxDeltaT(Phi = Phi)
#'
#' # If you would use the drift matrix Drift as input, then use:
#' ##if (!require("expm")) install.packages("expm") # Use expm package for function logm()
#' ##library(expm)
#' ##Drift <- logm(Phi)/DeltaT
#' #calc.MaxDeltaT(Drift = Drift)
#' #calc.MaxDeltaT(Drift)
#'
#'
#' # Note that the function 'PhiPlot' can help to see (per element) whether a minimum or maximum is reached.
#' ##if (!require("expm")) install.packages("expm") # Use expm package for function logm()
#' ##library(expm)
#' ##Drift <- logm(Phi)/DeltaT
#' #PhiPlot(DeltaT = 1, Drift)
#'

calc.MaxDeltaT <- function(Drift = NULL, Phi = NULL) {

  #if (!require("expm")) install.packages("expm")
  #library(expm)
  #if (!require("nleqslv")) install.packages("nleqslv")
  #library(nleqslv)

  # Drift = A = -B
  # B is drift matrix that is pos def, so Phi(DeltaT) = expm(-B*DeltaT)

    if(is.null(Drift)){
    if(!is.null(Phi)){
      B <- -logm(Phi)/1
    }else{ # is.null(Phi)
      ("Either the drift matrix Drift or the autoregressive matrix Phi should be input to the function.")
      #("Note that Phi(DeltaT) = expm(-B*DeltaT).")
      stop()
    }
  }else{ # !is.null(Drift)
    B <- -Drift
    if(all(eigen(B)$val < 0)){
      #("All the eigenvalues of the drift matrix B are negative; therefore. I assume the input was A=-B instead of B. I will use -A=B in the calculation.")
      #("Note that Phi(DeltaT) = expm(-B*DeltaT).")
      ("All the eigenvalues of the drift matrix Drift are positive. Therefore. I assume the input was B=-A instead of A. I will use -B=A in the calculation.")
      B = -B
    }
  }
  # Check on B
  if(any(eigen(B)$val <= 0)){
    #("The function stopped, since some of the eigenvalues of the drift matrix B are negative or zero.")
    ("The function stopped, since some of the eigenvalues of the drift matrix Drift are positive or zero.")
    stop()
  }
  if(dim(B)[1] != dim(B)[2]){
    print(paste("The matrix (Drift or Phi) should be a square (q times q) matrix."))
    stop()
  }


  message = "There is a DeltaT such that the Phi(DeltaT) functions reach a minimum or maximum."

  if(length(B) == 1){
    q <- 1
  }else{
    q <- dim(B)[1]
  }


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
    message = "There is no DeltaT such that the Phi(DeltaT) functions reach a minimum or maximum."
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
    ("The nleqslv-function terminated.")
    ("Hence, there is no DeltaT such that the Phi(DeltaT) functions reach a minimum or maximum.")
    message = "There is no DeltaT such that the Phi(DeltaT) functions reach a minimum or maximum."
    #stop()
  }
  #
  #
  MaxDeltaT_mx[i,j] <- MaxDeltaT_ij
  }} #end loop i and j

  ## For a single i and j
  #MaxDeltaT_ij
  #print("Phi(Max DeltaT)_ij = ")
  #print((expm(-B*MaxDeltaT_mx[i,j]))[i,j])
  #print("Check on zero for 1st order deriv = ")
  #print((B %*% expm(-B*MaxDeltaT_mx[i,j]))[i,j])
  #print("")

  # If loop is in place
  MaxDeltaT_mx
  Phi_MaxDeltaT_mx <- matrix(NA, nrow = q, ncol = q)
  for(i in 1:q){
    for(j in 1:q){
      Phi_MaxDeltaT_mx[i,j] <- (expm(-B*MaxDeltaT_mx[i,j]))[i,j]
      #print("Phi(Max DeltaT)_ij = ")
      #print((expm(-B*MaxDeltaT_mx[i,j]))[i,j])
      #print("Check on zero for 1st order deriv = ")
      #print((B %*% expm(-B*MaxDeltaT_mx[i,j]))[i,j])
      #print("")
    }
  }

  ############################################################################################################

  final <- list(DeltaT_MinOrMaxPhi = MaxDeltaT_mx,
                MinOrMaxPhi = Phi_MaxDeltaT_mx)

  return(final)

}
