#' Time-interval (DeltaT) for which the (discrete-time) residual covariance matrix is diagonal
#'
#' Time-interval (DeltaT) for which the (discrete-time) residual covariance matrix (i.e., SigmaVAR(DeltaT)) is diagonal (together with that diagonal SigmaVAR and corresponding lagged-effects matrix Phi). The interactive web application 'Phi-and-Psi-Plots and Find DeltaT' also contains this functionality, you can find it on my website: \url{https://www.uu.nl/staff/RMKuiper/Websites\%20\%2F\%20Shiny\%20apps}.
#'
#' @param Phi Matrix of size q times q of (un)standardized lagged effects of the first-order discrete-time vector autoregressive (DT-VAR(1)) model.
#' It also takes a fitted object from the classes "varest" (from the VAR() function in vars package) and "ctsemFit" (from the ctFit() function in the ctsem package); see example below. From such an object, the (standardized) Drift matrix is calculated/extracted.
#' @param SigmaVAR Residual covariance matrix of the first-order discrete-time vector autoregressive (DT-VAR(1)) model.
#' @param Drift Optional (either Phi or Drift). Underling first-order continuous-time lagged effects matrix (i.e., Drift matrix) of the discrete-time lagged effects matrix Phi(DeltaT). By default, input for Phi is used; only when Phi = NULL, Drift will be used.
#' @param Sigma Optional (either SigmaVAR, Sigma or Gamma). Residual covariance matrix of the first-order continuous-time (CT-VAR(1)) model, that is, the diffusion matrix. By default, input for SigmaVAR is used; only when SigmaVAR = NULL, Sigma will be used.
#' @param Gamma Optional (either SigmaVAR, Sigma or Gamma). Stationary covariance matrix, that is, the contemporaneous covariance matrix of the data. By default, input for SigmaVAR is used; only when SigmaVAR = NULL, Gamma will be used.
#' Note that if Phi and SigmaVAR (or Drift and Sigma) are known, Gamma can be calculated; hence, only one out of SigmaVAR, Sigma, and Gamma is needed as input.
#' @param xstart_DeltaT Optional. Starting value for DeltaT. If you see in the SigmaVAR-plot a DeltaT for which SigmaVAR is diagonal (i.e., the covariances are zero) and the function renders DeltaT_diag = 0 as a solution, then change this start value accordingly. By default, xstart_DeltaT = 1
#'
#' @return The output renders the time-interval for which the (discrete-time) residual covariance matrix (SigmaVAR) is diagonal (DeltaT_diag), together with that diagonal SigmaVAR and corresponding lagged-effects matrix Phi (i.e., SigmaVAR(DeltaT_diag) and Phi(DeltaT_diag)).
#' @importFrom expm expm
#' @importFrom nleqslv nleqslv
#' @export
#' @examples
#'
#' # library(CTmeta)
#'
#' ## Example 1 ##
#' # Here, DeltaT_diag = 0
#'
#' ##################################################################################################
#' # Input needed in examples below with q=2 variables.
#' # I will use the example matrices stored in the package:
#' Phi <- myPhi[1:2, 1:2] # = Phi(1), that is, Phi when DeltaT = 1
#' q <- dim(Phi)[1]
#' #SigmaVAR <- diag(q) # Then, DeltaT_diag = 1 (since, implicitly, DeltaT = 1)
#' Gamma <- matrix(c(1, 0.5, 0.4, 1), byrow = T, nrow = q, ncol = q)
#' SigmaVAR <- Gamma - Phi %*% Gamma %*% t(Phi)
#' # or:
#' Drift <- myDrift
#' Sigma <- diag(2) # for ease. Note that this is not the CT-equivalent of SigmaVAR.
#' ##################################################################################################
#'
#' DiagDeltaT(Phi, SigmaVAR)
#'
#' # If you would use the drift matrix Drift as input, then use:
#' DiagDeltaT(Drift = Drift, Sigma = Sigma)
#'
#'
#' # Note that the function 'SigmaVARPlot' can help to see whether there is a DeltaT for which SigmaVAr(DeltaT) is diagonal.
#' SigmaVARPlot(DeltaT, Phi, SigmaVAR)
#'
#'
#'
#' ## Example 2: input from fitted object of class "varest" ##
#' # Here, there exists a DeltaT_diag unequal to 0
#' #
#' data <- myData
#' if (!require("vars")) install.packages("vars")
#' library(vars)
#' out_VAR <- VAR(data, p = 1)
#'
#' DiagDeltaT(out_VAR)
#'
#' # Note that the function 'SigmaVARPlot' can help to see whether there is a DeltaT for which SigmaVAr(DeltaT) is diagonal.
#' SigmaVARPlot(DeltaT, out_VAR)
#' SigmaVARPlot(DeltaT, out_VAR, Min = 0, Max = 1)
#'
#'
#'
#' ## Example 3 ##
#' Here, there exists a DeltaT_diag unequal to 0 (DeltaT_diag = 0.3565)
#'
#' q <- 2
#' Phi_2 <- matrix(c(0.5543025, 0.2547274,
#'                   0.1477329, 0.6612970),
#'                 nrow = q, ncol = q)
#' SigmaVAR_2 <- matrix(c(1.1990749, 0.1990749,
#'                        0.1990749, 1.1990749),
#'                      nrow = q, ncol = q)
#'
#' DiagDeltaT(Phi_2, SigmaVAR_2)
#'
#' # Note that the function 'SigmaVARPlot' can help to see whether there is a DeltaT for which SigmaVAr(DeltaT) is diagonal.
#' DeltaT <- 1
#' SigmaVARPlot(DeltaT, Phi_2, SigmaVAR_2)
#' SigmaVARPlot(DeltaT, Phi_2, SigmaVAR_2, Min = 0, Max = 1)
#'


DiagDeltaT <- function(Phi = NULL, SigmaVAR = NULL, Drift = NULL, Sigma = NULL, Gamma = NULL, xstart_DeltaT = 1) {
  #xstart_DeltaT <- 1

  DeltaT <- 1 # Needed for determining B/Drift.
  # Btw Cannot use it as option, yet, then I have to change code such that DeltaT_diag is adjusted accordingly (by transfarming Phi).

  # Needed: check on CTM param, Phi, and Gamma

  # Check on Phi
  if(any(class(Phi) == "varest")){
    SigmaVAR <- cov(resid(Phi))
    Phi <- Acoef(Phi)[[1]]
    #
    Gamma <- Gamma.fromVAR(Phi, SigmaVAR)
    #
    #CTparam <- CTMparam (DeltaT, Phi, SigmaVAR)
    #Drift <- CTparam$Drift
    #Sigma <- CTparam$Sigma
    CTMp <- CTMparam(DeltaT, Phi)
    if(is.null(CTMp$ErrorMessage)){
      Drift <- CTMp$Drift  # Drift <- logm(Phi)/DeltaT  # Phi <- expm(Drift * DeltaT)
    }else{
      ErrorMessage <- CTMp$ErrorMessage
      return(ErrorMessage)
      stop(ErrorMessage)
    }
  } else if(any(class(Phi) == "ctsemFit")){
    Drift <- summary(Phi)$DRIFT
    Sigma <- summary(Phi)$DIFFUSION
    #
    Gamma <- Gamma.fromCTM(Drift, Sigma)
    #
    VarEst <- VARparam(DeltaT, Drift, Sigma)
    Phi <- VarEst$Phi
    #SigmaVAR <- VarEst$SigmaVAR
  } else{

    # Drift = A = -B
    # B is drift matrix that is pos def, so Phi(DeltaT) = expm(-B*DeltaT)
    if(is.null(Phi)){
      if(!is.null(Drift)){

        if(length(Drift) == 1){
          Phi <- exp(Drift*DeltaT)
        }else{
          Phi <- expm(Drift*DeltaT)
        }

        B <- -Drift

        # Check on B
        if(length(B) > 1){
          Check_B_or_Phi(B)
          if(all(Re(eigen(B)$val) < 0)){
            cat("All (the real parts of) the eigenvalues of the drift matrix Drift are positive. Therefore. I assume the input for Drift was B = -A instead of A (or -Phi instead of Phi). I will use Drift = -B = A.")
            cat("Note that Phi(DeltaT) = expm(-B*DeltaT) = expm(A*DeltaT) = expm(Drift*DeltaT).")
            Drift <- -B
            #
            if(length(Drift) == 1){
              Phi <- exp(Drift*DeltaT)
            }else{
              Phi <- expm(Drift*DeltaT)
            }
          }
          if(any(Re(eigen(B)$val) < 0)){
            #ErrorMessage <- ("The function stopped, since some of (the real parts of) the eigenvalues of the drift matrix Drift are positive.")
            #stop(ErrorMessage)
            cat("If the function stopped, this is because some of (the real parts of) the eigenvalues of the drift matrix Drift are positive.")
          }
        }
      }else{ # is.null(Drift)
        ErrorMessage <- ("Either the drift matrix Drift or the autoregressive matrix Phi should be input to the function.")
        #("Note that Phi(DeltaT) = expm(-B*DeltaT).")
        return(ErrorMessage)
        stop(ErrorMessage)
      }
    }else{ # Phi not NULL
      if(length(Phi) != 1){
        Check_Phi_or_B(Phi)
      }
    }
    #
    if(length(Phi) == 1){
      q <- 1
    }else{
      q <- dim(Phi)[1]
    }
    #
    if(is.null(Drift)){
      CTMp <- CTMparam(DeltaT, Phi)
      if(is.null(CTMp$ErrorMessage)){
        Drift <- CTMp$Drift  # Drift <- logm(Phi)/DeltaT  # Phi <- expm(Drift * DeltaT)
      }else{
        ErrorMessage <- CTMp$ErrorMessage
        return(ErrorMessage)
        stop(ErrorMessage)
      }
    }

    # Check on SigmaVAR, Sigma, and Gamma
    if(is.null(SigmaVAR) & is.null(Gamma) & is.null(Sigma)){ # All three unknown
      ErrorMessage <- (paste0("The arguments SigmaVAR, Sigma, or Gamma are not found: one should be part of the input. Notably, in case of the first matrix, specify 'SigmaVAR = SigmaVAR'."))
      return(ErrorMessage)
      stop(ErrorMessage)
    }else if(is.null(Gamma)){ # Gamma unknown

      if(!is.null(SigmaVAR)){ # SigmaVAR known, use SigmaVAR and Phi or Drift

        # Check on SigmaVAR
        Check_SigmaVAR(SigmaVAR, q)

        # Calculate Gamma
        Gamma <- Gamma.fromVAR(Phi, SigmaVAR)

      }else if(!is.null(Sigma)){ # Sigma known

        # Check on Sigma
        Check_Sigma(Sigma, q)

        # Calculate Gamma
        Gamma <- Gamma.fromCTM(Drift, Sigma)

      }

    }else if(!is.null(Gamma)){ # Gamma known

      # Checks on Gamma
      Check_Gamma(Gamma, q)

    }

  }
  #
  if(length(Phi) == 1){
    q <- 1
  }else{
    q <- dim(Phi)[1]
  }

  # CHECKs and calculate corresponding Gamma (Check all matrices pos def and B=-Drift also not complex)
  outcome.GammaAndChecks <- ChecksCTM(Drift, Gamma = Gamma)
  #Gamma <- outcome.GammaAndChecks$Gamma
  ChecksAreFine <- outcome.GammaAndChecks$ChecksAreFine
  errorMatrices <- outcome.GammaAndChecks$error


  message_start <- "No message / warning / error. Hence, there is positive DeltaT for which the diagonals/variances in SigmaVAR are positive (and off-diagonals 0)."
  message <- message_start

  message_startvalues <- "In case the Psi-plot/SigmaVAR-plot does show a solution (or another solution) for DeltaT such that Psi is diagonal (i.e., the covariances are 0), \n
  alter the starting value for 'DeltaT_diag'. Notably, by default, the value 1 is used."
  # Note that in theory the starting values for the q variances can be a problem as well:
  # I may want to adjust that first 9depending on the solution and message obtained).
  if(is.null(xstart_DeltaT)){
    xstart_DeltaT <- 1
  }else{
    if(length(xstart_DeltaT) != 1){
      message_startvalues <- "The starting value for 'DeltaT_diag' should be 1 number. \n Since it is not the case, the value 1 is used."
      #cat(message_startvalues)
      xstart_DeltaT <- 1
    }
  }
  #
  xstart <- c(rep(1,q), xstart_DeltaT)


  EigenPhi <- eigen(Phi)
  eigenvaluesPhi <- EigenPhi$val
  D_Phi <- diag( eigenvaluesPhi )
  V <- EigenPhi$vectors
  invV <- solve(V)

  Q = invV %*% Gamma %*% t(invV)
  DD = eigenvaluesPhi %*% t(eigenvaluesPhi)


  # General set of equations/restrictions
  SolveForSigmaAndDelta_fie <- function(q, Q, invV, DD) {
    teller = 0
    function(x) {
      y <- numeric((q+1)*q/2)
      for(i in 1:q){
        for(j in i:q){
          teller = teller + 1
          sum = 0
          for(h in 1:q){
            sum = sum + x[h] * (invV[j,h] * invV[i,h])
          }
          y[teller] <- Q[i,j] - sum - Q[i,j] * (DD[i,j]^x[q+1])
        }
      }
      y
    }
  }
  SolveForSigmaAndDelta <- SolveForSigmaAndDelta_fie(q, Q, invV, DD)
  #
  # Note that there are (q+1)*q/2 equations/restrictions and q+1 unknowns, so if q > 2 delete (the last) (q+1)*(q-2)/2 lines!
  # If q >= 5, then, some equations are not inspected...
  #
  # Next, the first q+1 equations
  SolveForSigmaAndDelta_first_fie <- function(q, Q, invV, DD) {
    teller = 0
    function(x) {
      y <- numeric(q+1)
      for(i in 1:q){
        if (teller == q+1+1){ break }
        for(j in i:q){
          teller = teller + 1
          if (teller == q+1+1){ break }
          sum = 0
          for(h in 1:q){
            sum = sum + x[h] * (invV[j,h] * invV[i,h])
          }
          y[teller] <- Q[i,j] - sum - Q[i,j] * (DD[i,j]^x[q+1])
        }
      }
      y
    }
  }
  SolveForSigmaAndDelta_first <- SolveForSigmaAndDelta_first_fie(q, Q, invV, DD)
  #
  #xstart <- rep(1,q+1)
  #xstart <- c(diag(SigmaVAR), 1)
  #SolveForSigmaAndDelta_first(xstart)
  fstart_first <- SolveForSigmaAndDelta_first(xstart)
  if(all(is.nan(fstart_first) == FALSE) == FALSE){
    message_nan <- "For the first (q+1) equations, the starting values for the q variances and for DeltaT_diag do not render values from the function."
    #cat(message_nan)
  }


  # Note that there are (q+1)*q/2 equations and q+1 unknowns, so if q > 2 delete (the last) (q+1)*(q-2)/2 lines!
  # Next, the last q+1 equations
  SolveForSigmaAndDelta_last_fie <- function(q, Q, invV, DD) {
    teller = 0
    function(x) {
      y <- numeric(q+1)
      for(i in q:1){
        if (teller == q+1+1){ break }
        for(j in q:i){
          teller = teller + 1
          if (teller == q+1+1){ break }
          sum = 0
          for(h in 1:q){
            sum = sum + x[h] * (invV[j,h] * invV[i,h])
          }
          y[teller] <- Q[i,j] - sum - Q[i,j] * (DD[i,j]^x[q+1])
        }
      }
      y
    }
  }
  SolveForSigmaAndDelta_last <- SolveForSigmaAndDelta_last_fie(q, Q, invV, DD)
  #
  #xstart <- rep(1,q+1)
  #xstart <- rep(0.2,q+1)
  #xstart <- rep(0.4,q+1)
  #xstart <- c(diag(SigmaVAR), 1)
  #SolveForSigmaAndDelta_last(xstart)
  fstart_last <- SolveForSigmaAndDelta_last(xstart)
  if(all(is.nan(fstart_last) == FALSE) == FALSE){
    if(all(is.nan(fstart_first) == FALSE) == FALSE){
      message_nan <- "For the first and last (q+1) equations, the starting values for the q variances and for DeltaT_diag do not render values from the function."
      #cat(message_nan)
    }else{
      message_nan <- "For the last (q+1) equations, the starting values for the q variances and for DeltaT_diag do not render values from the function."
      #cat(message_nan)
    }
  }

  # Note: If there is a solution to all equations, then this is also a solution to the subset.
  #However, the subset can have more solutions, not only the solution for all equations.
  # In that sense, multiple starting values should be used. But, I did not do that here.
  # I do check whether the found solution is in agreement with all equations.


  # Find solutions based on subsets of equations
  #
  #if (!require("nleqslv")) install.packages("nleqslv")
  #library(nleqslv)
  #
  #xstart <- rep(0.2,q+1)
  #
  sol_first <- nleqslv(xstart, SolveForSigmaAndDelta_first, control=list(btol=.001, allowSingular = TRUE)) #, method="Newton")
  DiagAndDelta_first <- sol_first$x
  #
  sol_last <- nleqslv(xstart, SolveForSigmaAndDelta_last, control=list(btol=.01))
  DiagAndDelta_last <- sol_last$x
  #
  #sol_first$x
  #sol_last$x
  #
  # Function terminated if sol_first$termcd does not equal 1 (and 1 is "Function criterion is near zero. Convergence of function values has been achieved.")
  if(sol_first$termcd != 1 & sol_last$termcd != 1){ # Both terminated
    message_termcd <- "The nleqslv-function terminated (for the first and last q+1 equations)."
    #cat(message_termcd)
    Terminated = 3 # = both terminated
  }else{ # Not both terminated, but one or none
    Terminated = 0 # = none terminated
    if(sol_first$termcd != 1){
      message_termcd <- "The nleqslv-function terminated (for the first q+1 equations)."
      #cat(message_termcd)
      Terminated = 1 # = first terminated
    }
    if(sol_last$termcd != 1){
      message_termcd <- "The nleqslv-function terminated (for the last q+1 equations)."
      #cat(message_termcd)
      Terminated = 2 # = last terminated
    }
  }

  # Check whether the solution from part of the equations are in agreement with all equations.
  sol_All_first = 0
  sol_All_last = 0
  if(Terminated == 0 | Terminated == 1){
    if(all(abs(SolveForSigmaAndDelta(sol_first$x)) < 0.000001)){
      sol_All_first = 1
    }
  }
  if(Terminated == 0 | Terminated == 2){
    if(all(abs(SolveForSigmaAndDelta(sol_last$x)) < 0.000001)){
      sol_All_last = 1
    }
  }
  #
  # If q > 4, then some equations not inspected.
  # Hence, if q > 4 and both do not render a solution to all equations, then look at other equations as wel!
  #if(q > 4 & sol_All_first == 0 & sol_All_last == 0){
  #  cat("'q > 4 & sol_All_first == 0 & sol_All_last == 0' happened. \n
  #        Note that there is now also not a zero solution - which should always exist.")
  #  # Then, inspect the (q+1)*(q-4)/2 not seen equations. See how many sets of q+1 needed.
  #  # Notably, using other subsets or other starting values can help as well.
  #}



  if(sol_All_first == 1 & sol_All_last == 1){ # Then, both have a solution

    if(DiagAndDelta_first[q+1] > 0.0001 | DiagAndDelta_last[q+1] > 0.0001){ # Then, at least one positive DeltaT
      if(abs(DiagAndDelta_first[q+1] - DiagAndDelta_last[q+1]) < 0.0001){ # Then, same solution for DeltaT
        DiagAndDelta <- DiagAndDelta_first
        #
        message_DiagAndDelta <- "There are two solutions and they are the same."
        #cat(message_DiagAndDelta)
      }else{ # Note same solution
        if(DiagAndDelta_first[q+1] > 0.0001 & DiagAndDelta_last[q+1] > 0.0001){ # Then, both positive solution and highest may be infinity, so take lowest one.
          DiagAndDelta <- DiagAndDelta_last
          if(DiagAndDelta_first[q+1] < DiagAndDelta_last[q+1]){
            DiagAndDelta <- DiagAndDelta_first
          }
          #
          message_DiagAndDelta <- "There are two positive solutions, the lowest DeltaT is used."
          #cat(message_DiagAndDelta)
        }else{ # Then, only one positive and thus take highest one.
          DiagAndDelta <- DiagAndDelta_last
          if(DiagAndDelta_first[q+1] > DiagAndDelta_last[q+1]){
            DiagAndDelta <- DiagAndDelta_first
          }
          #
          message_DiagAndDelta <- "There are two solutions, but only one positive DeltaT."
          #cat(message_DiagAndDelta)
        }
      }
      if(any(DiagAndDelta[1:q] < 0)){ # All variances should be positive
        message <- "There is no positive DeltaT such that SigmaVAR is a 'positive diagonal' matrix. \n That is, the solution contains one or more negative variances."
        #cat(message)
      }
    }else{ # Then, both < 0.0001 and probably (at least one) near 0.
      message <- "There is no non-negative DeltaT such that SigmaVAR is a diagonal matrix. \n Only for DeltaT approximately 0, it is (approximately) a 0-matrix."
      #cat(message)
      #
      message_DiagAndDelta <- "There are two solutions and both are negative. At least one of them is probaly near zero."
      #cat(message_DiagAndDelta)
      #
      # If q > 4, inspect the (q+1)*(q-4)/2 not seen equations. See how many sets of q+1 needed.
      # Notably, using other subsets or other starting values can help as well.
    }

  }else{

    #Only one of them has a solultion, obtain that one:
    if(sol_All_first == 1){DiagAndDelta <- DiagAndDelta_first}
    if(sol_All_last == 1){DiagAndDelta <- DiagAndDelta_last}

    # Check positive DeltaT and positive variances/diagonals:
    if(DiagAndDelta[q+1] < 0.0001 & any(DiagAndDelta[1:q] < 0)){
      message <- "There is no non-negative DeltaT such that SigmaVAR is a 'positive diagonal' matrix. \n Only for DeltaT approximately 0, it is (approximately) a 0-matrix."
      #cat(message)
      #
      # If q > 4, inspect the (q+1)*(q-4)/2 not seen equations. See how many sets of q+1 needed.
      # Notably, using other subsets or other starting values can help as well.
    }else if(DiagAndDelta[q+1] < 0.0001){
      message <- "There is no non-negative DeltaT such that SigmaVAR is a diagonal matrix. \n Only for DeltaT approximately 0, it is (approximately) a 0-matrix."
      #cat(message)
      #
      # If q > 4, inspect the (q+1)*(q-4)/2 not seen equations. See how many sets of q+1 needed.
      # Notably, using other subsets or other starting values can help as well.
    }else if(any(DiagAndDelta[1:q] < 0)){ # All variants should be positive
      message <- "There is no positive DeltaT such that SigmaVAR is a 'positive diagonal' matrix. \n That is, the solution contains one or more negative variances."
      #cat(message)
    }

    }


  if(message == message_start){

    S <- matrix(0,q,q)
    for(i in 1:q){
    S[i,i] <- DiagAndDelta_first[i]
    }
    DeltaT <- DiagAndDelta_first[q+1]

    Phi_DeltaT <- V %*% (D_Phi^DeltaT) %*% invV

    #Gamma <- Gamma.fromVAR <- function(Phi_DeltaT, S)
    Sxy <- sqrt(diag(diag(Gamma)))
    Gamma_s <- solve(Sxy) %*% Gamma %*% solve(Sxy)
    Phi_DeltaT_s <- solve(Sxy) %*% Phi_DeltaT %*% Sxy
    S_s <- solve(Sxy) %*% S %*% solve(Sxy)


    final <- list(DeltaT_diag = DeltaT,
                  message = message,
                  Phi_DeltaT_diag = Phi_DeltaT,
                  SigmaVAR_DeltaT_diag = S,
                  Gamma = Gamma,
                  StandPhi_DeltaT_diag = Phi_DeltaT_s,
                  StandSigmaVAR_DeltaT_diag = S_s,
                  Gamma_s = Gamma_s,
                  #MaxDiff_IfNonzeroIndicationMultipleSolutions = max(abs(DiagAndDelta_first - DiagAndDelta_last)),
                  errorMatrices=errorMatrices,
                  message_startvalues=message_startvalues)

    if(q > 4){
      final <- list(final,
                    Warning <- "Since q > 4, some equations in calculating DeltaT_diag were not used. \n
                    When from the Psi-plot/SigmaVAR-plot it is clear that there exist a positive (non-zero) 'DeltaT_diag' solution and \n
                    adjusting the starting value accordingly does not help, please contact me (r.m.kuiper@uu.nl). \n
                    Then, I will add a part to the code where the currently un-used equations are inspected."
                    # Note that in theory the starting values for the q variances can be a problem as well:
                    # I may want to adjust that first 9depending on the solution and message obtained).
      )
    }

  }else{
    final <- list(DeltaT_diag = 0,
                  message = message,
                  Phi_DeltaT_diag = diag(q),
                  SigmaVAR_DeltaT_diag = matrix(0,q,q),
                  errorMatrices = errorMatrices,
                  message_startvalues = message_startvalues)
  }

  return(final)
}
