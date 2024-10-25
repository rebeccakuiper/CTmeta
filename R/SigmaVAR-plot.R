#' Psi-Plot: Plot of Psi / SigmaVAR
#'
#' This function makes a plot of Psi(DeltaT) / SigmaVAR(DeltaT), the residual covariance matrix of the discrete-time model, for a range of time intervals based on its underlying drift matrix. There is also an interactive web application on my website to create a Phi-plot: Phi-and-Psi-Plots and Find DeltaT (\url{https://www.uu.nl/staff/RMKuiper/Websites\%20\%2F\%20Shiny\%20apps}).
#'
#' @param DeltaT Optional. The time interval used. By default, DeltaT = 1.
#' @param Phi Matrix of size q times q of (un)standardized lagged effects of the first-order discrete-time vector autoregressive (DT-VAR(1)) model.
#' It also takes a fitted object from the classes "varest" (from the VAR() function in vars package) and "ctsemFit" (from the ctFit() function in the ctsem package); see example below. From such an object, the (standardized) Phi/Drift and SigmaVAR/Sigma matrices are calculated/extracted.
#' @param SigmaVAR Residual covariance matrix of the first-order discrete-time vector autoregressive (DT-VAR(1)) model.
#' @param Drift Optional (either Phi or Drift). Underling first-order continuous-time lagged effects matrix (i.e., Drift matrix) of the discrete-time lagged effects matrix Phi(DeltaT).
#' @param Sigma Optional (either SigmaVAR, Sigma, or Gamma). Residual covariance matrix of the first-order continuous-time (CT-VAR(1)) model, that is, the diffusion matrix.
#' @param Gamma Optional (either SigmaVAR, Sigma, or Gamma). Stationary covariance matrix, that is, the contemporaneous covariance matrix of the data.
#' Note that if Phi and SigmaVAR (or Drift and Sigma) are known, Gamma can be calculated; hence, only one out of SigmaVAR, Sigma, and Gamma is needed as input.
#' @param AddGamma Optional. Indicator (0/1) for including horizontal lines at the values for Gamma in the plot. By default, AddGamma = 1.
#' Note that SigmaVAR converges to Gamma, so the time-interval dependent curves of SigmaVAR will converge for large time-intervals to the Gamma-lines.
#' @param Stand Optional. Indicator for whether Phi (or Drift) and SigmaVAR (or Sigma) should be standardized (1) or not (0). By default, Stand = 0.
#' @param Min Optional. Minimum time interval used in the Phi-plot. By default, Min = 0.
#' @param Max Optional. Maximum time interval used in the Phi-plot. By default, Max = 10.
#' @param Step Optional. The step-size taken in the time intervals. By default, Step = 0.05. Hence, using the defaults, the Phi-plots is based on the values of Phi(DeltaT) for DeltaT = 0, 0.05, 0.10, ..., 10. Note: Especially in case of complex eigenvalues, this step size should be very small (then, the oscillating behavior can be seen best).
#' @param WhichElements Optional. Matrix of same size as Drift denoting which element/line should be plotted (1) or not (0). By default, WhichElements  selects the unique elements (that is, all variances and the unique covariances). Note that even though not all lines have to be plotted, the full Phi/Drift and Sigma(VAR)/Gamma matrices are needed to determine the selected lines.
#' @param Labels Optional. Vector with (character) labels of the lines to be plotted. The length of this vector equals the number of 1s in WhichElements (or equals q*(q+1)/2). Note, if AddGamma = 1, then twice this number is needed. By default, Labels = NULL, which renders labels with Greek letter of SigmaVAR (as a function of the time-interval); and, if AddGamma = 1, also for Gamma.
#' @param Col Optional. Vector with color values (integers) of the lines to be plotted. The length of this vector equals the number of 1s in WhichElements (or equals q*(q+1)/2, the unique elements in the symmetric matrix SigmaVAR). By default, Col = NULL, which renders the same color for effects that belong to the same outcome variable (i.e. a row in the SigmaVAR matrix). See \url{https://www.statmethods.net/advgraphs/parameters.html} for more information about the values.
#' @param Lty Optional. Vector with line type values (integers) of the lines to be plotted. The length of this vector equals the number of 1s in WhichElements (or equals q*(q+1)/2). By default, Lty = NULL, which renders solid lines for the variances and the same type of dashed line for the covariances. See \url{https://www.statmethods.net/advgraphs/parameters.html} for more information about the values.
#' @param Title Optional. A character or a list consisting of maximum 3 character-strings or 'expression' class objects that together represent the title of the Phi-plot. By default, Title = NULL, then the following code will be used for the title: as.list(expression(Sigma[VAR](Delta[t])~plot), "How do the VAR(1) (co)variance parameters vary", "as a function of the time-interval").
#' @param Diag Optional. An indicator (TRUE/FALSE) whether a vertical line for the DeltaT for which the residual covariance matrix is diagonal should be added to the plot (TRUE) or not (FALSE). This value is obtained by the function DiagDeltaT(). By default, Diag = FALSE; hence, by default, no vertical line is added.
#'
#' @return This function returns a Psi/SigmaVAR-plot for a range of time intervals.
#' @importFrom expm expm
#' @export
#' @examples
#'
#' # library(CTmeta)
#'
#' ### Make Psi-plot/SigmaVAR-plot ###
#'
#' ## Example 1 ##
#'
#' # Phi(DeltaT)
#' DeltaT <- 1
#' Phi <- myPhi[1:2,1:2]
#' q <- dim(Phi)[1]
#' SigmaVAR <- diag(q) # for ease
#' #
#' # or
#' Drift <- myDrift
#' q <- dim(Drift)[1]
#' Sigma <- diag(q) # for ease. Note that this is not the CT-equivalent of SigmaVAR.
#'
#' # Example 1.1: unstandardized Phi&SigmaVAR #
#' #
#' # Make plot of SigmaVAR (3 examples):
#' SigmaVARPlot(DeltaT, Phi, SigmaVAR)
#' SigmaVARPlot(DeltaT, Phi, SigmaVAR, Min = 0, Max = 6, Step = 0.01)                # Specifying range x-axis and precision
#' SigmaVARPlot(DeltaT, Drift = Drift, Sigma = Sigma, Min = 0, Max = 6, Step = 0.01) # Using Drift&Sigma instead of Phi&SigmaVAR
#'
#'
#' # Example 1.2: standardized Phi&SigmaVAR #
#' SigmaVARPlot(DeltaT, Phi, SigmaVAR, Stand = 1)
#'
#'
#' ## Example 2: input from fitted object of class "varest" ##
#'
#' DeltaT <- 1
#' data <- myData
#' if (!require("vars")) install.packages("vars")
#' library(vars)
#' out_VAR <- VAR(data, p = 1)
#'
#' # Example 2.1: unstandardized Phi #
#' SigmaVARPlot(DeltaT, out_VAR)
#'
#' # Example 2.2: standardized Phi #
#' SigmaVARPlot(DeltaT, out_VAR, Stand = 1)
#'
#'
#' ## Example 3: Change plot options ##
#' DeltaT <- 1
#' Phi <- myPhi[1:2,1:2]
#' q <- dim(Phi)[1]
#' SigmaVAR <- diag(q) # for ease
#' #
#' WhichElements <- matrix(1, ncol = q, nrow = q) # Now, all elements are 1
#' diag(WhichElements) <- 0 # Now, the covariances are excluded by setting the diagonals to 0.
#' Lab <- c("12", "21")
#' LabelsS <- NULL
#' LabelsG <- NULL
#' for(i in 1:length(Lab)){
#'  e <- bquote(expression(Sigma[VAR](Delta[t])[.(Lab[i])]))
#'  LabelsS <- c(LabelsS, eval(e))
#'  e <- bquote(expression(Gamma[.(Lab[i])]))
#'  LabelsG <- c(LabelsG, eval(e))
#' }
#' Labels <- c(LabelsS, LabelsG)
#' Col <- c(1,2,1,2)
#' Lty <- c(1,2,1,2)
#' # Standardized Phi and SigmaVAR
#' SigmaVARPlot(DeltaT, Phi, SigmaVAR, Stand = 1, Min = 0, Max = 10, Step = 0.05, WhichElements = WhichElements, Labels = Labels, Col = Col, Lty = Lty)
#'


SigmaVARPlot <- function(DeltaT = 1, Phi = NULL, SigmaVAR = NULL, Drift = NULL, Sigma = NULL, Gamma = NULL, AddGamma = 1, Stand = 0, Min = 0, Max = 10, Step = 0.05, WhichElements = NULL, Labels = NULL, Col = NULL, Lty = NULL, Title = NULL, Diag = FALSE) {
#DeltaT = 1; Phi = NULL; SigmaVAR = NULL;
#Drift = NULL; Sigma = NULL; Gamma = NULL; AddGamma = 1; Stand = 0; Min = 0; Max = 10; Step = 0.05; WhichElements = NULL; Labels = NULL; Col = NULL; Lty = NULL; Title = NULL; Diag = FALSE
#library(CTmeta); Phi <- myPhi[1:2,1:2]; SigmaVAR <- diag(2)

  #q <- dim(Phi)[1]
  #WhichElements <- matrix(0, ncol = q, nrow = q)
  #WhichElements[upper.tri(WhichElements)] <- 1
  # Drift = NULL; Sigma = NULL; Gamma = NULL; AddGamma = 1; Stand = 0; Min = 0; Max = 10; Step = 0.05; Labels = NULL; Col = NULL; Lty = NULL; Title = NULL; Diag = FALSE
  # Drift = NULL; Sigma = NULL; Gamma = NULL; AddGamma = 0; Stand = 0; Min = 0; Max = 10; Step = 0.05; Labels = NULL; Col = NULL; Lty = NULL; Title = NULL; Diag = FALSE
  # Diag = TRUE

  #  #######################################################################################################################
  #
  #  #if (!require("expm")) install.packages("expm")
  #  library(expm)
  #
  #  #######################################################################################################################

  # Checks:
  if(length(DeltaT) != 1){
    ErrorMessage <- (paste0("The argument DeltaT should be a scalar, that is, one number, that is, a vector with one element. Currently, DeltaT = ", DeltaT))
    return(ErrorMessage)
    stop(ErrorMessage)
  }
  if(Stand != 0 & Stand != 1){
    ErrorMessage <- (paste0("The argument Stand should be a 0 or a 1, not ", Stand))
    return(ErrorMessage)
    stop(ErrorMessage)
  }
  if(length(Min) != 1){
    ErrorMessage <- (paste0("The argument Min should be a scalar, that is, one number, that is, a vector with one element. Currently, Min = ", Min))
    return(ErrorMessage)
    stop(ErrorMessage)
  }
  if(length(Max) != 1){
    ErrorMessage <- (paste0("The argument Max should be a scalar, that is, one number, that is, a vector with one element. Currently, Max = ", Max))
    return(ErrorMessage)
    stop(ErrorMessage)
  }
  if(length(Step) != 1){
    ErrorMessage <- (paste0("The argument Step should be a scalar, that is, one number, that is, a vector with one element. Currently, Step = ", Step))
    return(ErrorMessage)
    stop(ErrorMessage)
  }
  if(!is.logical(Diag) & Diag != FALSE & Diag != TRUE){
    ErrorMessage <- (paste0("The argument 'Diag' should be T(RUE) or F(ALSE) (or 1 or 0), not ", Diag))
    stop(ErrorMessage)
  }
  #
  # Check on Phi
  if(any(class(Phi) == "varest")){
    SigmaVAR_VARest <- cov(resid(Phi))
    Phi <- Acoef(Phi)[[1]]
    CTMp <- CTMparam(DeltaT, Phi)
    if(is.null(CTMp$ErrorMessage)){
      Drift <- CTMp$Drift  # Drift <- logm(Phi)/DeltaT  # Phi <- expm(Drift * DeltaT)
    }else{
      ErrorMessage <- CTMp$ErrorMessage
      return(ErrorMessage)
      stop(ErrorMessage)
    }
    Gamma <- Gamma.fromVAR(Phi, SigmaVAR_VARest)
  } else if(any(class(Phi) == "ctsemFit")){
    Drift <- summary(Phi)$DRIFT
    Sigma_ctsem <- summary(Phi)$DIFFUSION
    Gamma <- Gamma.fromCTM(Drift, Sigma_ctsem)
  } else{

    if(is.null(Drift)){
      if(!is.null(Phi)){
        CTMp <- CTMparam(DeltaT, Phi)
        if(is.null(CTMp$ErrorMessage)){
          Drift <- CTMp$Drift  # Drift <- logm(Phi)/DeltaT  # Phi <- expm(Drift * DeltaT)
        }else{
          ErrorMessage <- CTMp$ErrorMessage
          return(ErrorMessage)
          stop(ErrorMessage)
        }
      }else{ # is.null(Phi)
        ErrorMessage <- ("Either the drift matrix Drift or the autoregressive matrix Phi should be input to the function.")
        #("Note that Phi(DeltaT) = expm(Drift*DeltaT).")
        return(ErrorMessage)
        stop(ErrorMessage)
      }
    }
    #
    # Check on B=-Drift
    if(length(Drift) > 1){
      Check_B_or_Phi(B=-Drift)
      if(all(Re(eigen(Drift)$val) > 0)){
        cat("All (the real parts of) the eigenvalues of the drift matrix Drift are positive. Therefore. I assume the input for Drift was B = -A instead of A (or -Phi instead of Phi). I will use Drift = -B = A.")
        cat("Note that Phi(DeltaT) = expm(-B*DeltaT) = expm(A*DeltaT) = expm(Drift*DeltaT).")
        Drift = -Drift
      }
      if(any(Re(eigen(Drift)$val) > 0)){
        #ErrorMessage <- ("The function stopped, since some of (the real parts of) the eigenvalues of the drift matrix Drift are positive.")
        #return(ErrorMessage)
        #stop(ErrorMessage)
        cat("If the function stopped, this is because some of (the real parts of) the eigenvalues of the drift matrix Drift are positive.")
      }
    }
  }
  #
  if(length(Drift) == 1){
    q <- 1
  }else{
    q <- dim(Drift)[1]
  }
  #
  #
  # Check on SigmaVAR, Sigma, and Gamma
  if(is.null(SigmaVAR) & is.null(Gamma) & is.null(Sigma)){ # All three unknown
    ErrorMessage <- (paste0("The arguments SigmaVAR, Sigma, or Gamma are not found: one should be part of the input. Notably, in case of the first matrix, specify 'SigmaVAR = SigmaVAR'."))
    return(ErrorMessage)
    stop(ErrorMessage)
  }else if(is.null(Gamma)){ # Gamma unknown, calculate Gamma from Phi & SigmaVAR or Drift & Sigma

    if(!is.null(SigmaVAR)){ # SigmaVAR known, calculate Gamma from Phi & SigmaVAR

      # Check on SigmaVAR
      Check_SigmaVAR(SigmaVAR, q)

      # Calculate Gamma
      if(is.null(Phi)){
        if(q == 1){
          Phi <- exp(Drift*DeltaT)
        }else{
          Phi <- expm(Drift*DeltaT)
        }
      }
      Gamma <- Gamma.fromVAR(Phi, SigmaVAR)

    }else if(!is.null(Sigma)){ # Sigma known, calculate Gamma from Drift & Sigma

      # Check on Sigma
      Check_Sigma(Sigma, q)

      # Calculate Gamma
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
      Gamma <- Gamma.fromCTM(Drift, Sigma)

    }

  }else if(!is.null(Gamma)){ # Gamma known, only check on Gamma needed

    # Checks on Gamma
    Check_Gamma(Gamma, q)

  }
  #
  #
  if(Stand == 1){
    # Standardize Drift and Gamma
    Sxy <- sqrt(diag(diag(Gamma)))
    Gamma <- solve(Sxy) %*% Gamma %*% solve(Sxy)
    Drift <- solve(Sxy) %*% Drift %*% Sxy
    #Sigma_s <- solve(Sxy) %*% Sigma %*% solve(Sxy)
  }
  #
  #
  if(!is.null(WhichElements)){
    # Check on WhichElements
    Check_WhichElts(WhichElements, q)
    nrLines <- sum(WhichElements)
  } else{
    WhichElements <- matrix(1, ncol = q, nrow = q)
    WhichElements[lower.tri(WhichElements)] <- 0
    nrLines <- sum(WhichElements)
  }
  WhichTF <- matrix(as.logical(WhichElements), q, q)
  #
  if(!is.null(Labels)){
    if(AddGamma == 1){
      if(length(Labels) != 2*nrLines){
        ErrorMessage <- (paste0("The argument Labels should contain ", 2*nrLines, " elements, that is, q*(q+1) or twice the number of 1s in WhichElements (or WhichElements is incorrectly specified); not ", length(Labels), ". Note that Labels are needed for both SigmaVAR and Gamma (since AddGamma = 1)."))
        return(ErrorMessage)
        stop(ErrorMessage)
      }
    }else{
      if(length(Labels) != nrLines){
        ErrorMessage <- (paste0("The argument Labels should contain ", nrLines, " elements, that is, q*(q+1)/2 or the number of 1s in WhichElements (or WhichElements is incorrectly specified); not ", length(Labels)))
        return(ErrorMessage)
        stop(ErrorMessage)
      }
    }
    #if(any(!is.character(Labels))){ # Note: This does not suffice, since it could also be an expression
    #  ErrorMessage <- (paste0("The argument Labels should consist of solely characters."))
    #  return(ErrorMessage)
    #  stop(ErrorMessage)
    #}
  }
  if(!is.null(Col)){
    if(AddGamma == 1){
      if(length(Col) != 2*nrLines){
        ErrorMessage <- (paste0("The argument Col should contain ", 2*nrLines, " elements, that is, q*(q+1) or twice the number of 1s in WhichElements (or WhichElements is incorrectly specified); not ", length(Col), ". Note that values (integers) are needed for both SigmaVAR and Gamma."))
        return(ErrorMessage)
        stop(ErrorMessage)
      }
    }else{
      if(length(Col) != nrLines){
        ErrorMessage <- (paste0("The argument Col should contain ", nrLines, " elements, that is, q*(q+1)/2 or the number of 1s in WhichElements (or WhichElements is incorrectly specified); not ", length(Col)))
        return(ErrorMessage)
        stop(ErrorMessage)
      }
    }
    if(any(Col %% 1 != 0)){
      ErrorMessage <- (paste0("The argument Col should consist of solely integers."))
      return(ErrorMessage)
      stop(ErrorMessage)
    }
  }
  if(!is.null(Lty)){
    if(AddGamma == 1){
      if(length(Lty) != 2*nrLines){
        ErrorMessage <- (paste0("The argument Lty should contain ", 2*nrLines, " elements, that is, q*(q+1) or twice the number of 1s in WhichElements (or WhichElements is incorrectly specified); not ", length(Lty), ". Note that values (integers) are needed for both SigmaVAR and Gamma."))
        return(ErrorMessage)
        stop(ErrorMessage)
      }
    }else{
      if(length(Lty) != nrLines){
        ErrorMessage <- (paste0("The argument Lty should contain ", nrLines, " elements, that is, q*(q+1)/2 or the number of 1s in WhichElements (or WhichElements is incorrectly specified); not ", length(Lty)))
        return(ErrorMessage)
        stop(ErrorMessage)
      }
    }
    if(any(Lty %% 1 != 0)){
      ErrorMessage <- (paste0("The argument Lty should consist of solely integers."))
      return(ErrorMessage)
      stop(ErrorMessage)
    }
  }
  if(!is.null(Title)){
    if(length(Title) != 1 & !is.list(Title)){
      ErrorMessage <- (paste0("The argument Title should be a character or a list (containing at max 3 items)."))
      return(ErrorMessage)
      stop(ErrorMessage)
    }
    if(length(Title) > 3){
      ErrorMessage <- (paste0("The list Title should at max contain 3 items. Currently, it consists of ", length(Title), " items."))
      return(ErrorMessage)
      stop(ErrorMessage)
    }
    # Option: Also check whether each element in list either a "call" or a 'character' is...
  }


  if(is.null(Labels)){
    #subscripts = NULL
    #for(i in 1:q){
    #  subscripts = c(subscripts, paste0(i, 1:q, sep=""))
    #}
    subscripts = NULL
    for(j in 1:q){
      for(i in 1:q){
        if(WhichElements[j,i] == 1){
          subscripts = c(subscripts, paste0(j, i, sep=""))
        }
      }
    }
    legendT = NULL
    for(i in 1:length(subscripts)){
      e <- bquote(expression(Sigma[VAR](Delta[t])[.(subscripts[i])]))
      legendT = c(legendT, eval(e))
    }
    if(AddGamma == 1){
      legendG = NULL
      for(i in 1:length(subscripts)){
        e <- bquote(expression(Gamma[.(subscripts[i])]))
        legendG = c(legendG, eval(e))
      }
      legendT = c(legendT, legendG)
    }
  } else{
    legendT <- as.vector(Labels)
  }

  if(is.null(Title)){
    Title_1 <- expression(Sigma[VAR](Delta[t])~plot)
    Title_2 <- "How do the VAR(1) (co)variance parameters vary"
    Title_3 <- "as a function of the time-interval"
  }else{
    Title_1 <- NULL
    Title_2 <- NULL
    if(length(Title) == 1){
      if(is.list(Title)){
        Title_3 <- Title[[1]]
      }else{
        Title_3 <- Title
      }
    }else if(length(Title) == 2){
      Title_2 <- Title[[1]]
      Title_3 <- Title[[2]]
    }else if(length(Title) == 3){
      Title_1 <- Title[[1]]
      Title_2 <- Title[[2]]
      Title_3 <- Title[[3]]
    }
  }

  if(is.null(Col)){
    Col_mx <- matrix(NA, ncol = q, nrow = q)
    for(i in 1:q){
      Col_mx[i, 1:q] <- i
    }
    Col <- t(Col_mx)[t(WhichTF)]
    #Col <- as.vector(t(Col_mx))
  }
  #
  if(is.null(Lty)){
    Lty_mx <- matrix(NA, ncol = q, nrow = q)
    diag(Lty_mx) <- 1
    Lty_mx[upper.tri(Lty_mx, diag = FALSE)] <- 2:(1+length(Lty_mx[upper.tri(Lty_mx, diag = FALSE)]))
    Lty_mx[lower.tri(Lty_mx, diag = FALSE)] <- Lty_mx[upper.tri(Lty_mx, diag = FALSE)]
    Lty <- t(Lty_mx)[t(WhichTF)]
    #Lty <- as.vector(t(Lty_mx))
  }
  #
  LWD_S <- 2.5
  LWD_G <- 1
  LWD_0 <- 1.5
  if(AddGamma == 0){
    LWD <- rep(LWD_S, length(Lty))
  }else{
    LWD <- c(rep(LWD_S, length(Lty)), rep(LWD_G, length(Lty)))
  }


  ####


  if(any(is.complex(eigen(Drift)$val))){
    while (!is.null(dev.list()))  dev.off()  # to reset the graphics pars to defaults
    # Multiple solutions, then 2x2 plots
    op <- par(mfrow=c(2,2))
    complex <- TRUE
    #nf <- layout(matrix(c(1,2,5,3,4,6),2,3,byrow = TRUE), c(3,3,1), c(2,2,1), TRUE)
    nf <- layout(matrix(c(1,2),1,2,byrow = TRUE), c(6), c(4,1), TRUE)
    #layout.show(nf)
  } else{
    op <- par(mfrow=c(1,1))
    complex <- FALSE
    #
    while (!is.null(dev.list()))  dev.off()  # to reset the graphics pars to defaults
    par(mar=c(par('mar')[1:3], 0)) # optional, removes extraneous right inner margin space
    plot.new()
    l <- legend(0, 0,
                legend = legendT, #cex=CEX,
                bty = "n",
                lty=Lty, # gives the legend appropriate symbols (lines)
                lwd=LWD,
                col=Col # gives the legend lines the correct color and width
    )
    # calculate right margin width in ndc
    w <- 1.5 *( grconvertX(l$rect$w, to='ndc') - grconvertX(0, to='ndc') )
    par(omd=c(0, 1-w, 0, 1))
    #
  }


  q <- dim(Drift)[1]

  DeltaTs<-seq(Min,Max,by=Step)

  SigmaVARDeltaTs<-array(data=NA,dim=c(q,q,length(DeltaTs)))
  if(length(Drift) == 1){
    for(i in 1:length(DeltaTs)){
      SigmaVARDeltaTs[,,i] <- Gamma - exp(Drift*DeltaTs[i]) * Gamma * t(exp(Drift*DeltaTs[i]))
    }
  }else{
    for(i in 1:length(DeltaTs)){
      SigmaVARDeltaTs[,,i] <- Gamma - expm(Drift*DeltaTs[i]) %*% Gamma %*% t(expm(Drift*DeltaTs[i]))
    }
  }

  Xlab <- expression(Time-interval~(Delta[t]))
  Ylab <- expression(Sigma[VAR](Delta[t])~values)
  #
  #wd <- getwd()
  #dev.copy(png, filename = paste0(wd, "/www/SigmaVARPlot.png"))
  teller <- 1
  # Determine YLIM based on what to be plotted (and making sure 0 is in it)
  # Elements to be plotted:
  #WhichTF <- matrix(as.logical(WhichElements), q, q)
  WhichTF_array <- array(WhichTF, dim = dim(SigmaVARDeltaTs))
  EltsInPlot <- SigmaVARDeltaTs[WhichTF_array]
  if(AddGamma == 1){
    EltsGammaInPlot <- Gamma[WhichTF]
    YLIM=c(min(EltsInPlot, EltsGammaInPlot, 0), max(EltsInPlot, 0, EltsGammaInPlot, 0))
  }else{
    YLIM=c(min(EltsInPlot, 0), max(EltsInPlot, 0))
  }
  #YLIM=c(min(SigmaVARDeltaTs, Gamma), max(SigmaVARDeltaTs, Gamma))
  plot(y=rep(0, length(DeltaTs)), x=DeltaTs, type="l", ylim=YLIM,
       ylab = Ylab,
       xlab = Xlab,
       col=1000, lwd=LWD_0, lty=1,
       main = mtext(side=3, line=2, adj=0, as.expression(Title_1)),
       sub = mtext(side=3, line=c(1,0), adj=0, c(as.expression(Title_2), as.expression(Title_3)))
  )
  #
  teller <- 0
  for(j in 1:q){
    #for(i in j:q){
    for(i in 1:q){
      #j = 1; i = 1
      if(WhichElements[j,i] == 1){
        teller <- teller + 1
        lines(y=SigmaVARDeltaTs[j,i,], x=DeltaTs, col=Col[teller], lwd=LWD_S, lty=Lty[teller])
        if(AddGamma == 1){
          lines(y=rep(Gamma[j,i], length(DeltaTs)), x = DeltaTs, col=Col[teller], lwd=LWD_G, lty=Lty[1])
        }
      }
    }
  }

  #Add lines for DiagDeltaT (if Diag == TRUE)
  if(Diag == TRUE){
    if(is.null(Phi)){
      #Phi <- expm(Drift*DeltaT)
      VarEst <- VARparam(DeltaT, Drift, Sigma)
      Phi <- VarEst$Phi
      SigmaVAR <- VarEst$SigmaVAR
    }
    DiagDt <- DiagDeltaT(Phi = Phi, SigmaVAR = SigmaVAR)
    #DiagDt <- DiagDeltaT(Drift = Drift, Sigma = Sigma)
    if(is.null(DiagDt$ErrorMessage)){
      DeltaT_DiagDt <- round(DiagDt$DeltaT_diag, 2)
      if(DeltaT_DiagDt >= Min & DeltaT_DiagDt <= Max){
        abline(v=DeltaT_DiagDt,
               col="gray", lwd=LWD_0)
        axis(side = 1, pos = (offset = YLIM[1]), DeltaT_DiagDt, cex.axis = .7, col.axis = "darkgray", col = "darkgray", lwd=LWD_0) # 3 = Add axis on top
      }
    }else{
      ErrorMessage <- DiagDt$ErrorMessage
      return(ErrorMessage)
      # TO DO is dit nu informatief?
      # TO DO zegt nu:
      # Error in DiagDt$ErrorMessage : $ operator is invalid for atomic vectors
    }
  }


  if(complex == TRUE){
    if(q<4){CEX = 1}else{CEX = (1.4-q/10)} # Check for optimal values!
    par(mar = c(0,0,0,0))
    plot.new()
    legend(par('usr')[2], par('usr')[4], xpd=NA,
           legend = legendT, cex=CEX,
           bty = "n",
           lty=Lty, # gives the legend appropriate symbols (lines)
           lwd=LWD,
           col=Col # gives the legend lines the correct color and width
    )
  }else{
    if(q<4){CEX = 1}else{CEX = (1.4-q/10)} # Check for optimal values!
    #
    #legend("topright",
    legend(par('usr')[2], par('usr')[4], xpd=NA,
           legend = legendT, cex=CEX,
           bty = "n",
           lty=Lty, # gives the legend appropriate symbols (lines)
           lwd=LWD,
           col=Col # gives the legend lines the correct color and width
    )
  }



#dev.off()

par(op)


############################################################################################################

#final <- list(.. = ...)
#return(final)

}


