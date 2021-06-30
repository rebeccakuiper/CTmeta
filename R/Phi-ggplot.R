#' Phi-plot of Phi based on its underlying drift matrix
#'
#' This function makes a Phi-plot of Phi(DeltaT) for a range of time intervals based on its underlying drift matrix. There is also an interactive web application on my website to create a Phi-plot: Phi-and-Psi-Plots and Find DeltaT (\url{https://www.uu.nl/staff/RMKuiper/Websites\%20\%2F\%20Shiny\%20apps}).
#'
#' @param DeltaT Optional. The time interval used. By default, DeltaT = 1.
#' @param Phi Matrix of size q times q of (un)standardized lagged effects. Note that the Phi (or Drift) matrix should be standardized to make a fair comparison between cross-lagged effects.
#' It also takes a fitted object from the classes "varest" (from the VAR() function in vars package) and "ctsemFit" (from the ctFit() function in the ctsem package); see example below. From such an object, the (standardized) Phi/Drift matrix is calculated/extracted.
#' @param Drift Optional (either Phi or Drift). Underling continuous-time lagged effects matrix (i.e., Drift matrix) of the discrete-time lagged effects matrix Phi(DeltaT). By default, input for Phi is used; only when Phi = NULL, Drift will be used.
#' @param Stand Optional. Indicator for whether Phi (or Drift) should be standardized (1) or not (0). In case Stand = 1, one of the following matrices should be input as well: SigmaVAR, Sigma, or Gamma (or it is subtracted from a varest or ctsemFit object). By default, Stand = 0.
#' @param SigmaVAR Optional (if Stand = 1, then either SigmaVAR, Sigma, or Gamma needed). Residual covariance matrix of the first-order discrete-time vector autoregressive (DT-VAR(1)) model.
#' @param Sigma Optional (if Stand = 1, then either SigmaVAR, Sigma, or Gamma needed). Residual covariance matrix of the first-order continuous-time (CT-VAR(1)) model, that is, the diffusion matrix.
#' @param Gamma Optional (either SigmaVAR, Sigma or Gamma). Stationary covariance matrix, that is, the contemporaneous covariance matrix of the data.
#' @param Min Optional. Minimum time interval used in the Phi-plot. By default, Min = 0.
#' @param Max Optional. Maximum time interval used in the Phi-plot. By default, Max = 10.
#' @param Step Optional. The step-size taken in the time intervals. By default, Step = 0.05. Hence, using the defaults, the Phi-plots is based on the values of Phi(DeltaT) for DeltaT = 0, 0.05, 0.10, ..., 10. Note: Especially in case of complex eigenvalues, this step size should be very small (then, the oscillating behavior can be seen best).
#' @param WhichElements Optional. Matrix of same size as Drift denoting which element/line should be plotted (1) or not (0). By default, WhichElements = NULL. Note that even though not all lines have to be plotted, the full Drift matrix is needed to determine the selected lines.
#' @param Labels Optional. Vector with (character) labels of the lines to be plotted. The length of this vector equals the number of 1s in WhichElements (or equals q*q). By default, Labels = NULL, which renders labels with Greek letter of Phi (as a function of the time-interval) together with the indices (of outcome and predictor variables).
#' @param Col Optional. Vector with color values (integers) of the lines to be plotted. The length of this vector equals the number of 1s in WhichElements (or equals q*q). By default, Col = NULL, which renders the same color for effects that belong to the same outcome variable (i.e. a row in the Drift matrix). See \url{https://www.statmethods.net/advgraphs/parameters.html} for more information about the values.
#' @param Lty Optional. Vector with line type values (integers) of the lines to be plotted. The length of this vector equals the number of 1s in WhichElements (or equals q*q). By default, Lty = NULL, which renders solid lines for the autoregressive effects and the same type of dashed line for reciprocal effects (i.e., same type for Phi_ij as for Phi_ji). See \url{https://www.statmethods.net/advgraphs/parameters.html} for more information about the values.
#' @param Title Optional. A character or a list consisting of maximum 2 character-strings or 'expression' class objects that together represent the title of the Phi-plot. By default, Title = NULL, then the following code will be used for the title: as.list(expression(Phi(Delta[t])~plot), "How do the lagged parameters vary \n as a function of the time-interval")). Note that the default 2-items list will result in 3 lines because of the use of '\n'.
#'
#' @return This function returns a Phi-plot for a range of time intervals.
#' @importFrom expm expm
#' @importFrom expm logm
#' @importFrom purrr map
#' @importFrom ggplot2 ggplot
#' @importFrom ggpubr ggarrange
#' @import dplyr
#' @export
#' @examples
#'
#' # library(CTmeta)
#'
#' ### Make Phi-plot ###
#'
#' ## Example 1 ##
#'
#' # Phi(DeltaT)
#' DeltaT <- 1
#' Phi <- myPhi[1:2,1:2]
#' # or: Drift
#' Drift <- myDrift
#'
#' # Example 1.1: unstandardized Phi #
#' #
#' # Make plot of Phi
#' ggPhiPlot(DeltaT, Phi)
#' ggPhiPlot(DeltaT, Phi, Min = 0, Max = 10, Step = 0.01)           # Specifying range x-axis and precision
#' ggPhiPlot(DeltaT, Drift = Drift, Min = 0, Max = 10, Step = 0.01) # Using Drift instead of Phi
#'
#'
#' # Example 1.2: standardized Phi #
#' q <- dim(Phi)[1]
#' SigmaVAR <- diag(q) # for ease
#' ggPhiPlot(DeltaT, Phi, Stand = 1, SigmaVAR = SigmaVAR)
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
#' ggPhiPlot(DeltaT, out_VAR)
#'
#' # Example 2.2: standardized Phi #
#' ggPhiPlot(DeltaT, out_VAR, Stand = 1)
#'
#'
#' ## Example 3: Change plot options ##
#' DeltaT <- 1
#' Phi <- myPhi[1:2,1:2]
#' q <- dim(Phi)[1]
#' SigmaVAR <- diag(q) # for ease
#' #
#' WhichElements <- matrix(1, ncol = q, nrow = q) # Now, all elements are 1
#' diag(WhichElements) <- 0 # Now, the autoregressive parameters are excluded by setting the diagonals to 0.
#' Lab <- c("12", "21")
#' Labels <- NULL
#' for(i in 1:length(Lab)){
#'  e <- bquote(expression(Phi(Delta[t])[.(Lab[i])]))
#'  Labels <- c(Labels, eval(e))
#' }
#' Col <- c(1,2)
#' Lty <- c(1,2)
#' # Standardized Phi
#' ggPhiPlot(DeltaT = 1, Phi, Stand = 1, SigmaVAR = SigmaVAR, Min = 0, Max = 10, Step = 0.05, WhichElements = WhichElements, Labels = Labels, Col = Col, Lty = Lty)
#'


ggPhiPlot <- function(DeltaT = 1, Phi = NULL, Drift = NULL, Stand = 0, SigmaVAR = NULL, Sigma = NULL, Gamma = NULL, Min = 0, Max = 10, Step = 0.05, WhichElements = NULL, Labels = NULL, Col = NULL, Lty = NULL, Title = NULL) {
# DeltaT = 1; Drift = NULL; Stand = 0; SigmaVAR = NULL; Sigma = NULL; Gamma = NULL; Min = 0; Max = 10; Step = 0.05; WhichElements = NULL; Labels = NULL; Col = NULL; Lty = NULL; Title = NULL
# library(expm); library(purrr); library(ggplot2); library(dplyr); library(ggpubr) # library(tidyverse)
# library(CTmeta); Phi <- myPhi[1:2,1:2]

  # Note needed:
  #@import tidyverse
  #@import ggpubr


  #  #######################################################################################################################
  #
  #  #if (!require("expm")) install.packages("expm")
  #  library(expm)
  #
  #  #######################################################################################################################

  # Checks:
  if(length(DeltaT) != 1){
    print(paste0("The argument DeltaT should be a scalar, that is, one number, that is, a vector with one element. Currently, DeltaT = ", DeltaT))
    stop()
  }
  if(Stand != 0 & Stand != 1){
    print(paste0("The argument Stand should be a 0 or a 1, not ", Stand))
    stop()
  }
  if(length(Min) != 1){
    print(paste0("The argument Min should be a scalar, that is, one number, that is, a vector with one element. Currently, Min = ", Min))
    stop()
  }
  if(length(Max) != 1){
    print(paste0("The argument Max should be a scalar, that is, one number, that is, a vector with one element. Currently, Max = ", Max))
    stop()
  }
  if(length(Step) != 1){
    print(paste0("The argument Step should be a scalar, that is, one number, that is, a vector with one element. Currently, Step = ", Step))
    stop()
  }
  #
  # Check on Phi
  if(any(class(Phi) == "varest")){
    Phi_VARest <- Acoef(Phi)[[1]]
    Drift <- logm(Phi_VARest)/DeltaT # Phi = expm(Drift * DeltaT)
  } else if(any(class(Phi) == "ctsemFit")){
    Drift <- summary(Phi)$DRIFT
  } else{

    if(is.null(Drift)){
      if(!is.null(Phi)){
        if(length(Phi) == 1){
          Drift <- log(Phi)/DeltaT
        }else{
          Drift <- logm(Phi)/DeltaT
        }
      }else{ # is.null(Phi)
        ("Either the drift matrix Drift or the autoregressive matrix Phi should be input to the function.")
        #("Note that Phi(DeltaT) = expm(Drift*DeltaT).")
        stop()
      }
    }
    #
    # Check on B
    if(length(Drift) > 1){
      Check_B_or_Phi(B=-Drift)
      if(all(Re(eigen(Drift)$val) > 0)){
        ("All (the real parts of) the eigenvalues of the drift matrix Drift are positive. Therefore. I assume the input for Drift was B = -A instead of A. I will use Drift = -B = A.")
        ("Note that Phi(DeltaT) = expm(-B*DeltaT) = expm(A*DeltaT) = expm(Drift*DeltaT).")
        Drift = -Drift
      }
      if(any(Re(eigen(Drift)$val) >= 0)){
        ("The function stopped, since some of (the real parts of) the eigenvalues of the drift matrix Drift are positive or zero.")
        stop()
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
  if(Stand == 1){
    # Check on SigmaVAR, Sigma, and Gamma
    if(any(class(Phi) == "varest")){
      SigmaVAR <- cov(resid(Phi))
      Phi <- Phi_VARest
      Gamma <- Gamma.fromVAR(Phi, SigmaVAR)
    }else if(any(class(Phi) == "ctsemFit")){
      Sigma <- summary(Phi)$DIFFUSION
      Gamma <- Gamma.fromCTM(Drift, Sigma)
    }else if(is.null(SigmaVAR) & is.null(Gamma) & is.null(Sigma)){ # All three unknown
      print(paste0("The arguments SigmaVAR, Sigma, or Gamma are not found: one should be part of the input (when Stand = 1). Notably, in case of the first matrix, specify 'SigmaVAR = SigmaVAR'."))
      stop()
    }else if(is.null(Gamma)){ # Gamma unknown, calculate Gamma from Phi & SigmaVAR or Drift & Sigma

      if(!is.null(SigmaVAR)){ # SigmaVAR known, calculate Gamma from Phi & SigmaVAR

        # Check on SigmaVAR
        Check_SigmaVAR(SigmaVAR, q)

        # Calculate Gamma
        if(is.null(Phi)){
          if(q == 1){
            Phi <- exp(-B*DeltaT)
          }else{
            Phi <- expm(-B*DeltaT)
          }
        }
        Gamma <- Gamma.fromVAR(Phi, SigmaVAR)


      }else if(!is.null(Sigma)){ # Sigma known, calculate Gamma from Drift & Sigma

        # Check on Sigma
        Check_Sigma(Sigma, q)

        # Calculate Gamma
        if(is.null(Drift)){
          if(q == 1){
            Drift <- log(Phi)/DeltaT
          }else{
            Drift <- logm(Phi)/DeltaT # Phi = expm(Drift * DeltaT)
          }
        }
        Gamma <- Gamma.fromCTM(Drift, Sigma)

      }

    }else if(!is.null(Gamma)){ # Gamma known, only check on Gamma needed

      # Checks on Gamma
      Check_Gamma(Gamma, q)

    }

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
    nrLines <- q*q #<- sum(WhichElements)
  }
  #
  if(!is.null(Labels)){
    if(length(Labels) != nrLines){
      print(paste0("The argument Labels should contain ", nrLines, " elements, that is, q*q or the number of 1s in WhichElements (or WhichElements is incorrectly specified); not ", length(Labels)))
      stop()
    }
    #if(any(!is.character(Labels))){ # Note: This does not suffice, since it could also be an expression
    #  print(paste0("The argument Labels should consist of solely characters."))
    #  stop()
    #}
  }
  if(!is.null(Col)){
    if(length(Col) != nrLines){
      print(paste0("The argument Col should contain ", nrLines, " elements, that is, q*q or the number of 1s in WhichElements (or WhichElements is incorrectly specified); not ", length(Col)))
      stop()
    }
    if(any(Col %% 1 != 0)){
      print(paste0("The argument Col should consist of solely integers."))
      stop()
    }
  }
  if(!is.null(Lty)){
    if(length(Lty) != nrLines){
      print(paste0("The argument Lty should contain ", nrLines, " elements, that is, q*q or the number of 1s in WhichElements (or WhichElements is incorrectly specified); not ", length(Lty)))
      stop()
    }
    if(any(Lty %% 1 != 0)){
      print(paste0("The argument Lty should consist of solely integers."))
      stop()
    }
  }
  if(!is.null(Title)){
    if(length(Title) != 1 & !is.list(Title)){
      print(paste0("The argument Title should be a character or a list (containing at max 2 items)."))
      stop()
    }
    if(length(Title) > 2){
      print(paste0("The list Title should at max contain 2 items. Currently, it consists of ", length(Title), " items."))
      stop()
    }
    # Option: Also check whether each element in list either a "call" or a 'character' is...
  }


  #def.par <- par(no.readonly = TRUE) # save default, for resetting...
  #par(def.par)  #- reset to default

  if(is.null(Labels)){
    subscripts = NULL
    for(i in 1:q){
      subscripts = c(subscripts, paste0(i, 1:q, sep=""))
    }
    legendT = NULL
    for(i in 1:(q*q)){
      e <- bquote(expression(Phi(Delta[t])[.(subscripts[i])]))
      legendT <- c(legendT, eval(e))
    }
  } else{
    legendT <- as.vector(Labels)
  }

  if(is.null(Col)){
    Col <- matrix(NA, ncol = q, nrow = q)
    for(i in 1:q){
      Col[i, 1:q] <- i
    }
    Col <- as.vector(t(Col))
  }

  if(is.null(Lty)){
    Lty <- matrix(NA, ncol = q, nrow = q)
    diag(Lty) <- 1
    Lty[upper.tri(Lty, diag = FALSE)] <- 2:(1+length(Lty[upper.tri(Lty, diag = FALSE)]))
    Lty[lower.tri(Lty, diag = FALSE)] <- Lty[upper.tri(Lty, diag = FALSE)]
    Lty <- as.vector(t(Lty))
  }

  if(is.null(Title)){
    Title_1 <- expression(Phi(Delta[t])~plot)
    Title_2 <- "How do the lagged parameters vary \n as a function of the time-interval"
  }else{
    if(length(Title) == 1){
      if(is.list(Title)){
        Title_1 <- Title[[1]]
      }else{
        Title_1 <- Title
      }
      Title_2 <- NULL
    }else if(length(Title) == 2){
      Title_1 <- Title[[1]]
      Title_2 <- Title[[2]]
    }
  }

###############################################################################################


  while (!is.null(dev.list()))  dev.off()  # to reset the graphics pars to defaults


  if(any(is.complex(eigen(Drift)$val))){
    complex <- TRUE
  } else{
    complex <- FALSE
  }


  DeltaTs<-seq(Min,Max,by=Step)

  PhiDeltaTsDF <- map(DeltaTs, function(x) {
    if(length(Drift) == 1) {exp(Drift * x)}
    else {expm(Drift * x)}
  }) %>%
    map(function(x) data.frame(Values = as.vector(t(x)))) %>%
    bind_rows %>%
    bind_cols(WhichElements = rep(as.vector(WhichElements), length(DeltaTs))) %>%
    filter(WhichElements == 1) %>%
    bind_cols(DeltaTs = rep(DeltaTs, each = sum(WhichElements)),
              Color = rep(as.character(Col), length(DeltaTs)),
              LineType = rep(as.character(Lty), length(DeltaTs)),
              Labels = rep(as.character(legendT), length(DeltaTs)))

  Xlab <- expression(Time-interval (Delta[t]))
  Ylab <- expression(Phi(Delta[t])~values)
  #
  phi_plot <- ggplot(PhiDeltaTsDF, aes(DeltaTs, Values, color = Labels, linetype = Labels)) +
    geom_line(lwd = 0.75) +
    geom_abline(intercept = 0, slope = 0, alpha = .5) +
    scale_linetype_manual(name = " ", values = Lty, labels = legendT) +
    scale_color_manual(name = " ", values = Col, labels = legendT) +
    ylab(Ylab) +
    xlab(Xlab) +
    #labs(title = as.expression(Title_1),
    #     subtitle = as.expression(Title_2)) +
    ggtitle(as.expression(Title_1), subtitle = Title_2) +
    theme_classic() +
    theme(plot.title = element_text(margin = margin(t = 20))) +
    ylim(0,1) +
    theme(
      legend.key.width = unit(2, "lines"),
      legend.spacing.x = unit(1.5, "lines"),
      legend.text = element_text(size = 12)
    ) #; phi_plot


  if(complex){
    # Multiple solutions and add 3 plots (2 for 2 different solutions and one scatter plot)
    Title_2_N <- "using an 'aliasing' matrix \n (i.e., another solution for A)"
    Title_1_N2 <- expression(Phi(Delta[t])~scatter~plot~'for'~multiples~of~Delta[t])
    Title_2_N2 <- expression(Note~that~'for'~multiples~of~Delta[t]~Phi(Delta[t])~is~unique)
    #Title_2_N2 <- expression(paste("Note that for multiples of ", Delta[t], "\n", Phi(Delta[t]), "is unique"))

    EigenDrift <- eigen(Drift)
    V <- EigenDrift$vector

    PhiDeltaTsDF_L <- list(NULL)

    for(N in 1:2){ # Note: last plot is scatter plot
      im <- complex(real = 0, imaginary = 1)
      diagN <- matrix(0, ncol = q, nrow = q)
      # Note: ordering eigenvalues is based on Mod(eigenvalues): so, if find one complex then next is its conjugate.
      W_complex <- which(Im(EigenDrift$val) != 0)
      NrComplexPairs <- length(W_complex)/2
      tellerComplex = -1
      for(i in 1:NrComplexPairs){
        tellerComplex = tellerComplex + 2
        index <- W_complex[tellerComplex]
        diagN[index,index] <- 1
        diagN[index+1,index+1] <- -diagN[index,index]
        # Note if nr of complex pairs > 1: 'diagN' should always be x and -x within a conjugate pair, but over the complex pairs x does not have to be the same...
      }
      diagN <- N * diagN
      A_N = Drift + (2 * base::pi * im / DeltaT) * V %*% diagN %*% solve(V)
      #A_N
      #print(A_N)
      Drift_N <- Re(A_N)
      #
      PhiDeltaTsDF_N <- map(DeltaTs, function(x) {
        expm(Drift_N * x)
        }) %>%
        map(function(x) data.frame(Values = as.vector(t(x)))) %>%
        bind_rows %>%
        bind_cols(WhichElements = rep(as.vector(WhichElements), length(DeltaTs))) %>%
        filter(WhichElements == 1) %>%
        bind_cols(DeltaTs = rep(DeltaTs, each = sum(WhichElements)),
                  Color = rep(as.character(Col), length(DeltaTs)),
                  LineType = rep(as.character(Lty), length(DeltaTs)),
                  Labels = rep(as.character(legendT), length(DeltaTs)))

      PhiDeltaTsDF_L[[N]] <- PhiDeltaTsDF_N
      #
    }
    # In case last plot is scatter plot
    # In last plot a scatter plot, for multiples of DeltaT, from Min to Max.
    Min_ <- Min + Min%%DeltaT # last part is remainder after integer division
    Max_ <- Max - Max%%DeltaT # last part is remainder after integer division
    DeltaTs <- seq(Min_, Max_, by=DeltaT)
    #
    PhiDeltaTsDF_4 <- map(DeltaTs, function(x) {
      expm(Drift_N * x)
    }) %>%
      map(function(x) data.frame(Values = as.vector(t(x)))) %>%
      bind_rows %>%
      bind_cols(WhichElements = rep(as.vector(WhichElements), length(DeltaTs))) %>%
      filter(WhichElements == 1) %>%
      bind_cols(DeltaTs = rep(DeltaTs, each = sum(WhichElements)),
                Color = rep(as.character(Col), length(DeltaTs)),
                LineType = rep(as.character(Lty), length(DeltaTs)),
                Labels = rep(as.character(legendT), length(DeltaTs)))

    for (i in 1:2) {
      p.plot <- ggplot(PhiDeltaTsDF_L[[i]], aes(DeltaTs, Values, color = Labels, linetype = Labels)) +
        geom_line(lwd = 0.75) +
        geom_abline(intercept = 0, slope = 0, alpha = .5) +
        scale_linetype_manual(name = " ", values = Lty, labels = legendT) +
        scale_color_manual(name = " ", values = Col, labels = legendT) +
        ylab(Ylab) +
        xlab(Xlab) +
        #labs(title = as.expression(Title_1),
        #     subtitle = Title_2_N) +
        ggtitle(as.expression(Title_1), subtitle = Title_2_N) +
        theme_classic() +
        theme(plot.title = element_text(margin = margin(t = 20))) +
        ylim(0,1) +
        theme(
          legend.key.width = unit(2, "lines"),
          legend.spacing.x = unit(1.5, "lines"),
          legend.text = element_text(size = 12)
        )

      PlotName <- paste0("Plot_", i)
      assign(PlotName, p.plot)
      #Plot_1
      #Plot_2
    }

    Plot_3 <- ggplot(PhiDeltaTsDF_4, aes(DeltaTs, Values, color = Labels, shape = Labels)) +
    geom_point(show.legend = ) +
    geom_abline(intercept = 0, slope = 0, alpha = .5) +
    scale_shape_manual(name = " ", values = Lty, labels = legendT) +
    scale_color_manual(name = " ", values = Col, labels = legendT) +
    ylab(Ylab) +
    xlab(Xlab) +
    #labs(title = as.expression(Title_1_N2),
    #     subtitle = as.expression(Title_2_N2)) +
    ggtitle(as.expression(Title_1_N2), subtitle = Title_2_N2) +
    theme_classic() +
    theme(plot.title = element_text(margin = margin(t = 20))) +
    ylim(0,1) +
    theme(
      legend.key.width = unit(2, "lines"),
      legend.spacing.x = unit(1, "lines"),
      legend.text = element_text(size = 12)
    )
    # Plot_3

    Plot <- phi_plot + theme(legend.position = "none")
    Plot_2_ <- Plot_2 + theme(legend.position = "none")
    Plot_complex <- ggarrange(plotlist = list(Plot, Plot_1, Plot_2_, Plot_3), ncol = 2, nrow = 2,
              widths = c(3,4)) %>% show

    final <- list(PhiPlot = phi_plot,
                  complex = complex,
                  PhiPlot_aliasing_1 = Plot_1,
                  PhiPlot_aliasing_2 = Plot_2,
                  PhiPlot_scatter = Plot_3,
                  PhiPlot_all = Plot_complex)

  }else{ # if not complex, then only one plot
    final <- list(PhiPlot = phi_plot,
                  complex = complex)
    print(phi_plot)
  }


  ############################################################################################################

  #final <- list(.. = ...)
  return(invisible(final))

}


