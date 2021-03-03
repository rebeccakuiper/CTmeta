#' Phi-plot of Phi based on its underling drift matrix
#'
#' This function makes a Phi-plot of Phi(DeltaT) for a range of time intervals based on its underling drift matrix. There is also an interactive web application on my website to create a Phi-plot: Phi-and-Psi-Plots and Find DeltaT (\url{https://www.uu.nl/staff/RMKuiper/Websites\%20\%2F\%20Shiny\%20apps}).
#'
#' @param DeltaT Optional. The time interval used. By default, DeltaT = 1.
#' @param Phi Matrix of size q times q of (un)standardized lagged effects. Note that the Phi (or Drift) matrix should be standardized to make a fair comparison between cross-lagged effects.
#' It also takes a fitted object from the classes "varest" (from the VAR() function in vars package) and "ctsemFit" (from the ctFit() function in the ctsem package); see example below. From such an object, the (standardized) Drift matrix is calculated/extracted.
#' @param Drift Optional (either Phi or Drift). Underling continuous-time lagged effects matrix (i.e., Drift matrix) of the discrete-time lagged effects matrix Phi(DeltaT). By default, input for Phi is used; only when Phi = NULL, Drift will be used.
#' @param Min Optional. Minimum time interval used in the Phi-plot. By default, Min = 0.
#' @param Max Optional. Maximum time interval used in the Phi-plot. By default, Max = 10.
#' @param Step Optional. The step-size taken in the time intervals. By default, Step = 0.05. Hence, using the defaults, the Phi-plots is based on the values of Phi(DeltaT) for DeltaT = 0, 0.05, 0.10, ..., 10. Note: Especially in case of complex eigenvalues, this step size should be very small (then, the oscillating behaviour can be seen best).
#' @param WhichElements Optional. Matrix of same size as Drift denoting which element/line should be plotted (1) or not (0). By default, WhichElements = NULL. Note that even though not all lines have to be plotted, the full Drift matrix is needed to determine the selected lines.
#' @param Labels Optional. Vector with (character) labels of the lines to be plotted. The length of this vector equals the number of 1s in WhichElements (or equals q*q). By default, Labels = NULL, which renders labels with Greek letter of Phi (as a function of the time-interval) together with the indices (of outcome and predictor variables).
#' @param Col Optional. Vector with color values (integers) of the lines to be plotted. The length of this vector equals the number of 1s in WhichElements (or equals q*q). By default, Col = NULL, which renders the same color for effects that belong to the same outcome variable (i.e. a row in the Drift matrix). See \url{https://www.statmethods.net/advgraphs/parameters.html} for more information about the values.
#' @param Lty Optional. Vector with line type values (integers) of the lines to be plotted. The length of this vector equals the number of 1s in WhichElements (or equals q*q). By default, Lty = NULL, which renders solid lines for the autoregressive effects and the same type of dashed line for reciprocal effects (i.e., same type for Phi_ij as for Phi_ji). See \url{https://www.statmethods.net/advgraphs/parameters.html} for more information about the values.
#' @param Title Optional. A character or a list consisting of maximum 3 characters or 'call' class objects, like from the paste() function, that together represent the title of the Phi-plot. By default, Title = NULL, then the following code will be used for the title: as.list(expression(paste(Phi(Delta[t]), " plot:"), "How do the lagged parameters vary", "as a function of the time-interval"))).
#'
#' @return This function returns a Phi-plot for a range of time intervals.
#' @importFrom expm expm
#' @import plyr
#' @import tidyverse
#' @import ggpubr
#' @export
#' @examples
#' ### Make Phi-plot ###
#'
#' ## Example 1 ##
#' #
#' # Phi(DeltaT)
#' DeltaT <- 1
#' Phi <- myPhi[1:2,1:2] # For simplicity, it is assumed that this is a standardized Phi matrix.
#' #
#' # Determine the continuous-time equivalent, that is, the drift matrix
#' if (!require("expm")) install.packages("expm") # Use expm package for function logm()
#' library(expm)
#' Drift <- logm(Phi)/DeltaT
#' #
#' # Make plot of Phi
#' ggPhiPlot(DeltaT, Phi)
#' ggPhiPlot(DeltaT, Phi, Min = 0, Max = 10, Step = 0.01)
#' ggPhiPlot(DeltaT, Drift = Drift, Min = 0, Max = 10, Step = 0.01)
#'
#'
#' ## Example 2: input from fitted object of class "varest" ##
#' #
#' DeltaT <- 1
#' data <- myData
#' if (!require("vars")) install.packages("vars")
#' library(vars)
#' out_VAR <- VAR(data, p = 1)
#' ggPhiPlot(DeltaT, out_VAR)
#'
#'
#' ## Example 3: Change plot options ##
#' # Note: use Phi or Drift from Example 1
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
#' ggPhiPlot(DeltaT = 1, Phi, Min = 0, Max = 10, Step = 0.5, WhichElements, Labels, Col, Lty)
#'


ggPhiPlot <- function(DeltaT = 1, Phi = NULL, Drift = NULL, Min = 0, Max = 10, Step = 0.05, WhichElements = NULL, Labels = NULL, Col = NULL, Lty = NULL, Title = NULL) {

  #  #######################################################################################################################
  #
  #  #if (!require("expm")) install.packages("expm")
  #  library(expm)
  #
  #  #######################################################################################################################

  # Checks:
  if(length(DeltaT) != 1){
    print(paste("The argument DeltaT should be a scalar, that is, one number, that is, a vector with one element."))
    stop()
  }
  if(length(Min) != 1){
    print(paste("The argument Min should be a scalar, that is, one number, that is, a vector with one element."))
    stop()
  }
  if(length(Max) != 1){
    print(paste("The argument Max should be a scalar, that is, one number, that is, a vector with one element."))
    stop()
  }
  if(length(Step) != 1){
    print(paste("The argument Step should be a scalar, that is, one number, that is, a vector with one element."))
    stop()
  }
  #
  # Check on Phi
  if(any(class(Phi) == "varest")){
    Phi_VARest <- Acoef(Phi)[[1]]
    Drift <- logm(Phi_VARest)/DeltaT # Phi = expm(Drift * deltaT)
    if(length(Drift) == 1){
      q <- 1
    }else{
      q <- dim(Drift)[1]
    }
    # TO DO bepaal standardized Phi en dus Drift!
  } else if(any(class(Phi) == "ctsemFit")){
    Drift <- summary(Phi)$DRIFT
    if(length(Drift) == 1){
      q <- 1
    }else{
      q <- dim(Drift)[1]
    }
    # TO DO bepaal standardized Drift!
  } else{

    if(is.null(Drift)){
      if(!is.null(Phi)){
        Drift <- logm(Phi)/1
      }else{ # is.null(Phi)
        ("Either the drift matrix Drift or the autoregressive matrix Phi should be input to the function.")
        ("Note that Phi(DeltaT) = expm(-B*DeltaT).")
        stop()
      }
    }else{ # !is.null(Drift)
      if(all(eigen(Drift)$val > 0)){
        ("All the eigenvalues of the drift matrix Drift are positive; therefore. I assume the input was B=-A instead of A. I will use -Drift in the calculation.")
        ("Note that Phi(DeltaT) = expm(-B*DeltaT) = expm(A*DeltaT).")
        Drift = -Drift
      }
    }

    if(length(Drift) == 1){
      q <- 1
    }else{
      #
      if(is.null(dim(Drift))){
        if(!is.null(length(Drift))){
          print(paste("The argument Drift (or Phi) is not a matrix of size q times q."))
          stop()
        }else{
          print(paste("The argument Drift (or Phi) is not found: The lagged effects matrix Drift is unknown, but should be part of the input."))
          stop()
        }
      }else{
        if(dim(Drift)[1] != dim(Drift)[2] | length(dim(Drift)) != 2){
          print(paste("The argument Drift (or Phi) is not a matrix of size q times q."))
          stop()
        }
        q <- dim(Drift)[1]
      }
    }
  }
  #
  #
  if(!is.null(WhichElements)){
    if(length(WhichElements) == 1){
      if(q != 1){
        print(paste("The argument WhichElements is one element and not a matrix of size q times q, with q = ", q, "."))
        stop()
      }
    } else if(dim(WhichElements)[1] != dim(WhichElements)[2] | length(dim(WhichElements)) != 2){
      print(paste("The argument WhichElements is not a (square) matrix. It should be a matrix of size q times q, with q = ", q, "."))
      stop()
    } else if(dim(WhichElements)[1] != q){
      print(paste("The argument WhichElements is not a matrix of size q times q, with q = ", q, "."))
      stop()
    }
    if(any(WhichElements != 0 & WhichElements != 1)){
      print(paste("The argument WhichElements should consist of solely 1s and 0s."))
      stop()
    }
    nrLines <- sum(WhichElements)
  } else{
    WhichElements <- matrix(1, ncol = q, nrow = q)
    nrLines <- q*q #<- sum(WhichElements)
  }
  if(!is.null(Labels)){
    if(length(Labels) != nrLines){
      print(paste("The argument Labels should contain ", nrLines, " elements, that is, q*q or the number of 1s in WhichElements (or WhichElements is incorrectly specified)."))
      stop()
    }
    #if(any(!is.character(Labels))){ # TO DO could also be an expression
    #  print(paste("The argument Labels should consist of solely characters."))
    #  stop()
    #}
  }
  if(!is.null(Col)){
    if(length(Col) != nrLines){
      print(paste("The argument Col should contain ", nrLines, " elements, that is, q*q or the number of 1s in WhichElements (or WhichElements is incorrectly specified)."))
      stop()
    }
    if(any(Col %% 1 != 0)){
      print(paste("The argument Col should consist of solely integers."))
      stop()
    }
  }
  if(!is.null(Lty)){
    if(length(Lty) != nrLines){
      print(paste("The argument Lty should contain ", nrLines, " elements, that is, q*q or the number of 1s in WhichElements (or WhichElements is incorrectly specified)."))
      stop()
    }
    if(any(Lty %% 1 != 0)){
      print(paste("The argument Lty should consist of solely integers."))
      stop()
    }
  }
  if(!is.null(Title)){
    if(length(Title) != 1 & !is.list(Title)){
      print(paste("The argument Title should be a character or a list (containing at max 3 items)."))
      stop()
    }
    if(length(Title) > 3){
      print(paste("The argument (list) Title should at max contain 3 items. Currently, it consists of ", length(Title), " items."))
      stop()
    }
    # TO DO check of elk element in list een "call" of een 'character' is...
  }

  # TO DO bepaal standardized Drift! Dus dan voor een VAR(1) de Sigma gebruiken of Gamma!
  # TO DO Geldt voor andere modellen ook dat ik Gamma kan gebruiken??


  #def.par <- par(no.readonly = TRUE) # save default, for resetting...
  #par(def.par)  #- reset to default

  if(is.null(Labels)){
    subscripts = NULL
    for(i in 1:q){
      subscripts = c(subscripts, paste(i, 1:q, sep=""))
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
    gg_title <- as.list(expression(paste(Phi(Delta[t]), " plot:"), "How do the lagged parameters vary as a function of the time-interval"))
  }else{
    if(length(Title) == 1){
      if(is.list(Title)){gg_title <- Title[[1]]}
      else {gg_title <- list(Title, " ")}
    }
    if(length(Title) == 2){
      title1 <- Title[[1]]
      title2 <- Title[[2]]
      gg_title <- list(title1, title2)
    }
  }
  if(is.null(Title)){
    #Title <- as.list(expression(paste(Phi(Delta[t]), " plot:"), "How do the overall lagged parameters vary", "as a function of the time-interval"))
    Title <- as.list(expression(paste(Phi(Delta[t]), " plot:"), "How do the lagged parameters vary", "as a function of the time-interval"))
  }else{
    if(length(Title) == 1){
      if(is.list(Title)){Title <- Title[[1]]}
      Title <- list(" ", " ", Title)
    }
    if(length(Title) == 2){
      title1 <- Title[[1]]
      title2 <- Title[[2]]
      Title <- list("", title1, title2)
    }
  }



  if(any(is.complex(eigen(Drift)$val))){
    while (!is.null(dev.list()))  dev.off()  # to reset the graphics pars to defaults
    # Multiple solutions, then 2x2 plots
    op <- par(mfrow=c(2,2))
    complex <- TRUE
    nf <- layout(matrix(c(1,2,5,3,4,6),2,3,byrow = TRUE), c(3,3,1), c(2,2,1), TRUE)
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
                lwd=rep(2, q*q),
                col=Col # gives the legend lines the correct color and width
    )
    # calculate right margin width in ndc
    w <- 1.5 *( grconvertX(l$rect$w, to='ndc') - grconvertX(0, to='ndc') )
    par(omd=c(0, 1-w, 0, 1))
    #
  }


  DeltaTs<-seq(Min,Max,by=Step)

  PhiDeltaTs<-array(data=NA,dim=c(q,q,length(DeltaTs)))
  if(length(Drift) == 1){
    for(i in 1:length(DeltaTs)){
      PhiDeltaTs[,,i]<-exp(Drift*DeltaTs[i])
    }
  }else{
    for(i in 1:length(DeltaTs)){
      PhiDeltaTs[,,i]<-expm(Drift*DeltaTs[i])
    }
  }

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

  phi_plot <- ggplot(PhiDeltaTsDF, aes(DeltaTs, Values, color = Labels, linetype = Labels)) +
    geom_line(lwd = 0.75) +
    geom_abline(intercept = 0, slope = 0, alpha = .5) +
    scale_linetype_manual(name = " ", values = Lty, labels = legendT) +
    scale_color_manual(name = " ", values = Col, labels = legendT) +
    ylab(expression(paste(Phi(Delta[t]), " values"))) +
    xlab(expression(paste("Time-interval (", Delta[t], ")", sep = ""))) +
    labs(title = as.expression(gg_title[1]),
         subtitle = as.expression(gg_title[-1])) +
    theme_classic() +
    theme(plot.title = element_text(margin = margin(t = 20))) +
    ylim(0,1)

  plot_list <- list(phi_plot + theme(legend.position = "none"))

  show(phi_plot)
  #wd <- getwd()
  #dev.copy(png, filename = paste0(wd, "/www/PhiPlot.png"))
  teller <- 1
  plot(y=rep(0, length(DeltaTs)), x=DeltaTs, type="l", ylim=c(min(PhiDeltaTs), max(PhiDeltaTs)),
       #ylab = expression(paste("Overall ", Phi(Delta[t]), " values")),
       ylab = expression(paste(Phi(Delta[t]), " values")),
       xlab = expression(paste("Time-interval (", Delta[t], ")", sep="")),
       col=1000, lwd=2, lty=1,
       main=mtext(do.call(expression, Title), side=3, line = c(2,1,0), cex = 1 )
       #"Effect lag curve: \n How do the VAR(1) parameters Phi vary \n as a function of the time-interval"
  )


  teller <- 0
  for(j in 1:q){
    for(i in 1:q){
      if(WhichElements[j,i] == 1){
        teller <- teller + 1
        lines(y=PhiDeltaTs[j,i,], x=DeltaTs, col=Col[teller], lwd=2, lty=Lty[teller])
      }
    }
  }

  if(complex == FALSE){
    if(q<4){CEX = 1}else{CEX = (1.4-q/10)} # Check for optimal values!
    #
    #legend("topright",
    legend(par('usr')[2], par('usr')[4], xpd=NA,
           legend = legendT, cex=CEX,
           bty = "n",
           lty=Lty, # gives the legend appropriate symbols (lines)
           lwd=rep(2, q*q),
           col=Col # gives the legend lines the correct color and width
    )
  }

  if(is.complex(eigen(Drift)$val)){
    # Multiple solutions and add 3 plots (2 for 2 different solutions and one scatter plot)
    EigenDrift <- eigen(Drift)
    V <- EigenDrift$vector
    title_N <- as.list(expression(paste(Phi(Delta[t]), " plot"), "using an 'aliasing' matrix", "(i.e., another solution for A)"))
    gg_title_N <- as.list(expression(paste(Phi(Delta[t]), " plot"), "\nusing an 'aliasing' matrix \n(i.e., another solution for A)"))

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
      A_N = Drift + (2 * base::pi * im / 1) * V %*% diagN %*% solve(V) # Here DeltaT=1, because A is input and thus DeltaT does not effect this
      #A_N
      #print(A_N)
      Drift_N <- Re(A_N)
      PhiDeltaTs_N<-array(data=NA,dim=c(q,q,length(DeltaTs)))
      for(i in 1:length(DeltaTs)){
        PhiDeltaTs_N[,,i]<-expm(Drift_N*DeltaTs[i])
      }

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
      plot(y=rep(0, length(DeltaTs)), x=DeltaTs, type="l", ylim=c(min(PhiDeltaTs_N), max(PhiDeltaTs_N)),
           ylab = expression(paste(Phi(Delta[t]), " values")), xlab = expression(paste("Time-interval (", Delta[t], ")", sep="")),
           col=1000, lwd=2, lty=1,
           main=mtext(do.call(expression, title_N), side=3, line = c(2,1,0), cex = 1 )
           #"Effect lag curve: \n How do the VAR(1) parameters Phi vary \n as a function of the time-interval"
      )
      #
      teller <- 0
      for(j in 1:q){
        for(i in 1:q){
          if(WhichElements[j,i] == 1){
            teller <- teller + 1
            lines(y=PhiDeltaTs_N[j,i,], x=DeltaTs, col=Col[teller], lwd=2, lty=Lty[teller])
          }
        }
      }
    }
    # In case last plot is scatter plot
    # In last plot a scatter plot, for multiples of DeltaT, from Min to Max.
    Min_ <- Min + Min%%DeltaT # last part is remainder after integer division
    Max_ <- Max - Max%%DeltaT # last part is remainder after integer division
    DeltaTs <- seq(Min_, Max_, by=DeltaT)
    PhiDeltaTs_N<-array(data=NA,dim=c(q,q,length(DeltaTs)))
    for(i in 1:length(DeltaTs)){
      PhiDeltaTs_N[,,i]<-expm(Drift_N*DeltaTs[i])
    }

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
    #
    title_N <- as.list(
      expression(paste(Phi(Delta[t]), " scatter plot for multiples of ", Delta ["t"]),
                 paste("(Note: For multiples of ", Delta ["t"], ", "),
                 paste(Phi(Delta[t]), " is unique)")
      )
    )
    gg_title_N2 <- as.list(
      expression(paste(Phi(Delta[t]), " scatter plot for multiples of ", Delta ["t"]),
                 paste("\n(Note: For multiples of ", Delta ["t"], ", "),
                 paste(Phi(Delta[t]), " is unique)")
      )
    )

    plot(y=rep(0, length(DeltaTs)), x=DeltaTs, type="l", ylim=c(min(PhiDeltaTs_N), max(PhiDeltaTs_N)),
         ylab = expression(paste(Phi(Delta[t]), " values")), xlab = expression(paste("Time-interval (", Delta[t], ")", sep="")),
         col=1000, lwd=2, lty=1,
         main=mtext(do.call(expression, title_N), side=3, line = c(2,1,0), cex = 1 )
         #"Effect lag curve: \n How do the VAR(1) parameters Phi vary \n as a function of the time-interval"
    )
    #
    teller <- 0
    for(j in 1:q){
      for(i in 1:q){
        if(WhichElements[j,i] == 1){
          teller <- teller + 1
          points(y=PhiDeltaTs_N[j,i,], x=DeltaTs, col=Col[teller], lwd=2, pch=Lty[teller])
        }
      }
    }

  for (i in 1:2) {
    plot_list[[length(plot_list) + 1]] <- ggplot(PhiDeltaTsDF_L[[i]], aes(DeltaTs, Values, color = Labels, linetype = Labels)) +
      geom_line(lwd = 0.75) +
      geom_abline(intercept = 0, slope = 0, alpha = .5) +
      scale_linetype_manual(name = " ", values = Lty, labels = legendT) +
      scale_color_manual(name = " ", values = Col, labels = legendT) +
      ylab(expression(paste(Phi(Delta[t]), " values"))) +
      xlab(expression(paste("Time-interval (", Delta[t], ")", sep = ""))) +
      labs(title = as.expression(gg_title_N[1]),
           subtitle = as.expression(gg_title_N[-1])) +
      theme_classic() +
      theme(plot.title = element_text(margin = margin(t = 20))) +
      ylim(0,1)

    if (i == 2) {
      plot_list[[length(plot_list)]] <- plot_list[[length(plot_list)]] + theme(legend.position = "none")
    }
  }

  plot_list[[length(plot_list) + 1]] <- ggplot(PhiDeltaTsDF_4, aes(DeltaTs, Values, color = Labels, shape = Labels)) +
    geom_point(show.legend = ) +
    geom_abline(intercept = 0, slope = 0, alpha = .5) +
    scale_shape_manual(name = " ", values = Lty, labels = legendT) +
    scale_color_manual(name = " ", values = Col, labels = legendT) +
    ylab(expression(paste(Phi(Delta[t]), " values"))) +
    xlab(expression(paste("Time-interval (", Delta[t], ")", sep = ""))) +
    labs(title = as.expression(gg_title_N2[1]),
         subtitle = as.expression(gg_title_N2[-1])) +
    theme_classic() +
    theme(plot.title = element_text(margin = margin(t = 20))) +
    ylim(0,1)

  ggarrange(plotlist = plot_list, ncol = 2, nrow = 2,
            widths = c(3,4)) %>% show
  }

  if(complex == TRUE){
    if(q<4){CEX = 1}else{CEX = (1.4-q/10)} # Check for optimal values!
    par(mar = c(0,0,0,0))
    plot.new()
    legend("center",
           legend = legendT, cex=CEX,
           bty = "n",
           lty=Lty, # gives the legend appropriate symbols (lines)
           lwd=rep(2, q*q),
           col=Col # gives the legend lines the correct color and width
    )
    plot.new()
    legend("center",
           legend = legendT, cex=CEX,
           bty = "n",
           pch=Lty, # gives the legend appropriate symbols (lines)
           #lwd=rep(2, q*q),
           col=Col # gives the legend lines the correct color and width
    )
  }

  #dev.off()

  par(op)


  ############################################################################################################

  #final <- list(.. = ...)
  #return(final)

}


