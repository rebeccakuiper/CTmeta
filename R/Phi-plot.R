#' Phi-plot of Phi based on its underling drift matrix
#'
#' This function makes a Phi-plot of Phi(DeltaT) for a range of time intervals based on its underling drift matrix. There is also an interactive web application on my website to create a Phi-plot: Phi-and-Psi-Plots and Find DeltaT (\url{https://www.uu.nl/staff/RMKuiper/Websites\%20\%2F\%20Shiny\%20apps}).
#'
#' @param DeltaT The time interval used. By default, DeltaT = 1.
#' @param Drift Underling continuous-time lagged effects matrix (i.e., Drift matrix) of the discrete-time lagged effects matrix Phi(DeltaT). Note that the Drift (Phi) matrix should be standardized to make a fair comparison between cross-legged effects.
#' It also takes a fitted object from the classes "varest" (from the VAR() function in vars package) and "ctsemFit" (from the ctFit() function in the ctsem package); see example below. From such an object, the (standardized) Drift matrix is calculated/extracted.
#' @param Min Minimum time interval used in the Phi-plot. By default, Min = 0.
#' @param Max Maximum time interval used in the Phi-plot. By default, Max = 10.
#' @param Step The step-size taken in the time intervals. By default, Step = 0.05. Hence, using the defaults, the Phi-plots is based on the values of Phi(DeltaT) for DeltaT = 0, 0.05, 0.10, ..., 10. Note: Especially in case of complex eigenvalues, this step size should be very small (then, the oscillating behaviour can be seen best).
#' @param WhichElements Optional. Matrix of same size as Drift denoting which element/line should be plotted (1) or not (0). By default, WhichElements = NULL. Note that even though not all lines have to be plotted, the full Drift matrix is needed to determine the selected lines.
#' @param Labels Optional. Vector with (character) labels of the lines to be plotted. The length of this vector equals the number of 1s in WhichElements (or equals q*q). By default, Labels = NULL, which renders labels with Greek letter of Phi (as a function of the time-interval) together with the indices (of outcome and predictor variables).
#' @param Col Optional. Vector with color values (integers) of the lines to be plotted. The length of this vector equals the number of 1s in WhichElements (or equals q*q). By default, Col = NULL, which renders the same color for effects that belong to the same outcome variable (i.e. a row in the Drift matrix). See \url{https://www.statmethods.net/advgraphs/parameters.html} for more information about the values.
#' @param Lty Optional. Vector with line type values (integers) of the lines to be plotted. The length of this vector equals the number of 1s in WhichElements (or equals q*q). By default, Lty = NULL, which renders solid lines for the autoregressive effects and the same type of dashed line for reciprocal effects (i.e., same type for Phi_ij as for Phi_ji). See \url{https://www.statmethods.net/advgraphs/parameters.html} for more information about the values.
#' @param Title Optional. A character or a list consisting of maximum 3 characters (or 'call' class objects, like from the paste() function) that together represent the title of the Phi-plot. By default, Title = NULL, then the following code will be used for the title: as.list(expression(paste(Phi(Delta[t]), " plot:"), "How do the lagged parameters vary", "as a function of the time-interval")))
#'
#' @return This function returns a Phi-plot for a range of time intervals.
#' @importFrom expm expm
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
#' PhiPlot(DeltaT, Drift)
#' PhiPlot(DeltaT, Drift, Min = 0, Max = 10, Step = 0.01)
#'
#'
#' ## Example 2: input from fitted object of class "varest" ##
#' #
#' DeltaT <- 1
#' data <- myData
#' if (!require("vars")) install.packages("vars")
#' library(vars)
#' out_VAR <- VAR(data, p = 1)
#' PhiPlot(DeltaT, out_VAR)
#'
#'
#' ## Example 3: Change plot options ##
#' # Note: use Drift from Example 1
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
#' PhiPlot(DeltaT = 1, Drift, Min = 0, Max = 10, Step = 0.5, WhichElements, Labels, Col, Lty)
#'


PhiPlot <- function(DeltaT = 1, Drift, Min = 0, Max = 10, Step = 0.05, WhichElements = NULL, Labels = NULL, Col = NULL, Lty = NULL, Title = NULL) {

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
  # Check on Drift
  if(class(Drift) == "varest"){
    Phi_VARest <- Acoef(Drift)[[1]]
    Drift <- logm(Phi_VARest)/DeltaT # Phi = expm(Drift * deltaT)
    q <- dim(Drift)[1]
    # TO DO bepaal standardized Phi en dus Drift!
  } else if(class(Drift) == "ctsemFit"){
    Drift <- summary(Drift)$DRIFT
    q <- dim(Drift)[1]
    # TO DO bepaal standardized Drift!
  } else{
    if(length(Drift) == 1){
      q <- 1
    }else{
      #
      if(is.null(dim(Drift))){
        if(!is.null(length(Drift))){
          print(paste("The argument Drift is not a matrix of size q times q."))
          stop()
        }else{
          print(paste("The argument Drift is not found: The continuous-time lagged effects matrix Drift is unknown, but should be part of the input."))
          stop()
        }
      }else{
        if(dim(Drift)[1] != dim(Drift)[2] | length(dim(Drift)) != 2){
          print(paste("The argument Drift is not a matrix of size q times q."))
          stop()
        }
        q <- dim(Drift)[1]
      }
    }
  }
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
    #title <- as.list(expression(paste(Phi(Delta[t]), " plot:"), "How do the overall lagged parameters vary", "as a function of the time-interval"))
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
           legend = legendT, cex=CEX,
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
  #
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
    #
    title_N <- as.list(
      expression(paste(Phi(Delta[t]), " scatter plot for multiples of ", Delta ["t"]),
                 paste("(Note: For multiples of ", Delta ["t"], ", "),
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
  } # end if complex


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


