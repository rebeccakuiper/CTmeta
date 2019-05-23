
#' Makes Phi-plot of Phi based on its underling drift matrix
#'
#' @param DeltaT The time interval used.
#' @param Drift Underling drift matrix of the Phi; e.g., the overall Phi obtained from CT meta-analysis.
#' @param Min Minimum time interval used in the plot. By default, Min = 0.
#' @param Max Maximum time interval used in the plot. By default, Max = 100.
#' @param Step The step-size taking in the time intervals. By default, Step = 0.5. Hence, using the defaults, the values of Phi(DeltaT) are determined for 0, 0.5, 1, 1.5, ..., 100.
#'
#' @return Phi-plot.
#' @importFrom expm expm
#' @export
#' @examples
#'
#' # Make PhiPlot of overallPhi obtained from CTmeta (with the CTMA function)
#' # Input for CTMA to obtain overallPhi
#' Phi <- matrix(c(0.25, 0.10,
#'                 0.20, 0.36,
#'                 0.35, 0.20,
#'                 0.30, 0.46,
#'                 0.15, 0.00,
#'                 0.10, 0.26), byrow=T, ncol = q)
#' SigmaVAR_s <- diag(q) # for ease
#' SigmaVAR <- rbind(SigmaVAR_s, SigmaVAR_s, SigmaVAR_s)
#' # If Phi and SigmaVAR are known, one can calculate Gamma:
#' Gamma <- array(data=NA, dim=c(S*q,q))
#' teller <- 1
#' for(s in 1:S){
#'   Gamma[teller:(teller+1),] <- calc.Gamma.fromVAR(Phi[teller:(teller+1),], SigmaVAR[teller:(teller+1),])
#'   teller <- teller + q
#' }
#' DeltaT <- c(2, 3, 1)
#' N <- c(643, 651, 473)
#' DeltaTStar <- 1
#' out_CTmeta <- CTMA(N, 0, Phi, SigmaVAR, Gamma, DeltaTStar, DeltaT, Moderators = 0, Mod = NULL, FEorRE = 1, alpha=0.05)
#' overallPhi <- matrix(out_CTmeta$Overall_standPhi_DeltaTStar, byrow = T, ncol = q)
#' overallDrift <- logm(overallPhi)/DeltaTStar # Use expm package
#'
#' # Make plot of above obtained overallPhi
#' PhiPlot(DeltaTStar, overallDrift, Min = 0, Max = 40, Step = 0.5)




PhiPlot <- function(DeltaT, Drift, Min = 0, Max = 100, Step = 0.5) {

#  #######################################################################################################################
#
#  #if (!require("expm")) install.packages("expm")
#  library(expm)
#
#  #######################################################################################################################


  #def.par <- par(no.readonly = TRUE) # save default, for resetting...
  #par(def.par)  #- reset to default

  q <- dim(Drift)[1]

  DeltaTs<-seq(Min,Max,by=Step)

  if(any(is.complex(eigen(Drift)$val))){
    # Multiple solutions, then 2x2 plots
    op <- par(mfrow=c(2,2))
    complex = TRUE
    nf <- layout(matrix(c(1,2,5,3,4,6),2,3,byrow = TRUE), c(3,3,1), c(2,2,1), TRUE)
    #layout.show(nf)
  } else{
    op <- par(mfrow=c(1,1))
    complex = FALSE
  }


  PhiDeltaTs<-array(data=NA,dim=c(q,q,length(DeltaTs)))
  for(i in 1:length(DeltaTs)){
    PhiDeltaTs[,,i]<-expm(Drift*DeltaTs[i])
  }
  tellerCol = 0
  tellerLTY = 1
  j = 1
  teller = 0
  Col = matrix(NA, ncol = q, nrow = q)
  Lty = matrix(NA, ncol = q, nrow = q)
  title <- as.list(expression(paste(Phi(Delta[t]), " plot:"), "How do the overall lagged parameters vary", "as a function of the time-interval"))

  #wd <- getwd()
  #dev.copy(png, filename = paste0(wd, "/www/PhiPlot.png"))
  plot(y=PhiDeltaTs[1,1,], x=DeltaTs, type="l", ylim=c(min(PhiDeltaTs), max(PhiDeltaTs)),
       ylab = expression(paste("Overall", Phi(Delta[t]), " values")), xlab = expression(paste("Time-interval (", Delta[t], ")", sep="")),
       col=(tellerCol+1), lwd=2, lty=tellerLTY,
       main=mtext(do.call(expression, title), side=3, line = c(2,1,0), cex = 1 )
       #"Effect lag curve: \n How do the VAR(1) parameters Phi vary \n as a function of the time-interval"
  )
  #
  Col[1,1] = (tellerCol+1)
  Lty[1,1] = tellerLTY
  #
  for(i in 2:q){
    teller = teller + 1
    lines(y=PhiDeltaTs[j,i,],x=DeltaTs, col=(tellerCol+j), lwd=2, lty=(tellerLTY+teller))
    lines(y=PhiDeltaTs[i,j,],x=DeltaTs, col=(tellerCol+i), lwd=2, lty=(tellerLTY+teller))
    #
    Col[j,i] = (tellerCol+j)
    Lty[j,i] = (tellerLTY+teller)
    Col[i,j] = (tellerCol+i)
    Lty[i,j] = (tellerLTY+teller)
  }
  for(j in 2:q){
    for(i in j:q){
      if(i == j){
        lines(y=PhiDeltaTs[j,i,],x=DeltaTs, col=(tellerCol+j), lwd=2, lty=tellerLTY)
        #
        Col[j,i] = (tellerCol+j)
        Lty[j,i] = tellerLTY
      } else{
        teller = teller + 1
        lines(y=PhiDeltaTs[j,i,],x=DeltaTs, col=(tellerCol+j), lwd=2, lty=(tellerLTY+teller))
        lines(y=PhiDeltaTs[i,j,],x=DeltaTs, col=(tellerCol+i), lwd=2, lty=(tellerLTY+teller))
        #
        Col[j,i] = (tellerCol+j)
        Lty[j,i] = (tellerLTY+teller)
        Col[i,j] = (tellerCol+i)
        Lty[i,j] = (tellerLTY+teller)
      }
    }
  }

  #lines(y=rep(0, length(DeltaTs)),x=DeltaTs, col=1000, lwd=2, lty=1)


  subscripts = NULL
  for(i in 1:q){
    subscripts = c(subscripts, paste(i, 1:q, sep=""))
  }
  legendT = NULL
  for(i in 1:(q*q)){
    e <- bquote(expression(Phi(Delta[t])[.(subscripts[i])]))
    legendT = c(legendT, eval(e))
  }



  if(complex == FALSE){
    if(q<4){CEX = 1}else{CEX = (1.4-q/10)} # Check for optimal values!
    legend("topright",
           legend = legendT, cex=CEX,
           bty = "n",
           lty=as.vector(t(Lty)), # gives the legend appropriate symbols (lines)
           lwd=rep(2, q*q),
           col=as.vector(t(Col)) # gives the legend lines the correct color and width
    )
  }


  if(is.complex(eigen(Drift)$val)){
    # Multiple solutions and add 3 plots for 3 different solutions
    EigenDrift <- eigen(Drift)
    V <- EigenDrift$vector
    #tellerCol_N <- 0 # if in one plot
    title_N <- as.list(expression(paste(Phi(Delta[t]), " plot"), "using an 'aliasing' matrix", "(i.e., another solution for A)"))
    #for(N in 1:3){
    for(N in 1:2){ # In case last plot is scatter plot
      tellerCol_N <- 0 # if seperate plots
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
      Col_N = matrix(NA, ncol = q, nrow = q)
      Lty_N = matrix(NA, ncol = q, nrow = q)
      #tellerCol_N <- tellerCol_N + q*q # if in one plot
      #
      # In case of seperate plots:
      plot(y=PhiDeltaTs_N[1,1,], x=DeltaTs, type="l", ylim=c(min(PhiDeltaTs_N), max(PhiDeltaTs_N)),
           ylab = expression(paste(Phi(Delta[t]), " values")), xlab = expression(paste("Time-interval (", Delta[t], ")", sep="")),
           col=(tellerCol_N+1), lwd=2, lty=tellerLTY,
           main=mtext(do.call(expression, title_N), side=3, line = c(2,1,0), cex = 1 )
           #"Effect lag curve: \n How do the VAR(1) parameters Phi vary \n as a function of the time-interval"
      )
      #
      Col_N[1,1] = (tellerCol_N+1)
      Lty_N[1,1] = tellerLTY
      #
      j = 1
      teller = 0
      for(i in 2:q){
        teller = teller + 1
        lines(y=PhiDeltaTs_N[j,i,],x=DeltaTs, col=(tellerCol_N+j), lwd=2, lty=(tellerLTY+teller))
        lines(y=PhiDeltaTs_N[i,j,],x=DeltaTs, col=(tellerCol_N+i), lwd=2, lty=(tellerLTY+teller))
        #
        Col_N[j,i] = (tellerCol_N+j)
        Lty_N[j,i] = (tellerLTY+teller)
        Col_N[i,j] = (tellerCol_N+i)
        Lty_N[i,j] = (tellerLTY+teller)
      }
      for(j in 2:q){
        #
        # In case of adding to one plot
        #for(j in 1:q){
        #
        #
        for(i in j:q){
          if(i == j){
            lines(y=PhiDeltaTs_N[j,i,],x=DeltaTs, col=(tellerCol_N+j), lwd=2, lty=tellerLTY)
            #
            Col_N[j,i] = (tellerCol_N+j)
            Lty_N[j,i] = tellerLTY
          } else{
            teller = teller + 1
            lines(y=PhiDeltaTs_N[j,i,],x=DeltaTs, col=(tellerCol_N+j), lwd=2, lty=(tellerLTY+teller))
            lines(y=PhiDeltaTs_N[i,j,],x=DeltaTs, col=(tellerCol_N+i), lwd=2, lty=(tellerLTY+teller))
            #
            Col_N[j,i] = (tellerCol_N+j)
            Lty_N[j,i] = (tellerLTY+teller)
            Col_N[i,j] = (tellerCol_N+i)
            Lty_N[i,j] = (tellerLTY+teller)
          }
        }
      }
      #
      #
      # In case of adding to one plot:
      #Col <- c(Col, Col_N)
      #Lty <- c(Lty, Lty_N)
      #legendT = c(legendT, expression(Phi_N(Delta[t])[11]), expression(Phi_N(Delta[t])[12]), expression(Phi_N(Delta[t])[21]), expression(Phi_N(Delta[t])[22]))
      #
      ## In case of seperate plots:
      ##legendT = c(expression(Phi_N(Delta[t])[11]), expression(Phi_N(Delta[t])[12]), expression(Phi_N(Delta[t])[21]), expression(Phi_N(Delta[t])[22]))
      #legend("topright",
      #       legend = legendT,
      #       bty = "n",
      #       lty=as.vector(t(Lty_N)), # gives the legend appropriate symbols (lines)
      #       lwd=rep(2, q*q),
      #       col=as.vector(t(Col_N)) # gives the legend lines the correct color and width
      #)
      #
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
    tellerCol_N <- 0 # if seperate plots
    Col_N = matrix(NA, ncol = q, nrow = q)
    Lty_N = matrix(NA, ncol = q, nrow = q)
    plot(y=PhiDeltaTs_N[1,1,], x=DeltaTs, type="p", ylim=c(min(PhiDeltaTs_N), max(PhiDeltaTs_N)),
         ylab = expression(paste(Phi(Delta[t]), " values")), xlab = expression(paste("Time-interval (", Delta[t], ")", sep="")),
         col=(tellerCol_N+1), lwd=2, pch=tellerLTY,
         main=mtext(do.call(expression, title_N), side=3, line = c(2,1,0), cex = 1 )
         #"Effect lag curve: \n How do the VAR(1) parameters Phi vary \n as a function of the time-interval"
    )
    #
    Col_N[1,1] = (tellerCol_N+1)
    Lty_N[1,1] = tellerLTY
    #
    j = 1
    teller = 0
    for(i in 2:q){
      teller = teller + 1
      points(y=PhiDeltaTs_N[j,i,],x=DeltaTs, col=(tellerCol_N+j), lwd=2, pch=(tellerLTY+teller))
      Col_N[j,i] = (tellerCol_N+j)
      Lty_N[j,i] = (tellerLTY+teller)
      #
      teller = teller + 1
      points(y=PhiDeltaTs_N[i,j,],x=DeltaTs, col=(tellerCol_N+i), lwd=2, pch=(tellerLTY+teller))
      Col_N[i,j] = (tellerCol_N+i)
      Lty_N[i,j] = (tellerLTY+teller)
    }
    for(j in 2:q){
      #
      # In case of adding to one plot
      #for(j in 1:q){
      #
      #
      for(i in j:q){
        if(i == j){
          points(y=PhiDeltaTs_N[j,i,],x=DeltaTs, col=(tellerCol_N+j), lwd=2, pch=tellerLTY)
          #
          Col_N[j,i] = (tellerCol_N+j)
          Lty_N[j,i] = tellerLTY
        } else{
          teller = teller + 1
          points(y=PhiDeltaTs_N[j,i,],x=DeltaTs, col=(tellerCol_N+j), lwd=2, pch=(tellerLTY+teller))
          Col_N[j,i] = (tellerCol_N+j)
          Lty_N[j,i] = (tellerLTY+teller)
          #
          teller = teller + 1
          points(y=PhiDeltaTs_N[i,j,],x=DeltaTs, col=(tellerCol_N+i), lwd=2, pch=(tellerLTY+teller))
          Col_N[i,j] = (tellerCol_N+i)
          Lty_N[i,j] = (tellerLTY+teller)
        }
      }
    }
    #legend("topright",
    #       legend = legendT,
    #       bty = "n",
    #       pch=as.vector(t(Lty_N)), # gives the legend appropriate symbols (lines)
    #       #lwd=rep(2, q*q),
    #       col=as.vector(t(Col_N)) # gives the legend lines the correct color and width
    #)
  } # end if complex



  if(complex == TRUE){
    if(q<4){CEX = 1}else{CEX = (1.4-q/10)} # Check for optimal values!
    par(mar = c(0,0,0,0))
    plot.new()
    legend("center",
           legend = legendT, cex=CEX,
           bty = "n",
           lty=as.vector(t(Lty)), # gives the legend appropriate symbols (lines)
           lwd=rep(2, q*q),
           col=as.vector(t(Col)) # gives the legend lines the correct color and width
    )
    plot.new()
    legend("center",
           legend = legendT, cex=CEX,
           bty = "n",
           pch=as.vector(t(Lty_N)), # gives the legend appropriate symbols (lines)
           #lwd=rep(2, q*q),
           col=as.vector(t(Col_N)) # gives the legend lines the correct color and width
    )
  }

  #dev.off()

  par(op)





  ############################################################################################################

  #final <- list(.. = ...)
  #return(final)

}


