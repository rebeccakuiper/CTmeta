% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Calc VARparam from CTMparam.r
\name{VARparam}
\alias{VARparam}
\title{Discrete-time estimates from continuous-time estimates}
\usage{
VARparam(DeltaT = 1, Drift, Sigma = NULL, Gamma = NULL)
}
\arguments{
\item{DeltaT}{Optional. The time interval used. By default, DeltaT = 1.}

\item{Drift}{Matrix of size q times q of (un)standardized continuous-time lagged effects, called drift matrix. Note that Phi(DeltaT) = expm(Drift*DeltaT). By default, input for Phi is used; only when Phi = NULL, Drift will be used.
It also takes a fitted object from the class "ctsemFit" (from the ctFit() function in the ctsem package); see example below. From such an object, the (standardized) Drift and Sigma matrices are extracted.}

\item{Sigma}{Residual covariance matrix of the first-order continuous-time (CT-VAR(1)) model, that is, the diffusion matrix. By default, input for Sigma is used; only when Sigma = NULL, Gamma will be used.}

\item{Gamma}{Optional (either Sigma or Gamma). Stationary covariance matrix, that is, the contemporaneous covariance matrix of the data. By default, input for Sigma is used; only when Sigma = NULL, Gamma will be used.
Note that if Drift and Sigma are known, Gamma can be calculated; hence, either Sigma or Gamma is needed as input.}
}
\value{
The output renders the discrete-times equivalent matrices of the continuous-time ones.
}
\description{
The discrete-time lagged-effects model matrices corresponding to the continuous-time ones. The interactive web application 'Phi-and-Psi-Plots and Find DeltaT' also contains this functionality, you can find it on my website: \url{https://www.uu.nl/staff/RMKuiper/Websites\%20\%2F\%20Shiny\%20apps}.
}
\examples{

# library(CTmeta)

## Example 1 ##

##################################################################################################
# Input needed in examples below with q=2 variables.
# I will use the example matrix stored in the package:
DeltaT <- 1
Drift <- myDrift
q <- dim(Drift)[1]
Sigma <- diag(q) # for ease
#
Gamma <- Gamma.fromCTM(Drift, Sigma)
##################################################################################################

DeltaT <- 1
VARparam(DeltaT, Drift, Sigma)
# or
VARparam(DeltaT, Drift, Gamma = Gamma)


## Example 2: input from fitted object of class "ctsemFit" ##
#
#data <- myData
#if (!require("ctsemFit")) install.packages("ctsemFit")
#library(ctsemFit)
#out_CTM <- ctFit(...)
#
#DeltaT <- 1
#VARparam(DeltaT, out_CTM)

}
