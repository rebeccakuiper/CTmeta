
#DeltaT <- 1
#kronsum <- kronecker(diag(q),B) + kronecker(B,diag(q))
#SigmaVAR <- matrix((solve(kronsum) %*% (diag(q*q) - expm(-kronsum * DeltaT)) %*% as.vector(Sigma)), ncol=q, nrow=q)
# So, we do not need all the stuff below!


calc.VARparam <- function(DeltaT, B, Sigma) {
# DeltaT = time interval
# q = Number of dimensions, that is, number of dependent variables


#######################################################################################################################

if (!require("expm")) install.packages("expm")
library(expm)

#######################################################################################################################

  if(any(class(B) == "ctsemFit")){
    B <- -1 * summary(B)$DRIFT
    Sigma <- summary(B)$DIFFUSION
  }

  if(length(B) == 1){
    q <- 1
  }else{
    q <- dim(B)[1]
  }

# Phi = exp{-B * DeltaT} = V exp{-Eigenvalues * DeltaT} V^-1 = V D_Phi V^-1
ParamVAR <- expm(-B*DeltaT)
#if(all(abs(Im(ParamVAR)) < 0.0001) == TRUE){ParamVAR <- Re(ParamVAR)}

kronsum <- kronecker(diag(q),B) + kronecker(B,diag(q))
Sigma_VAR <- matrix((solve(kronsum) %*% (diag(q*q) - expm(-kronsum * DeltaT)) %*% as.vector(Sigma)), ncol=q, nrow=q)
#if(all(abs(Im(Sigma_VAR)) < 0.0001) == TRUE){Sigma_VAR <- Re(Sigma_VAR)}


############################################################################################################

final <- list(Phi = ParamVAR, SigmaVAR = Sigma_VAR)
return(final)

}

