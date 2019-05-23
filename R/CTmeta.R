
#' Transforms the lagged effect estimates for a given interval to the ones corresponding to a different time-interval
#'
#' @param N Number of persons (panel data) or measurement occasions - 1 (time series data) used in the S studies. Matrix of size S times 1.
#' @param TypeMx Indicator of type of matrix for Phi, SigmaVAR, and Gamma: 0) Stacked matrix of size S q times q; 2) Array with dimensions q times q times S; with S the number of studies in the meta-analysis and q the number of variables (leading to an q times q lagged effects matrix Phi).
#' @param Phi Stacked matrix (TypeMx == 0) or array (TypeMx == 1) of (un)standardized lagged effects matrices for all S studies in the meta-analysis.
#' @param SigmaVAR Stacked matrix (TypeMx == 0) or array (TypeMx == 1) of residual covariance matrices of the first-order discrete-time vector autoregressive (DT-VAR(1)) model.
#' @param Gamma Stacked matrix (TypeMx == 0) or array (TypeMx == 1) of stationary covariance matrices, that is, the contemporaneous covariance matrices of the data sets.
#' @param DeltaTStar The time interval (scalar) to which the standardized lagged effects matrix should be transformed to.
#' @param DeltaT The time intervals used in the S studies. Matrix of size S times 1.
#' @param Moderators Indicator whether there are moderators to be included (1) or not (0; default).
#' @param Mod Matrix of moderators to be included in the analysis when Moderators == 1. By default, Mod = NULL.
#' @param FEorRE Indicator whether continuous-time meta-analysis (CTmeta) should use a fixed-effects model (1; default) or random-effects model (2).
#' @param alpha The alpha level in determining the (1-alpha)*100 percent confidence interval (CI). By default, alpha is set to 0.05, resulting in a 95 percent CI.
#'
#' @return vectorized transformed standardized lagged effects, their covariance matrix, and the corresponding elleptical/multivariate 95 percent CI.
#' @importFrom expm expm
#' @importFrom fastDummies dummy_cols
#' @importFrom metafor rma.mv
#' @importFrom metafor rma.uni
#' @export
#' @examples
#'
#' ###########################################################
#' ### Input for examples ###
#' q = 2 # Number of variables in the process/system
#' S = 3 # Number of studies included in the meta-analysis
#'
#' ## Phi as S q times q matrix (TypeMx == 0) ##
#'
#' # Based on discrete-time model #
#' Phi <- matrix(c(0.25, 0.10,
#'                 0.20, 0.36,
#'                 0.35, 0.20,
#'                 0.30, 0.46,
#'                 0.15, 0.00,
#'                 0.10, 0.26), byrow = T, ncol = q)
#'
#' # Based on correlation matrix #
#' corr_YXYX <- matrix(c(1.00, 0.40, 0.63, 0.34,
#'                       0.40, 1.00, 0.31, 0.63,
#'                       0.63, 0.31, 1.00, 0.41,
#'                       0.34, 0.63, 0.41, 1.00), byrow = T, ncol = 2*q)
#' out <- calc.TransPhi_Corr(12, 24, 2235, corr_YXYX)
#' Phi_1 <- matrix(out$vecStandPhi_DeltaTStar, byrow = T, ncol = q)
#' Phi_2 <- matrix(out$vecStandPhi_DeltaTStar, byrow = T, ncol = q)
#' Phi_3 <- matrix(out$vecStandPhi_DeltaTStar, byrow = T, ncol = q)
#' Phi <- rbind(Phi_1, Phi_2, Phi_3)
#'
#'
#' ## SigmaVAR and Gamma as S q times q matrix (TypeMx == 0) ##
#' # Based on discrete-time model, using SigmaVAR #
#' SigmaVAR_s <- diag(q) # for ease
#' SigmaVAR <- rbind(SigmaVAR_s, SigmaVAR_s, SigmaVAR_s)
#' # If Phi and SigmaVAR are known, one can calculate Gamma:
#' Gamma <- array(data=NA, dim=c(S*q,q))
#' teller <- 1
#' for(s in 1:S){
#'   Gamma[teller:(teller+1),] <- calc.Gamma.fromVAR(Phi[teller:(teller+1),], SigmaVAR[teller:(teller+1),])
#'   teller <- teller + q
#' }
#'
#' # Based on discrete-time model, using Gamma #
#' # Use Gamma from above
#' # If Phi and SigmaVAR are known, one can calculate Gamma:
#' SigmaVAR <- array(data=NA, dim=c(S*q,q))
#' teller <- 1
#' for(s in 1:S){
#'   SigmaVAR[teller:(teller+1),] <- Gamma[teller:(teller+1),] - Phi[teller:(teller+1),] %*% Gamma[teller:(teller+1),] %*% t(Phi[teller:(teller+1),])
#'   teller <- teller + q
#' }
#'
#' # Based on correlation matrix #
#' # Use 'out' from above. # out <- calc.TransPhi_Corr(12, 24, 2235, corr_YXYX)
#' SigmaVAR_1 <- out$SigmaVAR_DeltaTStar
#' Gamma_1 <- out$Gamma
#' # Do this for each study to make a stacked matrix out of this
#'
#' ## Phi as array with dimensions q times q times S (TypeMx == 1) ##
#'
#' # Based on discrete-time model #
#' # Use Phi from above
#' Phi_studies <- array(data=NA, dim=c(q,q,S))
#' teller = 0
#' for(s in 1:S){
#'   Phi_studies[1:q,1:q,s] <- matrix(t(Phi)[(teller+1):(teller+q*q)], byrow = T, ncol=q*q)
#'   Phi_studies[1:q,1:q,s] <- t(Phi_studies[1:q,1:q,s])
#'   teller = teller + q*q
#' }
#' #Phi <- Phi_studies
#'
#' # Use SigmaVAR from above
#' SigmaVAR_studies <- array(data=NA, dim=c(q,q,S))
#' teller = 0
#' for(s in 1:S){
#'   SigmaVAR_studies[1:q,1:q,s] <- SigmaVAR[(teller+1):(teller+q),1:q]
#'   teller = teller + q
#' }
#' #SigmaVAR <- SigmaVAR_studies
#' Gamma_studies <- array(data=NA, dim=c(q,q,S))
#' for(s in 1:S){
#'   Gamma_studies[,,s] <- calc.Gamma.fromVAR(Phi_studies[,,s], SigmaVAR_studies[,,s])
#' }
#' #Gamma <- CovMx_studies
#'
#' # Use Gamma from above
#' Gamma_studies <- array(data=NA, dim=c(q,q,S))
#' teller = 0
#' for(s in 1:S){
#'   Gamma_studies[1:q,1:q,s] <- Gamma[(teller+1):(teller+q),1:q]
#'   teller = teller + q
#' }
#' #Gamma <- CovMx_studies
#' SigmaVAR_studies <- array(data=NA, dim=c(q,q,S))
#' for(s in 1:S){
#'   SigmaVAR_studies[,,s] <- Gamma_studies[,,s] - Phi_studies[,,s] %*% Gamma_studies[,,s] %*% t(Phi_studies[,,s])
#' }
#' #SigmaVAR <- SigmaVAR_studies
#'
#'
#' # The study-specific sampel size (N) and time interval (DeltaT)
#' DeltaT <- matrix(c(2, 3, 1))
#' N <- matrix(c(643, 651, 473))
#' ###########################################################
#'
#'
#' # Based on input above, one can run one of these two CTmeta analyses
#' DeltaTStar <- 1
#' CTMA(N, 1, Phi_studies, SigmaVAR_studies, Gamma_studies, DeltaTStar, DeltaT)
#' CTMA(N, 0, Phi, SigmaVAR, Gamma, DeltaTStar, DeltaT)
#'
#'
#' # In case of moderators
#' Mod <- matrix(c(64,65,47)) # 1 moderator
#' #Mod <- matrix(cbind(c(64,65,47), c(78,89,34)), ncol = q); colnames(Mod) <- c("Mod1", "Mod2") # two moderators, in each column 1
#' CTMA(N, 1, Phi_studies, SigmaVAR_studies, Gamma_studies, DeltaTStar, DeltaT, Moderators = 1, Mod)
#'
#'
#' # Afterwards: Make PhiPlot of overallPhi
#' DeltaTStar <- 1
#' out_CTmeta <- CTMA(N, 0, Phi, SigmaVAR, Gamma, DeltaTStar, DeltaT, Moderators = 0, Mod = NULL, FEorRE = 1, alpha=0.05)
#' q <- sqrt(length(out_CTmeta$Overall_standPhi_DeltaTStar))
#' overallPhi <- matrix(out_CTmeta$Overall_standPhi_DeltaTStar, byrow = T, ncol = q)
#' overallDrift <- logm(overallPhi)/DeltaTStar # Use expm package
#' PhiPlot(DeltaTStar, overallDrift, Min = 0, Max = 40, Step = 0.5)


CTMA <- function(N, TypeMx, Phi, SigmaVAR, Gamma, DeltaTStar, DeltaT, Moderators = 0, Mod = NULL, FEorRE = 1, alpha=0.05) {

#  #######################################################################################################################
#  #if (!require("expm")) install.packages("expm")
#  library(expm)
#  if (!require("fastDummies")) install.packages("fastDummies")
#  library(fastDummies)
#  if (!require("metafor")) install.packages("metafor")
#  library(metafor)
#  #######################################################################################################################


  S <- length(N) #dim(N)[1]
  q <- dim(Phi)[2]

      if(TypeMx  == 0){
        if(dim(Phi)[1] != S*q | dim(SigmaVAR)[1] != S*q | dim(Gamma)[1] != S*q){
          print("When you state 'TypeMx  = 0', the dimensions of Phi, SigmaVAR, and Gamma should be S*q times q.")
        }
      }else{
        if(dim(Phi)[3] != S | dim(SigmaVAR)[3] != S | dim(Gamma)[3] != S){
          print("When you state 'TypeMx  = 1', the dimensions of Phi, SigmaVAR, and Gamma should be q times q times S.")
        }
      }

      if(TypeMx  == 0){
        Phi_studies <- array(data=NA, dim=c(q,q,S))
        teller <- 0
        for(s in 1:S){
          Phi_studies[1:q,1:q,s] <- matrix(t(Phi)[(teller+1):(teller+q*q)], byrow = T, ncol=q*q)
          Phi_studies[1:q,1:q,s] <- t(Phi_studies[1:q,1:q,s])
          teller <- teller + q*q
        }
        Phi <- Phi_studies
        #
        SigmaVAR_studies <- array(data=NA, dim=c(q,q,S))
        teller <- 0
        for(s in 1:S){
          SigmaVAR_studies[1:q,1:q,s] <- SigmaVAR[(teller+1):(teller+q),1:q]
          teller <- teller + q
        }
        SigmaVAR <- SigmaVAR_studies
        #
        Gamma_studies <- array(data=NA, dim=c(q,q,S))
        teller <- 0
        for(s in 1:S){
          Gamma_studies[1:q,1:q,s] <- Gamma[(teller+1):(teller+q),1:q]
          teller <- teller + q
        }
        Gamma <- Gamma_studies
      }


      # Start with multivariate meta-an (GLS) on transformed variables
      Multivar <- 1
      Trans <- 1
      #
      StudiesComplexEV <- NULL
      StudiesNegEV <- NULL
      #
      StudiesCovMxNotPosDef <- NULL
      #
      #
      vecStandPhi <- array(data=NA, dim=c(q*q,S))
      CovMxPhi <- array(data=NA, dim=c(q*q,q*q,S))
      #Warning <- array(data=NA, dim=c(S))
      #
      # Studies where Phi has complex eigenvalues or at least one negative eigenvalue
      messageTrans <- "All eigenvalues are positive and real; hence, the Phi's are transformed to Phi(DeltaT*) to account for the time-interval dependency."
      for(s in 1:S){
        EV <- eigen(Phi[,,s])$values
        if(any(is.complex(EV))){ StudiesComplexEV <- c(StudiesComplexEV, s) }
        if(any(Re(EV) < 0)){ StudiesNegEV <- c(StudiesNegEV, s) }
      }
      # If ratioDeltaT is integer then also complex or neg eigenvalues work...
      ratioDeltaT <- DeltaTStar / DeltaT
      # If the studies for which the EV are neg or complex have integer ratioDeltaT then it is fine, otherwise it is not
      if(!is.null(StudiesComplexEV) | !is.null(StudiesNegEV)){
        messageTrans <- "Some studies have eigenvalues of Phi that are negative and/or real, but for those studies the ratio DeltaT*/DeltaT is integer; therefore, the Phi's are transformed to Phi(DeltaT*) to account for the time-interval dependency."
        #if(any(!is.integer(ratioDeltaT[StudiesComplexEV]))){ # NB ik moet checken of als ik deel door hoogste geheel getal er dan geen rest waarde is!!
        if( any( (ratioDeltaT[StudiesComplexEV] %% 1) == 0 ) ){
          Trans <- 0
          messageTrans <- "There is at least one study for which the eigenvalues of Phi are complex and for which the ratio DeltaT*/DeltaT is not integer; therefore, dummy variables are used to account for the time-interval dependency."
        }
        if( any( (ratioDeltaT[StudiesNegEV] %% 1) == 0 ) ){
          Trans <- 0
          messageTrans <- "There is at least one study for which the eigenvalues of Phi are negative and for which the ratio DeltaT*/DeltaT is not integer; therefore, dummy variables are used to account for the time-interval dependency."
        }
        if( any( (ratioDeltaT[StudiesComplexEV] %% 1) == 0 ) & any( (ratioDeltaT[StudiesNegEV] %% 1) == 0 ) ){
          Trans <- 0
          messageTrans <- "There is at least one study for which the eigenvalues of Phi are complex and/or negative and for which the ratio DeltaT*/DeltaT is not integer; therefore, dummy variables are used to account for the time-interval dependency."
        }
      }
      #
      #
      messageMultivar <- "For each study, the covariance matrix is positive definite; hence, a multivariate approach is used."
      if(Trans == 1){
        # Calculate standardized transformed Phi
        #s = 0
        #while(s < S & Trans == 1){ # go through all, unless a warning of not unique solution (which I already cover above...)
        #  s = s +1
        for(s in 1:S){
          out <- calc.CovMxStandTransPhi(DeltaTStar, DeltaT[s], N[s], Phi[,,s], Gamma[,,s], SigmaVAR[,,s])
          vecStandPhi[,s] <- out$vecStandPhi_DeltaTStar
          CovMxPhi[,,s] <- out$CovMx_vecStandPhi_DeltaTStar
          #Warning[s] <- out$warning
          #Trans <- 1 - as.numeric(out$Warning)
          #
          # If covariance matrix is not positive definite then do univariate approach (and use only variances)
          if(any(eigen(CovMxPhi[,,s])$val < 0)){
            StudiesCovMxNotPosDef <- c(StudiesCovMxNotPosDef, s) # for which studies is this the case
            Multivar <- 0 # we need to conduct a univariate meta-analysis now
          }
        }
        if(!is.null(StudiesCovMxNotPosDef)){
          messageMultivar <- "There is at least one study for which the covariance matrix is not positive definite; therefore, a univariate approach is used."
        }
      }
      #
      if(Trans == 0){ # then we cannot transform all the Phi's (uniquely), so use dummy methode
        # Calculate standardized Phi
        for(s in 1:S){
          out <- calc.CovMxStandPhi(N[s], Phi[,,s], Gamma[,,s], SigmaVAR[,,s])
          vecStandPhi[,s] <- out$vecStandPhi_DeltaT
          CovMxPhi[,,s] <- out$CovMx_vecStandPhi_DeltaT
          #
          # If covariance matrix is not positive definite then do univariate approach (and use only variances)
          if(any(eigen(CovMxPhi[,,s])$val < 0)){
            StudiesCovMxNotPosDef <- c(StudiesCovMxNotPosDef, s) # for which studies is this the case
            Multivar <- 0 # we need to conduct a univariate meta-analysis now
          }
          if(!is.null(StudiesCovMxNotPosDef)){
            messageMultivar <- "There is at least one study for which the covariance matrix is not positive definite; therefore, a univariate approach is used."
          }
        }
        # Calculate dummy variables
        if(Multivar == 0){ # if not multivar
          TI <- DeltaT
        } else{ # multivariate
          TI <- matrix(rep(as.matrix(DeltaT), each = q*q), ncol = dim(DeltaT)[2])
        }
        G <- length(unique(TI))
        D <- dummy_cols(TI)
        # First column contain data itself, so take that out:
        if(dim(D)[2] > G){
          D <- D[,-1]
        }
        D. <- as.matrix(D)
        D. <- D.[,order(unique(TI))]
        if(G > 1){
          colnames(D.) <- unique(TI)[order(unique(TI))]
        }
      }
      #
      #
      if(Multivar == 0){ # if not multivar, then we need only the variances and the original moderators and grouping variables
        varPhi <- apply(CovMxPhi, 3, diag) # Gives q times S matrix
        #
        if(Moderators == 1){
          Mod. <- Mod
        }
      }else{ # if multivar, then we need to add an outcome variable and 'copy' the moderators and grouping variables
        sub = NULL
        for(i in 1:q){
          sub = c(sub, paste(i, 1:q, sep=""))
        }
        outcome <- rep(sub, S)
        #
        #
        if(Moderators == 1){
          Mod. <- matrix(rep(as.matrix(Mod), each = q*q), ncol = dim(Mod)[2])
          colnames(Mod.) <- colnames(Mod)
        }
        #
        if(FEorRE == 2){
          RandomPart <- matrix(rep(1:S, each = q*q), ncol = 1)
        }
      }
      #
      #
      #
      # Meta-analysis
      CovMx <- array(0, dim = c(S*q*q, S*q*q))
      for(s in 1:S){
        CovMx[((s-1)*q^2+1):(s*q^2),((s-1)*q^2+1):(s*q^2)] = CovMxPhi[,,s]
      }
      #
      # Multivar and Trans
      if(Multivar == 1 & Trans == 1){
        vecVecStandPhi <- as.vector(vecStandPhi) # vecStandPhi[q:q*q,1:S]
        if(FEorRE == 2){ # RE
          if(Moderators == 1){ # Moderators
            metaan <- rma.mv(yi=vecVecStandPhi, V=CovMx, mods = ~ outcome + outcome:Mod. - 1, random = ~ outcome | RandomPart,
                           struct = "UN", method = "ML")  # With struct="UN", the random effects are allowed to have different variances for each outcome and are allowed to be correlated.
          }else{ # No Moderators
            metaan <- rma.mv(yi=vecVecStandPhi, V=CovMx, mods = ~ outcome - 1, random = ~ outcome | RandomPart,
                             struct = "UN", method = "ML")  # With struct="UN", the random effects are allowed to have different variances for each outcome and are allowed to be correlated.
          }
          tau2_metaan_MV <- metaan$tau2
        }else{# FE
          if(Moderators == 1){ # Moderators
            metaan <- rma.mv(yi=vecVecStandPhi, V=CovMx, mods = ~ outcome - 1 + outcome:Mod., method = "FE")
          }else{ # No Moderators
            metaan <- rma.mv(yi=vecVecStandPhi, V=CovMx, mods = ~ outcome - 1, method = "FE")
          }
        }
        Phi_metaan_MV <- coef(metaan)[1:q^2]
        CovMxPhi_metaan <- metaan$vb[1:q^2,1:q^2]
        # elliptical 95%CI
        # Determine points on 95% LL contour - Has to be done per unique time interval
        #minPhi <- matrix(NA, length(unique(TI)), q^2)
        #maxPhi <- matrix(NA, length(unique(TI)), q^2)
        #tellerCovMx <- 1
        #for(i in 1:length(unique(TI))){
          vecPhi_metaan <- Phi_metaan_MV #[i,]
          CovMx_metaan <- CovMxPhi_metaan #[i,,]
          #tellerCovMx <- tellerCovMx + 4
          eigenCovMx <- eigen(CovMx_metaan)
          lambda <- eigenCovMx$val
          E <- eigenCovMx$vec
          Chi2 <- qchisq(p=0.05, df=(q*q), lower.tail=FALSE)
          LB_vecPhi <- matrix(NA, nrow=q*q, ncol =q*q)
          UB_vecPhi <- matrix(NA, nrow=q*q, ncol =q*q)
          LL <- matrix(NA, nrow=q*q, ncol=2)
          teller <- 0
          for(row in 1:q){
            for(column in 1:q){
              teller <- teller + 1
              LB_vecPhi[teller,] <- vecPhi_metaan - sqrt(Chi2 * lambda[teller]) * E[,teller]
              UB_vecPhi[teller,] <- vecPhi_metaan + sqrt(Chi2 * lambda[teller]) * E[,teller]
              LL[teller,1] <- t(LB_vecPhi[teller,]-vecPhi_metaan) %*% solve(CovMx_metaan) %*% (LB_vecPhi[teller,]-vecPhi_metaan)
              LL[teller,2] <- t(UB_vecPhi[teller,]-vecPhi_metaan) %*% solve(CovMx_metaan) %*% (UB_vecPhi[teller,]-vecPhi_metaan)
              Phi_LB_t <- matrix(LB_vecPhi[teller,], ncol=q, nrow=q)
              Phi_UB_t <- matrix(UB_vecPhi[teller,], ncol=q, nrow=q)
            }
          }
          #
          # Note that UB can be smaller than LB and LB larger than UB!
          #minPhi_Trans[i,] <- apply(rbind(LB_vecPhi, UB_vecPhi), 2, min)
          #maxPhi_Trans[i,] <- apply(rbind(LB_vecPhi, UB_vecPhi), 2, max)
          minPhi <- apply(rbind(LB_vecPhi, UB_vecPhi), 2, min)
          maxPhi <- apply(rbind(LB_vecPhi, UB_vecPhi), 2, max)
          multiCI <- rbind(minPhi, maxPhi)
          rownames(multiCI) <- c("LB", "UB")
          sub = NULL
          for(i in 1:q){
            sub = c(sub, paste0("overallPhi", i, 1:q, sep=""))
          }
          colnames(multiCI) <- sub
        #}
      #
      # Multivar and Dummies
      }else if(Multivar == 1 & Trans == 0){
        vecVecStandPhi <- as.vector(vecStandPhi) # vecStandPhi[q:q*q,1:S]
        if(FEorRE == 2){ # RE
          if(G > 1){
            if(Moderators == 1){ # Moderators
              metaan <- rma.mv(yi=vecVecStandPhi, V=CovMx, mods = ~ -1 + outcome:D. + outcome:Mod.,
                               random = ~ outcome | RandomPart, struct = "UN", method = "ML")  # With struct="UN", the random effects are allowed to have different variances for each outcome and are allowed to be correlated.
            }else{ # No Moderators
              metaan <- rma.mv(yi=vecVecStandPhi, V=CovMx, mods = ~ -1 + outcome:D.,
                               random = ~ outcome | RandomPart, struct = "UN", method = "ML")  # With struct="UN", the random effects are allowed to have different variances for each outcome and are allowed to be correlated.
            }
            tau2_metaan_MV <- metaan$tau2
          } else{ # G = 1, no use for dummies (only one group)
            if(Moderators == 1){ # Moderators
              metaan <- rma.mv(yi=vecVecStandPhi, V=CovMx, mods = ~ outcome - 1 + outcome:Mod.,
                               random = ~ outcome | RandomPart, struct = "UN", method = "ML")
            }else{ # No Moderators
              metaan <- rma.mv(yi=vecVecStandPhi, V=CovMx, mods = ~ outcome - 1,
                               random = ~ outcome | RandomPart, struct = "UN", method = "ML")
            }
          }
        }else{# FE
          if(G > 1){
            if(Moderators == 1){ # Moderators
              metaan <- rma.mv(yi=vecVecStandPhi, V=CovMx, mods = ~ -1 + outcome:D. + outcome:Mod., method = "FE")
            }else{ # No Moderators
              metaan <- rma.mv(yi=vecVecStandPhi, V=CovMx, mods = ~ -1 + outcome:D., method = "FE")
            }
          } else{ # G = 1, no use for dummies (only one group)
            if(Moderators == 1){ # Moderators
              metaan <- rma.mv(yi=vecVecStandPhi, V=CovMx, mods = ~ outcome - 1 + outcome:Mod., method = "FE")
            }else{ # No Moderators
              metaan <- rma.mv(yi=vecVecStandPhi, V=CovMx, mods = ~ outcome - 1, method = "FE")
            }
          }
        }
        Phi_metaan_MV <- matrix(coef(metaan), length(unique(TI)), q^2, byrow = T)
        CovMxPhi_metaan <- metaan$vb
        # elliptical 95%CI
        # Determine points on 95% LL contour - Has to be done per unique time interval
        minPhi <- matrix(NA, length(unique(TI)), q^2)
        maxPhi <- matrix(NA, length(unique(TI)), q^2)
        sub = NULL
        for(i in 1:q){
          sub = c(sub, paste0("overallPhi", i, 1:q, sep=""))
        }
        dimnames(minPhi) <- dimnames(maxPhi) <- list(paste0("DeltaT=", unique(TI)[order(unique(TI))]), sub)
        tellerCovMx <- 1
        for(i in 1:length(unique(TI))){
          vecPhi_metaan <- Phi_metaan_MV[i,]
          CovMx_metaan <- CovMxPhi_metaan[tellerCovMx:(tellerCovMx+3),tellerCovMx:(tellerCovMx+3)]
          tellerCovMx <- tellerCovMx + 4
          eigenCovMx <- eigen(CovMx_metaan)
          lambda <- eigenCovMx$val
          E <- eigenCovMx$vec
          Chi2 <- qchisq(p=0.05, df=(q*q), lower.tail=FALSE)
          LB_vecPhi <- matrix(NA, nrow=q*q, ncol =q*q)
          UB_vecPhi <- matrix(NA, nrow=q*q, ncol =q*q)
          LL <- matrix(NA, nrow=q*q, ncol=2)
          teller = 0
          for(row in 1:q){
            for(column in 1:q){
              teller = teller + 1
              LB_vecPhi[teller,] <- vecPhi_metaan - sqrt(Chi2 * lambda[teller]) * E[,teller]
              UB_vecPhi[teller,] <- vecPhi_metaan + sqrt(Chi2 * lambda[teller]) * E[,teller]
              LL[teller,1] <- t(LB_vecPhi[teller,]-vecPhi_metaan) %*% solve(CovMx_metaan) %*% (LB_vecPhi[teller,]-vecPhi_metaan)
              LL[teller,2] <- t(UB_vecPhi[teller,]-vecPhi_metaan) %*% solve(CovMx_metaan) %*% (UB_vecPhi[teller,]-vecPhi_metaan)
              Phi_LB_t <- matrix(LB_vecPhi[teller,], ncol=q, nrow=q)
              Phi_UB_t <- matrix(UB_vecPhi[teller,], ncol=q, nrow=q)
            }
          }
          # Note that UB can be smaller than LB and LB larger than UB!
          minPhi[i,] <- apply(rbind(LB_vecPhi, UB_vecPhi), 2, min)
          maxPhi[i,] <- apply(rbind(LB_vecPhi, UB_vecPhi), 2, max)
        }
      #
      # Univar and Trans
      }else if(Multivar == 0 & Trans == 1){
        Phi_metaan <- array(data=NA, dim=c(q,q))
        sePhi_metaan <- array(data=NA, dim=c(q,q))
        minPhi <- array(data=NA, dim=c(q,q))
        maxPhi <- array(data=NA, dim=c(q,q))
        tau2_metaan <- array(data=NA, dim=c(q,q))
        QEp_metaan <- array(data=NA, dim=c(q,q))
        QMp_metaan <- array(data=NA, dim=c(q,q))
        for(teller in 1:(q^2)){ # Do meta-analysis (over S studies) for each of the q*q SLEs
          j = (teller+q-1)%/%q
          k = teller-(j-1)*q
          Phi_jk <- vecStandPhi[teller,]
          var <- varPhi[teller,]
          if(FEorRE == 2){ # RE
            if(Moderators == 1){ # Moderators
              metaan <- rma.uni(yi=Phi_jk, vi=var, mods = Mod., method = "ML")
            }else{ # No Moderators
              metaan <- rma.uni(yi=Phi_jk, vi=var, method = "ML")
            }
            tau2_metaan[j,k] <- metaan$tau2
          }else{# FE
            if(Moderators == 1){ # Moderators
              metaan <- rma.uni(yi=Phi_jk, vi=var, mods = Mod., method = "FE")
            }else{ # No Moderators
              metaan <- rma.uni(yi=Phi_jk, vi=var, method = "FE")
            }
          }
          Phi_metaan[j,k] <- coef(metaan)[1]
          sePhi_metaan[j,k] <- metaan$se[1] # This is se for DeltaT^g
          minPhi[j,k] <- metaan$ci.lb[1]
          maxPhi[j,k] <- metaan$ci.ub[1]
          QEp_metaan[j,k] <- metaan$QEp
          if(Moderators == 1){ # Moderators
            QMp_metaan[j,k] <- metaan$QMp
          }
        }
      #
      # Univar and Dummies
      }else if(Multivar == 0 & Trans == 0){
        Phi_metaan_dum <- array(data=NA, dim=c(q,q,G))
        sePhi_metaan_dum <- array(data=NA, dim=c(q,q,G))
        minPhi <- array(data=NA, dim=c(q,q,G))
        maxPhi <- array(data=NA, dim=c(q,q,G))
        dimnames(minPhi)[3] <- dimnames(maxPhi)[3] <- list(paste0("DeltaT=",unique(TI)[order(unique(TI))]))
        tau2_metaan_dum <- array(data=NA, dim=c(q,q,G))
        QEp_metaan <- array(data=NA, dim=c(q,q))
        QMp_metaan <- array(data=NA, dim=c(q,q))
        for(teller in 1:(q^2)){ # Do meta-analysis (over S studies) for each of the q*q SLEs
          j = (teller+q-1)%/%q
          k = teller-(j-1)*q
          Phi_jk <- vecStandPhi[teller,]
          var <- varPhi[teller,]
          if(FEorRE == 2){ # RE
            if(G > 1){
              if(Moderators == 1){ # Moderators
                metaan_g <- rma.uni(yi=Phi_jk, vi=var, mods = ~ -1 + D. + Mod., method = "ML")
              }else{ # No Moderators
                metaan_g <- rma.uni(yi=Phi_jk, vi=var, mods = ~ -1 + D., method = "ML")
              }
            } else{ # G == 1: no dummies needed
              if(Moderators == 1){ # Moderators
                metaan_g <- rma.uni(yi=Phi_jk, vi=var, mods = ~ 1 + Mod., method = "ML")
              }else{ # No Moderators
                metaan_g <- rma.uni(yi=Phi_jk, vi=var, method = "ML")
              }
            }
            tau2_metaan_dum[j,k,] <- metaan_g$tau2
          }else{# FE
            if(G > 1){
              if(Moderators == 1){ # Moderators
                metaan_g <- rma.uni(yi=Phi_jk, vi=var, mods = ~ -1 + D. + Mod., method = "FE")
              }else{ # No Moderators
                metaan_g <- rma.uni(yi=Phi_jk, vi=var, mods = D., intercept=FALSE, method = "FE")
                #metaan_g <- rma.uni(yi=Phi_jk, vi=var, mods = ~ -1 + D., method = "FE")
              }
            } else{ # G == 1: no dummies needed
              if(Moderators == 1){ # Moderators
                metaan_g <- rma.uni(yi=Phi_jk, vi=var, mods = ~ 1 + Mod., method = "FE")
              }else{ # No Moderators
                #metaan_g <- rma.uni(yi=Phi_jk, vi=var, method = "FE")
                metaan_g <- rma.uni(yi=Phi_jk, vi=var, intercept=TRUE, method = "FE")
              }
            }
          }
          if(Moderators == 1){ # Moderators
            Phi_metaan_dum[j,k,] <- coef(metaan_g)[1:G]
            sePhi_metaan_dum[j,k,] <- metaan_g$se[1:G] # This is se for DeltaT^g
            # These are the values when Mod. = 0.
          }else{ # No Moderators
            Phi_metaan_dum[j,k,] <- coef(metaan_g)
            sePhi_metaan_dum[j,k,] <- metaan_g$se # This is se for DeltaT^g
          }
          minPhi[j,k,] <- metaan_g$ci.lb[1:G]
          maxPhi[j,k,] <- metaan_g$ci.ub[1:G]
          QEp_metaan[j,k] <- metaan_g$QEp
          if(Moderators == 1){ # Moderators
            QMp_metaan[j,k] <- metaan_g$QMp
          }
        }
      }
      #
      #
      #
      # Output
      dT <- unique(DeltaT)[order(unique(DeltaT))]
      #names(dT) <- "Unique DeltaT values used in the studies"
      dT_star <- DeltaTStar
      #names(dT_start) <- "DeltaT_star value"
      if(!is.null(StudiesComplexEV)) names(StudiesComplexEV) <- "Study indices for which (some of) the eigenvalues are complex"
      if(!is.null(StudiesNegEV)) names(StudiesNegEV) <- "Study indices for which (some of) the eigenvalues are negative"
      if(!is.null(StudiesCovMxNotPosDef)) names(StudiesCovMxNotPosDef) <- "Study indices for which the covariance matrix is not positive definite"
      ratioDeltaT <- matrix(ratioDeltaT)
      colnames(ratioDeltaT) <- "Ratio DeltaT*/DeltaT"
      rownames(ratioDeltaT) <- paste("Study ", rep(1:S))
      # Multivar and Trans
      if(Multivar == 1 & Trans == 1){
        if(FEorRE == 1){ # FE
          final <- list(DeltaTStar = dT_star,
                        Overall_standPhi_DeltaTStar = Phi_metaan_MV,
                        elliptical_CI = multiCI,
                        alpha_CI = alpha,
                        CovMx_OverallPhi_DeltaTStar = CovMxPhi_metaan,
                        messageTrans = messageTrans, messageMultivar = messageMultivar,
                        StudiesComplexEV = StudiesComplexEV, StudiesNegEV = StudiesNegEV, StudiesCovMxNotPosDef = StudiesCovMxNotPosDef,
                        ratioDeltaT = ratioDeltaT,
                        summaryMetaAnalysis = summary(metaan))
        } else{ # RE
          final <- list(DeltaTStar = dT_star,
                        Overall_standPhi_DeltaTStar = Phi_metaan_MV,
                        elliptical_CI = multiCI,
                        alpha_CI = alpha,
                        CovMx_OverallPhi_DeltaTStar = CovMxPhi_metaan,
                        tau2 = tau2_metaan_MV,
                        messageTrans = messageTrans, messageMultivar = messageMultivar,
                        StudiesComplexEV = StudiesComplexEV, StudiesNegEV = StudiesNegEV, StudiesCovMxNotPosDef = StudiesCovMxNotPosDef,
                        ratioDeltaT = ratioDeltaT,
                        summaryMetaAnalysis = summary(metaan))
        }
      # Multivar and Dummies
      }else if(Multivar == 1 & Trans == 0){
        if(FEorRE == 1){ # FE
          final <- list(UniqueTimeIntervals = dT,
                        Overall_standPhi_PerUniqueTimeInterval = Phi_metaan_MV,
                        LB_elliptical_CI = minPhi,
                        UB_elliptical_CI = maxPhi,
                        alpha_CI = alpha,
                        CovMx_OverallPhi_PerUniqueTimeInterval = CovMxPhi_metaan,
                        messageTrans = messageTrans, messageMultivar = messageMultivar,
                        StudiesComplexEV = StudiesComplexEV, StudiesNegEV = StudiesNegEV, StudiesCovMxNotPosDef = StudiesCovMxNotPosDef,
                        ratioDeltaT = ratioDeltaT,
                        summaryMetaAnalysis = summary(metaan))
        } else{ # RE
          final <- list(UniqueTimeIntervals = dT,
                        Overall_standPhi_PerUniqueTimeInterval = Phi_metaan_MV_dum,
                        LB_elliptical_CI = minPhi,
                        UB_elliptical_CI = maxPhi,
                        alpha_CI = alpha,
                        CovMx_OverallPhi_PerUniqueTimeInterval = CovMxPhi_metaan,
                        tau2 = tau2_metaan_MV,
                        messageTrans = messageTrans, messageMultivar = messageMultivar,
                        StudiesComplexEV = StudiesComplexEV, StudiesNegEV = StudiesNegEV, StudiesCovMxNotPosDef = StudiesCovMxNotPosDef,
                        ratioDeltaT = ratioDeltaT,
                        summaryMetaAnalysis = summary(metaan))
        }
      # Univar and Trans
      }else if(Multivar == 0 & Trans == 1){
        if(FEorRE == 1){ # FE
          final <- list(DeltaTStar = dT_star,
                        Overall_standPhi_DeltaTStar = Phi_metaan,
                        LB_CI = minPhi,
                        UB_CI = maxPhi,
                        alpha_CI = alpha,
                        se_Overall_standPhi_DeltaTStar = sePhi_metaan,
                        pvalue_HeterogeneityTest = QEp_metaan,
                        pvalue_OmnibusCoefTest_IfModerator = QMp_metaan,
                        messageTrans = messageTrans, messageMultivar = messageMultivar,
                        StudiesComplexEV = StudiesComplexEV, StudiesNegEV = StudiesNegEV, StudiesCovMxNotPosDef = StudiesCovMxNotPosDef,
                        ratioDeltaT = ratioDeltaT)
        } else{ # RE
          final <- list(DeltaTStar = dT_star,
                        Overall_standPhi_DeltaTStar = Phi_metaan,
                        LB_CI = minPhi,
                        UB_CI = maxPhi,
                        alpha_CI = alpha,
                        se_Overall_standPhi_DeltaTStar = sePhi_metaan,
                        tau2 = tau2_metaan,
                        pvalue_HeterogeneityTest = QEp_metaan,
                        pvalue_OmnibusCoefTest_IfModerator = QMp_metaan,
                        messageTrans = messageTrans, messageMultivar = messageMultivar,
                        StudiesComplexEV = StudiesComplexEV, StudiesNegEV = StudiesNegEV, StudiesCovMxNotPosDef = StudiesCovMxNotPosDef,
                        ratioDeltaT = ratioDeltaT)
        }
      # Univar and Dummies
      }else if(Multivar == 0 & Trans == 0){
        if(FEorRE == 1){ # FE
          final <- list(UniqueTimeIntervals = dT,
                        Overall_standPhi_PerUniqueTimeInterval = Phi_metaan_dum,
                        LB_CI = minPhi,
                        UB_CI = maxPhi,
                        alpha_CI = alpha,
                        se_Overall_standPhi_PerUniqueTimeInterval = sePhi_metaan_dum,
                        pvalue_HeterogeneityTest = QEp_metaan,
                        pvalue_OmnibusCoefTest = QMp_metaan,
                        messageTrans = messageTrans, messageMultivar = messageMultivar,
                        StudiesComplexEV = StudiesComplexEV, StudiesNegEV = StudiesNegEV, StudiesCovMxNotPosDef = StudiesCovMxNotPosDef,
                        ratioDeltaT = ratioDeltaT)
        } else{ # RE
          final <- list(UniqueTimeIntervals = dT,
                        Overall_standPhi_PerUniqueTimeInterval = Phi_metaan_dum,
                        LB_CI = minPhi,
                        UB_CI = maxPhi,
                        alpha_CI = alpha,
                        se_Overall_standPhi_PerUniqueTimeInterval = sePhi_metaan_dum,
                        tau2 = tau2_metaan_dum,
                        pvalue_HeterogeneityTest = QEp_metaan,
                        pvalue_OmnibusCoefTest = QMp_metaan,
                        messageTrans = messageTrans, messageMultivar = messageMultivar,
                        StudiesComplexEV = StudiesComplexEV, StudiesNegEV = StudiesNegEV, StudiesCovMxNotPosDef = StudiesCovMxNotPosDef,
                        ratioDeltaT = ratioDeltaT)
        }
      }
      return(final)
    }
