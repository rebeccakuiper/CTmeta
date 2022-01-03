#' Continuous-time meta-analysis on standardized lagged effects
#'
#' Continuous-time meta-analysis (CTmeta) on standardized lagged effects (Phi) taking into account the various time-intervals used in the primary studies. There is also an interactive web application on my website to perform CTmeta: \url{https://www.uu.nl/staff/RMKuiper/Websites\%20\%2F\%20Shiny\%20apps}.
#'
#' @param N Number of persons (panel data) or number of measurement occasions - 1 (time series data) used in the S primary studies. Matrix of size S times 1.
#' @param DeltaT The time intervals used in the S primary studies. Matrix of size S times 1. Note that all the time intervals should be on the same scale (e.g., two time-intervals of 60 minutes and 2 hours, should be either 60 and 120 or 1 and 2).
#' @param DeltaTStar The time interval (scalar) to which the standardized lagged effects matrix should be transformed to.
#' @param Phi Stacked matrix of size S*q times q or array with dimensions q times q times S of (un)standardized lagged effects matrices for all S primary studies in the meta-analysis; with q the number of variables (leading to an q times q lagged effects matrix in a single primary study). Note: In case primary studies report (lagged) correlation matrices, one can use the function 'TransPhi_Corr' to transform those to the corresponding standardized lagged effects matrices (see ?TransPhi_Corr and examples below).
#' @param SigmaVAR Stacked matrix of size S*q times q or array with dimensions q times q times S of residual covariance matrices of the first-order discrete-time vector autoregressive (DT-VAR(1)) model.
#' @param Gamma Optional (either SigmaVAR or Gamma). Stacked matrix of size S*q times q or array with dimensions q times q times S of stationary covariance matrices, that is, the contemporaneous covariance matrices of the data sets.
#' Note that if Phi and Gamma are known, SigmaVAR can be calculated. Hence, only SigmaVAR or Gamma is needed as input (if only Gamma, then use 'Gamma = Gamma' or set SigmaVAR to NULL, see examples below).
#' @param Moderators Optional. Indicator (TRUE/FALSE or 1/0) whether there are moderators to be included (TRUE or 1) or not (FALSE or 0). By default, Moderators = 0.
#' @param Mod Optional. An S x m matrix of m moderators to be included in the analysis when Moderators == 1. By default, Mod = NULL.
#' @param FEorRE Optional. Indicator (1/2) whether continuous-time meta-analysis should use a fixed-effects model (1) or random-effects model (2). By default, FEorRE = 1.
#' @param BetweenLevel Optional. Needed in case of a 2-level multilevel meta-analysis. BetweenLevel should be an S-vector or S x 1 matrix (namely one value for each study). It will only be used, if FEorRE == 2 (i.e., in a random-effects model). Then, one can add a between-level to the random part (e.g., sample number if multiple studies use the sample such that there is dependency between those studies). Note that the within level is Study number. By default, BetweenLevel = NULL.
#' @param Label Optional. If one creates, for example, a funnel or forest plot the labeling used in the rma.mv function is used. Label should be an q*q*S-vector (namely one value for each of the elements in a study-specific Phi (of size q x q) and for each study). It will only be used, in case the multivariate approach can be used (in case of the univariate approach, it will always use the labeling Study 1 to Study S). By default, Label = NULL; then it will use the labeling Study 1 Phi11, Study Phi 12, ..., Study 1 Phi qq, ... Study S Phi11, ..., Study S Phiqq.
#' @param alpha Optional. The alpha level in determining the (1-alpha)*100\% confidence interval (CI). By default, alpha = 0.05; resulting in a 95\% CI.
#' @param PrintPlot Optional. Indicator (TRUE/FALSE or 1/0) for rendering a Phi-plot (TRUE or 1) or not (FALSE or 0). By default, PrintPlot = FALSE.
#'
#' @return The output comprises, among others, the overall vectorized transformed standardized lagged effects, their covariance matrix, and the corresponding elliptical/multivariate 95\% CI.
#' @importFrom fastDummies dummy_cols
#' @importFrom metafor rma.mv
#' @importFrom metafor rma.uni
#' @export print.CTmeta
#' @export summary.CTmeta
#' @export coef.CTmeta
#' @export vcov.CTmeta
#' @export
#' @examples
#'
#' # library(CTmeta)
#'
#' ##################################################################################################
#' # Input needed in examples below with q=2 variables and S=3 primary studies
#' #
#' N <- matrix(c(643, 651, 473))
#' DeltaT <- matrix(c(2, 3, 1))
#' DeltaTStar <- 1
#' #
#' # I will use the example matrices stored in the package:
#' Phi <- myPhi
#' SigmaVAR <- mySigmaVAR
#' Gamma <- myGamma # Note: CTmeta does not need both SigmaVAR and Gamma, as denomstrated below.
#' # These are all three stacked matrices of size S*q times q.
#' # The CTmeta function will standardize these matrices (to make comparison of effects meaningful).
#' #
#' Moderators = 0 # By default set to 0. Hence, not per se needed, as demonstrated below.
#' ##################################################################################################
#'
#'
#' # Below, you can find example code. Note that, here, only 3 primary studies are used.
#' # In practice, one would normally have (many) more, but the code stays (more or less) the same.
#'
#'
#' ### Examples without comments ###
#'
#'
#' ## Example without moderators ##
#'
#' # Fixed effects model #
#'
#' # Run CTmeta
#' CTma <- CTmeta(N, DeltaT, DeltaTStar, Phi, SigmaVAR)
#' CTma
#' summary(CTma)
#'
#'
#' # Random effects (RE) model #
#'
#' CTmeta(N, DeltaT, DeltaTStar, Phi, SigmaVAR, FEorRE = 2)
#'
#' BetweenLevel <- c(1, 1, 2) # Assuming the first two studies used the same sample/dataset
#' CTmeta(N, DeltaT, DeltaTStar, Phi, SigmaVAR, FEorRE = 2, BetweenLevel = BetweenLevel) # Two-level RE meta-analysis example
#'
#'
#' ## Example with moderators ##
#'
#' Mod <- matrix(c(64,65,47)) # 1 moderator
#' CTma <- CTmeta(N, DeltaT, DeltaTStar, Phi, SigmaVAR, Moderators = 1, Mod = Mod) # fixed effects model
#' summary(CTma)
#'
#' q <- dim(Phi)[2]; Mod <- matrix(cbind(c(64,65,47), c(78,89,34)), ncol = q); colnames(Mod) <- c("Mod1", "Mod2") # two moderators, in each column 1
#' CTma <- CTmeta(N, DeltaT, DeltaTStar, Phi, SigmaVAR, Moderators = 1, Mod = Mod, FEorRE = 2); CTma$tau2 # random effects model
#' summary(CTma)
#'
#' Mod <- matrix(c(64,65,47)) # 1 moderator
#' BetweenLevel <- c(1, 1, 2) # Assuming the first two studies used the same sample/dataset.
#' CTmeta(N, DeltaT, DeltaTStar, Phi, SigmaVAR, Moderators = 1, Mod = Mod, FEorRE = 2, BetweenLevel = BetweenLevel) # Two-level RE meta-analysis example
#'
#'
#' ## funnel and forest plots ##
#' CTma <- CTmeta(N, DeltaT, DeltaTStar, Phi, SigmaVAR, FEorRE = 2)
#' funnel(CTma$summaryMetaAnalysis, label = 'out')
#' forest(CTma$summaryMetaAnalysis)
#' # One can do the same for the two-level analysis.
#' #
#' # Note, in case a univariate approach had to be taken, leading to multiple analyses, then one should create a plot per analysis:
#' # lapply(CTma$summaryMetaAnalysis_jk, funnel, label = 'out')
#' #
#' # In case you want to create a plot per element of the study-specific Phi's:
#' #elt <- rep(F, (q*q))
#' #elt[1] <- T # First element out of q*q true, so referring to element Phi11.
#' #yi_Phi11 <- CTma$summaryMetaAnalysis$yi[elt]
#' #vi_Phi11 <- CTma$summaryMetaAnalysis$vi[elt]
#' #funnel(yi_Phi11, vi_Phi11, label = 'out')
#' #forest(yi_Phi11, vi_Phi11)
#' elt <- rep(F, (q*q))
#' elt[2] <- T # Second element out of q*q true, so referring to element Phi12.
#' yi_Phi12 <- CTma$summaryMetaAnalysis$yi[elt]
#' vi_Phi12 <- CTma$summaryMetaAnalysis$vi[elt]
#' funnel(yi_Phi12, vi_Phi12, label = 'out')
#' forest(yi_Phi12, vi_Phi12)
#'
#'
#' ## Make customized Phi-plot of resulting overall Phi ##
#'
#' CTma <- CTmeta(N, DeltaT, DeltaTStar, Phi, SigmaVAR, PrintPlot = TRUE)
#' CTma$PhiPlot
#'
#' out_CTmeta <- CTmeta(N, DeltaT, DeltaTStar, Phi, SigmaVAR)
#' overallPhi <- out_CTmeta$Overall_standPhi
#' Title <- as.list(expression(paste0(Phi(Delta[t]), " plot:"),
#'    "How do the overall lagged parameters vary as a function of the time-interval"))
#' PhiPlot(DeltaTStar, overallPhi, Min = 0, Max = 40, Step = 0.5, Title = Title)
#' # or
#' ggPhiPlot <- ggPhiPlot(DeltaTStar, overallPhi, Min = 0, Max = 40, Step = 0.5, Title = Title)
#' ggPhiPlot$PhiPlot
#'
#'
#' ## Evaluate dominance of (absolute values of) cross-lagged effects ##
#'
#' if (!require("restriktor")) install.packages("restriktor") # Use restriktor package for function goric().
#'                                                            # Authors of goric(): Vanbrabant and Kuiper.
#' library(restriktor)
#'
#' out_CTmeta <- CTmeta(N, DeltaT, DeltaTStar, Phi, SigmaVAR)
#' H1 <- "abs(overallPhi12) < abs(overallPhi21)"
#' goric(out_CTmeta, H1, type = "gorica", comparison = "complement")
#'
#' out_CTmeta <- CTmeta(N, DeltaT, DeltaTStar, Phi, SigmaVAR)
#' H1 <- "abs(overallPhi12) < abs(overallPhi21); abs(overallPhi11) < abs(overallPhi22)"
#' goric(out_CTmeta, H1, type = "gorica", comparison = "complement")
#'
#' out_CTmeta <- CTmeta(N, DeltaT, DeltaTStar, Phi, SigmaVAR)
#' est <- coef(out_CTmeta)  # or: est  <- out_CTmeta$Overall_vecStandPhi_DeltaTStar
#' VCOV <- vcov(out_CTmeta) # or: VCOV <- out_CTmeta$CovMx_OverallPhi_DeltaTStar
#' goric(est, VCOV = VCOV, H1, type = "gorica", comparison = "complement")
#'
#'
#' ## What if primary studies report a (lagged) correlation matrix ##
#'
#' q <- 2
#' corr_YXYX <- matrix(c(1.00, 0.40, 0.63, 0.34,
#'                       0.40, 1.00, 0.31, 0.63,
#'                       0.63, 0.31, 1.00, 0.41,
#'                       0.34, 0.63, 0.41, 1.00), byrow = T, ncol = 2*q)
#' N <- matrix(c(643, 651, 473))
#' DeltaT <- matrix(c(2, 3, 1))
#' DeltaTStar <- 1
#' #
#' # Create input:
#' out_1 <- TransPhi_Corr(DeltaTStar = DeltaT[1], DeltaT = 1, N = N[1], corr_YXYX)
#' Phi_1 <- out_1$standPhi_DeltaTStar
#' SigmaVAR_1 <- out_1$standSigmaVAR_DeltaTStar
#' out_2 <- TransPhi_Corr(DeltaTStar = DeltaT[2], DeltaT = 1, N = N[2], corr_YXYX)
#' Phi_2 <- out_2$standPhi_DeltaTStar
#' SigmaVAR_2 <- out_2$standSigmaVAR_DeltaTStar
#' out_3 <- TransPhi_Corr(DeltaTStar = DeltaT[3], DeltaT = 1, N = N[3], corr_YXYX)
#' Phi_3 <- out_3$standPhi_DeltaTStar
#' SigmaVAR_3 <- out_3$standSigmaVAR_DeltaTStar
#' #
#' Phi <- rbind(Phi_1, Phi_2, Phi_3) # This, returns a stacked matrix of size S q times q.
#' SigmaVAR <- rbind(SigmaVAR_1, SigmaVAR_2, SigmaVAR_3)
#' out_CTmeta <- CTmeta(N, DeltaT, DeltaTStar, Phi, SigmaVAR)
#' out_CTmeta$Overall_standPhi
#'
#'
#' ##############################
#'
#'
#' ### Examples with comments ###
#'
#'
#' ## Example without moderators ##
#'
#' # Fixed effects model #
#'
#' # Run CTmeta with, for instance,
#' CTmeta(N, DeltaT, DeltaTStar, Phi, SigmaVAR)
#'
#' # There are multiple options; use one of the following:
#' CTmeta(N, DeltaT, DeltaTStar, Phi, SigmaVAR, Gamma, Moderators, Mod, 1) # The 1, here, says FEorRE = 1
#' CTmeta(N, DeltaT, DeltaTStar, Phi, SigmaVAR, Gamma, Moderators, Mod)    # Notably, if Moderators = 0, 'Mod' is not inspected.
#' CTmeta(N, DeltaT, DeltaTStar, Phi, SigmaVAR, Gamma, Moderators)
#' CTmeta(N, DeltaT, DeltaTStar, Phi, SigmaVAR, Gamma)
#' CTmeta(N, DeltaT, DeltaTStar, Phi, SigmaVAR)        # Says, implicitly, Gamma = NULL
#' CTmeta(N, DeltaT, DeltaTStar, Phi, Gamma = Gamma)   # Says, implicitly, SigmaVAR = NULL
#' CTmeta(N, DeltaT, DeltaTStar, Phi, NULL, Gamma)     # Says SigmaVAR = NULL
#'
#' # Note: Do NOT use
#' #CTmeta(N, DeltaT, DeltaTStar, Phi, Gamma)
#' # Then, CTmeta incorrectly uses SigmaVAR = Gamma.
#'
#' # Different types of output options are possible:
#' CTma <- CTmeta(N, DeltaT, DeltaTStar, Phi, SigmaVAR)
#' CTma # gives same as print(CTma)
#' summary(CTma)
#' print(CTma, digits = 4)
#' summary(CTma, digits = 4)
#' # In Rstudio, use 'CTma$' to see what output there is. For example:
#' CTma$summaryMetaAnalysis
#'
#'
#' # Random effects model #
#'
#' # Add "FEorRE = 2"; e.g.,
#' CTmeta(N, DeltaT, DeltaTStar, Phi, SigmaVAR, FEorRE = 2)
#'
#' # Two-level RE meta-analysis example
#' BetweenLevel <- c(1, 1, 2) # Assuming the first two studies used the same sample/dataset
#' CTmeta(N, DeltaT, DeltaTStar, Phi, SigmaVAR, FEorRE = 2, BetweenLevel = BetweenLevel) # Two-level RE meta-analysis example
#' # Note, one can also use this in case there are moderators (as in the example below).
#'
#'
#' ## Example with moderators ##
#'
#' Mod <- matrix(c(64,65,47)) # 1 moderator
#' #q <- dim(Phi)[2]; Mod <- matrix(cbind(c(64,65,47), c(78,89,34)), ncol = q); colnames(Mod) <- c("Mod1", "Mod2") # two moderators, in each column 1
#' CTma <- CTmeta(N, DeltaT, DeltaTStar, Phi, SigmaVAR, Moderators = 1, Mod = Mod) # fixed effects model
#' #CTma <- CTmeta(N, DeltaT, DeltaTStar, Phi, SigmaVAR, Moderators = 1, Mod = Mod, FEorRE = 2); CTma$tau2 # random effects model
#' summary(CTma)
#'
#' # Two-level RE meta-analysis example
#' Mod <- matrix(c(64,65,47)) # 1 moderator
#' BetweenLevel <- c(1, 1, 2) # Assuming the first two studies used the same sample/dataset
#' CTmeta(N, DeltaT, DeltaTStar, Phi, SigmaVAR, Moderators = 1, Mod = Mod, FEorRE = 2, BetweenLevel = BetweenLevel) # Two-level RE meta-analysis example
#'
#'
#' ## funnel and forest plots ##
#' CTma <- CTmeta(N, DeltaT, DeltaTStar, Phi, SigmaVAR, FEorRE = 2)
#' funnel(CTma$summaryMetaAnalysis, label = 'out')
#' forest(CTma$summaryMetaAnalysis)
#' # One can do the same for the two-level analysis.
#' #
#' # Note, in case a univariate approach had to be taken, leading to multiple analyses, then one should create a plot per analysis:
#' # lapply(CTma$summaryMetaAnalysis_jk, funnel, label = 'out')
#' #
#' # Notes on funnel and forest:
#' # - These plots are now based on the q*q elements in the study-specific Phi's and the S studies.
#' # See below, how you can create these plots per element of Phi. This should then be done separately for all q*q elements.
#' # - In case label names are too long or not insightful, one can change them by using the Label argument.
#' # In case of the funnel plot, one can also decide to not use the "label = 'out'" part.
#'
#' #
#' # Notes on forest:
#' # - For random-effects models of class "rma.mv" (see rma.mv) with multiple  values, the addpred argument can be used to specify for which level of the inner factor the prediction interval should be provided (since the intervals differ depending on the value).
#' # - One can also look into the functionality of addpoly().
#' #
#' # In case you want to create a plot per element of the study-specific Phi's:
#' #elt <- rep(F, (q*q))
#' #elt[1] <- T # First element out of q*q true, so referring to element Phi11.
#' #yi_Phi11 <- CTma$summaryMetaAnalysis$yi[elt]
#' #vi_Phi11 <- CTma$summaryMetaAnalysis$vi[elt]
#' #funnel(yi_Phi11, vi_Phi11, label = 'out')
#' #forest(yi_Phi11, vi_Phi11)
#' elt <- rep(F, (q*q))
#' elt[2] <- T # Second element out of q*q true, so referring to element Phi12.
#' yi_Phi12 <- CTma$summaryMetaAnalysis$yi[elt]
#' vi_Phi12 <- CTma$summaryMetaAnalysis$vi[elt]
#' funnel(yi_Phi12, vi_Phi12, label = 'out')
#' forest(yi_Phi12, vi_Phi12)
#'
#'
#' ## Make customized Phi-plot of resulting overall Phi ##
#'
#' # Option 1: Using the plot option in the function:
#' CTmeta(N, DeltaT, DeltaTStar, Phi, SigmaVAR, PrintPlot = TRUE)
#' # The plot can be stored and retrieved as an object as follows:
#' # CTma <- CTmeta(N, DeltaT, DeltaTStar, Phi, SigmaVAR, PrintPlot = TRUE)
#' # CTma$PhiPlot
#'
#'
#' # Option 2: A customized Phi-plot can be made using the function 'PhiPlot' (see below) or by using the interactive web app from my website (\url{https://www.uu.nl/staff/RMKuiper/Websites\%20\%2F\%20Shiny\%20apps}).
#' # Alternatively, one can use the function 'ggPhiPlot' instead of 'PhiPlot'.
#'
#' # Extract the (q times q) overall Phi matrix
#' out_CTmeta <- CTmeta(N, DeltaT, DeltaTStar, Phi, SigmaVAR)
#' # resulting overall Phi:
#' overallPhi <- out_CTmeta$Overall_standPhi
#'
#' # Make Phi-plot:
#' Title <- as.list(expression(paste0(Phi(Delta[t]), " plot:"),
#'    "How do the overall lagged parameters vary as a function of the time-interval"))
#' PhiPlot(DeltaTStar, overallPhi, Min = 0, Max = 40, Step = 0.5, Title = Title)
#' # or
#' ggPhiPlot(DeltaTStar, overallPhi, Min = 0, Max = 40, Step = 0.5, Title = Title)
#' # The plot can be stored and retrieved as an object as follows:
#' # ggPhiPlot <- ggPhiPlot(DeltaTStar, overallPhi, Min = 0, Max = 40, Step = 0.5, Title = Title)
#' # ggPhiPlot$PhiPlot
#'
#'
#' ## Evaluate dominance of (absolute values of) cross-lagged effects ##
#'
#' if (!require("restriktor")) install.packages("restriktor") # Use restriktor package for function goric().
#'                                                            # Authors of goric(): Vanbrabant and Kuiper.
#' #
#' # If version from github needed:
#' #if (!require("devtools")) install.packages("devtools")
#' #library(devtools)
#' #install_github("LeonardV/restriktor")
#' #
#' library(restriktor)
#'
#' # Option 1
#' # Use CTmeta object
#' out_CTmeta <- CTmeta(N, DeltaT, DeltaTStar, Phi, SigmaVAR)
#' #
#' # Example 1: dominance of (absolute values of) cross-lagged effects
#' # Specify hypothesis
#' H1 <- "abs(overallPhi12) < abs(overallPhi21)"
#' #H2 <- "abs(overallPhi12) > abs(overallPhi21)" # = complement of H1 and does not need to be specified, see below.
#' # Evaluate dominance of cross-lagged effects via AIC-type criterion called the GORICA (Altinisik, Nederhof, Hoijtink, Oldehinkel, Kuiper, accepted 2021).
#' #goric(out_CTmeta, H1, H2, type = "gorica", comparison = "none")
#' # or, since H2 is complement of H1, equivalently:
#' goric(out_CTmeta, H1, type = "gorica", comparison = "complement")
#' #
#' # Example 2: dominance of (absolute values of) cross-lagged effects and autoregressive effects
#' H1 <- "abs(overallPhi12) < abs(overallPhi21); abs(overallPhi11) < abs(overallPhi22)"
#' # Note that, now, specification of complement 'by hand' harder.
#' goric(out_CTmeta, H1, type = "gorica", comparison = "complement")
#' #
#' #
#' # Option 2
#' # Extract the vectorized overall standardized overallPhi matrix and its covariance matrix
#' # using the functions coef() and vcov()
#' out_CTmeta <- CTmeta(N, DeltaT, DeltaTStar, Phi, SigmaVAR)
#' est <- coef(out_CTmeta)  # or: est  <- out_CTmeta$Overall_vecStandPhi_DeltaTStar
#' VCOV <- vcov(out_CTmeta) # or: VCOV <- out_CTmeta$CovMx_OverallPhi_DeltaTStar
#' goric(est, VCOV = VCOV, H1, type = "gorica", comparison = "complement")
#'
#'
#' ## What if primary studies report a (lagged) correlation matrix ##
#'
#' # Suppose, for ease, that all S=3 primary studies reported the following lagged correlation matrix:
#' q <- 2
#' corr_YXYX <- matrix(c(1.00, 0.40, 0.63, 0.34,
#'                       0.40, 1.00, 0.31, 0.63,
#'                       0.63, 0.31, 1.00, 0.41,
#'                       0.34, 0.63, 0.41, 1.00), byrow = T, ncol = 2*q)
#'
#' # In the example below, the same N and DeltaT(Star) values are used:
#' N <- matrix(c(643, 651, 473))
#' DeltaT <- matrix(c(2, 3, 1))
#' DeltaTStar <- 1
#'
#' # Use the function 'TransPhi_Corr' to calculate the corresponding standardized lagged effects matrix per primary study.
#' # Note that one can already make the time-intervals equal via the arguments DeltaTStar and DeltaT, but CTmeta can as well.
#' # In this example, I deliberately make the time-intervals unequal, such that the example is in line with the input (i.e., DeltaT <- matrix(c(2, 3, 1))) and such the resulting overall Phi should equal the Phi that underlies this lagged correlation matrix (which I check at the end).
#' out_1 <- TransPhi_Corr(DeltaTStar = DeltaT[1], DeltaT = 1, N = N[1], corr_YXYX)
#' Phi_1 <- out_1$standPhi_DeltaTStar
#' SigmaVAR_1 <- out_1$standSigmaVAR_DeltaTStar
#' out_2 <- TransPhi_Corr(DeltaTStar = DeltaT[2], DeltaT = 1, N = N[2], corr_YXYX)
#' Phi_2 <- out_2$standPhi_DeltaTStar
#' SigmaVAR_2 <- out_2$standSigmaVAR_DeltaTStar
#' out_3 <- TransPhi_Corr(DeltaTStar = DeltaT[3], DeltaT = 1, N = N[3], corr_YXYX)
#' Phi_3 <- out_3$standPhi_DeltaTStar
#' SigmaVAR_3 <- out_3$standSigmaVAR_DeltaTStar
#'
#' # Make Phi
#' Phi <- rbind(Phi_1, Phi_2, Phi_3) # This, returns a stacked matrix of size S q times q.
#' SigmaVAR <- rbind(SigmaVAR_1, SigmaVAR_2, SigmaVAR_3)
#' # For more details, see ?TransPhi_Corr
#'
#' # The example CTmeta() code above can be run using this Phi; e.g.,
#' CTmeta(N, DeltaT, DeltaTStar, Phi, SigmaVAR)
#'
#' # The overall q-times-q (here, 2x2) lagged effects matrix Phi
#' out_CTmeta <- CTmeta(N, DeltaT, DeltaTStar, Phi, SigmaVAR)
#' out_CTmeta$Overall_standPhi
#' #
#' # As a check, to see indeed that CTmeta works properly (where the resulting Phi is independent of the choise of N).
#' TransPhi_Corr(DeltaTStar = 1, DeltaT = 1, N = 100, corr_YXYX)$standPhi_DeltaTStar
#' # Note that is normally not a check you would do.
#'


CTmeta <- function(N, DeltaT, DeltaTStar, Phi, SigmaVAR = NULL, Gamma = NULL, Moderators = 0, Mod = NULL, FEorRE = 1, BetweenLevel = NULL, Label = NULL, alpha=0.05, PrintPlot = FALSE) {
  #Gamma = NULL; Moderators = 0; Mod = NULL; FEorRE = 1; BetweenLevel = NULL; Label = NULL; alpha=0.05; PrintPlot = FALSE

#  #######################################################################################################################
#  if (!require("fastDummies")) install.packages("fastDummies")
#  library(fastDummies)
#  if (!require("metafor")) install.packages("metafor")
#  library(metafor)
#  #######################################################################################################################

  S <- length(N) #dim(N)[1]
  
  # Return an error if there are any NAs in N
  if (anyNA(N)) {stop("There are NA values in N.")}
  
  # Return an error if there are any NAs in DeltaT
  if (anyNA(DeltaT)) {stop("There are NA values in DeltaT.")}
  
  # Return an error if there are any NAs in SigmaVAR
  if (anyNA(SigmaVAR)) {stop("There are NA values in SigmaVAR.")}
  
  # Check
  if(S != length(DeltaT)){
    ErrorMessage <- (paste0("The length of the arguments N and DeltaT are not the same. They should both equate to S, the number of primary studies included in the meta-analysis. \n Here, the length of N is ", S, " and the length of DeltaT is ", length(DeltaT), "."))
    stop(ErrorMessage)
  }
  #
  # Make sure N and DeltaT are matrices
  if(length(dim(N)) < 2){
    N <- matrix(N, nrow = S, ncol = 1)
  }
  if(length(dim(DeltaT)) < 2){
    DeltaT <- matrix(DeltaT, nrow = S, ncol = 1)
  }
  #
  if(length(dim(N)) == 2){
    if(dim(N)[1] == 1 & dim(N)[2] != 1){
      N <- matrix(N, nrow = S, ncol = 1)
    }
    if(dim(N)[1] != 1 & dim(N)[2] != 1){
      ErrorMessage <- (paste0("The argument N should be a S x 1 matrix or an S-vector. Currently, it is of size ", dim(N)))
      return(ErrorMessage)
      stop(ErrorMessage)
    }
  }
  if(length(dim(DeltaT)) == 2){
    if(dim(DeltaT)[1] == 1 & dim(DeltaT)[2] != 1){
      DeltaT <- matrix(DeltaT, nrow = S, ncol = 1)
    }
    if(dim(DeltaT)[1] != 1 & dim(DeltaT)[2] != 1){
      ErrorMessage <- (paste0("The argument DeltaT should be a S x 1 matrix or an S-vector. Currently, it is of size ", dim(DeltaT)))
      return(ErrorMessage)
      stop(ErrorMessage)
    }
  }
  #
  if(length(dim(N)) > 3){
    ErrorMessage <- (paste0("The argument N should be a S x 1 matrix or an S-vector. Currently, it is of size ", dim(N)))
    return(ErrorMessage)
    stop(ErrorMessage)
  }
  if(length(dim(DeltaT)) > 3){
    ErrorMessage <- (paste0("The argument DeltaT should be a S x 1 matrix or an S-vector. Currently, it is of size ", dim(DeltaT)))
    return(ErrorMessage)
    stop(ErrorMessage)
  }

  # Check DeltaTStar
  
  if(length(DeltaTStar) == 0){
    stop("DeltaTStar is not specified.")
  }
  
  if(length(DeltaTStar) != 1){
    ErrorMessage <- (paste0("The argument DeltaTStar should be a scalar (i.e. one number or a vector with one element). To inspect multiple DeltaTStar values, the analysis should be done for each value separately."))
    stop(ErrorMessage)
  }

  # Check other arguments
  #
  if(length(alpha) != 1){
    ErrorMessage <- (paste0("The argument alpha should be a scalar, that is, one number, that is, a vector with one element. Currently, alpha = ", alpha))
    return(ErrorMessage)
    stop(ErrorMessage)
  }
  #
  if(!is.logical(Moderators) & Moderators != FALSE & Moderators != TRUE){
    ErrorMessage <- (paste0("The argument Moderators should be logical, that is, have the value T(RUE) or F(ALSE); or 1 or 0; not ", Moderators))
    return(ErrorMessage)
    stop(ErrorMessage)
  }
  if(Moderators == TRUE){
    if(dim(Mod)[1] != S){
      ErrorMessage <- (paste0("The argument Mod should be a S times m matrix, with m the number of moderators to be included in the model.
                   Thus, the number of rows of Mod should equal S = ", S, " not ", dim(Mod)[1], "."))
      return(ErrorMessage)
      stop(ErrorMessage)
    }
  }
  #
  if(FEorRE != 1 & FEorRE != 2){
    ErrorMessage <- (paste0("The argument FEorRE should be 1 or 2; not ", FEorRE))
    return(ErrorMessage)
    stop(ErrorMessage)
  }
  if(FEorRE == 2 & !is.null(BetweenLevel)){
    if(length(BetweenLevel) != S){
      ErrorMessage <- (paste0("The argument BetweenLevel should be a S vector or S x 1 matrix.
                   Thus, the number of elements in BetweenLevel should equal S = ", S, " not ", length(BetweenLevel), "."))
      return(ErrorMessage)
      stop(ErrorMessage)
    }
  }
  #
  if(!is.null(Label)){
    if(length(Label) != q*q*S){
      ErrorMessage <- (paste0("The argument Label should be a q*q*S vector.
                     Thus, the number of elements in Label should equal q*q*S = ", q*q*S, " not ", length(Label), "."))
        return(ErrorMessage)
        stop(ErrorMessage)
    }
    if(!is.character(Label)){
      ErrorMessage <- (paste0("The argument Label should be a character (containing q*q*S labels/names)."))
      return(ErrorMessage)
      stop(ErrorMessage)
    }
  }
  #
  if(!is.logical(PrintPlot) & PrintPlot != FALSE & PrintPlot != TRUE){
    ErrorMessage <- (paste0("The argument 'PrintPlot' should be T(RUE) or F(ALSE); or 1 or 0; not ", PrintPlot))
    return(ErrorMessage)
    stop(ErrorMessage)
  }


  # Check on Phi
  
  # Return an error message is Phi is not numeric
  if (!is.numeric(Phi)) {
    stop("Phi contains non-numerical values.")
  }
  
  # Add a warning when Phi is an array with too many panels
  if (length(dim(Phi)) > 2 && dim(Phi)[3] > S) {
    warning(paste0("Phi is an array with ", dim(Phi)[3] - S, " more panels than S (the number of studies). The last ", dim(Phi)[3] - S, " panels were ignored."))
  }
  # Add a warning when Phi is a matrix with too many rows
  if (is.matrix(Phi) && dim(Phi)[1] > S*dim(Phi)[2]) {
    warning(paste0("Phi is a stacked matrix with ", dim(Phi)[1] - S*dim(Phi)[2], " more rows than S*q. The last ", dim(Phi)[1] - S*dim(Phi)[2], " rows were ignored."))
  }
  
  # Give a proper error when Phi is an array with the right number of elements but the wrong number of panels
  if (length(dim(Phi)) > 2 && dim(Phi)[3] < S) {
    stop(paste0("Phi is an array with ", S - dim(Phi)[3], " less panels than S (the number of studies). Phi should have S panels."))
  }
  
  if(length(Phi) == 1){
    ErrorMessage <- (paste0("The argument Phi should not consist of one element: a meta-analysis on one single element is not meaningfull. Notably, it should be a stacked matrix of size S*q times q or array with dimensions q times q times S, with S = ", S, " the number of primary studies and q = ", q, " the number of variables."))
    return(ErrorMessage)
    stop(ErrorMessage)
  }
  #
  if(is.null(dim(Phi))){
    if(!is.null(length(Phi))){
      ErrorMessage <- (paste0("The argument Phi is not a stacked matrix of size S*q times q or array with dimensions q times q times S, with q = ", q, " and S = ", S, ". Currently, it is a vector with ", length(Phi), "elements."))
      return(ErrorMessage)
      stop(ErrorMessage)
    }else{
      ErrorMessage <- (paste0("The argument Phi is not found: The lagged effects matrix Phi is unknown, but should be part of the input."))
      return(ErrorMessage)
      stop(ErrorMessage)
    }
  }

  if(length(Phi) == S){
    q <- 1
  }else{
    q <- dim(Phi)[2]
  }

  # Phi
  if(length(dim(Phi)) < 2){
    ErrorMessage <- (paste0("The lagged effects matrix Phi should be an S*q times q matrix or a q times q times S array, with S = ", S, " and q = ", q))
    return(ErrorMessage)
    stop(ErrorMessage)
  }
  if(length(dim(Phi)) == 2){
    Phi_studies <- array(data=NA, dim=c(q,q,S))
    teller <- 0
    for(s in 1:S){
      Phi_studies[1:q,1:q,s] <- matrix(t(Phi)[(teller+1):(teller+q*q)], byrow = T, ncol=q*q)
      Phi_studies[1:q,1:q,s] <- t(Phi_studies[1:q,1:q,s])
      teller <- teller + q*q
    }
    Phi <- Phi_studies
  }else if(length(dim(Phi)) > 3){
    ErrorMessage <- (paste0("The lagged effects matrix Phi should be an S*q times q matrix or a q times q times S array, with S = ", S, " and q = ", q, ". Currently, it is of size ", dim(Phi)))
    return(ErrorMessage)
    stop(ErrorMessage)
  }

  if(is.null(SigmaVAR) & is.null(Gamma)){ # Both SigmaVAR and Gamma unknown
    ErrorMessage <- (paste0("The arguments SigmaVAR and Gamma are not found: Both SigmaVAR and Gamma are unknown; either one (or both) should be part of the input. In case of first matrix, specify 'SigmaVAR = SigmaVAR'."))
    return(ErrorMessage)
    stop(ErrorMessage)
  }else if(is.null(Gamma)){ # Gamma unknown, calculate Gamma from SigmaVAR and Phi

    # Check on SigmaVAR
    if(length(SigmaVAR) == 1){
      ErrorMessage <- (paste0("The argument SigmaVAR should not consist of one element. Notably, it should be a stacked matrix of size S*q times q or array with dimensions q times q times S, with S = ", S, " the number of primary studies and q = ", q, " the number of variables."))
      return(ErrorMessage)
      stop(ErrorMessage)
    }
    
    # Add a warning when SigmaVAR is a matrix with too many rows
    if (is.matrix(SigmaVAR) && dim(SigmaVAR)[1] > S*q) {
      warning(paste0("SigmaVAR is a stacked matrix with ", dim(SigmaVAR)[1] - S*q, " more rows than S*q. The last ", dim(SigmaVAR)[1] - S*q, " rows were ignored."))
    }
    # Add a warning when SigmaVAR is a matrix with too many columns
    if (is.matrix(SigmaVAR) && dim(SigmaVAR)[2] > q) {
      warning(paste0("SigmaVAR is a stacked matrix with ", dim(SigmaVAR)[2] - q, " more columns than q. The last ", dim(SigmaVAR)[2] - q, " columns were ignored."))
    }

    # SigmaVAR
    if(length(dim(SigmaVAR)) < 2){
      ErrorMessage <- (paste0("The residual covariance matrix SigmaVAR should be an (S*q)*q matrix or a q*q*S array, with S = ", S, " and q = ", q))
      stop(ErrorMessage)
    }
    if(length(dim(SigmaVAR)) == 2){
      SigmaVAR_studies <- array(data=NA, dim=c(q,q,S))
      teller <- 0
      for(s in 1:S){
        SigmaVAR_studies[1:q,1:q,s] <- SigmaVAR[(teller+1):(teller+q),1:q]
        teller <- teller + q
      }
      SigmaVAR <- SigmaVAR_studies
    }else if(length(dim(SigmaVAR)) > 3){
      ErrorMessage <- (paste0("The residual covariance matrix SigmaVAR should be an (S*q)*q matrix or a q*q*S array, with S = ", S, " and q = ", q, ". Currently, it is of size ", dim(SigmaVAR)))
      return(ErrorMessage)
      stop(ErrorMessage)
    }

    # Calculate Gamma
    # If Phi and SigmaVAR are known, one can calculate Gamma:
    #Gamma <- array(data=NA, dim=c(S*q,q))
    #teller <- 1
    #for(s in 1:S){
    #  Gamma[teller:(teller+1),] <- Gamma.fromVAR(Phi[teller:(teller+1),], SigmaVAR[teller:(teller+1),])
    #  teller <- teller + q
    #}
    
    # Check SigmaVAR is numerical otherwise we won't be able to get Gamma
    if (!is.null(SigmaVAR) && !is.numeric(SigmaVAR)) {
      stop("SigmaVAR contains non-numerical values.")
    }
    
    Gamma_studies <- array(data=NA, dim=c(q,q,S))
    for(s in 1:S){
      Gamma_studies[1:q,1:q,s] <- Gamma.fromVAR(Phi[1:q,1:q,s], SigmaVAR[1:q,1:q,s])
    }
    Gamma <- Gamma_studies
    
    # Add a warning when SigmaVAR is an array with too many panels
    if (length(dim(SigmaVAR)) > 2 && dim(SigmaVAR)[3] > S) {
      warning(paste0("SigmaVAR is an array with ", dim(SigmaVAR)[3] - S, " more panels than S (the number of studies). The last ", dim(SigmaVAR)[3] - S, " panels were ignored."))
    }

  }else if(is.null(SigmaVAR)){ # SigmaVAR unknown, calculate SigmaVAR from Gamma and Phi

    # Check on Gamma
    if(length(Gamma) == 1){
      ErrorMessage <- (paste0("The argument Gamma should not consist of one element. Notably, it should be a stacked matrix of size S*q times q or array with dimensions q times q times S, with S = ", S, " the number of primary studies and q = ", q, " the number of variables."))
      stop(ErrorMessage)
    }

    # Gamma
    if(length(dim(Gamma)) < 2){
      ErrorMessage <- (paste0("The covariance matrix Gamma should be an S*q times q matrix or a q times q times S array, with S = ", S, " and q = ", q))
      return(ErrorMessage)
      stop(ErrorMessage)
    }
    if(length(dim(Gamma)) == 2){
      Gamma_studies <- array(data=NA, dim=c(q,q,S))
      teller <- 0
      for(s in 1:S){
        Gamma_studies[1:q,1:q,s] <- Gamma[(teller+1):(teller+q),1:q]
        teller <- teller + q
      }
      Gamma <- Gamma_studies
    }else if(length(dim(Gamma)) > 3){
      ErrorMessage <- (paste0("The covariance matrix Gamma should be an S*q times q matrix or a q times q times S array, with S = ", S, " and q = ", q, ". Currently, it is of size ", dim(Gamma)))
      return(ErrorMessage)
      stop(ErrorMessage)
    }

    # Calculate SigmaVAR
    # If Phi and SigmaVAR are known, one can calculate Gamma:
    #SigmaVAR <- array(data=NA, dim=c(S*q,q))
    #teller <- 1
    #for(s in 1:S){
    #  SigmaVAR[teller:(teller+1),] <- Gamma[teller:(teller+1),] - Phi[teller:(teller+1),] %*% Gamma[teller:(teller+1),] %*% t(Phi[teller:(teller+1),])
    #  teller <- teller + q
    #}
    SigmaVAR_studies <- array(data=NA, dim=c(q,q,S))
    for(s in 1:S){
      Gam <- Gamma[1:q,1:q,s]
      Ph <- Phi[1:q,1:q,s]
      SigmaVAR_studies[1:q,1:q,s] <- Gam - Ph %*% Gam %*% t(Ph)
    }
    SigmaVAR <- SigmaVAR_studies

  }else{ # Both SigmaVAR and Gamma known

    # Check on SigmaVAR
    if(length(SigmaVAR) == 1){
      ErrorMessage <- (paste0("The argument SigmaVAR should not consist of one element. It should be a stacked matrix of size S*q times q or array with dimensions q times q times S, with S = ", S, " the number of primary studies and q = ", q, " the number of variables."))
      return(ErrorMessage)
      stop(ErrorMessage)
    }
    # Check on Gamma
    if(length(Gamma) == 1){
      ErrorMessage <- (paste0("The argument Gamma should not consist of one element. It should be a stacked matrix of size S*q times q or array with dimensions q times q times S, with S = ", S, " the number of primary studies and q = ", q, " the number of variables."))
      return(ErrorMessage)
      stop(ErrorMessage)
    }

    # SigmaVAR
    if(length(dim(SigmaVAR)) < 2){
      ErrorMessage <- (paste0("The residual covariance matrix SigmaVAR should be an S*q times q matrix or a q times q times S array, with S = ", S, " and q = ", q))
      return(ErrorMessage)
      stop(ErrorMessage)
    }
    if(length(dim(SigmaVAR)) == 2){
      SigmaVAR_studies <- array(data=NA, dim=c(q,q,S))
      teller <- 0
      for(s in 1:S){
        SigmaVAR_studies[1:q,1:q,s] <- SigmaVAR[(teller+1):(teller+q),1:q]
        teller <- teller + q
      }
      SigmaVAR <- SigmaVAR_studies
    }else if(length(dim(SigmaVAR)) > 3){
      ErrorMessage <- (paste0("The residual covariance matrix SigmaVAR should be an S*q times q matrix or a q times q times S array, with S = ", S, " and q = ", q, ". Currently, it is of size ", dim(SigmaVAR)))
      return(ErrorMessage)
      stop(ErrorMessage)
    }

    # Gamma
    if(length(dim(Gamma)) < 2){
      ErrorMessage <- (paste0("The covariance matrix Gamma should be an S*q times q matrix or a q times q times S array, with S = ", S, " and q = ", q))
      return(ErrorMessage)
      stop(ErrorMessage)
    }
    if(length(dim(Gamma)) == 2){
      Gamma_studies <- array(data=NA, dim=c(q,q,S))
      teller <- 0
      for(s in 1:S){
        Gamma_studies[1:q,1:q,s] <- Gamma[(teller+1):(teller+q),1:q]
        teller <- teller + q
      }
      Gamma <- Gamma_studies
    }else if(length(dim(Gamma)) > 3){
      ErrorMessage <- (paste0("The covariance matrix Gamma should be an S*q times q matrix or a q times q times S array, with S = ", S, " and q = ", q, ". Currently, it is of size ", dim(Gamma)))
      return(ErrorMessage)
      return(ErrorMessage)
      stop(ErrorMessage)
    }

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
      messageTrans <- "All eigenvalues are positive and real. Hence, the Phi's are transformed to Phi(DeltaT*) to account for the time-interval dependency and standardized (to make comparison of effects meaningful)."
      for(s in 1:S){
        EV <- eigen(Phi[,,s])$values
        if(any(is.complex(EV))){ StudiesComplexEV <- c(StudiesComplexEV, s) }
        if(any(Re(EV) < 0)){ StudiesNegEV <- c(StudiesNegEV, s) }
      }
      # If ratioDeltaT is integer then also complex or neg eigenvalues work...
      ratioDeltaT <- DeltaTStar / DeltaT
      # If the studies for which the EV are neg or complex have integer ratioDeltaT then it is fine, otherwise it is not
      if(!is.null(StudiesComplexEV) | !is.null(StudiesNegEV)){
        messageTrans <- "Some studies have eigenvalues of Phi that are negative and/or real, but for those studies the ratio DeltaT*/DeltaT is integer. Therefore, the Phi's are transformed to Phi(DeltaT*) to account for the time-interval dependency."
        #if(any(!is.integer(ratioDeltaT[StudiesComplexEV]))){ # NB ik moet checken of als ik deel door hoogste geheel getal er dan geen rest waarde is!!
        if( any( (ratioDeltaT[StudiesComplexEV] %% 1) == 0 ) ){
          Trans <- 0
          messageTrans <- "There is at least one study for which the eigenvalues of Phi are complex and for which the ratio DeltaT*/DeltaT is not integer. Therefore, dummy variables are used to account for the time-interval dependency."
        }
        if( any( (ratioDeltaT[StudiesNegEV] %% 1) == 0 ) ){
          Trans <- 0
          messageTrans <- "There is at least one study for which the eigenvalues of Phi are negative and for which the ratio DeltaT*/DeltaT is not integer. Therefore, dummy variables are used to account for the time-interval dependency."
        }
        if( any( (ratioDeltaT[StudiesComplexEV] %% 1) == 0 ) & any( (ratioDeltaT[StudiesNegEV] %% 1) == 0 ) ){
          Trans <- 0
          messageTrans <- "There is at least one study for which the eigenvalues of Phi are complex and/or negative and for which the ratio DeltaT*/DeltaT is not integer. Therefore, dummy variables are used to account for the time-interval dependency."
        }
      }
      #
      #
      messageMultivar <- "For each study, the covariance matrix is positive definite. Hence, a multivariate approach is used."
      if(Trans == 1){
        # Calculate standardized transformed Phi
        messagePhi <- matrix(NA, ncol = 1, nrow = S)
        #s = 0
        #while(s < S & Trans == 1){ # go through all, unless a warning of not unique solution (which I already cover above...)
        #  s = s +1
        for(s in 1:S){
          out <- StandTransPhi(DeltaTStar, DeltaT[s], N[s], Phi[,,s], Gamma = Gamma[,,s], SigmaVAR = SigmaVAR[,,s])
          vecStandPhi[,s] <- out$vecStandPhi_DeltaTStar
          CovMxPhi[,,s] <- out$CovMx_vecStandPhi_DeltaTStar
          messagePhi[s] <- out$warning
          #Trans <- 1 - as.numeric(out$Warning)
          #
          # If covariance matrix is not positive definite then do univariate approach (and use only variances)
          if(any(eigen(CovMxPhi[,,s])$val < 0)){
            StudiesCovMxNotPosDef <- c(StudiesCovMxNotPosDef, s) # for which studies is this the case
            Multivar <- 0 # we need to conduct a univariate meta-analysis now
          }
        }
        if(!is.null(StudiesCovMxNotPosDef)){
          messageMultivar <- "There is at least one study for which the covariance matrix is not positive definite. Therefore, a univariate approach is used."
        }
      }
      #
      if(Trans == 0){ # then we cannot transform all the Phi's (uniquely), so use dummy method
        # Calculate standardized Phi
        for(s in 1:S){
          out <- StandPhi(N[s], Phi[,,s], Gamma[,,s], SigmaVAR[,,s])
          vecStandPhi[,s] <- out$vecStandPhi_DeltaT
          CovMxPhi[,,s] <- out$CovMx_vecStandPhi_DeltaT
          #
          # If covariance matrix is not positive definite then do univariate approach (and use only variances)
          if(any(eigen(CovMxPhi[,,s])$val < 0)){
            StudiesCovMxNotPosDef <- c(StudiesCovMxNotPosDef, s) # for which studies is this the case
            Multivar <- 0 # we need to conduct a univariate meta-analysis now
          }
          if(!is.null(StudiesCovMxNotPosDef)){
            messageMultivar <- "There is at least one study for which the covariance matrix is not positive definite. Therefore, a univariate approach is used."
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
      }else{ # if multivar, then we need to add an overallPhi variable and 'copy' the moderators and grouping variables
        sub = NULL
        for(i in 1:q){
          sub = c(sub, paste0(i, 1:q, sep=""))
        }
        overallPhi <- rep(sub, S)
        #
        Label_S <- paste0("Study ", rep(1:S, each = q*q))
        Label_S
        Label_Phi <- paste0(Label_S, " Phi", overallPhi)
        #
        #
        if(Moderators == 1){
          Mod. <- matrix(rep(as.matrix(Mod), each = q*q), ncol = dim(Mod)[2])
          colnames(Mod.) <- colnames(Mod)
        }
        #
        if(FEorRE == 2){
          if(Multivar == 1){
            Study <- matrix(rep(1:S, each = q*q), ncol = 1)
            if(!is.null(BetweenLevel)){
              BetweenLevel <- rep(BetweenLevel, each = q*q)
            }
          }else{
            Study <- matrix(rep(1:S), ncol = 1)
          }
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
            if(is.null(BetweenLevel)){
              metaan <- rma.mv(yi=vecVecStandPhi, V=CovMx, mods = ~ overallPhi + overallPhi:Mod. - 1,
                               random = ~ overallPhi | Study,
                               slab = Label_Phi,
                               struct = "UN", method = "ML")  # With struct="UN", the random effects are allowed to have different variances for each overallPhi and are allowed to be correlated.
            }else{
              metaan <- rma.mv(yi=vecVecStandPhi, V=CovMx, mods = ~ overallPhi + overallPhi:Mod. - 1,
                               random = ~ 1 | BetweenLevel / Study,
                               slab = Label_Phi,
                               method = "ML")
            }
          }else{ # No Moderators
            if(is.null(BetweenLevel)){
              metaan <- rma.mv(yi=vecVecStandPhi, V=CovMx, mods = ~ overallPhi - 1,
                               random = ~ overallPhi | Study,
                               slab = Label_Phi,
                               struct = "UN", method = "ML")  # With struct="UN", the random effects are allowed to have different variances for each overallPhi and are allowed to be correlated.
            }else{
              metaan <- rma.mv(yi=vecVecStandPhi, V=CovMx, mods = ~ overallPhi - 1,
                               random = ~ 1 | BetweenLevel / Study,
                               slab = Label_Phi,
                               method = "ML")
            }
          }
          tau2_metaan_MV <- metaan$tau2
        }else{# FE
          if(Moderators == 1){ # Moderators
            metaan <- rma.mv(yi=vecVecStandPhi, V=CovMx, mods = ~ overallPhi - 1 + overallPhi:Mod.,
                             slab = Label_Phi,
                             method = "FE")
          }else{ # No Moderators
            metaan <- rma.mv(yi=vecVecStandPhi, V=CovMx, mods = ~ overallPhi - 1,
                             slab = Label_Phi,
                             method = "FE")
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
          vecPhi <- Phi_metaan_MV
      #
      # Multivar and Dummies
      }else if(Multivar == 1 & Trans == 0){
        vecVecStandPhi <- as.vector(vecStandPhi) # vecStandPhi[q:q*q,1:S]
        if(FEorRE == 2){ # RE
          if(G > 1){
            if(Moderators == 1){ # Moderators
              if(is.null(BetweenLevel)){
                metaan <- rma.mv(yi=vecVecStandPhi, V=CovMx, mods = ~ -1 + overallPhi:D. + overallPhi:Mod.,
                                 random = ~ overallPhi | Study,
                                 slab = Label_Phi,
                                 struct = "UN", method = "ML")  # With struct="UN", the random effects are allowed to have different variances for each overallPhi and are allowed to be correlated.
              }else{
                metaan <- rma.mv(yi=vecVecStandPhi, V=CovMx, mods = ~ -1 + overallPhi:D. + overallPhi:Mod.,
                                 random = ~ 1 | BetweenLevel / Study,
                                 slab = Label_Phi,
                                 method = "ML")
              }
            }else{ # No Moderators
              if(is.null(BetweenLevel)){
                metaan <- rma.mv(yi=vecVecStandPhi, V=CovMx, mods = ~ -1 + overallPhi:D.,
                                 random = ~ overallPhi | Study,
                                 slab = Label_Phi,
                                 struct = "UN", method = "ML")  # With struct="UN", the random effects are allowed to have different variances for each overallPhi and are allowed to be correlated.
              }else{
                metaan <- rma.mv(yi=vecVecStandPhi, V=CovMx, mods = ~ -1 + overallPhi:D.,
                                 random = ~ 1 | BetweenLevel / Study,
                                 slab = Label_Phi,
                                 method = "ML")
              }
            }
            tau2_metaan_MV <- metaan$tau2
          } else{ # G = 1, no use for dummies (only one group)
            if(Moderators == 1){ # Moderators
              if(is.null(BetweenLevel)){
                metaan <- rma.mv(yi=vecVecStandPhi, V=CovMx, mods = ~ overallPhi - 1 + overallPhi:Mod.,
                                 random = ~ overallPhi | Study,
                                 slab = Label_Phi,
                                 struct = "UN", method = "ML")  # With struct="UN", the random effects are allowed to have different variances for each overallPhi and are allowed to be correlated.
              }else{
                metaan <- rma.mv(yi=vecVecStandPhi, V=CovMx, mods = ~ overallPhi - 1 + overallPhi:Mod.,
                                 random = ~ 1 | BetweenLevel / Study,
                                 slab = Label_Phi,
                                 method = "ML")
              }
            }else{ # No Moderators
              if(is.null(BetweenLevel)){
                metaan <- rma.mv(yi=vecVecStandPhi, V=CovMx, mods = ~ overallPhi - 1,
                                 random = ~ overallPhi | Study,
                                 slab = Label_Phi,
                                 struct = "UN", method = "ML")  # With struct="UN", the random effects are allowed to have different variances for each overallPhi and are allowed to be correlated.
              }else{
                metaan <- rma.mv(yi=vecVecStandPhi, V=CovMx, mods = ~ overallPhi - 1,
                                 random = ~ 1 | BetweenLevel / Study,
                                 slab = Label_Phi,
                                 method = "ML")
              }
            }
          }
        }else{# FE
          if(G > 1){
            if(Moderators == 1){ # Moderators
              metaan <- rma.mv(yi=vecVecStandPhi, V=CovMx, mods = ~ -1 + overallPhi:D. + overallPhi:Mod.,
                               slab = Label_Phi,
                               method = "FE")
            }else{ # No Moderators
              metaan <- rma.mv(yi=vecVecStandPhi, V=CovMx, mods = ~ -1 + overallPhi:D.,
                               slab = Label_Phi,
                               method = "FE")
            }
          } else{ # G = 1, no use for dummies (only one group)
            if(Moderators == 1){ # Moderators
              metaan <- rma.mv(yi=vecVecStandPhi, V=CovMx, mods = ~ overallPhi - 1 + overallPhi:Mod.,
                               slab = Label_Phi,
                               method = "FE")
            }else{ # No Moderators
              metaan <- rma.mv(yi=vecVecStandPhi, V=CovMx, mods = ~ overallPhi - 1,
                               slab = Label_Phi,
                               method = "FE")
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
        vecPhi <- Phi_metaan_MV
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
        summaryMetaAnalysis_jk <- list()
        for(teller in 1:(q^2)){ # Do meta-analysis (over S studies) for each of the q*q SLEs
          j = (teller+q-1)%/%q
          k = teller-(j-1)*q
          Phi_jk <- vecStandPhi[teller,]
          var <- varPhi[teller,]
          if(FEorRE == 2){ # RE
            if(Moderators == 1){ # Moderators
              if(is.null(BetweenLevel)){
                metaan <- rma.uni(yi=Phi_jk, vi=var, mods = Mod., method = "ML")
              }else{
                metaan <- rma.mv(yi=Phi_jk, V=diag(var), mods = Mod., method = "ML",
                                 random = ~ 1 | BetweenLevel/Study)
              }
            }else{ # No Moderators
              if(is.null(BetweenLevel)){
                metaan <- rma.uni(yi=Phi_jk, vi=var, method = "ML")
              }else{
                metaan <- rma.mv(yi=Phi_jk, V=diag(var), method = "ML",
                                 random = ~ 1 | BetweenLevel/Study)
              }
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
          NameList <- paste0("j", j, "_k", k)
          summaryMetaAnalysis_jk[[NameList]] <- summary(metaan)
        }
        vecPhi <- Phi_metaan
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
        summaryMetaAnalysis_jk <- list()
        for(teller in 1:(q^2)){ # Do meta-analysis (over S studies) for each of the q*q SLEs
          j = (teller+q-1)%/%q
          k = teller-(j-1)*q
          Phi_jk <- vecStandPhi[teller,]
          var <- varPhi[teller,]
          if(FEorRE == 2){ # RE
            if(G > 1){
              if(Moderators == 1){ # Moderators
                if(is.null(BetweenLevel)){
                  metaan_g <- rma.uni(yi=Phi_jk, vi=var, mods = ~ -1 + D. + Mod., method = "ML")
                }else{
                  metaan_g <- rma.mv(yi=Phi_jk, V=diag(var), mods = ~ -1 + D. + Mod., method = "ML",
                                   random = ~ 1 | BetweenLevel/Study)
                }
              }else{ # No Moderators
                if(is.null(BetweenLevel)){
                  metaan_g <- rma.uni(yi=Phi_jk, vi=var, mods = ~ -1 + D., method = "ML")
                }else{
                  metaan_g <- rma.mv(yi=Phi_jk, V=diag(var), mods = ~ -1 + D., method = "ML",
                                   random = ~ 1 | BetweenLevel/Study)
                }
              }
            } else{ # G == 1: no dummies needed
              if(Moderators == 1){ # Moderators
                if(is.null(BetweenLevel)){
                  metaan_g <- rma.uni(yi=Phi_jk, vi=var, mods = ~ 1 + Mod., method = "ML")
                }else{
                  metaan_g <- rma.mv(yi=Phi_jk, V=diag(var), mods = ~ 1 + Mod., method = "ML",
                                     random = ~ 1 | BetweenLevel/Study)
                }
              }else{ # No Moderators
                if(is.null(BetweenLevel)){
                  metaan_g <- rma.uni(yi=Phi_jk, vi=var, method = "ML")
                }else{
                  metaan_g <- rma.mv(yi=Phi_jk, V=diag(var), method = "ML",
                                     random = ~ 1 | BetweenLevel/Study)
                }
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
          NameList <- paste0("j", j, "_k", k)
          summaryMetaAnalysis_jk[[NameList]] <- summary(metaan_g)
        }
        vecPhi <- Phi_metaan_dum
      }
      #
      #
      #
      # Output
      #
      # Plot
      if(PrintPlot == TRUE){
        # Extract the q times q overall Phi matrix
        q <- sqrt(length(vecPhi))
        overallPhi <- matrix(vecPhi, byrow = T, ncol = q) # resulting overall Phi
        # Make Phi-plot:
        Title <- as.list(expression(Phi(Delta[t])~plot), "How do the overall lagged parameters vary as a function of the time-interval")
        min <- min(DeltaT)
        max <- max(DeltaT)
        step <- (max - min + 1)/50 #(max - min + 1)/10
        #PhiPlot(DeltaTStar, overallPhi, Min = min, Max = max, Step = step, Title = Title)
        phi_plot <- ggPhiPlot(DeltaTStar, overallPhi, Min = min, Max = max, Step = step, Title = Title)
        PhiPlot <- phi_plot$PhiPlot
        if(is.null(phi_plot$ErrorMessage)){
          PhiPlot <- phi_plot$PhiPlot
        }else{
          return(ErrorMessage)
          stop(ErrorMessage)
        }
      }
      #
      #
      # Other output
      dT <- unique(DeltaT)[order(unique(DeltaT))]
      #names(dT) <- "Unique DeltaT values used in the studies"
      dT_star <- DeltaTStar
      #names(dT_start) <- "DeltaT_star value"
      if(!is.null(StudiesComplexEV)) names(StudiesComplexEV) <- "Study indices for which (some of) the eigenvalues are complex"
      if(!is.null(StudiesNegEV)) names(StudiesNegEV) <- "Study indices for which (some of) the eigenvalues are negative"
      if(!is.null(StudiesCovMxNotPosDef)) names(StudiesCovMxNotPosDef) <- "Study indices for which the covariance matrix is not positive definite"
      ratioDeltaT <- matrix(ratioDeltaT)
      colnames(ratioDeltaT) <- "Ratio DeltaT*/DeltaT"
      rownames(ratioDeltaT) <- paste0("Study ", rep(1:S))
      # Multivar and Trans
      if(Multivar == 1 & Trans == 1){
        if(FEorRE == 1){ # FE
          final <- list(DeltaTStar = dT_star,
                        Overall_standPhi_DeltaTStar = matrix(Phi_metaan_MV, byrow = T, ncol = q),
                        Overall_vecStandPhi_DeltaTStar = Phi_metaan_MV,
                        #elliptical_CI = multiCI,
                        LB_elliptical_CI = multiCI[1,],
                        UB_elliptical_CI = multiCI[2,],
                        alpha_CI = alpha,
                        CovMx_OverallPhi_DeltaTStar = CovMxPhi_metaan,
                        messageTrans = messageTrans, messageMultivar = messageMultivar,
                        StudiesComplexEV = StudiesComplexEV, StudiesNegEV = StudiesNegEV, StudiesCovMxNotPosDef = StudiesCovMxNotPosDef,
                        ratioDeltaT = ratioDeltaT,
                        messagePhi = messagePhi,
                        summaryMetaAnalysis = summary(metaan))
        } else{ # RE
          final <- list(DeltaTStar = dT_star,
                        Overall_standPhi_DeltaTStar = matrix(Phi_metaan_MV, byrow = T, ncol = q),
                        Overall_vecStandPhi_DeltaTStar = Phi_metaan_MV,
                        #elliptical_CI = multiCI,
                        LB_elliptical_CI = multiCI[1,],
                        UB_elliptical_CI = multiCI[2,],
                        alpha_CI = alpha,
                        CovMx_OverallPhi_DeltaTStar = CovMxPhi_metaan,
                        tau2 = tau2_metaan_MV,
                        messageTrans = messageTrans, messageMultivar = messageMultivar,
                        StudiesComplexEV = StudiesComplexEV, StudiesNegEV = StudiesNegEV, StudiesCovMxNotPosDef = StudiesCovMxNotPosDef,
                        ratioDeltaT = ratioDeltaT,
                        messagePhi = messagePhi,
                        summaryMetaAnalysis = summary(metaan))
        }
      # Multivar and Dummies
      }else if(Multivar == 1 & Trans == 0){
        if(FEorRE == 1){ # FE
          final <- list(UniqueTimeIntervals = dT,
                        Overall_standPhi_PerUniqueTimeInterval = matrix(Phi_metaan_MV, byrow = T, ncol = q),
                        Overall_vecStandPhi_PerUniqueTimeInterval = Phi_metaan_MV,
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
                        Overall_standPhi_PerUniqueTimeInterval = matrix(Phi_metaan_MV, byrow = T, ncol = q),
                        Overall_vecStandPhi_PerUniqueTimeInterval = Phi_metaan_MV,
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
                        Overall_standPhi_DeltaTStar = matrix(Phi_metaan, byrow = T, ncol = q),
                        Overall_vecStandPhi_DeltaTStar = Phi_metaan,
                        LB_CI = minPhi,
                        UB_CI = maxPhi,
                        alpha_CI = alpha,
                        se_Overall_standPhi_DeltaTStar = sePhi_metaan,
                        pvalue_HeterogeneityTest = QEp_metaan,
                        pvalue_OmnibusCoefTest_IfModerator = QMp_metaan,
                        messageTrans = messageTrans, messageMultivar = messageMultivar,
                        StudiesComplexEV = StudiesComplexEV, StudiesNegEV = StudiesNegEV, StudiesCovMxNotPosDef = StudiesCovMxNotPosDef,
                        ratioDeltaT = ratioDeltaT,
                        summaryMetaAnalysis_jk = summaryMetaAnalysis_jk)
        } else{ # RE
          final <- list(DeltaTStar = dT_star,
                        Overall_standPhi_DeltaTStar = matrix(Phi_metaan, byrow = T, ncol = q),
                        Overall_vecStandPhi_DeltaTStar = Phi_metaan,
                        LB_CI = minPhi,
                        UB_CI = maxPhi,
                        alpha_CI = alpha,
                        se_Overall_standPhi_DeltaTStar = sePhi_metaan,
                        tau2 = tau2_metaan,
                        pvalue_HeterogeneityTest = QEp_metaan,
                        pvalue_OmnibusCoefTest_IfModerator = QMp_metaan,
                        messageTrans = messageTrans, messageMultivar = messageMultivar,
                        StudiesComplexEV = StudiesComplexEV, StudiesNegEV = StudiesNegEV, StudiesCovMxNotPosDef = StudiesCovMxNotPosDef,
                        ratioDeltaT = ratioDeltaT,
                        summaryMetaAnalysis_jk = summaryMetaAnalysis_jk)
        }
      # Univar and Dummies
      }else if(Multivar == 0 & Trans == 0){
        if(FEorRE == 1){ # FE
          final <- list(UniqueTimeIntervals = dT,
                        Overall_standPhi_PerUniqueTimeInterval = matrix(Phi_metaan_dum, byrow = T, ncol = q),
                        Overall_vecStandPhi_PerUniqueTimeInterval = Phi_metaan_dum,
                        LB_CI = minPhi,
                        UB_CI = maxPhi,
                        alpha_CI = alpha,
                        se_Overall_standPhi_PerUniqueTimeInterval = sePhi_metaan_dum,
                        pvalue_HeterogeneityTest = QEp_metaan,
                        pvalue_OmnibusCoefTest = QMp_metaan,
                        messageTrans = messageTrans, messageMultivar = messageMultivar,
                        StudiesComplexEV = StudiesComplexEV, StudiesNegEV = StudiesNegEV, StudiesCovMxNotPosDef = StudiesCovMxNotPosDef,
                        ratioDeltaT = ratioDeltaT,
                        summaryMetaAnalysis_jk = summaryMetaAnalysis_jk)
        } else{ # RE
          final <- list(UniqueTimeIntervals = dT,
                        Overall_standPhi_PerUniqueTimeInterval = matrix(Phi_metaan_dum, byrow = T, ncol = q),
                        Overall_vecStandPhi_PerUniqueTimeInterval = Phi_metaan_dum,
                        LB_CI = minPhi,
                        UB_CI = maxPhi,
                        alpha_CI = alpha,
                        se_Overall_standPhi_PerUniqueTimeInterval = sePhi_metaan_dum,
                        tau2 = tau2_metaan_dum,
                        pvalue_HeterogeneityTest = QEp_metaan,
                        pvalue_OmnibusCoefTest = QMp_metaan,
                        messageTrans = messageTrans, messageMultivar = messageMultivar,
                        StudiesComplexEV = StudiesComplexEV, StudiesNegEV = StudiesNegEV, StudiesCovMxNotPosDef = StudiesCovMxNotPosDef,
                        ratioDeltaT = ratioDeltaT,
                        summaryMetaAnalysis_jk = summaryMetaAnalysis_jk)
        }
      }

      if(PrintPlot == T){
        final[["PhiPlot"]] <- PhiPlot
        print(PhiPlot)
      }

      class(final) <- c("CTmeta", "list")
      final
    }
