#' Simulate values for initializing integrated model for a single area
#'
#' @param nim.data list of input objects representing data
#' @param nim.constants list of input objects representing constants
#' @param R_perF logical. If TRUE, treats recruitment rate as juvenile per adult female.
#' If FALSE, treats recruitment rate as juvenile per adult (sum of both sexes).
#' @param survVarT logical. If TRUE, survival is simulated including annual variation.
#' @param fitRodentCov logical. If TRUE, initial values are generated for rodent 
#' covariate (effect) and covariate is included in data simulation.
#' 
#' @return A list containing one complete set of initial values for the model.
#' @export
#'
#' @examples

simulateInits_singleArea <- function(nim.data, nim.constants, R_perF, survVarT, fitRodentCov, initVals.seed){
  
  set.seed(initVals.seed)
  
  # Limits and constants #
  #----------------------#
  
  if(nim.constants$N_areas > 1){
    stop("Input data contains > 1 area. Use function simulateInits() instead. ")
  }
  
  N_areas <- nim.constants$N_areas
  N_ageC <- nim.constants$N_ageC
  N_years <- nim.constants$N_years
  N_sites <- nim.constants$N_sites
  
  L <- nim.data$L
  W <- nim.constants$W
  pi <- 3.141593
  A <- nim.data$A
  
  
  # Missing covariate values #
  #--------------------------#
  
  if(!is.null(nim.data$RodentOcc)){
    RodentOcc <- nim.data$RodentOcc
  }else{
    RodentOcc <- rep(0, N_years)
  }
  
  if(NA %in% RodentOcc){
    RodentOcc[which(is.na(RodentOcc))] <- runif(length(which(is.na(RodentOcc))), 0, 1)
  }
  
  Inits_RodentOcc <- RodentOcc
  Inits_RodentOcc[which(!is.na(nim.data$RodentOcc))] <- NA
  
  
  # Vital rates #
  #-------------#
  
  ## Area-specific survival parameters
  Mu.S1 <- runif(1, 0.6, 0.7)
  Mu.S <- runif(1, 0, Mu.S1)
  sigmaT.S <- runif(1, 0, 0.05)
  
  eps.S1.prop <- runif(1, 0.3, 0.8)
  epsT.S <- rep(0, N_years)

  S <- plogis(qlogis(Mu.S) + epsT.S[1:N_years])
  S1 <- plogis(qlogis(Mu.S1) + eps.S1.prop*epsT.S[1:N_years])
  S2 <- S/S1
  
  ## Area-specific reproductive parameters
  Mu.R <- runif(1, 1, 4)
  
  if(fitRodentCov){
    betaR.R <- runif(1, 0, 0.1)
  }else{
    betaR.R <- 0
  }
  
  sigmaT.R <- runif(1, 0, 2)
  
  epsT.R <- rep(0, N_years)
  #epsT.R <- rnorm(N_year, 0, sigmaT.R)
    
  R_year <- exp(log(Mu.R) + betaR.R*RodentOcc[1:N_years] + epsT.R[1:N_years])

  
  # Detection parameters #
  #----------------------#
  
  ## Area-specific detection parameters
  mu.dd <- runif(1, 2, 5)
  sigmaT.dd <- runif(1, 1, 10)
  
  epsT.dd <- rep(0, N_years)
  #epsT.dd <- rnorm(N_years, 0, sd = sigmaT.dd)
    
  sigma <- exp(mu.dd + epsT.dd[1:N_years])
  esw <- sqrt(pi * sigma^2 / 2) 
  
  p <- rep(NA, N_years)
  for(t in 1:N_years){
    p[t] <- min(esw[t], W) / W
  }

  
  # Population model #
  #------------------#
  
  ## Initial densities / population sizes
  sigma.D <- runif(1, 0.05, 2)
  ratio.JA1 <- runif(1, 0.2, 0.6)
  
  N_exp <- Density <- array(0, dim = c(N_ageC, max(N_sites), N_years))
  
  D_x_sum <- apply(nim.data$N_a_line_year, c(2,3), sum) / (L*W*2)
  D_data <- D_x_sum[which(!is.na(D_x_sum) & D_x_sum > 0)]
  Mu.D1 <- runif(1, quantile(D_data, 0.25), quantile(D_data, 0.75))  
  
  for(j in 1:N_sites){
    
    Density[1, j, 1] <- Mu.D1*ratio.JA1
    Density[2, j, 1] <- Mu.D1*(1-ratio.JA1)
    
    N_exp[1, j, 1] <- extraDistr::rtpois(1, lambda = Density[1, j, 1]*L[j, 1]*W*2, a = 1)
    N_exp[2, j, 1] <- extraDistr::rtpois(1, lambda = Density[2, j, 1]*L[j, 1]*W*2, a = 1)
  }
  
  ## Population projection over time
  for(j in 1:N_sites){
    for(t in 2:N_years){
      
      Density[2, j, t] <- sum(Density[1:N_ageC, j, t-1])*S[t-1] # Adults

      if(R_perF){
        Density[1, j, t] <- (Density[2, j, t]/2)*R_year[t] 
      }else{
        Density[1, j, t] <- Density[2, j, t]*R_year[t]
      }
      
      N_exp[1:N_ageC, j, t] <- Density[1:N_ageC, j, t]*L[j, t]*W*2
    }
  }
  
  ## Area-specific population size and density
  N_tot_exp <- rep(NA, N_years)
  
  for (t in 1:N_years){
    N_tot_exp[t] <- sum(N_exp[1, 1:N_sites, t] + N_exp[2, 1:N_sites, t])    ## Summing up expected number of birds in covered area; 
  }
  
  # Assembly #
  #----------#
  
  InitVals <- list(
    
    Mu.D1 = Mu.D1, 
    sigma.D = sigma.D,
    eps.D1 = rep(0, max(N_sites)),
    ratio.JA1 = ratio.JA1,
    
    Mu.R = Mu.R,
    sigmaT.R = sigmaT.R,
    sigmaR.R = 0,
    epsT.R = epsT.R,
    epsR.R = rep(0, max(N_sites)),
    R_year = R_year,
    
    mu.dd = mu.dd,
    sigmaT.dd = sigmaT.dd,
    sigmaR.dd = 0,
    epsT.dd = epsT.dd,
    epsR.dd = rep(0, max(N_sites)),
    sigma = sigma, sigma2 = sigma^2,
    esw = esw,
    p = p,
    
    Mu.S = Mu.S, 
    sigmaT.S = sigmaT.S,
    sigmaR.S = 0,
    epsT.S = epsT.S,
    epsR.S = rep(0, max(N_sites)),
    Mu.S1 = Mu.S1,
    eps.S1.prop = eps.S1.prop,
    S1 = S1, S2 = S2, S = S,
    
    Density = Density,
    N_exp = N_exp,
    N_tot_exp = N_tot_exp
  )
  
  if(fitRodentCov){
    InitVals$betaR.R <- betaR.R
    InitVals$RodentOcc <- Inits_RodentOcc
  }
  
  return(InitVals)
}
