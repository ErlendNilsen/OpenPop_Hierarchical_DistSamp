#' Write integrated distance sampling model code
#'
#' @param survVarT logical. If TRUE, writes code for a model including random
#' year variation in survival probability. If FALSE, assumed constant survival
#' probability across time. 
#' @param telemetryData logical. If TRUE, uses information from telemetry data
#' from Lierne. If FALSE, only line transect data is used. 
#' @return an R call object specifying the model structure for integrated 
#' distance sampling model. 
#' @export
#'
#' @examples

writeModelCode_singleArea <- function(survVarT, telemetryData){
  
  IDSM.code <- nimble::nimbleCode({
    
    # N_sites = number of sites
    # N_ageC = number of age classes
    # N_years = number of years
    # N_sumR_obs = number of data points in juvenile:adult ratio counts
    
    # N_exp[a, j, t] = Number of age class a individuals in site j in year t
    # Density[a, j, t] = Density of age class a individuals in site j in year t
    # L[j, t] = length of transect line in site j in year t
    # W = truncation distance for line transect surveys
    
    # Mu.D1 = average initial density
    
    # S[t] = annual survival from year t to t+1
    # R_year[t] = recruitment rate in year t
    # p[t] = average distance sampling detection rate in year t
    # sigma[t] = average distance sampling detection decay rate in year t
    
    # eps.D1[j] = random site effect on initial density (site j)
    
    
    
    ####################
    # POPULATION MODEL #
    ####################
    
    #-----------------------------------------#
    # Initial population size/density (t = 1) #
    #-----------------------------------------#
    
    for (j in 1:N_sites){
      
      ## Adult densities
      Density[2, j, 1] <- exp(log(Mu.D1) + eps.D1[j])
      
      ## Juvenile densities
      if(R_perF){
        Density[1, j, 1] <- (Density[2, j, 1]/2)*R_year[1] 
      }else{
        Density[1, j, 1] <- Density[2, j, 1]*R_year[1]
      }
      
      ## Adult and juvenile numbers
      N_exp[1:N_ageC, j, 1] <- Density[1:N_ageC, j, 1]*L[j, 1]*W*2      
    }
    
    #-------------------------------#
    # Population dynamics for t > 1 #
    #-------------------------------#
    
    for(j in 1:N_sites){
      for(t in 2:N_years){
        
        ## Adult densities
        Density[2, j, t] <- sum(Density[1:N_ageC, j, t-1])*S[t-1] 
        
        ## Juvenile densities
        if(R_perF){
          Density[1, j, t] <- (Density[2, j, t]/2)*R_year[t] 
        }else{
          Density[1, j, t] <- Density[2, j, t]*R_year[t]
        }
        
        ## Adult and juvenile numbers
        N_exp[1:N_ageC, j, t] <- Density[1:N_ageC, j, t]*L[j, t]*W*2
      }
    }
    
    #--------------------#
    # Derived parameters #
    #--------------------#
    
    ## Area- and year-specific total densities
    for (t in 1:N_years){
      N_tot_exp[t] <- sum(N_exp[1, 1:N_sites, t] + N_exp[2, 1:N_sites, t])
    }
    
    ## Area-, year-, and age-class specific density (for monitoring)
    for(a in 1:N_ageC){
      for(t in 1:N_years){
        meanDens[a, t] <- mean(Density[a, 1:N_sites, t])
      } # t
    } # a
    
    
    ####################
    # DATA LIKELIHOODS #
    ####################
    
    ## Age-specific line transect counts
    # N_a_line_year[a, j, t] = number of age class a individuals detected in site j of in year t
    for(j in 1:N_sites){
      for(t in 1:N_years){
        for(a in 1:N_ageC){
          
          N_a_line_year[a, j, t] ~ dpois(p[t]*N_exp[a, j, t])
        }
      }
    }
    
    ## Juvenile:adult ratios from line transect observations
    # N_sumR_obs = number of observations in juvenile:adult count data
    # sumR_obs[i] = i'th entry in juvenile count data
    # sumAd_obs[i] = i'th entry in adult count data
    for (i in 1:N_sumR_obs){
      
      sumR_obs[i] ~ dpois(R_year[sumR_obs_year[i]]*sumAd_obs[i])
    }
    
    
    ## Line transect observation distances (likelihood using nimbleDistance::dHN)
    # N_obs = number of observations in detection distance data
    # y[i] = i'th entry in detection distance data
    for (i in 1:N_obs){ 
      
      y[i] ~ dHN(sigma = sigma[Year_obs[i]], Xmax = W, point = 0)
    }
    
    
    if(telemetryData){
      ## Known-fate telemetry data
      # N_years_RT = number of yeard for which radio-telemetry data is available
      # Survs1[t, z] = numbers of collared individuals released (z = 1) and 
      # recovered alive (z = 2) during season 1 of year t
      # Survs2[t, z] = numbers of collared individuals released (z = 1) and 
      # recovered alive (z = 2) during season 2 of year t
      # S1[t] = survival probability through season 1 of year t
      # S2[t] = survival probability through season 2 of year t
      
      for (t in 1:N_years_RT){
        
        Survs1[t, 2] ~ dbinom(S1[year_Survs[t]], Survs1[t, 1])
        Survs2[t, 2] ~ dbinom(S2[year_Survs[t]], Survs2[t, 1])
      }
    }
    
    
    ################################
    # PARAMETER MODELS/CONSTRAINTS #
    ################################
    
    ## Distance sampling detection parameters
    
    for(t in 1:N_years){
      
      # Detection decay
      log(sigma[t]) <- mu.dd + epsT.dd[t]
      sigma2[t] <- sigma[t] * sigma[t]
      
      # Effective strip width
      esw[t] <- sqrt(pi * sigma2[t] / 2) 
      
      # Average detection rate 
      p[t] <- min(esw[t], W) / W
    }
    
    
    
    ## Annual recruitment rates
    
    if(fitRodentCov){
      R_year[1:N_years] <- exp(log(Mu.R) + betaR.R*RodentOcc[1:N_years] + epsT.R[1:N_years])
    }else{
      R_year[1:N_years] <- exp(log(Mu.R) + epsT.R[1:N_years])
    }
    
    
    
    ## Annual survival probabilities
    
    if(survVarT){
      logit(S[1:(N_years-1)]) <- logit(Mu.S) + epsT.S[1:(N_years-1)]
    }else{
      S[1:(N_years-1)] <- Mu.S
    }
    
    
    ## Seasonal survival probabilities in area with radiotelemetry data
    # Season 1
    if(survVarT){
      logit(S1[1:(N_years-1)]) <- logit(Mu.S1) + eps.S1.prop*(epsT.S[1:(N_years-1)])
    }else{
      S1[1:(N_years-1)] <- Mu.S1
    }
    
    # Season 2
    S2[1:(N_years-1)] <- S[1:(N_years-1)] / S1[1:(N_years-1)]
    
    
    ###########
    # PRIORS  #
    ###########
    
    #-----------------------#
    # Intercepts / averages #
    #-----------------------#
    
    Mu.R  ~ dunif(0, 20) # Recruitment
    Mu.S ~ dunif(0, 1) # Survival
    mu.dd ~ dunif(-10, 100) # Detection
    
    ## Initial density
    Mu.D1 ~ dunif(0, 10)

    
    #----------------#
    # Random effects #
    #----------------#
    
    ## Standard deviations
    
    # Recruitment
    sigmaT.R ~ dunif(0, 5)

    # Survival 
    Mu.S1 ~ dunif(0, 1)
    
    if(survVarT){
      sigmaT.S ~ dunif(0, 5)
      eps.S1.prop ~ dunif(0, 1)
    }
    
    # Detection
    sigmaT.dd ~ dunif(0, 20)

    # Initial density
    sigma.D ~ dunif(0, 20)
    
    
    ## Random effect levels
    
    # Shared year variation
    for(t in 1:N_years){
      
      epsT.R[t] ~ dnorm(0, sd = sigmaT.R) # Recruitment
      epsT.dd[t] ~ dnorm(0, sd = sigmaT.dd) # Detection
    }
    
    for(t in 1:(N_years-1)){
      
      if(survVarT){
        epsT.S[t] ~ dnorm(0, sd = sigmaT.S) # Survival
      }
    }
    
    # Site/transect variation
    for(j in 1:N_sites){
      eps.D1[j] ~ dnorm(0, sd = sigma.D)
    }
    
    #-------------------#
    # Covariate effects #
    #-------------------#
    
    ## Rodent effect on reproduction
    if(fitRodentCov){
      betaR.R ~ dunif(-5, 5)
    }
    
    
    #------------------#
    # Other parameters #
    #------------------#
    
    pi <- 3.141593
    
    
    ###############################
    # COVARIATE IMPUTATION MODELS #
    ###############################
    
    if(fitRodentCov){
      for (t in 1:N_years){
        
        RodentOcc[t] ~ dnorm(mean = 0, sd = 1)
      }
    }
    
    
  })
  
  return(IDSM.code)
}