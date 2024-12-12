#' Set up model, data, and initial values for running MCMC
#'
#' @param modelCode an R call object specifying the model structure for integrated 
#' distance sampling model
#' @param nim.data list of input objects representing data
#' @param nim.constants list of input objects representing constants
#' @param R_perF logical. If TRUE, treats recruitment rate as juvenile per adult female.
#' If FALSE, treats recruitment rate as juvenile per adult (sum of both sexes).
#' @param survVarT logical. If TRUE, survival is simulated including annual variation.
#' @param fitRodentCov logical. If TRUE, rodent covariate on reproduction is included.
#' @param addDummyDim logical. If TRUE (default) adds a dummy "area" dimension when 
#' simulating initial values for a single area implementation. This is necessary 
#' for the multi-area setup/model to run with data from only one area. 
#' @param niter integer. Number of MCMC iterations (default = 25000)
#' @param nthin integer. Thinning factor (default = 5)
#' @param nburn integer. Number of iterations to discard as burn-in (default = 5000)
#' @param nchains integer. Number of chains to run.
#' @param testRun logical. If TRUE, sets up for a test run with 10 iterations,
#' no thinning, and no burn-in (default = FALSE)
#' @param initVals.seed integer. Seed to use for inital value simulation.
#'
#' @return list of list containing all components necessary for running model 
#' with `nimble::nimbleMCMC()`
#' @export
#'
#' @examples

setupModel <- function(modelCode, customDist,
                       nim.data, nim.constants,
                       R_perF, survVarT, fitRodentCov, addDummyDim = TRUE,
                       niter = 200000, nthin = 30, nburn = 111000, nchains = 3,
                       testRun = FALSE, initVals.seed){
  
  require('nimble')
  require('nimbleDistance')
  
  ## Set parameters to monitor
  params <- c("esw", "p", #"D",
              "R_year", "Mu.R", "h.Mu.R", "h.sigma.R", "sigmaT.R", "sigmaR.R",
              "sigma", "mu.dd", "h.mu.dd", "h.sigma.dd", "sigmaT.dd", "sigmaR.dd",
              "meanDens", 
              "Mu.D1", "sigma.D",
              "S", "Mu.S", "h.Mu.S", "h.sigma.S",
              "Mu.S1")
  
  if(survVarT){
    params <- c(params, "sigmaT.S", "sigmaR.S", "eps.S1.prop")
  }
  
  if(fitRodentCov){
    params <- c(params, "betaR.R", "h.Mu.betaR.R", "h.sigma.betaR.R", "RodentOcc")
  }
  
  if(nim.constants$N_areas == 1 & !addDummyDim){
    hyperparam.idx <- which(params %in% c("h.Mu.R", "h.sigma.R", "h.Mu.S", "h.sigma.S", "h.Mu.betaR.R", "h.sigma.betaR.R", "h.mu.dd", "h.sigma.dd"))
    resvar.idx <- which(params %in% c("sigmaR.R", "sigmaR.S", "sigmaR.dd"))
    params <- params[-c(hyperparam.idx, resvar.idx)]
  }
  
  ## Simulate initial values
  #set.seed(initVals.seed)
  initVals <- list()
  for(c in 1:nchains){
    
    if(nim.constants$N_areas == 1 & !addDummyDim){
      
      initVals[[c]] <- simulateInits_singleArea(nim.data = nim.data, 
                                                nim.constants = nim.constants, 
                                                R_perF = R_perF,
                                                survVarT = survVarT,
                                                fitRodentCov = fitRodentCov,
                                                initVals.seed = initVals.seed[c])
    }else{
      
      initVals[[c]] <- simulateInits(nim.data = nim.data, 
                                     nim.constants = nim.constants, 
                                     R_perF = R_perF, 
                                     survVarT = survVarT,
                                     fitRodentCov = fitRodentCov,
                                     initVals.seed = initVals.seed[c])
    }

  }
  
  ## Adjust MCMC parameters if doing a test run
  if(testRun){
    niter <- 50
    nthin <- 1
    nburn <- 0
  }
  
  ## Collate model setup variables in a list
  setup <- list(
    modelCode = modelCode,
    modelParams = params,
    initVals = initVals,
    mcmcParams = list(niter = niter, nthin = nthin, 
                      nburn = nburn, nchains = nchains)
  )
   
  return(setup) 
}
