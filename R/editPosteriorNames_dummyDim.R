#' Helper function that reformats a single area model fit without a dummy dimension to work with plotting code
#'
#' @param IDSM.out.tidy an mcmc list containing posterior samples from a model run.
#' @param N_areas integer. The number of areas used in the model. 
#' Function only relevant when N_areas = 1.
#' @param addDummyDim logical. No default. Indicator of whether a dummy 
#' dimension was added for analysis (TRUE) or not (FALSE). Function only 
#' relevant when addDummyDim = FALSE. 
#'
#' @return
#' @export
#'
#' @examples
#' 
editPosteriorNames_dummyDim <- function(IDSM.out.tidy, N_areas, addDummyDim){
  
  
  if(N_areas > 1 | addDummyDim){
    stop("This function is only relevant for plotting outputs of model runs with one area/locality AND addDummyDim = FALSE")
  }
  
  orig_names <- unlist(dimnames(IDSM.out.tidy[[1]])[2])
  new_names <- orig_names
  
  for(i in 1:length(orig_names)){
    
    if(orig_names[i] == "Mu.S1" | grepl("sigmaT", orig_names[i])){
      next()
    }
    
    if(grepl("Mu", orig_names[i]) | 
       grepl("mu", orig_names[i]) |
       orig_names[i] %in% c("betaR.R", "sigma.D")){
      new_names[i] <- paste0(orig_names[i], "[1]")
      
    }else{
      fragments <- stringr::str_split_fixed(orig_names[i], pattern = "\\[", n = 2)
      new_names[i] <- paste0(fragments[1, 1], "[1, ", fragments[1, 2])
    }
  }
  
  for(c in 1:length(IDSM.out.tidy)){
    dimnames(IDSM.out.tidy[[c]])[2] <- list(new_names)
  }
  
  return(IDSM.out.tidy)
}
