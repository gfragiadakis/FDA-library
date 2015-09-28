#' Get annotations for a set of FCS files
#' @param FCS_files Object containing FCS files
#' @export



get_annotations <- function(FCS_files){

  annotations <- t(sapply(FCS_files$annotations, "[[", "value"))
  colnames(annotations) <- FCS_files$annotations[[1]]$type

  data.frame(file_ID = FCS_files$`_id`, filename = FCS_files$filename, annotations)
}

#' Display populations and reagents
#' Provides a list of populations and reagents that can be used to select which are desired features
#' for subsequent analysis
#' @param FCS_files Object containing FCS files
#' @param experimentID the experiment ID as a string
#' @param access_key An object containing server URL and authentication info: see "authenticate"
#' @export

display_parameters <- function(FCS_files, experimentID, access_key){

  populations <- get_populations(experimentID, access_key)$name

  get_reagents <- function(FCS_files){
    FCS_files$panel[[1]]$reagent
  }
  reagents <- get_reagents(FCS_files)

  parameters <- list(populations, reagents)
  names(parameters) <- c("populations", "reagents")
  return(parameters)
}



