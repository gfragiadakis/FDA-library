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

#' Get set of statistics when you input populations and reagents
#' @export
get_statistics_set <- function(experimentID, FCS_files, access_key, populations, reagents, statistic_type, k = NULL){

  # make feature set data frame (populations x reagents)
  feature_set <- matrix(0, nrow = length(populations) * length(reagents), ncol = 4)
  colnames(feature_set) <- c("population", "reagent", "channel", "feature")
  feature_set <- as.data.frame(feature_set)
  r <- 1
  x <- FCS_files$panel[[1]]
  for (i in populations){
    for (j in reagents){
      feature_set$population[r] <- i
      feature_set$reagent[r] <- j
      feature_set$channel[r] <- x[(x$reagent == j) %in% TRUE, "channel"]
      feature_set$feature[r] <- paste(i, j, sep = "_")
      r <- r + 1
    }
  }

  # outer(X, Y, FUN)
  stat_wrapper <- function(FCS_fileID, feature){

    pop_object <- get_populations(experimentID, access_key)
    pop_name <- feature_set$population[feature_set$feature == feature][1]
    print(pop_object$name)
    print(pop_name)
    populationID <- pop_object$`_id`[pop_object$name == pop_name]

    channel_name <- feature_set$channel[feature_set$feature == feature][1]
    print(channel_name)

    stat <- get_statistic(experimentID, FCS_fileID, access_key, channel_name, statistic_type, k, populationID)
    return(stat)
  }

  X <- FCS_files$`_id`
  Y <- feature_set$feature
  outer(X, Y, FUN = stat_wrapper)

}



