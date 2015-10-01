#' Get annotations for a set of FCS files
#'
#' @param FCS_files Object containing FCS files
#' @export


get_annotations <- function(FCS_files){

  annotations <- t(sapply(FCS_files$annotations, "[[", "value"))
  colnames(annotations) <- FCS_files$annotations[[1]]$type

  return(data.frame(file_ID = FCS_files$`_id`, filename = FCS_files$filename, annotations))
}

#' Display populations and reagents
#'
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

#' Convert names
#'
#' Convert names for compatibility in R
#' @param x A vector of names
#' @export

convert_names <- function(x){
  new_names <- make.names(x)
}

#' Get set of statistics for chosen features
#'
#' Specify populations and reagents to generate a feature set of each pair and retrieve statistics
#' Annotations are optional
#' @param experimentID the experiment ID as a string
#' @param FCS_fileID the FCS file ID as a string
#' @param access_key An object containing server URL and authentication info: see "authenticate"
#' @param populations A character vector of population names
#' @param reagents A character vector of reagent names
#' @param statistic_type Statistic desired, "mean", "median", "quantile", "eventcount"
#' @param k required for statistic "quantile", number from 0.0 to 1.0
#' @param annotate Logical operator to indicate whether annotations should be included in returned object
#' @export
#'
get_statistics_set <- function(experimentID, FCS_files, access_key, populations, reagents,
                               statistic_type, k = NULL, annotate = FALSE){
  # TO DO: make colnames of stats_frame match contents in feature set (+ and space issue)
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

  # retrieve statistics and organize into a data frame

  stats_frame <- matrix(0, nrow = length(FCS_files$`_id`), ncol = length(feature_set$feature))
  rownames(stats_frame) <- FCS_files$`_id`
  colnames(stats_frame) <- convert_names(feature_set$feature)
  stats_frame <- as.data.frame(stats_frame)

  for (X in 1:length(FCS_files$`_id`)){

    for (Y in 1:length(feature_set$feature)){

      FCS_fileID <- FCS_files$`_id`[X]
      feature <- feature_set$feature[Y]
      pop_object <- get_populations(experimentID, access_key)
      pop_name <- feature_set$population[feature_set$feature == feature][1]
      populationID <- pop_object$`_id`[pop_object$name == pop_name]

      channel_name <- feature_set$channel[feature_set$feature == feature][1]

      stat <- get_statistic(experimentID, FCS_fileID, access_key, channel_name, statistic_type, k, populationID)
      stats_frame[X, Y] <- stat[[1]]
    }
  }

  # convert names in feature object
  feature_set <- as.data.frame(apply(feature_set, 2, convert_names))

 # inculstion of annotations

  if (annotate == TRUE){
    annotation_frame <- get_annotations(FCS_files)

    if (all(rownames(stats_frame) ==  annotation_frame$file_ID) == TRUE){
      annotated_stats_frame <- data.frame(annotation_frame, stats_frame, row.names = NULL)

      return(list(statistics = annotated_stats_frame, feature_object = feature_set))

    } else {
      print("ERROR: file IDs do not match in annotations and statistics")
    }
  } else {
    return(list(statistics = stats_frame, feature_object = feature_set))
  }
}



