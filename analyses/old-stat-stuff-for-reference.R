
#' Get set of statistics for chosen features
#'
#' Specify populations and reagents to generate a feature set of each pair and retrieve statistics
#' Annotations are optional
#' @param experimentID the experiment ID as a string
#' @param FCS_fileID the FCS file ID as a string
#' @param access_key an object containing server URL and authentication info: see "authenticate"
#' @param populations character vector of population names
#' @param reagents character vector of reagent names
#' @param statistic_type statistic desired, "mean", "median", "quantile", "eventcount"
#' @param k required for statistic "quantile", number from 0.0 to 1.0
#' @param annotate logical operator to indicate whether annotation00s should be included in returned object
#' @examples
#' access_key <- authenticate("my_username", "my_password",
#'                            baseURL =  "http://52.27.144.218/api/v1")
#' experiments <- get_experiments(access_key)
#' experimentID <- experiments[experiments$name == "experiment_name", '_id']
#' FCS_fileID <- get_FCS_files(experimentID, access_key)[1, "_id"]
#' populations <- c("CD4+T cells", "CD8+T cells", "CD14+ Monocytes")
#' reagents <- c("pSTAT1", "pSTAT3", "pMAPKAPK2")
#' get_statistic_set(experimentID, FCS_files, access_key, populations, reagents,
#'                    statistic_type = "median", k = NULL, annotate = TRUE)
#' @export
#'
get_statistics_set <- function(experimentID, FCS_files, access_key, populations, reagents,
                               statistic_type, k = NULL, annotate = FALSE){
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

      stat <- get_statistic(experimentID, FCS_fileID, access_key, channel_name, statistic_type,
                            k, populationID)
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

#' Get set of statistics for chosen features
#'
#' Specify populations and reagents to generate a feature set of each pair and retrieve statistics
#' Annotations are optional
#' @param experimentID the experiment ID as a string
#' @param FCS_fileID the FCS file ID as a string
#' @param access_key an object containing server URL and authentication info: see "authenticate"
#' @param populations character vector of population names
#' @param reagents character vector of reagent names
#' @param statistic_type statistic desired, "mean", "median", "quantile", "eventcount"
#' @param k required for statistic "quantile", number from 0.0 to 1.0
#' @param annotate logical operator to indicate whether annotation00s should be included in returned object
#' @examples
#' access_key <- authenticate("my_username", "my_password",
#'                            baseURL =  "http://52.27.144.218/api/v1")
#' experiments <- get_experiments(access_key)
#' experimentID <- experiments[experiments$name == "experiment_name", '_id']
#' FCS_fileID <- get_FCS_files(experimentID, access_key)[1, "_id"]
#' populations <- c("CD4+T cells", "CD8+T cells", "CD14+ Monocytes")
#' reagents <- c("pSTAT1", "pSTAT3", "pMAPKAPK2")
#' get_statistic_set(experimentID, FCS_files, access_key, populations, reagents,
#'                    statistic_type = "median", k = NULL, annotate = TRUE)
#' @export
#'
get_statistics_set_parallel <- function(experimentID, FCS_files, access_key, populations, reagents,
                                        statistic_type, k = NULL, b = batch_size, annotate = FALSE){
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

  k <- 1
  stat_urls <- vector(, length = length(FCS_files$`_id`)*length(feature_set$feature))

  for (X in 1:length(FCS_files$`_id`)){

    for (Y in 1:length(feature_set$feature)){

      FCS_fileID <- FCS_files$`_id`[X]
      feature <- feature_set$feature[Y]
      pop_object <- get_populations(experimentID, access_key)
      pop_name <- feature_set$population[feature_set$feature == feature][1]
      populationID <- pop_object$`_id`[pop_object$name == pop_name]

      channel_name <- feature_set$channel[feature_set$feature == feature][1]

      stat_url <- get_statistic_url(experimentID, FCS_fileID, access_key, channel_name, statistic_type,
                                    k, populationID)
      stat_urls[k] <- stat_url

      k <- k + 1

    }
  }

  #d <- 1:length(stat_urls)
  #ind_list <- split(d, ceiling(seq_along(d)/b))
  #r <- c()
  #for (i in 1:length(ind_list)){
  #  ind <- ind_list[[i]]
  #  r0 <- getURL(url = stat_urls[ind], .opts = access_key$opts)
  #  r <- c(r, r0)
  #}
  r <- getURIAsynchronous(url = stat_urls, .opts = access_key$opts)
  r2 <- lapply(r, fromJSON)
  stats_vec <- sapply(r2, function(x) x[[1]])

  stats_frame <- matrix(stats_vec, byrow = TRUE, ncol = length(feature_set$feature))
  rownames(stats_frame) <- FCS_files$`_id`
  colnames(stats_frame) <- convert_names(feature_set$feature)
  stats_frame <- as.data.frame(stats_frame)


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
    return(list(statistics = stats_frame, feature_object = feature_set, stat_urls = stat_urls))
  }
}
