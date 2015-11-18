# These are functions pertaining to creating the building block structures of the data
# including annotations, features, and queried statistics

#' Get annotations for a set of FCS files
#'
#' If some files are not annotated, it will return a warning and only return the fcs files and annotations only
#' @param FCS_files object containing FCS files from get_FCS_files()
#' @export


get_annotations <- function(FCS_files){

  if (sum(sapply(FCS_files$annotations, nrow) == 0) != 0){
    FCS_files <- FCS_files[sapply(FCS_files$annotations, nrow) != 0, ]
    print("Warning: missing annotations in experiment, returning truncated file frame")
  } else {}

  annotations <- t(sapply(FCS_files$annotations, "[[", "value"))
  colnames(annotations) <- FCS_files$annotations[[1]]$type

  return(data.frame(file_ID = FCS_files$`_id`, filename = FCS_files$filename, annotations))
}

#' Display populations and reagents
#'
#' Provides a list of population names, channel names, and reagent names that can be used to select which are desired features
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

  get_channels <- function(FCS_files){
    FCS_files$panel[[1]]$channel
  }
  channels <- get_channels(FCS_files)

  parameters <- list(populations, channels, reagents)
  names(parameters) <- c("populations", "channels", "reagents")
  return(parameters)
}

#' Get channel names from reagent names
#'
#' Returns a vector of channel names that corresponds to the reagent names from an FCS_files set (see get FCS_files)
#' @param reagent_names a vector of reagent names
#' @param FCS_files a set of FCS files from get_FCS_files()
#' @examples
#' reagent_names <- c("pSTAT1", "pSTAT3", "pMAPKAPK2")
#' FCS_files <- get_FCS_files(experimentID, access_key)
#' channel_names <- convert_reagents_to_channels(reagent_names, FCS_files)
#' @export
#'

convert_reagents_to_channels <- function(reagent_names, FCS_files){

  df <- FCS_files$panel[[1]]
  channel_names <- df[df$reagent %in% reagent_names, "channel"]
  return(channel_names)

}

#' Convert names
#'
#' Convert names for compatibility in R, removing promblematic characters such as "+" of spaces
#' @param x A vector of names
#' @export

convert_names <- function(x){
  new_names <- make.names(x)
}

#' Make a bulk statistics request
#'
#' Makes a bulk statistics request for a set of FCS files, populations, channels, and statistics types
#' Note limitations may exist on the number of queries in a bulk request (queries = # files x # populations x # channels)
#' @param experimentID the experiment ID as a string
#' @param FCS_fileIDs vector of FCS file IDs as strings
#' @param access_key An object containing server URL and authentication info: see "authenticate"
#' @param channel_names A vector of channel names as strings (note: not reagent names; for conversion of reagents to channels, see function convert_reagents_to_channels())
#' @param statistic_types A vector of statistic types returned including "mean", "median", "quantile", "eventcount"
#' @param q specified quantile (numeric)
#' @param populationIDs vector of population IDs as strings
#' @export

get_bulk_statistics <- function(experimentID, FCS_fileIDs, access_key, channel_names, statistic_types, q=NULL, populationIDs = NULL) {
  opts <- access_key$opts
  baseURL <- access_key$baseURL

  url <- paste(
    paste(baseURL, "experiments", experimentID, "bulkstatistics", sep="/"),
    paste(
      paste("fcsFileIds", toJSON(FCS_fileIDs), sep="="),
      paste("channels", toJSON(channel_names), sep="="),
      paste("statistics", toJSON(statistic_types), sep="="),
      paste("q", q, sep="="),
      paste("populationIds", toJSON(populationIDs), sep="="), sep="&"), sep="?")

  return(fromJSON(getURL(url, .opts = opts)))
}

#' Get set of statistics using bulk requests
#'
#' @param experimentID the experiment ID as a string
#' @param FCS_files FCS files to take statistics from that contains IDs; can be either a data frame (e.g. output from get_FCS_files) or a vector of IDs as strings
#' @param access_key An object containing server URL and authentication info: see "authenticate"
#' @param channel_names A vector of channel names as strings (can specify either reagents or channels; leave the other as NULL)
#' @param reagent_names A vector of reagent names as strings (can specify either reagents or channels; leave the other as NULL; if using reagent_names, FCS_files must be an output of get_FCS_files())
#' @param statistic_types A vector of statistic types returned including "mean", "median", "quantile", "eventcount"
#' @param q specified quantile (numeric)
#' @param population_names vector of population names as strings
#' @param query_limit the number of queries the server can handle at once (if exceeds query limit, it will batch requests)
#' @export
#' @examples
#' population_names <- c("CD4+T cells", "CD8+T cells", "CD14+ Monocytes")
#' reagent_names <- c("pSTAT1", "pSTAT3", "pMAPKAPK2")
#' statistics_frame <- get_statistics_set_bulk(experimentID, FCS_files, access_key, reagent_names = reagent_names,
#'                      statistic_types = c("median"), q=NULL, population_names = population_names, query_limit = 500)

get_statistics_set_bulk <- function(experimentID, FCS_files, access_key, channel_names = NULL, reagent_names = NULL,
                                    statistic_types, q=NULL, population_names = NULL, query_limit = 2500){

  # FCS_fileIDs
  if (is.data.frame(FCS_files)){
    FCS_fileIDs <- FCS_files$`_id`
  } else {
    FCS_fileIDs <- FCS_files
  }
  # generate channel_names
  if (length(reagent_names) != 0 & length(channel_names) == 0){
    channel_names <- convert_reagents_to_channels(reagent_names, FCS_files)
  } else if (length(reagent_names) == 0 & length(channel_names) == 0){
    print("Error: Must specify channel_names or reagent_names")
  }

  # generate population IDs
  population_frame <- get_populations(experimentID, acess_key)
  populationIDs <- population_frame[population_frame$name %in% population_names, "_id"]

  # subdivide into queries of appropriate length
  query_size <- length(FCS_fileIDs)*length(populationIDs)*length(channel_names)

  if (query_size <= query_limit){

    stat_frame <- get_bulk_statistics(experimentID, FCS_fileIDs, access_key, channel_names,
                                      statistic_types, q=q, populationIDs = populationIDs)

    } else if (query_size > query_limit){

    FCS_batch_size <- floor(query_limit/(length(populationIDs)*length(channel_names)))
    d <- 1:length(FCS_fileIDs)
    ind_list <- split(d, ceiling(seq_along(d)/FCS_batch_size))

    stat_frame = c()
    for (i in 1:length(ind_list)){
      FCS_fileID_set <- FCS_fileIDs[ind_list[[i]]]
      print(i)
      new_frame <- get_bulk_statistics(experimentID, FCS_fileID_set, access_key, channel_names,
                                        statistic_types, q=NULL, populationIDs = populationIDs)
      stat_frame <- rbind(stat_frame, new_frame)
    }
  }

return(stat_frame)


# reformat everything post query for output; can keep in same format but should subsitute reagents, populations, and filenames
# fuse into features?? (probably somewhere else...)

}

#' Generate fold change statistics object
#'
#' Takes a statistics object and presents statistics relative to designated baseline.
#' Options exist for asinh ratio or fold change on raw counts.
#' If using a statistics object that is an output of a get_statistics function, use statistics_object$statistics
#' Must contain "Condition" column
#' Unique identifier columns such as filename, File ID, etc must be specified as arguments for removal
#'
#'@param statistics_frame data frame of statistics
#'@param basal_name name of the basal condition as string e.g. "Basal"
#'@param fold_type either "asinh", which calculates asinh ratio, or "raw", which calculates the fold change in raw counts
#'@param unique_IDs a vector of strings that are the names of the unique identifier columns (e.g. file_ID or filename)
#'@param annotation_columns vector of strings that are the names of the annotation columns (e.g. "Donor", "Species", "Gender", "Condition")
#'@examples
#' get_folds(statistics_frame = statistics_object$statistics, basal_name = "Basal", fold_type = "asinh", unique_IDs = c("file_ID", "filename"),
#'              annotation_columns = c("Donor","Species", "Gender","Condition"))
#'@export
#'@import reshape2
#'@import dplyr

get_folds <- function(statistics_frame, basal_name, fold_type, unique_IDs = c("file_ID", "filename"),
                      annotation_columns = c("Donor","Species", "Gender","Condition")){

  df <- statistics_frame
  df <- df[, !(colnames(df) %in% unique_IDs)]

  formula <- as.formula(paste(paste(annotation_columns, collapse = " + "), "variable", sep = " ~ "))

  base_df <- df %>%
    dplyr::filter(Condition == basal_name) %>%
    dplyr::rename(Base_Condition = Condition) %>%
    reshape2::melt() %>%
    dplyr::rename(Base_Value = value)


  stim_df <- df %>%
    dplyr::filter(Condition != basal_name) %>%
    dplyr::rename(Stim_Condition = Condition) %>%
    reshape2::melt() %>%
    dplyr::rename(Stim_Value = value)

  if (fold_type == "asinh"){

    base_df$Base_Value <- asinh(base_df$Base_Value/5)
    stim_df$Stim_Value <- asinh(stim_df$Stim_Value/5)

    new_df <- merge(stim_df, base_df) %>%
      dplyr::mutate(fold_value = Stim_Value - Base_Value) %>%
      dplyr::select(-Base_Condition, -Base_Value, -Stim_Value) %>%
      dplyr::rename(Condition = Stim_Condition) %>%
      reshape2::dcast(formula)

  } else if (fold_type == "raw"){

    new_df <- merge(stim_df, base_df) %>%
      dplyr::mutate(fold_value = Stim_Value / Base_Value) %>%
      dplyr::select(-Base_Condition, -Base_Value, -Stim_Value) %>%
      dplyr::rename(Condition = Stim_Condition) %>%
      reshape2::dcast(formula)
  }
  return(new_df)

}

#' Filter features at a given threshold
#'
#' Filters features per condition based on whether the absolute value of the median across
#' samples is greater than the specified threshold
#' @param data data frame of fold features (see get_folds function)
#' @param threshold threshold for absolute value of the median fold feature per condition
#' @export
#' @import dplyr
#' @import reshape2
#' @examples
#' fold_feature_data <- get_folds(statistics_frame = statistics_object$statistics, basal_name = "Basal", fold_type = "asinh", unique_IDs = c("file_ID", "filename"),
#'              annotation_columns = c("Donor","Species", "Gender","Condition"))
#' threshold_features(data = fold_feature_data, threshold = 0.2)

threshold_features <- function(data, threshold){

  compiled_df <- data.frame()

  for (i in levels(data$Condition)){

    if (i != "Basal"){

      df <- data %>%
        dplyr::filter(Condition == i) %>%
        reshape2::melt(variable.name = "Feature")

      df_grouped <- dplyr::group_by(df, Feature)
      df_stats <- dplyr::summarise(df_grouped, mean = mean(value))

      thresholded_features <- df_stats$Feature[abs(df_stats$mean) >= threshold]

      new_df <- dplyr::filter(df, Feature %in% thresholded_features)
      compiled_df <- rbind(compiled_df, new_df)
    }

  }

  compiled_df <- dplyr::mutate(compiled_df, Condition_Feature = paste(Condition, Feature, sep = "_"))
  compiled_df$Condition_Feature <- as.factor(compiled_df$Condition_Feature)
  return(compiled_df)
}
