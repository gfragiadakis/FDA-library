#' Cluster rows and columns of a matrix
#' @param data a data matrix
#' @export
#'
cluster_matrix <- function(data){
  rd<-dist(data)
  rc<-hclust(rd)
  cd<-dist(t(data))
  cc<-hclust(cd)
  ordered_data <- data[rc$order, cc$order]
  cluster_object <- list(data = ordered_data, rc = rc, cc =cc)
  return(cluster_object)
}

#' Make an adjacency matrix
#'
#' @param correlation_matrix a matrix of correlated values
#' @param cor_threshold the significance cutoff for the adjacency matrix
#' @export

make_adjacency <- function(correlation_matrix, cor_threshold = 0.5){
  adj <- correlation_matrix
  adj[adj >= cor_threshold] <- 1
  adj[adj < -cor_threshold] <- -1
  adj[(adj > -cor_threshold) & (adj < cor_threshold)] <- 0
  return(adj)
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

