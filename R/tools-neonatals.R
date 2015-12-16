#' Neonatal data: separate Mom and Baby into two separate data frames
#'
#'@param df data frame of data
#'@param annotation_columns names of non-numeric columns in data frame
#'@param rename TRUE or FALSE if you want to append _M and _B on cohort-speciific features
#'@export


separate_mvb <- function(df, annotation_columns, rename = TRUE){

  baby_df <- dplyr::filter(df, SampleType == "B")
  mother_df <- dplyr::filter(df, SampleType == "M")

  # test pairs line up
  if (identical(baby_df$Pair, mother_df$Pair)){

    # restrict to numeric columns
    baby_df <- baby_df[, !(colnames(baby_df) %in% annotation_columns)]
    mother_df <- mother_df[, !(colnames(mother_df) %in% annotation_columns)]

    if (rename == TRUE){
      colnames(baby_df) <- paste(colnames(baby_df), "B", sep = "_")
      colnames(mother_df) <- paste(colnames(mother_df), "M", sep = "_")
    }
    df_list <- list(Baby = baby_df, Mother = mother_df)
    return(df_list)

  } else {
    print("Error: ordering doesn't match")
  }

}
