library(dplyr)
df <- results$statistics

# remove file ID, file name and donor

df <- dplyr::select(df, -file_ID, -filename, -Donor)

# choose one species and remove species column

df<- df %>% dplyr::filter(Species == "Cyno") %>% dplyr::select(-Species)

df <- dplyr::group_by(df, Condition, Gender)

test <- summarise_each(df, funs(mean(., na.rm = TRUE)))


# things I want to make with the data:
# fold change over baseline
# mapping of cell types or phosphos when I use melt (so feature can correspond)
# be able to plot summary statistics in a histogram
# filter? 3000 data points per donor
# need to see the distributions of the data
# think I need to start by filtering, once I have the fold change data

# to introduce your mapping (i.e. have cell types and phospho correspondence, for plotting)

melted_df <- melt(df, variable.name = "feature_names", value.name = "counts")
merge(select(fdf, -channel), melted_df, by.x = "feature", by.y = "feature_names")

## Here's the idea for the asynchronous:
# write a get_statistic_url function (just gets url, doesn't do the retrieval part)
# make a get_statistic_from_url function (takes url, returns statistic)
# in loop, get all the url's, then give them as a vector to the asynchronous funcion from RCurl
# put it into matrix form

# trying it out
FCS_fileID <- FCS_files$`_id`[1]
FCS_fileID_2 <- FCS_files$`_id`[2]
urls[1]<- get_statistic_url(experimentID, FCS_fileID, access_key, channel_name = "La139Di",
                            statistic_type = "median", k=NULL, populationID=NULL)
urls[2]<- get_statistic_url(experimentID, FCS_fileID_2, access_key, channel_name = "La139Di",
                            statistic_type = "median", k=NULL, populationID=NULL)
test <- getURIAsynchronous(url = urls, .opts = access_key$opts)
test2<- lapply(test, fromJSON)
sapply(test2, function(x) x[[1]])






