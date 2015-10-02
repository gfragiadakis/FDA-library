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








