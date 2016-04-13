# Read in frame of features

compiled_df <- read.csv("~/Documents/FDAlibrary/saved_data_structures/human_signaling_fold_asinh_0.2.csv", row.names = 1)
df <- compiled_df

# correlate all the features and then look at the distribution of the correlations

df <- dplyr::select(df, Donor, Condition_Feature, value)
df <- reshape2::dcast(df, Donor ~ Condition_Feature)
df <- data.frame(df[,-1], row.names = df$Donor)

correlations <- cor(df)
sorted <- sort(correlations)
hist(sorted)


plot_correlations(df, cor_threshold = 0.5, output_directory = "~/Documents/FDAlibrary/analyses/plots/thresholded_asinh_0.2/correlation/")
