library(FDAlibrary)
library(dplyr)
library(reshape2)

# Notes:
# 7605 donated twice: R20W3_7605 and R19W3_7605.
# 7962 donated twice: R20W11_7962 and R25W5_7962
# R23 has incorrect metadata

# set data for structure output
dataDirectory <- "~/Documents/FDAlibrary/saved_data_structures/"

# Read in data from External query

human_table <- read.table("~/Documents/FDAlibrary/saved_data_structures/signaling-for-gabi-v2.3.tsv", sep="\t")
colnames(human_table) <- c("filename", "Population", "Reagent", "Raw_median")
human_df <- data.frame(human_table)
human_df$Raw_median <- as.numeric(as.character(human_df$Raw_median))

# Convert names and store each
population_key <- data.frame(Original = levels(human_df$Population),
                             Converted =convert_names(levels(human_df$Population)))

reagent_key <- data.frame(Original = levels(human_df$Reagent),
                          Converted =convert_names(levels(human_df$Reagent)))

# Make features (Population x Reagent)
human_df$Population <- as.factor(convert_names(human_df$Population))
human_df$Reagent <- as.factor(convert_names(human_df$Reagent))

# Filter out extra populations

new_df <- human_df %>%
  dplyr::filter(Population != "CD4..CD8.") %>%
  dplyr::filter(Population != "HLA.DR..CD11c.") %>%
  dplyr::filter(Population != "singlets_basophils.") %>%
  dplyr::filter(Population != "pDCs") %>%
  dplyr::mutate(Feature = paste(Population, Reagent, sep = "_"))

new_df$Feature <- as.factor(new_df$Feature)
feature_matrix <- dplyr::select(new_df, Population, Reagent, Feature)

new_df <- dplyr::select(new_df, filename, Feature, Raw_median)


# Change to wide format

wide_df <- dcast(new_df, filename ~ Feature, value.var = "Raw_median")


# Get the annotations from the server
access_key <- authenticate("gfragiadakis", "password", baseURL =  "http://54.186.1.226/api/v1")

experiments <- get_experiments(access_key)

experimentID <- experiments[experiments$name == "All Species", '_id']
FCS_files <- get_FCS_files(experimentID, access_key)
annotations <- get_annotations(FCS_files)
human_annotations <- dplyr::filter(annotations,Species == "Human")

annotated_df <- merge(human_annotations, wide_df)

# 7605 donated twice: R20W3_7605 and R19W3_7605.
# 7962 donated twice: R20W11_7962 and R25W5_7962
# make donors unique

levels(annotated_df$Donor) <- c(levels(annotated_df$Donor), "7605_2", "7962_2")
annotated_df[grep("R20W3", annotated_df$filename), "Donor"] <- "7605_2"
annotated_df[grep("R25W5", annotated_df$filename), "Donor"] <- "7962_2"

write.csv(annotated_df, file = paste(dataDirectory, "human_signaling_raw.csv", sep = ""))

fold_df <- get_folds(annotated_df, basal_name = "Basal", fold_type = "raw")
asinh_df <- get_folds(annotated_df, basal_name = "Basal", fold_type = "asinh")

write.csv(fold_df, file = paste(dataDirectory, "human_signaling_fold.csv", sep = ""))
write.csv(asinh_df,file = paste(dataDirectory, "human_signaling_fold_asinh.csv", sep = ""))

compiled_df <- threshold_features(data = asinh_df, threshold = 0.2)
write.csv(compiled_df, file = paste(dataDirectory, "human_signaling_fold_asinh_0.2.csv", sep = ""))
