library(dplyr)
library(reshape2)

# Read in data from External query

human_table <- read.table("~/Documents/FDAlibrary/saved_data_structures/signaling-for-gabi.tsv", sep="\t")
colnames(human_table) <- c("Condition", "ID_unformatted", "Population", "Reagent", "Raw_median")
human_df <- data.frame(human_table)

# Convert names and store each
population_key <- cbind(levels(human_df$Population), convert_names(levels(human_df$Population)))
colnames(population_key) <- c("Original", "Converted")

reagent_key <- cbind(levels(human_df$Reagent), convert_names(levels(human_df$Reagent)))
colnames(reagent_key) <- c("Original", "Converted")

# Make features (Population x Reagent)
human_df$Population <- as.factor(convert_names(human_df$Population))
human_df$Reagent <- as.factor(convert_names(human_df$Reagent))
# Filter out those dumb pops

Feature <- paste(human_df$Population, human_df$Reagent, sep = "_")
human_df <- cbind(human_df, Feature)
feature_matrix <- dplyr::select(human_df, Population, Reagent, Feature)
human_df <- dplyr::select(human_df, Condition, ID_unformatted, Feature, Raw_median)

# Change to wide format
wide_human_df <- dcast(human_df, Condition + ID_unformatted ~ Feature)









## Random garbage may use later

access_key <- authenticate("gfragiadakis", "password", baseURL =  "http://52.27.144.218/api/v1")

experiments <- get_experiments(access_key)

experimentID <- experiments[experiments$name == "All Species", '_id']


# get_experiment(experimentID, access_key)

FCS_files <- get_FCS_files(experimentID, access_key)
FCS_files <- FCS_files[1:2624,]
annotations <- get_annotations(FCS_files)
Humans <- dplyr::filter(annotations,Species == "Human")
