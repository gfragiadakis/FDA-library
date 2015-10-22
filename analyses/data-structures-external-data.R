library(dplyr)
library(reshape2)
library(magrittr)

# Read in data from External query

human_table <- read.table("~/Documents/FDAlibrary/saved_data_structures/signaling-for-gabi.tsv", sep="\t")
colnames(human_table) <- c("Condition", "Donor", "Population", "Reagent", "Raw_median")
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
  dplyr::mutate(Feature = paste(Population, Reagent, sep = "_"))

new_df$Feature <- as.factor(new_df$Feature)
feature_matrix <- dplyr::select(new_df, Population, Reagent, Feature)

new_df <- dplyr::select(new_df, Donor, Condition, Feature, Raw_median)


# Change to wide format

# data clean up (Messups in the data that will be fixed with new sheets)
new_df <- unique(new_df)

wide_df <- dcast(new_df, Donor + Condition ~ Feature, value.var = "Raw_median")


# Get the annotations from the server

access_key <- authenticate("gfragiadakis", "password", baseURL =  "http://54.191.150.29/api/v1")

experiments <- get_experiments(access_key)

experimentID <- experiments[experiments$name == "All Species", '_id']
FCS_files <- get_FCS_files(experimentID, access_key)
FCS_files <- FCS_files[1:2624,] # should fix annotations function so don't have to do it
annotations <- get_annotations(FCS_files)
human_annotations <- dplyr::filter(annotations,Species == "Human")

annotated_df <- merge(human_annotations, wide_df)
