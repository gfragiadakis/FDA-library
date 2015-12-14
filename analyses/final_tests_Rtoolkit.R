# Final tests for R toolkit for blackbuck

access_key <- authenticate("gfragiadakis", "password", baseURL =  "http://52.33.207.241/api/v1")

experiments <- get_experiments(access_key)

experimentID <- experiments[experiments$name == "All Species", '_id']

experiment <- get_experiment(experimentID, access_key)

FCS_file_list <- get_FCS_files(experimentID, access_key)

download_fcs_file(experimentID, FCS_fileID = FCS_file_list$`_id`[1], local_file_path = "~/Desktop/", filename = "test.fcs", access_key)

gates <- get_gates(experimentID, access_key)

populations <- get_populations(experimentID, access_key)
populationID <- populations$`_id`[1]

scale_sets <- get_scale_sets(experimentID, access_key)

get_statistic(experimentID, FCS_fileID, access_key, channel_name = "La139Di", statistic_type = "median", k=NULL, populationID=NULL)
get_statistic(experimentID, FCS_fileID, access_key, channel_name = "La139Di", statistic_type = "median", k=NULL, populationID=populationID)

get_statistic_url(experimentID, FCS_fileID, access_key, channel_name = "La139Di",  statistic_type = "median", k=NULL, populationID=NULL)


event_table <- get_events(experimentID, FCS_fileID, access_key, populationID=populationID)

annotated_FCS_list <- get_annotations(FCS_file_list)

display_parameters(FCS_file_list, experimentID, access_key)

# Once have bulk statistics,
# test folds
# test thresholds

channel_names <- c("Dy161Di","Dy162Di","Dy163Di")
population_names <- c("pDCs","CD14+ Monocytes")
test <- get_statistics_set_bulk(experimentID, FCS_files, access_key, channel_names = channel_names,
  statistic_types = c("median"), q=NULL, population_names = population_names, query_limit = 900)

### Testing bulk statistics: Run this for testing

access_key <- authenticate("gfragiadakis", "Anatoli4%Riza", baseURL =  "http://52.33.207.241/api/v1")

experiments <- get_experiments(access_key)

experimentID <- experiments[experiments$name == "All Species", '_id']

experiment <- get_experiment(experimentID, access_key)

FCS_file_list <- get_FCS_files(experimentID, access_key)

annotations <- get_annotations(FCS_file_list)

FCS_files <- FCS_file_list[annotations$Species == "Human",]

disp <- display_parameters(FCS_files, experimentID, access_key)

# 10 populations
population_names <- disp$populations[c(2, 6, 9, 10, 14, 16, 18, 19, 20, 21)]

# 16 reagents
reagent_names <- disp$reagents[c(12, 17, 18, 19, 20, 22, 23, 25, 27, 33, 34, 36, 38, 42, 43, 45)]

# 1387 FCS_files; 221,920 requests
query_num <- 10*16*1387

statistics_frame <- get_statistics_set_bulk(experimentID, FCS_files, access_key, reagent_names = reagent_names,
                      statistic_types = c("median"), q=NULL, population_names = population_names, query_limit = 10000)
