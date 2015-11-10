# Final tests for R toolkit for blackbuck

access_key <- authenticate("gfragiadakis", "password", baseURL =  "http://54.186.1.226/api/v1")

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

