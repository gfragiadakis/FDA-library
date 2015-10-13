access_key <- authenticate("gfragiadakis", "password", baseURL =  "http://52.27.144.218/api/v1")

experiments <- get_experiments(access_key)

experimentID <- experiments[experiments$name == "All Species", '_id']

# get_experiment(experimentID, access_key)

FCS_files <- get_FCS_files(experimentID, access_key)
FCS_files <- FCS_files[1:10, ]

display_parameters(FCS_files, experimentID, access_key)

populations <- c("CD4+T cells", "CD8+T cells", "CD14+ Monocytes")
reagents <- c("pSTAT1", "pSTAT3", "pMAPKAPK2")

results <- get_statistics_set(experimentID, FCS_files, access_key, populations, reagents,
                              statistic_type = "median", k = NULL, annotate = TRUE)

results_parallel <- get_statistics_set_parallel(experimentID, FCS_files, access_key, populations, reagents,
                                    statistic_type = "median", k = NULL, annotate = TRUE)

# Run this for trials

time_original <- system.time((results = get_statistics_set(experimentID, FCS_files, access_key, populations, reagents,
                              statistic_type = "median", k = NULL, annotate = TRUE)))

time_parallel_14 <- system.time((results_parallel_14 = get_statistics_set_parallel(experimentID, FCS_files, access_key, populations, reagents,
                                                                             statistic_type = "median", k = NULL, b= 14, annotate = TRUE)))


time_parallel_2 <- system.time((get_statistics_set_parallel(experimentID, FCS_files, access_key, populations, reagents,
                                                                                   statistic_type = "median", k = NULL, b = 2, annotate = TRUE)))

time_parallel_4 <- system.time((get_statistics_set_parallel(experimentID, FCS_files, access_key, populations, reagents,
                                                            statistic_type = "median", k = NULL, b = 4, annotate = TRUE)))

time_parallel_6 <- system.time((get_statistics_set_parallel(experimentID, FCS_files, access_key, populations, reagents,
                                                            statistic_type = "median", k = NULL, b = 6, annotate = TRUE)))



