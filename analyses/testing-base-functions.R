access_key <- authenticate("gfragiadakis", "password", baseURL =  "http://52.27.144.218/api/v1")

experiments <- get_experiments(access_key)

experimentID <- experiments[experiments$name == "All Species", '_id']

# get_experiment(experimentID, access_key)

FCS_files <- get_FCS_files(experimentID, access_key)
FCS_files <- FCS_files[1:10, ]

display_parameters(FCS_files, experimentID, access_key)

populations <- c("CD4+T cells", "CD8+T cells", "CD14+ Monocytes")
reagents <- c("pSTAT1", "pSTAT3", "pMAPKAPK2")




FCS_files$annotations[1]


populationID <- get_populations(experimentID, access_key)[1, '_id']
scaleID = get_scale_sets(experimentID, access_key)[1, '_id']
