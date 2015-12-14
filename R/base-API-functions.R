#' Authenticate acess
#'
#' This function establishes the connection to your account on the remote server
#' It returns an object containing the server URL and your connection information
#' @param username username for the account
#' @param password password for the account
#' @param baseURL URL for the server
#' @examples
#' authenticate("my_username", "my_password",
#'              baseURL =  "https://www.immuneatlas.org/api/v1")
#' @import jsonlite
#' @import RCurl
#' @export
authenticate <- function(username, password, baseURL) {
  userinfo <- toJSON(list(username=username, password=password), auto_unbox = TRUE)
  header <- c(Accept="application/json; charset=UTF-8","Content-Type"="application/json")
  signinURL <- paste(baseURL, "signin", sep = "/")
  result <- postForm(signinURL,.opts=list(httpheader=header, postfields=userinfo))
  jwt <- fromJSON(result)$accessToken
  authToken <- paste("Authorization: Bearer", jwt)
  opts <- list(httpheader = c(authToken))
  return(list(baseURL = baseURL, opts = opts, jwt = jwt))
}

#' Retrieve all accessible experiments
#'
#' This function retrieves an object with the information about all acessible experiments
#' @param baseURL URL for the server
#' @param access_key an object containing server URL and authentication info: see "authenticate"
#' @examples
#' access_key <- authenticate("my_username", "my_password",
#'                            baseURL =  "http://54.186.1.226/api/v1")
#' experiments <- get_experiments(access_key)
#' @import jsonlite
#' @import RCurl
#' @export

get_experiments <- function(access_key) {
  opts <- access_key$opts
  baseURL <- access_key$baseURL
  url <- paste(baseURL, "experiments", sep = "/")
  return(fromJSON(getURL(url, .opts = opts)))
}

#' Retrieve one experiment
#'
#' This function retrieves an object with the information from one experiment specified by the experiment ID
#' @param experimentID the experiment ID as a string
#' @param access_key an object containing server URL and authentication info: see "authenticate"
#' @examples
#' access_key <- authenticate("my_username", "my_password",
#'                            baseURL =  "http://54.186.1.226/api/v1")
#' experiments <- get_experiments(access_key)
#' experimentID <- experiments[experiments$name == "experiment_name", '_id']
#' get_experiment(experimentID, access_key)
#' @export

get_experiment <- function(experimentID, access_key) {
  opts <- access_key$opts
  baseURL <- access_key$baseURL
  url <- paste(baseURL, "experiments", experimentID, sep = "/")
  return(fromJSON(getURL(url, .opts = opts)))
}

#' Get all FCS files for an experiment
#'
#' Returns a data frame containing all FCS files from an experiment based on the experiment ID
#' Lists information about event count, panel, and annotations
#' @param experimentID the experiment ID as a string
#' @param access_key an object containing server URL and authentication info: see "authenticate"
#' @examples
#' access_key <- authenticate("my_username", "my_password",
#'                            baseURL =  "http://54.186.1.226/api/v1")
#' experiments <- get_experiments(access_key)
#' experimentID <- experiments[experiments$name == "experiment_name", '_id']
#' get_FCS_files(experimentID, access_key)
#' @export

get_FCS_files <- function(experimentID, access_key) {
  opts <- access_key$opts
  baseURL <- access_key$baseURL
  url <- paste(baseURL, "experiments", experimentID, "fcsfiles", sep = "/")
  return(fromJSON(getURL(url, .opts = opts)))
}

#' Download an FCS file
#' @param experimentID the experiment ID as a string
#' @param FCS_fileID the FCS file ID as a string
#' @param local_file_path path to local directory for download
#' @param filename name of the file when downloaded (.fcs), may match the original or be newly defined in download
#' @param access_key an object containing server URL and authentication info: see "authenticate"
#' @examples
#' access_key <- authenticate("my_username", "my_password",
#'                            baseURL =  "http://54.186.1.226/api/v1")
#' experiments <- get_experiments(access_key)
#' experimentID <- experiments[experiments$name == "experiment_name", '_id']
#' FCS_file_list <- get_FCS_files(experimentID, access_key)
#' FCS_fileID <- FCS_file_list$`_id`[1]
#' download_fcs_file <- function(experimentID, FCS_fileID, local_file_path = "/Desktop/my_directory/", filename = "my_fcs_file.fcs", access_key)
#' filename <- FCS_file_list[FCS_file_list$`_id`==FCS_fileID, "filename"]
#' download_fcs_file <- function(experimentID, FCS_fileID, local_file_path = "/Desktop/my_directory/", filename = filename, access_key)
#' @export

download_fcs_file <- function(experimentID, FCS_fileID, local_file_path, filename, access_key) {
  opts <- access_key$opts
  baseURL <- access_key$baseURL
  full_file <- paste(local_file_path, filename, sep = "")
  f <- CFILE(full_file, mode="wb")
  url <- paste(paste(baseURL, "experiments", experimentID, "fcsfiles", FCS_fileID, sep = "/"), ".fcs", sep="")
  a <- curlPerform(url = url, .opts = opts, writedata = f@ref, noprogress=FALSE)
  close(f)
  return(a)
}

#' Get all gates for an experiment
#'
#' Returns a data frame with the information on all gates from an experiment including names, dimensions, tailored, etc.
#' @param experimentID the experiment ID as a string
#' @param access_key an object containing server URL and authentication info: see "authenticate"
#' @export

get_gates <- function(experimentID, access_key) {
  opts <- access_key$opts
  baseURL <- access_key$baseURL
  url <- paste(baseURL, "experiments", experimentID, "gates", sep = "/")
  return(fromJSON(getURL(url, .opts = opts)))
}

#' Get all populations for an experiment
#'
#' Returns a data frame that contains population IDs, names, and which gates define it
#' @param experimentID the experiment ID as a string
#' @param access_key an object containing server URL and authentication info: see "authenticate"
#' @examples
#' access_key <- authenticate("my_username", "my_password",
#'                            baseURL =  "http://52.27.144.218/api/v1")
#' experiments <- get_experiments(access_key)
#' experimentID <- experiments[experiments$name == "experiment_name", '_id']
#' get_populations(experimentID, access_key)
#' @export

get_populations <- function(experimentID, acess_key) {
  opts <- access_key$opts
  baseURL <- access_key$baseURL
  url <- paste(baseURL, "experiments", experimentID, "populations", sep = "/")
  return(fromJSON(getURL(url, .opts = opts)))
}

#' Get all scale sets for an experiment
#'
#' @param experimentID the experiment ID as a string
#' @param access_key an object containing server URL and authentication info: see "authenticate"
#' @examples
#' access_key <- authenticate("my_username", "my_password",
#'                            baseURL =  "http://52.27.144.218/api/v1")
#' experiments <- get_experiments(access_key)
#' experimentID <- experiments[experiments$name == "experiment_name", '_id']
#' get_scale_sets(experimentID, access_key)
#' @export

get_scale_sets <- function(experimentID, access_key) {
  opts <- access_key$opts
  baseURL <- access_key$baseURL
  url <- paste(baseURL, "experiments", experimentID, "scalesets", sep = "/")
  return(fromJSON(getURL(url, .opts = opts)))
}

#' Plot the data
#'
#' Creates a plot of specified channels and populations
#' This is in beta mode and needs to be refined for nicer output
#' @param experimentID the experiment ID as a string
#' @param FCS_fileID the FCS file ID as a string
#' @param access_key An object containing server URL and authentication info: see "authenticate"
#' @param x_channel_name The channel for the x axis named by metal
#' @param y_channel_name The channel for the y axis named by metal
#' @param populationID the population ID as a string
#' @export

get_plot <- function(experimentID, FCS_fileID, access_key, x_channel_name, y_channel_name, scale_setID, populationID = NULL, plot_type = "density",
                   ticks = "true", tick_labels = "false", axis_labels = "true") {
  opts <- access_key$opts
  baseURL <- access_key$baseURL
  url <- paste(
    paste(baseURL, "experiments", experimentID, "plot", sep="/"),
    paste(
      paste("fcsFileId", FCS_fileID, sep="="),
      paste("scaleSetId", scale_setID, sep="="),
      paste("xChannel", x_channel_name, sep="="),
      paste("yChannel", y_channel_name, sep="="),
      paste("plotType", plot_type, sep="="),
      paste("ticksQ", ticks, sep="="),
      paste("tickLabelsQ", tick_labels, sep="="),
      paste("axisLabelsQ", axis_labels, sep="="),
      paste("populationId", populationID, sep="="), sep="&"), sep="?")
  image <- png::readPNG(getURLContent(url, .opts = opts))
  # Windows, at least: can't get R to accept units="px". 96 is an example pixels-per-inch value.
  dev.new(width=228/96, height=228/96, units="in")
  par(mar = c(0,0,0,0))
  plot(0:1, 0:1, type="n", ann=FALSE, axes=FALSE, asp=1)
  rasterImage(image, 0, 0, 1, 1, interpolate=FALSE)
  dev.copy(pdf, "test.pdf")
  dev.off()
}

#' Retrieve statistic from a population
#'
#' Returns a desired statistic from a desired feature
#' @param experimentID the experiment ID as a string
#' @param FCS_fileID the FCS file ID as a string
#' @param access_key An object containing server URL and authentication info: see "authenticate"
#' @param channel_name Channel name
#' @param statistic_type Statistic desired, "mean", "median", "quantile", "eventcount
#' @param q required for statistic "quantile", number from 0.0 to 1.0
#' @param populationID Population ID of population of interest
#' @examples
#' access_key <- authenticate("my_username", "my_password",
#'                            baseURL =  "http://52.27.144.218/api/v1")
#' experiments <- get_experiments(access_key)
#' experimentID <- experiments[experiments$name == "experiment_name", '_id']
#' FCS_fileID <- get_FCS_files(experimentID, access_key)[1, "_id"]
#' get_statistic(experimentID, FCS_fileID, access_key, channel_name = "La139Di",
#'                statistic_type = "median", k=NULL, populationID=NULL)
#' @export

get_statistic <- function(experimentID, FCS_fileID, access_key, channel_name, statistic_type, q=NULL, populationID=NULL) {
  opts <- access_key$opts
  baseURL <- access_key$baseURL
  url <- paste(
    paste(baseURL, "experiments", experimentID, "statistics", sep="/"),
    paste(
      paste("fcsFileId", FCS_fileID, sep="="),
      paste("channel", channel_name, sep="="),
      paste("statistic", statistic_type, sep="="),
      paste("q", q, sep="="),
      paste("populationId", populationID, sep="="), sep="&"), sep="?")
  return(fromJSON(getURL(url, .opts = opts)))
}

#' Retrieve statistic URL from a population
#'
#' Returns a desired statistic from a desired feature
#' @param experimentID the experiment ID as a string
#' @param FCS_fileID the FCS file ID as a string
#' @param access_key An object containing server URL and authentication info: see "authenticate"
#' @param channel_name Channel name
#' @param statistic_type Statistic desired, "mean", "median", "quantile", "eventcount
#' @param k required for statistic "quantile", number from 0.0 to 1.0
#' @param populationID Population ID of population of interest
#' @examples
#' access_key <- authenticate("my_username", "my_password",
#'                            baseURL =  "http://52.27.144.218/api/v1")
#' experiments <- get_experiments(access_key)
#' experimentID <- experiments[experiments$name == "experiment_name", '_id']
#' FCS_fileID <- get_FCS_files(experimentID, access_key)[1, "_id"]
#' get_statistic_url(experimentID, FCS_fileID, access_key, channel_name = "La139Di",
#'                statistic_type = "median", k=NULL, populationID=NULL)
#' @export

get_statistic_url <- function(experimentID, FCS_fileID, access_key, channel_name, statistic_type, k=NULL, populationID=NULL) {
  opts <- access_key$opts
  baseURL <- access_key$baseURL
  url <- paste(
    paste(baseURL, "experiments", experimentID, "statistics", sep="/"),
    paste(
      paste("fcsFileId", FCS_fileID, sep="="),
      paste("channel", channel_name, sep="="),
      paste("statistic", statistic_type, sep="="),
      paste("k", k, sep="="),
      paste("populationId", populationID, sep="="), sep="&"), sep="?")
  return(url)
}

#' Retrieve events from a specified population
#'
#' Returns a data frame containing the cell events from a given experiment from a specified population
#'
#' @param experimentID the experiment ID as a string
#' @param FCS_fileID the FCS file ID as a string
#' @param access_key An object containing server URL and authentication info: see "authenticate"
#' @param populationID the population ID as a string
#' @export

get_events <- function(experimentID, FCS_fileID, access_key, populationID=NULL) {
  opts <- access_key$opts
  baseURL <- access_key$baseURL
  jwt <- access_key$jwt
  url <- paste(
    paste(paste(baseURL, "experiments", experimentID, "fcsfiles", FCS_fileID, sep="/"), ".tsv", sep=""),
    paste(
      paste("populationId", populationID, sep="="),
      paste("token", jwt, sep="="), sep="&"), sep="?")
  return(read.delim(url))
}
