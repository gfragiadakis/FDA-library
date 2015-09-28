#' Authenticate acess
#'
#' This function establishes the connection to your account on the remote server
#' It returns an object containing the server URL and your connection info: opts
#' @param username The username for the account
#' @param password The password for the account
#' @param baseURL The URL for the server
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
#' @param baseURL The URL for the server
#' @param access_key An object containing server URL and authentication info: see "authenticate"
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
#' @param access_key An object containing server URL and authentication info: see "authenticate"
#' @export

get_experiment <- function(experimentID, access_key) {
  opts <- access_key$opts
  baseURL <- access_key$baseURL
  url <- paste(baseURL, "experiments", experimentID, sep = "/")
  return(fromJSON(getURL(url, .opts = opts)))
}

#' Get all FCS files for an experiment
#'
#' Returns an object containing all FCS files from the experiments based on the experiment ID
#' @param experimentID the experiment ID as a string
#' @param access_key An object containing server URL and authentication info: see "authenticate"
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
#' @param local_file_path path where to download the data including file name
#' @param access_key An object containing server URL and authentication info: see "authenticate"
#' @export

download_fcs_file <- function(experimentID, FCS_fileID, local_file_path, access_key) {
  opts <- access_key$opts
  baseURL <- access_key$baseURL
  f <- CFILE(local_file_path, mode="wb")
  url <- paste(paste(baseURL, "experiments", experimentID, "fcsfiles", FCS_fileID, sep = "/"), ".fcs", sep="")
  a <- curlPerform(url = url, .opts = opts, writedata = f@ref, noprogress=FALSE)
  close(f)
  return(a)
}

#' Get all gates for an experiment
#' @param experimentID the experiment ID as a string
#' @param access_key An object containing server URL and authentication info: see "authenticate"
#' @export

get_gates <- function(experimentID, access_key) {
  opts <- access_key$opts
  baseURL <- access_key$baseURL
  url <- paste(baseURL, "experiments", experimentID, "gates", sep = "/")
  return(fromJSON(getURL(url, .opts = opts)))
}

#' Get all populations for an experiment
#' @param experimentID the experiment ID as a string
#' @param access_key An object containing server URL and authentication info: see "authenticate"
#' @export

get_populations <- function(experimentID, acess_key) {
  opts <- access_key$opts
  baseURL <- access_key$baseURL
  url <- paste(baseURL, "experiments", experimentID, "populations", sep = "/")
  return(fromJSON(getURL(url, .opts = opts)))
}

#' Get all scale sets for an experiment
#' @param experimentID the experiment ID as a string
#' @param access_key An object containing server URL and authentication info: see "authenticate"
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
#' @param k required for statistic "quantile", number from 0.0 to 1.0
#' @param access_key An object containing server URL and authentication info: see "authenticate"
#' @export

get_statistic <- function(experimentID, FCS_fileID, access_key, channel_name, statistic_type, k=NULL, populationID=NULL) {
  opts <- access_key$opts
  baseURL <- access_key$baseURL
  url = paste(
    paste(baseURL, "experiments", experimentID, "statistics", sep="/"),
    paste(
      paste("fcsFileId", FCS_fileID, sep="="),
      paste("channel", channel_name, sep="="),
      paste("statistic", statistic_type, sep="="),
      paste("k", k, sep="="),
      paste("populationId", populationID, sep="="), sep="&"), sep="?")
  print(url)
  return(fromJSON(getURL(url, .opts = opts)))
}

#' Retrieves events
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
