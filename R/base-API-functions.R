#' Authenticate acess
#' 
#' This function establishes the connection to your account on the remote server
#' @param username The username for the account
#' @param password The password for the account

authenticate = function(username, password) {
  userinfo = toJSON(list(username=username, password=password), auto_unbox = TRUE)
  header = c(Accept="application/json; charset=UTF-8","Content-Type"="application/json")
  signinURL = paste(baseURL, "signin", sep = "/")
  result = postForm(signinURL,.opts=list(httpheader=header, postfields=userinfo))
  jwt <<- fromJSON(result)$accessToken
  authToken = paste("Authorization: Bearer", jwt)
  opts <<- list(httpheader = c(authToken))
}

#' Retrieve all accessible experiments
#' 
#' This function retrieves an object with the information about all acessible experiments 

getExperiments = function() {
  url = paste(baseURL, "experiments", sep = "/")
  return(fromJSON(getURL(url, .opts = opts)))
}

#' Retrieve one experiment
#' 
#' This function retrieves an object with the information from one experiment specified by the experiment ID
#' @param experimentId the experiment ID as a string

getExperiment = function(experimentId) {
  url = paste(baseURL, "experiments", experimentId, sep = "/")
  return(fromJSON(getURL(url, .opts = opts)))
}

#' Get all FCS files for an experiment
#' 
#' Returns an object containing all FCS files from the experiments based on the experiment ID
#' @param experimentId the experiment ID as a string

getFcsFilesForExperiment = function(experimentId) {
  url = paste(baseURL, "experiments", experimentId, "fcsfiles", sep = "/")
  return(fromJSON(getURL(url, .opts = opts)))
}

#' Download an FCS file
#' @param experimentId the experiment ID as a string
#' @param fcsFileId the FCS file ID as a string
#' @param localFilePath path where to download the data

downloadFcsFile=function(experimentId, fcsFileId, localFilePath) {
  f = CFILE(localFilePath, mode="wb")
  url = paste(paste(baseURL, "experiments", experimentId, "fcsfiles", fcsFileId, sep = "/"), ".fcs", sep="")
  a = curlPerform(url = url, .opts = opts, writedata = f@ref, noprogress=FALSE)
  close(f)
  return(a)
}

#' Get all gates for an experiment
#' @param experimentId the experiment ID as a string

getGatesForExperiment = function(experimentId) {
  url = paste(baseURL, "experiments", experimentId, "gates", sep = "/")
  return(fromJSON(getURL(url, .opts = opts)))
}

#' Get all populations for an experiment
#' @param experimentId the experiment ID as a string

getPopulationsForExperiment = function(experimentId) {
  url = paste(baseURL, "experiments", experimentId, "populations", sep = "/")
  return(fromJSON(getURL(url, .opts = opts)))
}

#' Get all scale sets for an experiment
#' @param experimentId the experiment ID as a string
#' 
getScaleSetsForExperiment = function(experimentId) {
  url = paste(baseURL, "experiments", experimentId, "scalesets", sep = "/")
  return(fromJSON(getURL(url, .opts = opts)))
}

#' Plot the data
#' 
#' Creates a plot of specified channels and populations
#' @param experimentId the experiment ID as a string
#' @param fcsFileId the FCS file ID as a string
#' @param xChannelName The channel for the x axis named by metal
#' @param yChannelName The channel for the y axis named by metal
#' @param populationId the population ID as a string

getPlot = function(experimentId, fcsFileId, xChannelName, yChannelName, scaleSetId, populationId=NULL, plotType="density", ticks="true", tickLabels="false", axisLabels="true") {
  url = paste(
    paste(baseURL, "experiments", experimentId, "plot", sep="/"),
    paste(
      paste("fcsFileId", fcsFileId, sep="="),
      paste("scaleSetId", scaleSetId, sep="="),
      paste("xChannel", xChannelName, sep="="),
      paste("yChannel", yChannelName, sep="="),
      paste("plotType", plotType, sep="="),
      paste("ticksQ", ticks, sep="="),
      paste("tickLabelsQ", tickLabels, sep="="),
      paste("axisLabelsQ", axisLabels, sep="="),
      paste("populationId", populationId, sep="="), sep="&"), sep="?")
  image = readPNG(getURLContent(url, .opts = opts))
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
#' @param experimentId the experiment ID as a string
#' @param fcsFileId the FCS file ID as a string
#' @param channelName Channel name
#' @param statisticType Statistic desired, "mean", "median", "quantile", "eventcount
#' @param k required for statistic "quantile", number from 0.0 to 1.0

getStatistic = function(experimentId, fcsFileId, channelName, statisticType, k=NULL, populationId=NULL) {
  url = paste(
    paste(baseURL, "experiments", experimentId, "statistics", sep="/"),
    paste(
      paste("fcsFileId", fcsFileId, sep="="),
      paste("channel", channelName, sep="="),
      paste("statistic", statisticType, sep="="),
      paste("k", k, sep="="),
      paste("populationId", populationId, sep="="), sep="&"), sep="?")
  return(fromJSON(getURL(url, .opts = opts)))
}

#' Retrieves events
#' 
#' @param experimentId the experiment ID as a string
#' @param fcsFileId the FCS file ID as a string
#' @param populationId the population ID as a string

getEvents = function(experimentId, fcsFileId, populationId=NULL) {
  url = paste(
    paste(paste(baseURL, "experiments", experimentId, "fcsfiles", fcsFileId, sep="/"), ".tsv", sep=""),
    paste(
      paste("populationId", populationId, sep="="),
      paste("token", jwt, sep="="), sep="&"), sep="?")
  return(read.delim(url))
}
