# FDA-library

Tools for analyses on single-cell mass cytometry data, originally developed for 
the human FDA data set.  

## Functions to query the data server

This package contains functions part of an R toolkit to query a server of cytometry data
(note: base functions are all modifications of code originally written by ZBjornson).  
This includes getting a list of experiments, files, populations, gates, and statistics. NOTE: This has been moved to its own package in its own repo (blackbuck-R-toolkit)

## Functions to organize data into useful structures

This includes functions that will compile sets of queried statistics and annotations into
organized data frames that are easy to parse.  These are still being modified based on 
improvements to bulk statistics queries.  

## Analysis code

Contains markdowns chronicling the analysis of the datasets (analyses/analysis_markdowns)

