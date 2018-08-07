############################
# author: Roman Shopa
# Roman.Shopa[at]ncbj.gov.pl
############################

# Get system arguments
args <- commandArgs(trailingOnly = FALSE)
# detect directory for the srcipt
# solution from here: http://stackoverflow.com/a/32016824/538403
script_dir  <- dirname(normalizePath(sub("--file=", "", args[grep("--file=", args)])))
# will be used later
current_dir <- getwd()

# Check if the parameters (.json file) have been given
# reassign 'args' variable
if(length(args <- commandArgs(trailingOnly = TRUE)) == 0){
  cat("\nUsage: ")
  cat("Rscript [--vanilla] MapHitsToCentresOfStrips.R <parameters.json> \n\n")
  stop("No .json file given!")
}

# Load packages
require("methods", quietly = TRUE)    # essential for Rscript
require("jsonlite", quietly = TRUE)
require("data.table", quietly = TRUE) # for fread (data must fit in RAM)
# cd to source dir (== script_dir)
setwd(script_dir)
source("SetMappingEnvironment.R")
source("SetAnhPtsEnvironment.R")
source("InitParamsFromJSON.R")
# go back to initial dir
setwd(current_dir)

# read .json file name with parameters (first argument only)
json_params_file <- args[1]
# assign parameters
ini_params <- InitParamsFromJSON(json_file_path   = json_params_file,
                                 show_reader_info = TRUE)
# set the environment (fast functions) for processing the input data
ann_pts_environment <- SetAnhPtsEnvironment(ini_params[["json_reader"]],
                                            ini_params[["kde_controller"]])

# MAIN FUNCTION - remaps hits onto centres of strips and saves it to file
# validate_hits and progress_info could be changed here (the latter slows the calculation!)
ann_pts_environment$mapAndSave(validate_hits = FALSE,
                               progress_info = TRUE)
