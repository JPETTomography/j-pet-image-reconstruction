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
  cat("Rscript [--vanilla] EngageKDEReconstruction.R <parameters.json> \n\n")
  stop("No .json file given!")
}

# Load packages
require("methods", quietly = TRUE)    # essential for Rscript
require("jsonlite", quietly = TRUE)
require("data.table", quietly = TRUE) # for fread (data must fit in RAM)
# cd to source dir (== script_dir)
setwd(script_dir)
source("InitParamsFromJSON.R")
source("SetKDEEnvironment.R")
# go back to initial dir
setwd(current_dir)

# read .json file name with parameters (first argument only)
json_params_file <- args[1]
# assign parameters
ini_params <- InitParamsFromJSON(json_file_path   = json_params_file,
                                 show_reader_info = FALSE)

# --- IMPORT DATA ---
# assign correct filename for the data containing annihilation points (if existed)
ascii_tag <- ini_params[["json_reader"]][["mapper_options"]][["save_as_ascii"]][["annihilation_points"]]
if(ascii_tag){ input_file_name <- paste0(ini_params[["json_reader"]][["input_data"]],
                                        "_annihilation_points")
} else input_file_name <- paste0(ini_params[["json_reader"]][["input_data"]],
                                 "_annihilation_points.RData")
cat(paste0("\n----- File to proceed : -----\n", 
           input_file_name, "\n"))
# validate existence of a file
if(!file.exists(input_file_name)){
  cat("NO SUCH FILE! Execute DetectAnnihilationPoints.R first.\n\n")
  stop("No file for annihilation points detected! Process terminated.\n")
}

# assign data to the variable 'annihilation_points'
if(ascii_tag){
  annihilation_points <- fread(input_file_name, header = FALSE)
} else load(input_file_name) # .RData input format

# --- ESTIMATE KDE ---
require("ks", quietly = TRUE)         # for miltivariate KDE
# create KDE environment (fast functions) for processing the input data
kde_environment <- SetKDEEnvironment(ini_params[["kde_controller"]])
# estimate kde (normalise setting is set in reader)
cat("\n--- Starting KDE reconstruction ---\n")
image_object <- 
  kde_environment$reconstructImage(annihilation_points,
                                   ini_params[["json_reader"]][["kde_parameters"]][["normalise_by"]])
cat("-----------------------------------\n")

# --- SAVE THE IMAGE ---
# cd into working directory (might be dublicated)
setwd(ini_params[["directory"]])

if(ini_params[["json_reader"]][["kde_parameters"]][["save_as_ascii"]]){
  cat("Info: will be saved in ASCII format.\n")
  # save to separate directory
  system("mkdir -p ReconstructedImage")
  # save axis first
  axis <- image_object[["eval.points"]]
  fwrite(as.data.table(axis[[1]]), "ReconstructedImage/x_axis", 
         sep = "\t", row.names = FALSE, col.names = FALSE)
  fwrite(as.data.table(axis[[2]]), "ReconstructedImage/y_axis", 
         sep = "\t", row.names = FALSE, col.names = FALSE)
  fwrite(as.data.table(axis[[3]]), "ReconstructedImage/z_axis", 
         sep = "\t", row.names = FALSE, col.names = FALSE)
  # save image as melted data frame
  require("reshape2", warn.conflicts = FALSE, quietly = TRUE)
  image_melted <- as.data.table(melt(image_object[["estimate"]]))
  names(image_melted) <- c("X","Y","Z","Intensity")
  fwrite(image_melted, paste0("ReconstructedImage/",ini_params[["json_reader"]][["output_name"]]), 
         sep = "\t", row.names = FALSE, col.names = TRUE)
  cat(paste0("Successfully saved to directory:\n",
             ini_params[["directory"]],"ReconstructedImage/\n"))
} else {
  # save as a workspace .RData
  cat("Info: will be saved in .RData format.\n")
  image_name <- paste0(ini_params[["json_reader"]][["output_name"]],".RData")
  save(image_object, file = image_name)
  cat(paste0("Successfully saved to:\n", ini_params[["directory"]], image_name,"\n"))
}
