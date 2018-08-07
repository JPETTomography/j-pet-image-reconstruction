############################
# author: Roman Shopa
# Roman.Shopa[at]ncbj.gov.pl
############################

require("methods")    # essential for Rscript
require("jsonlite")

source("JSONReader.R")
source("KDEController.R")

# sets reader and KDE controller with the correct path to .json file with parameters
InitParamsFromJSON <- function(json_file_path, 
                               show_reader_info = TRUE){
  # create reader for parameters
  json_reader <- JSONReader(kde_json_params_file = json_file_path, 
                            show_info            = show_reader_info)
  # create KDE controller
  kde_controller <- KDEController(json_reader = json_reader)
  # full path to .json parameters
  json_file_dir_and_name <- json_reader$.splitFileAndDir(json_file_path)
  return(list(parameters_file = json_file_dir_and_name[2],
              directory       = paste0(json_file_dir_and_name[1],"/"),
              json_reader     = json_reader,
              kde_controller  = kde_controller))
}