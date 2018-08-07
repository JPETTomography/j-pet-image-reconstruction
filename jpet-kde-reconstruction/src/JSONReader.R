############################
# author: Roman Shopa
# Roman.Shopa[at]ncbj.gov.pl
############################

# A reader for .json files containing parameters
JSONReader <- setRefClass(
  "JSONReader",
  fields = list(
    input_data     = "character",
    output_name    = "character",
    scanner        = "list",
    mapper_options = "list",
    kde_parameters = "list",
    no_of_events   = "integer",
    # intended to be private:
    .initial_dir   = "character",
    .file_size     = "numeric",
    .available_RAM = "numeric",
    .total_RAM     = "numeric",
    .is_imported   = "logical"
  ),
  methods = list(
    # assign parameters from JSON file
    initialize = function(..., kde_json_params_file = NULL, show_info = TRUE){
      # initial directory
      .self[[".initial_dir"]] <- getwd()
      if(!is.null(kde_json_params_file)){
        callSuper(...)
        # split incomplete path to directory and filename
        directory_and_json_file <- .self$.splitFileAndDir(kde_json_params_file)
        # cd into dir where .json file is located
        cat("\n----- Working directory: -----\n")
        setwd(directory_and_json_file[1])
        cat(paste0(directory_and_json_file[1],"\n"))
        # read list
        json_list <- read_json(directory_and_json_file[2], simplifyVector = TRUE)
        # assign input/output file names and scanner
        .assignIO(json_list)
        .assignScanner(json_list)
        # detect RAM and the size of file
        .assignMemoryOptions(json_list)
        # mapper for hits and KDE parameters
        .assignMapperList(json_list)
        .assignKDEParameters(json_list)
        # set as successfully imported
        .self[[".is_imported"]] <- TRUE
      } else .self[[".is_imported"]] <- FALSE
      # show parameters or the usage if needed, with possible memory issues
      if(show_info) .self$.usageInfo()
      # get back to initial directory
      setwd(.self[[".initial_dir"]])
    },
    
    # ----- dot (.) denotes 'private' methods -----
    
    # ----- IO -----
    # assigns IO parameters
    .assignIO = function(json_list){
      IO_list <- json_list[["INPUT-OUTPUT"]]
      .self[["input_data"]]  <- .self$.fullPath(IO_list[["input-data-path"]]) # full path
      .self[["output_name"]] <- IO_list[["output-file-name"]]                 # only a name
    },
    # assign scanner
    .assignScanner = function(json_list){
      scanner_path <- .self$.fullPath(json_list[["INPUT-OUTPUT"]][["scanner-geometry-file"]])
      scanner_json <- read_json(scanner_path, simplifyVector = TRUE)
      .self[["scanner"]] <- .self$.generateScanner(scanner_json)
    },
    # ----- SCANNER -----
    # create scanner (as list)
    .generateScanner = function(scanner_json){
      if(.self$.validateScannerJSON(scanner_json)){
        # create an empty list
        output <- list(strip_dimensions = numeric(),
                       radius           = numeric(),
                       no_of_strips     = integer(),
                       angle_bias       = numeric())
        strip_size <- scanner_json[["scintillator-strip-size"]]
        output[["strip_dimensions"]] <- c(strip_size[["tangential-width"]],
                                          strip_size[["radial-width"]],
                                          strip_size[["length"]])
        # a loop through the names of layers
        layers <- scanner_json[["layers"]]
        for(n in names(layers)) {
          output[["radius"]] <- c(output[["radius"]], layers[[n]][["radius"]])
          n_strips <- layers[[n]][["number-of-strips"]]
          output[["no_of_strips"]] <- c(output[["no_of_strips"]], n_strips)
          output[["angle_bias"]]   <- c(output[["angle_bias"]],
                                        layers[[n]][["angle-bias-coefficient"]] * pi/n_strips)
        }
        return(output)
      } else {
        # ERROR!
        setwd(.self[[".initial_dir"]])
        stop("\nunsupported geometry of the scanner! Process terminated.\n")
        # return(list()) # returns an empty list
      }
    },
    # validate scanner (only cylindrical), not sensitive to .self (generic method)
    .validateScannerJSON = function(scanner_json){
      # another geometry to be added here in the future...
      return(scanner_json[["geometry"]] %in% c("cylindrical")) 
    },
    
    # ----- MEMORY -----
    # estimate RAM/filesize and assign big matrix parameters if needed
    .assignMemoryOptions = function(json_list){
      IO_list <- json_list[["INPUT-OUTPUT"]]
      input_data_path <- .self$.fullPath(IO_list[["input-data-path"]])
      .self[[".file_size"]]     <- as.numeric(file.info(input_data_path)[["size"]])
      # esimate RAM (in bytes) using shell command
      .self[[".available_RAM"]] <- as.numeric(system("awk '/MemFree/ {print $2}' /proc/meminfo ", 
                                                     intern = TRUE))*1e3
      .self[[".total_RAM"]] <- as.numeric(system("awk '/MemTotal/ {print $2}' /proc/meminfo ", 
                                                  intern = TRUE))*1e3
    },
    # --- MAPPER and KDE ---
    # assigns mapper_options (a list)
    .assignMapperList = function(json_list){
      # assign mapper (differs from initial JSON list)
      mapper_list <- json_list[["HIT-MAPPER-OPTIONS"]]
      map_options <- list(map_from    = .self$.getMapFrom(mapper_list[["map-from"]]),
                          filter_true = mapper_list[["filter-true-events"]])
      # set flags for save_as_ASCII
      map_options[["save_as_ascii"]] <- list(
        mapped_hits         = mapper_list[["save-as-ASCII"]][["mapped-hits"]],
        annihilation_points = mapper_list[["save-as-ASCII"]][["annihilation-points"]]
      )
      # final assignment
      .self[["mapper_options"]] <- map_options
    },
    # validates correct assignment for map_from field, not sensitive to .self
    .getMapFrom = function(map_tag){
      if(length(map_tag)!=0){
        map_tag <- tolower(map_tag)
        if(map_tag %in% c("strips_id","hits","none")) return(map_tag) 
        else {
          cat(paste0("\nWarning: unsuported map mode '", map_tag,"'!\n"))
          cat(paste0("Set back to 'strips_id' (default).\n"))
          return("strips_id")}
      } else return("strips_id")
    },
    
    # assigns parameters for KDE (not DRY, but easier)
    .assignKDEParameters = function(json_list){
      KDE_options <- json_list[["KS-PACKAGE-OPTIONS"]]
      KDE_list <- list(
        optimise_H_matrix = KDE_options[["optimise-H-matrix"]],
        normalise_by      = KDE_options[["normalise-by"]],
        resolution        = c(KDE_options[["resolution"]][["X"]],
                              KDE_options[["resolution"]][["Y"]],
                              KDE_options[["resolution"]][["Z"]]),
        # axis ranger will be set to NULL if auto-detect
        axis_ranges_list  = if(KDE_options[["axis-ranges"]][["auto-detect"]]) NULL
        else list(x = KDE_options[["axis-ranges"]][["X"]],
                  y = KDE_options[["axis-ranges"]][["Y"]],
                  z = KDE_options[["axis-ranges"]][["Z"]])
      )
      # set flags for save_as_ASCII
      KDE_list[["save_as_ascii"]] <- KDE_options[["save-as-ASCII"]]
      # final assignment
      .self[["kde_parameters"]] <- KDE_list
    },
    
    # ---- INFO ----
    # returns usage if incorrectly imported
    .usageInfo = function(){
      if(.self[[".is_imported"]]){
        cat("----- Input parameters -----\n")
        cat(paste("Input file path:",          .self[["input_data"]],"\n"))
        cat(paste("Output file name:",         .self[["output_name"]],"\n"))
        cat(paste("RAM available (bytes):",    .self[[".available_RAM"]], "\n"))
        cat(paste("Total RAM (bytes):"    ,    .self[[".total_RAM"]], "\n"))
        cat(paste("Size of file (bytes):",     .self[[".file_size"]], "\n"))
        cat(paste("Map mode:",                 .self[["mapper_options"]][["map_from"]], "\n"))
        cat(paste("Filter true coincidences:", .self[["mapper_options"]][["filter_true"]], "\n"))
        cat("\n----- Memory Warnings/Errors -----\n")
      }
      cat(.self$.memoryWarningsOrErrors()) # about memory
      cat("----------------------------------\n\n")
    },
    # mode warnings if needed
    # 0 (pass) - RAM preferrable, 
    # 1 (warning) - the size of file is comparable with RAM
    # 2 (error) - only disk might be used (RAM < size of file)
    .modeWarnings = function(){
      enormous_file <- .self[[".file_size"]] >= .self[[".available_RAM"]]
      # file is at least half size of memory available
      big_file <- .self[[".file_size"]] > .self[[".available_RAM"]]/2L
      if(!.self[[".is_imported"]]) return(as.integer(NA)) # if no args passed
      if(enormous_file) return(2L)
      else if(big_file) return(1L) else return(0L)
    },
    # messages about memory
    .memoryWarningsOrErrors = function(){
      if(!.self[[".is_imported"]]){
        # ERROR!
        stop("Incorrect import parameters! Process terminated.\n")
      }
      # 2 - file bigger than RAM
      if(.self$.modeWarnings()==2L){
        .self[[".is_imported"]] <- FALSE
        # ERROR!
        setwd(.self[[".initial_dir"]])
        stop(paste("Not sufficient memory for a large file!",
                     "Process terminated.\n"))
        # 1 - big file (> RAM/4)
      } else if(.self$.modeWarnings()==1L){
        return(paste("Warning: file might be too large for the RAM available.\n"))
      } else {
        return("No warnings/errors detected.\n")
      }
    },
    # --- functions related to handling full path to file ---
    # returns full file path
    .fullPath = function(incomplete_path){
      return(system(paste0("readlink -f ", incomplete_path), intern = TRUE))
    },
    # splits incomplete file path to full directory and filename
    .splitFileAndDir = function(file_path){
      full_path <- .self$.fullPath(file_path)
      directory <- dirname(full_path)
      # cut out directory path
      file_name <- unlist(strsplit(full_path, paste0(directory,"/")))[2]
      return(c(directory, file_name))
    }
  )
)
