############################
# author: Roman Shopa
# Roman.Shopa[at]ncbj.gov.pl
############################

# Sets the environment for the output as annihilation points 
# reader and KDE_controller (for range)
SetAnhPtsEnvironment <- function(reader, kde_controller){
  # validate input parameters
  if(!kde_controller[[".is_validated"]] | !reader[[".is_imported"]]){
    cat("ERROR: Bad parameters set in .json file and/or KDE controller.\n")
    return(NULL)
  }
  .input_file         <- reader[["input_data"]]
  .output_for_hits    <- paste0(.input_file,"_mapped")
  .output_for_points  <- paste0(.input_file,"_annihilation_points")
  .axis_ranges_matrix <- kde_controller[["axis_ranges_matrix"]]
  .mapper_type        <- reader[["mapper_options"]][["map_from"]]       # strip IDs or hits
  .filter_true        <- reader[["mapper_options"]][["filter_true"]]
  .save_as_ascii      <- reader[["mapper_options"]][["save_as_ascii"]]
  # environment for mapping hits or strip IDs onto Cartesians
  .mapperEnv <- SetMappingEnvironment(reader)
  # define strips map if needed
  if(.mapper_type == "strips_id") 
    .mapperEnv[[".strips_map"]] <- .mapperEnv$.defineStripsMap(reader[["scanner"]])
  
  # ----- IMPORT DATA -----
  # this would read the data and set transformation functions
  input_data <- fread(.input_file, header = FALSE)
  # filter only true coincidences if set (see GOJA documentation)
  if(.filter_true) input_data <- input_data[as.logical(input_data[,13] == 1),]
  
  # ----- FUNCTIONS -----
  
  
  # ----- SAVE DATA -----
  # function for saving mapped events 
  mapAndSave <- function(validate_hits = FALSE, 
                         progress_info = FALSE){
    # Transposition of events onto centres of strips
    mapped_hits <- .mapMultipleEvents(data          = input_data,
                                      mapper        = .mapperEvents,
                                      progress      = progress_info,
                                      progress_step = 7500L,
                                      validate_hits = validate_hits)
    cat("Info: remapped data will be saved in ")
    if(.save_as_ascii[["mapped_hits"]]){
      # save as ASCII
      cat("ASCII format.\n")
      # Transposition of events onto centres of strips
      fwrite(as.data.table(mapped_hits), .output_for_hits, 
                           sep = "\t", row.names = FALSE, col.names = FALSE)
      cat(paste("Saved to:\n", .output_for_hits, "\n"))
    } else {
      # save as workspace
      cat(".RData format.\n")
      file_name <- paste0(.output_for_hits,".RData")
      save(mapped_hits, file = file_name)
      cat(paste("Saved to:\n", file_name, "\n"))
    }
  }
  
  # function for mapping and saving annihilation points only
  estimateAndSaveAnnPoints <- function(validate_hits, progress_info){
    # Map hits first
    mapped_hits <- .mapMultipleEvents(data          = input_data,
                                      mapper        = .mapperEvents,
                                      progress      = progress_info,
                                      progress_step = 7500L,
                                      validate_hits = validate_hits)
    cat("Info: annihilation points will be saved in ")
    if(.save_as_ascii[["annihilation_points"]]){
      # save as ASCII
      cat("ASCII format.\n")
      file_name <- .output_for_points
      fwrite(as.data.table(composeAnhPts(mapped_hits)), file_name, 
             sep = "\t", row.names = FALSE, col.names = FALSE)
    } else {
      # save as workspace
      cat(".RData format.\n")
      annihilation_points <- composeAnhPts(mapped_hits)
      file_name <- paste0(.output_for_points,".RData")
      save(annihilation_points, file = file_name)
    }
    cat(paste0("Successfully saved to:\n", file_name,"\n"))
  }

  # --- Important function! ---
  # composes array of annihilation points in a given range (X,Y,Z)
  # does not use internal data! Any arbitrary data could be given
  composeAnhPts <- function(data, 
                            subrange           = TRUE, # logical
                            axis_ranges_matrix = .axis_ranges_matrix){
    # check if kde_controller is ok
    if(!subrange){
      axis_ranges_matrix <- NULL # for validation
    }
    cat("Estimation of annihilation points has started. Please, wait... ")
    # estimate annihilation points here (main loop)
    output <- t(apply(data, 1,
                      function(event){
                        event <- .mapperEnv[[".coerceToNumeric"]](event) # important for data.frame
                        xyz <- .detectAnnihilationPosition(event)        # points from event
                        if(.validateAnnihilationPoint(xyz,
                                                      subrange,
                                                      axis_ranges_matrix[1,],
                                                      axis_ranges_matrix[2,],
                                                      axis_ranges_matrix[3,]))
                          return(xyz)
                        else return(rep(NA,3))
                      }))
    # filter out NAs (by first column is enough)
    output <- output[!is.na(output[,1]),]
    cat("Done!\n")
    # add ranges if given (since KDE is sensitive) and return (REDUNDANT!)
    # if(!is.null(axis_ranges_matrix))  return(rbind(output, t(axis_ranges_matrix)))
    return(output)
  }  

  # ----- MAPPERS -----
  # defines mapper function for events (by hit or by strip_ID)
  .mapperEvents <- (
    if(.mapper_type == "none") NULL  # no return needed
    else function(event,
                  mapper_type   = .mapper_type,
                  validate_hits = TRUE){
      if(mapper_type == "hits") {
        # uses brackets notation for functions as it is faster than '$'
        # event <- .mapperEnv[[".coerceToNumeric"]](event)  # for data.frame (redundant!)
        return(round(c(.mapperEnv[["mapperFunction"]](event[1:2]), event[3:4],
                       .mapperEnv[["mapperFunction"]](event[5:6]), event[7:8]), 2))
      } else {
        event <- .mapperEnv[[".coerceToNumeric"]](event)  # for data.frame
        return(round(c(.mapperEnv[["mapperFunction"]](strip_ID = event[9],
                                                      hit = (if(validate_hits) event[1:2]
                                                             else NULL)),
                       event[3:4],
                       .mapperEnv[["mapperFunction"]](strip_ID = event[10],
                                                      hit = (if(validate_hits) event[5:6]
                                                             else NULL)),
                       event[7:8]),2))  
      }
    }
  )
  # returns remapped data according to the mapper given (additional parameters are optional)
  .mapMultipleEvents <- function(data, 
                                 mapper,
                                 progress      = TRUE,   # shows No of events processed
                                 progress_step = 7500L,  # info given every this step
                                 ...){                   # additional parameters allowed
    # add progress info if given (not very elegant solution)
    if(progress){
      iter <- 0L  # iterator
      cat("Remapping the data. Please, wait.\n")
      cat("Progress (events processed):\n")
      .mapper <- function(event, ...){
        iter <<- iter + 1L
        if(iter %% progress_step == 0){
          message(iter,"\r",appendLF = F)
          flush.console()
        }
        return(mapper(event, ...))
      }
    } else {
      # no progress but slightly faster
      cat("Remapping the data. Please, wait... ")
      .mapper <- mapper
      rm(mapper)
    }
    # -- main loop - apply mapper to every row --
    mapped_data <- t(apply(data, 1, .mapper, ...))
    # find indices with NA's and update data
    NA_factor <- apply(mapped_data, 1, function(event) sum(is.na(event)) > 0)
    mapped_data <- mapped_data[!NA_factor,]
    cat("Done!\n")
    return(mapped_data)
  }  

  # -------- routine functions ---------
  # detects centre of one LOR
  .detectCentre <- function(event) {
    (event[1:3] + event[5:7])/2
  }
  # detects Cartesians for the annihilation using TOF
  .detectAnnihilationPosition <- function(event, 
                                          speed_of_light = 299792458){
    LOR_length <- sqrt(sum((event[1:3]-event[5:7])^2))
    LOR_centre <- .detectCentre(event)
    # 1e-12 * 1e2 (into cm)
    delta_LOR <- as.numeric(speed_of_light*1e-10*(event[4]-event[8])/2)
    fraction <- 2*delta_LOR/LOR_length
    return(fraction*(event[5:7] - LOR_centre) + LOR_centre)
  }
  # checks whether the coordinates are non-NA and within the ranges if given
  .validateAnnihilationPoint <- function(coordinates,
                                         subrange = FALSE, # set TRUE if ranges are given
                                         x_range  = NULL,
                                         y_range  = NULL,
                                         z_range  = NULL){
    # non-NA
    NA_factor <- sum(is.na(coordinates)) != 0
    # not DRY but faster
    x <- coordinates[1]
    y <- coordinates[2]
    z <- coordinates[3]
    # check whether in range
    if(subrange){ 
      range_factor <- x >= x_range[1] & x <= x_range[2] &
                      y >= y_range[1] & y <= y_range[2] &
                      z >= z_range[1] & z <= z_range[2]
    } else range_factor <- TRUE
    # # prevent from being logical(0) (redundant!)
    # if(length(range_factor) == 0) range_factor <- TRUE 
    return(!NA_factor & range_factor)
  }
  # remove redundant fields
  rm(list=c("reader","kde_controller"))
  return(environment())
}
