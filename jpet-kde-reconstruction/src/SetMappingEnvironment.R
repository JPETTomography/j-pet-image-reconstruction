############################
# author: Roman Shopa
# Roman.Shopa[at]ncbj.gov.pl
############################

# Creates environment for mapping hits/events (faster than any type of class)
SetMappingEnvironment <- function(json_reader){
  # validate geometry of the scanner
  if(!json_reader[[".is_imported"]]){
    cat("ERROR: No scanners detected! No environment has been set.\n")
    return(NULL)
  }
  # these parameters will be separate (faster)
  .strip_dimensions <- json_reader[["scanner"]][["strip_dimensions"]]
  .radius           <- json_reader[["scanner"]][["radius"]]
  .no_of_strips     <- json_reader[["scanner"]][["no_of_strips"]]
  .angle_bias       <- json_reader[["scanner"]][["angle_bias"]]
  .mapper_type      <- json_reader[["mapper_options"]][["map_from"]]
  
  ##### Main functions #####
  mapperFunction <- (
    if(.mapper_type == "none") NULL # nothing in return
    else if(.mapper_type == "hits") 
      # maps 2D hit (X,Y) into discrete strip centres
      function(hit){
        # preprocessing
        special_output <- rep(as.numeric(NA),2)
        hit <- .coerceToNumeric(hit)                    # important for data frame
        if(!.validateHit(hit)) return(special_output)   # NAs if false
        # detect indices
        layer_ID <- .detectLayer(hit)
        strip_ID <- .getStripIndex(hit, 
                                   .no_of_strips[layer_ID], 
                                   .angle_bias[layer_ID])
        # key parameters
        hit_radius   <- .getRadius(hit)
        hit_angle    <- .getAngle(hit)
        strip_angle  <- .getDiscreteAngle(strip_ID, layer_ID)
        layer_radius <- .radius[layer_ID]
        # output (NAs or c(X,Y), remapped to the central axis of the strip)
        if(!.isInsideDetector(hit_radius, hit_angle, layer_radius, strip_angle) |
           !.isInLayer(hit_radius, layer_radius)) return(special_output)
        else return(round(layer_radius*c(cos(strip_angle),
                                         sin(strip_angle)),2))
      }
    else 
      # maps strip index onto pair of Cartesians (XY)
      function(strip_ID, hit = NULL){
        # preprocessing for the case if hit is given (additional validation)
        if(!is.null(hit)){
          hit <- .coerceToNumeric(hit)                         # important for data frame
          if(!.validateHit(hit)) return(rep(as.numeric(NA),2)) # special output - a pair of NAs
        }
        # detect layer index and reduced strip index
        layer_ID <- .strips_map[strip_ID, 2]
        strip_ID <- .strips_map[strip_ID, 1] - sum(.no_of_strips[0:(layer_ID-1)])
        # estimate angle
        strip_angle <- .getDiscreteAngle(strip_ID, layer_ID)
        return(round(.radius[layer_ID]*c(cos(strip_angle),
                                         sin(strip_angle)),2))
      }
  )

  # --- assumed to be PRIVATE: denoted by (.) dots ---

    # gets angle from X,Y pair (first 2 elevents of a hit as a vector)
  .getAngle <- function(hit){
    angle <- atan2(hit[2], hit[1])
    if(angle >= 0) return(angle)
    else return(angle + 2*pi)
  }
  # gets discrete angle from the pair of (strip ID, layer ID)
  .getDiscreteAngle <- function(detector_ID, 
                                layer_ID) {
    return(2*pi*(detector_ID - 1L)/.no_of_strips[layer_ID] + .angle_bias[layer_ID])
  }
  # prevents angle to be outside the range (0,2*pi)
  .weightedAngle <- function(alpha){
    return(alpha %% (2*pi))
  }
  # converts data.frame to humeric format
  .coerceToNumeric <- function(hit){
    if(class(hit)[1] != "numeric") return(as.numeric(hit)) # only first element (if many)
    else return(hit)
  }
  # radius from the centre (for X,Y pair)
  .getRadius <- function(hit){
    return(sqrt(hit[1]^2 + hit[2]^2))
  }
  # whether the hit is in the layer (depth of the strip validation)
  .isInLayer <- function(annihilation_radius,
                         layer_radius, 
                         strip_half_depth = .strip_dimensions[2]/2) {
    return(annihilation_radius >= layer_radius - strip_half_depth &
             annihilation_radius <= layer_radius + strip_half_depth)
  }
  # validate whether the hit is not outside the stip (width)
  .isInsideDetector <- function(radius,
                                angle,
                                layer_radius,
                                detector_angle){
    return(abs(radius*sin(angle-detector_angle)) <= .strip_dimensions[1]/2)
  }
  
  # validate if any elements are NA
  .hitIsNA <- function(hit){
    return(is.na(hit[1]) | is.na(hit[2]))
  }
  # transforms values to numeric if data frame
  # detects layer ID from hit (returns NA if outside)
  .detectLayer <- function(hit){
    r <- .getRadius(hit)
    layer_ID = 1L*.isInLayer(r, .radius[1]) + 
      2L*.isInLayer(r, .radius[2]) +
      3L*.isInLayer(r, .radius[3])
    if(layer_ID==0) return(NA) else return(layer_ID)
  }
  # validates whether the hit is inside detector 
  # and does not contain NAs
  .validateHit <- function(hit){
    if(.hitIsNA(hit)) return(FALSE)        # first check X and Y
    else return(!is.na(.detectLayer(hit))) # then whether it is within layer
  }
  # estimates ID of a strip
  # bias_for_first is the angular location of the first detector 
  # (such as n_i * 3.75deg + 1.875deg)
  .getStripIndex <- function(hit,
                             n_detectors,
                             bias_for_first = 0){
    # NA if no hit 
    if(.hitIsNA(hit)) return(as.integer(NA))
    # add half width of 'angle':
    # slight negative angle would still correspond to 1st detector
    angle <- .weightedAngle(.getAngle(hit) - bias_for_first + pi/n_detectors)
    # avoid zero index
    if(angle == 0) return(1L)
    else return(ceiling(n_detectors*(angle)/(2*pi)))
  }

  if(.mapper_type == "strips_id")
    # maps sequential IDs of multilayer scanner into table c(StripID,LayerID)
    # for cylindrical (ideal, big barrel) only! 
    .defineStripsMap <- function(scanner_obj){
      # create an empty map of the size (N_total_strips * 2)
      strips_in_total <- sum(scanner_obj[["no_of_strips"]])
      strMap <- array(rep(NA, strips_in_total * 2),
                      dim = c(strips_in_total, 2))
      for(i in 1:length(scanner_obj[["no_of_strips"]])) {
        IDs <- seq(1:scanner_obj[["no_of_strips"]][i]) + 
          sum(scanner_obj[["no_of_strips"]][0:(i-1)])
        strMap[IDs,] <- cbind(IDs, rep(i,length(IDs)))
      }
      return(strMap)
    }
  # -----------------------------------------------------
  # remove redundant fields
  rm(list = c("json_reader"))
  return(environment())
}
