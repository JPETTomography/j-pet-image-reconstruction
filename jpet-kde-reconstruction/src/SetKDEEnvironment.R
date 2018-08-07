############################
# author: Roman Shopa
# Roman.Shopa[at]ncbj.gov.pl
############################

# Sets the environment for KDE prodedure
SetKDEEnvironment <- function(kde_controller){
  # validate controller
  if(!kde_controller[[".is_validated"]]){
    stop("No parameters detected in KDE controller! No environment has been set.\n")
  }
  # Returns 'ks' object with image (optionally normalised), axis, H matrix etc.
  reconstructImage <- function(annihilation_points, 
                               normalise_by = NULL, ...){
    # get Hpi matrix
    cat("\t Estimating Hpi matrix... ")
    if(!kde_controller[["optimise_H_matrix"]]) H_matrix <- Hpi(annihilation_points, ...)
    else H_matrix <- 1.75e-1 * diag(kde_controller[["voxel_size"]]) # optimised value
    cat("Done!\n")
    # Execute multivariate KDE (expensive!)
    cat("\t Building fhat() function... ")
    f_hat <- kde(annihilation_points, 
                 H        = H_matrix,
                 xmin     = kde_controller[["axis_ranges_matrix"]][,1],
                 xmax     = kde_controller[["axis_ranges_matrix"]][,2],
                 gridsize = kde_controller[["resolution"]])
    cat("Done!\n")
    # normalise if needed
    if(is.null(normalise_by) ) return(f_hat)
    else {
      if(!.validate_normalise(normalise_by)){
        cat(paste("\t Warning: unknown normalise (must be 'peak' or 'sum').",
                  "The image has remained unchanged.\n"))
        return(f_hat)
      }
      f_hat[["estimate"]] <- .normalise(f_hat[["estimate"]], by = tolower(normalise_by))
      cat(paste0("\t Info: the output image has been normalised by ", normalise_by, ".\n"))
      return(f_hat)
    }
  }
  # ----- dot (.) denotes 'private' methods -----
  # 'sum' (or else) - normalise by sum, 'peak' - by peak
  .normalise <- function(x, by = "peak"){
    if(by == "peak") return((x-min(x))/(max(x)-min(x))) 
    else return(x/sum(x))
  }
  .validate_normalise <- function(normalise_by){
    return(tolower(normalise_by) %in% c("peak", "sum"))
  }
  return(environment())
}
