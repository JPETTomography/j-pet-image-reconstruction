############################
# author: Roman Shopa
# Roman.Shopa[at]ncbj.gov.pl
############################

# A controller for KDE procedure
KDEController <- setRefClass(
  "KDEController",
  fields = list(
    resolution          = "numeric",  # the size in pixels of the reconstructed image
    voxel_size          = "numeric",  # resulting size
    axis_ranges_matrix  = "matrix",   # sets boundaries for coordinates (or the full FOV)
    optimise_H_matrix   = "logical",  # if TRUE, no Hpi() will be executed (set from voxel size)
    save_as_ascii       = "logical",  # if FALSE, saves the output data as the workspace (.RData)
    .max_radius         = "numeric",  # a radius of the largest layer
    .scanner_length     = "numeric",  # PET scanner length
    .is_validated       = "logical"   # whether the input parameters are assigned
  ),
  methods = list(
    initialize = function(..., json_reader = NULL){
      if(!is.null(json_reader)){
        callSuper(..., .is_validated = TRUE)
        # assign variables from JSON reader
        kde_params <- json_reader[["kde_parameters"]]
        .self[["resolution"]]        <- kde_params[["resolution"]]
        .self[["optimise_H_matrix"]] <- kde_params[["optimise_H_matrix"]]
        .self[[".max_radius"]]       <- max(json_reader[["scanner"]][["radius"]])
        .self[[".scanner_length"]]   <- json_reader[["scanner"]][["strip_dimensions"]][3]
        # assign save as ascii
        .self[["save_as_ascii"]] <- kde_params[["save_as_ascii"]]
        # set ranges for subsetting annihilation points within the span
        axis_ranges_list <- kde_params[["axis_ranges_list"]]
        if(is.null(axis_ranges_list)){
          # set ranges as full FOV (R/sqrt(2)) for the largest radius
          max_xy_value <- .self[[".max_radius"]]/sqrt(2)
          max_z_value  <- .self[[".scanner_length"]]/2
          axis_ranges_list <- list(x = c(-max_xy_value, max_xy_value),
                                   y = c(-max_xy_value, max_xy_value),
                                   z = c(-max_z_value, max_z_value))
        } 
        # set axis ranges matrix
        .self[["axis_ranges_matrix"]] <- .self$.toAxisMatrix(axis_ranges_list)
        # estimate voxel size
        .self[["voxel_size"]] <- .self$.estimateVoxelSize(.self[["resolution"]], 
                                                          .self[["axis_ranges_matrix"]])
      } else .self[[".is_validated"]] <- FALSE
    },
    
    # ----- dot (.) denotes 'private' methods -----
    # returns axis matrix from axis list (input)
    .toAxisMatrix = function(axis_ranges_list){
      return(sapply(c(1,2), function(s) c(axis_ranges_list[["x"]][s],
                                          axis_ranges_list[["y"]][s],
                                          axis_ranges_list[["z"]][s])))
    },
    # estimate voxel size (vector pf 3 values)
    .estimateVoxelSize = function(resolution, axis_ranges_matrix){
      return(apply(axis_ranges_matrix, 1, diff)/resolution)
    }
  )
)
