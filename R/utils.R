#' Make 3D matrix of channel concentration values.
#' 
#' @param result_tbl A tibble of concentration values relative to the three dimensions of the eventual matrix (created using hydrogeom and flumeTracer)
#' @param target A character vector of length 1, default is "C_c" (channel concentration)
#' @return 3D concentration matrix with dimensionless concentration values 
#' (decaying from 1 to 0, e.g., scaled relative to initial concentration and 
#' final concentration); and relative to dimensionless time values 
#' (e.g., increasing from 0 to 1 as a fraction of maximum residence time, 
#' time 0 is when "C_c" = 1.0, time 1.0 is when "C_c" = 0) 
#' 
#' @description
#' Make 3D concentration matrix relative to the maximum residence time (tau_n), 
#' the ratio of volume between the hyporheic zone : total (V_ratio) 
#' and the curvature parameter in the residence time distribution 
#' (alpha for power law or sigma for exponential functional forms)
#'  
#' @rdname make_matrix
#' @export
make_matrix <- function(results_tbl, target = "C_c") {
  #extract the names of the axes of the tibble (aside from the target)
  axes_names <- names(results_tbl)[names(results_tbl) != target]
  #extract the range of values for each of the axes
  axis_values <- sapply(axes_names, \(x) sort(unique(results_tbl[[x]])), simplify = F)
  #create unique index for the axis values
  axis_idx <- sapply(axis_values, \(x) 1:length(x), simplify = F)
  #make tibbles of the axis values / indices
  axis_tbls <- Map(\(i,val) dplyr::tibble(idx = i, value = val), axis_idx, axis_values)
  #update names of axis tibbles
  axis_tbls <- Map(\(nm, tbl) {
    names(tbl)[1] <- paste0(nm, "_idx")
    names(tbl)[2] <- nm
    tbl
  }, 
  axes_names,
  axis_tbls)
  # join axis tibbles
  joined_result_tbl <- results_tbl
  for(join_tbl in axis_tbls) {
    joined_result_tbl <- inner_join(joined_result_tbl, join_tbl)
  }
  #query variables of axes
  vars <- paste0(rev(axes_names), "_idx")
  #arrange values 
  joined_result_tbl <- dplyr::arrange(joined_result_tbl, across(all_of(vars)))
  #create 3D matrix of values
  array(
    data = as.vector(unlist(joined_result_tbl[,target])),
    dim = sapply(axis_tbls, nrow),
    dimnames = lapply(axis_values, as.character))
}

#' Dedensify data using a log10 spacing
#'
#' Sample data using evenly spaced orders of magnitude.
#'
#' @param x vector of data values in x dimension or a data.frame
#' @param y vector of data values in y dimension (equal in length to x)
#' @param n number of samples to extract
#' @param x_offset value added to x before log10 is applied; large values
#'   provide more evenly spaced points.  Small values crowd points toward the
#'   left (smaller values of x).
#' @param x_col,y_col name or number of data.frame column containing x dimension
#'   and y dimension
#' @param ... other values passed on to dedensify.
#' @return a tibble of n {x,y} pairs
#' @export
#' @rdname dedensify
#' @export
dedensify <- function(x, y, n, x_offset = 1, decimals = 10) {
  
  #make sure x and y are ordered
  y = y[order(x)]
  x = sort(x)
  
  #add offset (see documentation - controls spacing of de-densified pts)
  x_off = x + x_offset
  
  #find log of offset x values
  log_x = log10(x_off)
  #de-densified x is n-values from min to max of logged values
  log_deden_x = seq(min(log_x), max(log_x), length.out = n)
  #exponentiate back to original x scale
  deden_x = round(10^log_deden_x - x_offset, decimals)
  #find y-values associated with de-densified x values using approx, and return tibble
  approx(x, y, deden_x) |>
    dplyr::as_tibble()
}

#' @rdname dedensify
#' @export
dedensify_df <- function(x, x_col, y_col, n, ...) {
  #if x is a dataframe (instead of individual vectors for x and y), use the column names provided to run dedensify()
  dedensify(x[[x_col]], x[[y_col]], n, ...)
}

