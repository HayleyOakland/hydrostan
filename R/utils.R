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
#' @param x_col,ycol name or number of data.frame column containing x dimension
#'   and y dimension
#' @param ... other values passed on to dedensify.
#' @return a tibble of n {x,y} pairs
#' @export
#' @rdname dedensify
#' @export
dedensify <- function(x, y, n, x_offset = 1, decimals = 10) {

  y = y[order(x)]
  x = sort(x)

  x_off = x + x_offset

  log_x = log10(x_off)
  log_deden_x = seq(min(log_x), max(log_x), length.out = n)
  deden_x = round(10^log_deden_x - x_offset, decimals)

  approx(x, y, deden_x) |>
    dplyr::as_tibble()
}

#' @rdname dedensify
#' @export
dedensify_df <- function(x, x_col, y_col, n, ...) {
  dedensify(x[[x_col]], x[[y_col]], n, ...)
}

