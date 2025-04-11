#' Dedensify data using a log10 spacing
#'
#' Sample data using evenly spaced orders of magnitude.
#'
#' @param x vector of data values in x dimension or a data.frame
#' @param x_col,ycol name or number of data.frame column containing x dimension and y
#'   dimension
#' @param y vector of data values in y dimension (equal in length to x)
#' @param n number of samples to extract
#' @return a tibble of n {x,y} pairs
#' @export
#' @rdname dedensify
#' @export
dedensify <- function(x, y, n) {
  y = y[order(x)]
  x = sort(x)

  n_zeros <- sum(x==0)
  if(sum(x<0) > 0 | n_zeros > 1) stop("Only one zero value and no negative values are allowed in x.")

  if(n_zeros) x = x[-1]
  log_x = log10(x)
  log_deden_x = seq(min(log_x), max(log_x), length.out = n)
  deden_x = 10^log_deden_x

  if(n_zeros) {
    x <- c(0, x)
    deden_x <- c(0, deden_x)
  }

  approx(x, y, deden_x) |>
    dplyr::as_tibble()
}

#' @rdname dedensify
#' @export
dedensify_df <- function(x, x_col, y_col, n) {
  dedensify(x[[x_col]], x[[y_col]], n = n)
}

