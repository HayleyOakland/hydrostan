% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{dedensify}
\alias{dedensify}
\alias{dedensify_df}
\title{Dedensify data using a log10 spacing}
\usage{
dedensify(x, y, n, x_offset = 1, decimals = 10, approx = T)

dedensify_df(x, x_col, y_col, n, ...)
}
\arguments{
\item{x}{vector of data values in x dimension or a data.frame}

\item{y}{vector of data values in y dimension (equal in length to x)}

\item{n}{number of samples to extract}

\item{x_offset}{value added to x before log10 is applied; large values
provide more evenly spaced points.  Small values crowd points toward the
left (smaller values of x).}

\item{approx}{logical (default=T) argument. If TRUE, dedensify() uses approx() 
to find y values associated with de-densified x values;
If FALSE, dedensify finds indices of original data x values closest to the de-densified
x values, then uses those indices to filter y values (from the original data, 
e.g., instead of approximating y values)}

\item{x_col, y_col}{name or number of data.frame column containing x dimension
and y dimension}

\item{...}{other values passed on to dedensify.}
}
\value{
a tibble of n {x,y} pairs
}
\description{
Sample data using evenly spaced orders of magnitude.
}
