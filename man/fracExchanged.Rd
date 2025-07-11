% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{fracExchanged}
\alias{fracExchanged}
\title{Compute fraction of hyporheic zone that has exchanged by time t}
\usage{
fracExchanged(shape, t, tau_0, tau_n, curvParVal)
}
\arguments{
\item{shape}{is a character vector of length 1, the shape of the residence time distribution (options: "exponent" or "powerLaw")}

\item{t}{is a numeric vector of length 1, the time of interest (same units as tau_0 and tau_n) to compute the fraction of the hyporheic zone that has exchanged with the channel}

\item{tau_0}{is a numeric vector of length 1, the minimum residence time (can be 0 if shape="exponent", must be >0 if shape="powerLaw") (same time units as t and tau_n)}

\item{tau_n}{is a numeric vector of length 1, the maximum residence time (same time units as t and tau_0)}

\item{curvParVal}{is a numeric vector of length 1, the curvature parameter in the residence time distribution function (sigma if shape = "exponent", alpha if shape = "powerLaw")}
}
\value{
Numeric vector of length 1, the value of the fraction of the 
hyporheic zone that has exchanged with the channel by given 
time, t, relative to residence time distribution parameters 
(tau_0, tau_n and the curvature parameter, sigma or alpha).
}
\description{
The fraction of the hyporheic zone that has exchanged by a given time `t`
can be computed as the integral from tau_0 to t of the washout function
divided by the integral of the washout function.
}
