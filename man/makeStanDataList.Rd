% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{makeStanDataList}
\alias{makeStanDataList}
\title{Make list of data required for stan run}
\usage{
makeStanDataList(
  concData,
  Cvals,
  flumeInfo,
  avgKprimeVals,
  fracExTau,
  model,
  resultsMatrix,
  include_drift,
  include_qChange
)
}
\arguments{
\item{concData}{Data frame containing concentration by time data per trial (trialIdx x C x time).}

\item{Cvals}{Data frame containing starting and ending concentration values per trial.}

\item{flumeInfo}{Data frame containing information about the flumes used in the trials, including `flumeIdx`, `phase` (1=beginnning of experiment, 2=end of experiment), and `finalCaddisN` (the final abundance of caddisflies in each flume)}

\item{avgKprimeVals}{Data frame containing avgKprime values per trial.}

\item{fracExTau}{Vector of length 1 containing numeric value of time to find the fraction of the hyporheic zone that has exchanged by this given time (e.g., "30", if concentration DF times are in minutes, will find the fraction of the hyporheic zone exchanged by 30 minutes); equivalent to I(tau)}

\item{model}{The model to run in stan. Options include "wellMixedModel" for a basic exponential decay model (assuming a well-mixed hyporheic zone), "powerLawRTD" for the power law functional form of an RTD model of the hyporheic zone, or "exponentialRTD" for the exponential function form of the RTD model of the hyporheic zone}

\item{resultsMatrix}{Needed if model != wellMixedModel. Matrix of concentrations from the hydrogeom model for interpolation, relative to values of the curvature parameter (alpha for PLRTD or sigma for ExpRTD), volume ratio, and time}

\item{include_drift}{Logical statement, whether to include drift as a parameter in the stan model.}

\item{include_qChange}{Logical statement, whether to include a change in q across trials (e.g., if there is more than one release per flume, are we testing the statistical model that q_down changed as a function of caddisfly density)}
}
\description{
Make list of data required for stan run
}
