% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CoxMDS.R
\name{medselect_surv}
\alias{medselect_surv}
\title{Data Splitting for Survival Response}
\usage{
medselect_surv(
  X,
  Y,
  MS,
  COV,
  penalty = "MCP",
  method = "ds",
  aggregation = "inclusion_rate",
  q = 0.1,
  cutoff = 0.5,
  gamma = 0.05,
  n_split = 25
)
}
\arguments{
\item{X}{a vector of the exposure.}

\item{Y}{a vector of the outcome.}

\item{MS}{a matrix of the mediators after the pre-liminary screening.}

\item{COV}{optional, a matrix of the potential covariates.}

\item{penalty}{"MCP", "lasso" or "SCAD",the penalty used in the cox regression.}

\item{method}{"ds" or "knockoff", corresponds to selection based on data splitting or knockoff.}

\item{aggregation}{"inclusion_rate" or "quantile", corresponds to aggregation method based on inclusion rate or quantile aggregation.}

\item{q}{the pre-defined FDR level in the data splitting procedure.}

\item{cutoff}{the proportion of the first split if we use data splitting method.}

\item{gamma}{the pre-defined quantile level if we use the quantile aggregation.}

\item{n_split}{the pre-defined number of multiple data splittings/knockoffs.}
}
\value{
A vector of mediator IDs selected by MDS/AKO.
}
\description{
Split the data for survival outcome and compute the corresponding statistics.
}
