% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wrapper_LbyL.R
\name{wrapper_LbyL}
\alias{wrapper_LbyL}
\title{wrapper for doing LbyL VI}
\usage{
wrapper_LbyL(spc, biochemphy, resp_fun, stat_fun, isSym = TRUE)
}
\arguments{
\item{spc}{spc}

\item{biochemphy}{colname of responser in SI}

\item{resp_fun}{function to calc VI, take param \code{b1, b2}, which are
reflectance vector corresponding to band1 and band2}

\item{stat_fun}{function to calc stat result, take param \code{x, y}, which are
numeric vector corresponding to predictor(VI) and responder(value of biochemphy)}

\item{isSym}{is Symmetry? TRUE/FALSE}
}
\value{
df
}
\description{
wrapper for doing LbyL VI
}
