% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ulits4tidymodels.R
\name{spc_add_inTrain}
\alias{spc_add_inTrain}
\title{add a inTrain column into spc's SI}
\usage{
spc_add_inTrain(spc, biochemphy, p = 2/3, groups = 7, seed = 324)
}
\arguments{
\item{spc}{the speclib obj}

\item{biochemphy}{the target var}

\item{p}{for createDataPartion}

\item{groups}{for createDataPartion}

\item{seed}{random seed}
}
\value{
new speclib lib
}
\description{
add a inTrain column into spc's SI
}
