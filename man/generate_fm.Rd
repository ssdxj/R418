% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ulits4caret.R
\name{generate_fm}
\alias{generate_fm}
\title{all bands model formula generator}
\usage{
generate_fm(df, yname, model = "B_")
}
\arguments{
\item{df}{the df with columns [yname, B_400, B_401, ...]}

\item{yname}{the target variable column name}
}
\value{
a formula object
}
\description{
all bands model formula generator
}
