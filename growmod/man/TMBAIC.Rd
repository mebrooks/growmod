% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{TMBAIC}
\alias{TMBAIC}
\title{Calculate the Akaike Information Criteria corrected for small sample size}
\usage{
TMBAIC(opt, n = NULL, correction = TRUE)
}
\arguments{
\item{opt}{an optimized TMB model}

\item{n}{number of observations}

\item{correction}{logical - Correct for small sample size? (recommended as it converges to other formula for large samples)}
}
\description{
Calculate the Akaike Information Criteria corrected for small sample size
}

