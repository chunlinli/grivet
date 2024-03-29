% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/discovery.R
\name{intdag.pmle.diff.aic}
\alias{intdag.pmle.diff.aic}
\title{Matrix V Estimation Function}
\usage{
intdag.pmle.diff.aic(X, Y, taus, gammas)
}
\arguments{
\item{X}{A n*q matrix for intervention variables.}

\item{Y}{A n*p matrix for primary variables.}

\item{taus}{Parameter used for TLP in each regression model where Y[,i] is the response with X the predictor matrix.}

\item{gammas}{Parameter used for TLP in each regression model where Y[,i] is the response with X the predictor matrix.}
}
\value{
An estimate of matrix V.
}
\description{
This function estimates matrix V.
}
\examples{
\dontrun{
to be added.
}
}
\references{
Discovery and inference of a causal network with hidden confounding.
}
