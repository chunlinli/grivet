% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mbtlp.R
\name{cv.MB_Union}
\alias{cv.MB_Union}
\title{Recovery of Precision Matrix's Support}
\usage{
cv.MB_Union(R, tau.list, gamma.list, n.fold)
}
\arguments{
\item{R}{Data matrix consisting of variables whose precision matrix is of interests.}

\item{tau.list}{A set of taus to be selected for TLP.}

\item{gamma.list}{A set pf gammas to be selected for TLP.}

\item{n.fold}{Number of folds in cross-validation.}
}
\value{
A list consits of an estimate of support, selected taus and gammas.
}
\description{
This function select best tunning parameters for TLP in recovering precision matrix's support.
}
\examples{
\dontrun{
to be added.
}
}
\references{
Discovery and inference of a causal network with hidden confounding.
}
