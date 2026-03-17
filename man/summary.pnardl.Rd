\name{summary.pnardl}
\alias{summary.pnardl}
\title{Summary Method for pnardl Objects}
\usage{
\method{summary}{pnardl}(object, ...)
}
\arguments{
\item{object}{An object of class \code{"pnardl"}.}

\item{...}{Further arguments passed to or from other methods (ignored).}
}
\value{
Invisibly returns a list with components:
  \item{coefficients}{Data frame with estimates, SE, t-stats, and p-values.}
  \item{asym_tests}{Asymmetry test results (may be \code{NULL}).}
  \item{ecm_summary}{Summary statistics of the ECT coefficients.}
}
\description{
Summary Method for pnardl Objects
}
