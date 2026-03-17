\name{pnardl}
\alias{pnardl}
\title{Panel Nonlinear ARDL Estimation}
\usage{
pnardl(
  data,
  depvar,
  lr,
  sr = character(0),
  asymmetric,
  id = "id",
  time = "time",
  estimator = c("pmg", "mg", "dfe"),
  p = 1L,
  asy_test = TRUE,
  multip = 0L,
  level = 0.95
)
}
\arguments{
\item{data}{A data frame containing all variables plus \code{id} and
\code{time} columns.}

\item{depvar}{Character. Name of the dependent variable.}

\item{lr}{Character vector. Names of the long-run regressors (levels).}

\item{sr}{Character vector. Names of the additional short-run regressors
(first-differences). These are included as \eqn{\Delta x} in the ECM.}

\item{asymmetric}{Character vector. Names of variables to decompose into
positive and negative partial sums. Must be a subset of \code{lr}.}

\item{id}{Character. Name of the panel identifier column in \code{data}.}

\item{time}{Character. Name of the time identifier column in \code{data}.}

\item{estimator}{Character. One of \code{"pmg"} (default), \code{"mg"}, or
\code{"dfe"}.}

\item{p}{Integer. Number of lags of \eqn{\Delta y} in the short-run
equation. Default is \code{1}.}

\item{asy_test}{Logical. If \code{TRUE} (default), perform Wald tests for
long-run and short-run asymmetry.}

\item{multip}{Integer. Number of periods for cumulative dynamic multipliers.
Set to \code{0} (default) to skip.}

\item{level}{Numeric. Confidence level for intervals (e.g., \code{0.95}).}
}
\value{
An object of class \code{"pnardl"}, which is a list with components:
  \item{coefficients}{Named numeric vector of pooled long-run and
    short-run coefficients.}
  \item{se}{Standard errors of \code{coefficients}.}
  \item{tstat}{t-statistics.}
  \item{pval}{p-values.}
  \item{panel_coefs}{List of per-panel estimated coefficient vectors.}
  \item{ecm_coefs}{Named vector of ECT (speed-of-adjustment) coefficients
    per panel.}
  \item{asym_tests}{Data frame of Wald test results for asymmetry (if
    \code{asy_test = TRUE}).}
  \item{multipliers}{Data frame of cumulative dynamic multipliers (if
    \code{multip > 0}).}
  \item{partial_sums}{Data frame with the generated partial-sum variables
    appended to the original data.}
  \item{estimator}{Character. Estimator used.}
  \item{depvar}{Character. Dependent variable name.}
  \item{asymmetric}{Character vector. Decomposed variable names.}
  \item{pos_vars}{Character vector. Names of positive partial-sum columns.}
  \item{neg_vars}{Character vector. Names of negative partial-sum columns.}
  \item{call}{The matched call.}
}
\description{
Estimates a Panel Nonlinear ARDL (PNARDL) model following Shin, Yu and
Greenwood-Nimmo (2014). Regressors listed in \code{asymmetric} are
decomposed into positive (\eqn{x^+}) and negative (\eqn{x^-}) partial
sums. The model is then estimated using within-group ordinary least squares
for each panel unit, and pooled estimates (PMG/MG/DFE-style) are computed.
}
\examples{
\donttest{
set.seed(42)
N <- 5; TT <- 30
df <- data.frame(
  id   = rep(1:N, each = TT),
  time = rep(1:TT, N),
  y    = rnorm(N * TT),
  x    = rnorm(N * TT)
)
res <- pnardl(data = df, depvar = "y", lr = c("y", "x"),
              sr = character(0), asymmetric = "x",
              id = "id", time = "time", p = 1, asy_test = FALSE)
print(res)
}
}
\references{
Shin, Y., Yu, B., & Greenwood-Nimmo, M. (2014). Modelling asymmetric
cointegration and dynamic multipliers in a nonlinear ARDL framework. In
R. C. Sickles & W. C. Horrace (Eds.), \emph{Festschrift in Honor of Peter
Schmidt} (pp. 281-314). Springer. \doi{10.1007/978-1-4899-8008-3_9}
}
