% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zipfpeVariance.R
\name{zipfpeVariance}
\alias{zipfpeVariance}
\title{Variance of the Zipf-PE distribution.}
\usage{
zipfpeVariance(alpha, beta, tolerance = 10^(-4))
}
\arguments{
\item{alpha}{Value of the \eqn{\alpha} parameter (\eqn{\alpha > 3}).}

\item{beta}{Value of the \eqn{\beta} parameter (\eqn{\beta \in (-\infty, +\infty)}).}

\item{tolerance}{Tolerance used in the calculations. (default = \eqn{10^{-4}})}
}
\value{
A positive real value corresponding to the variance of the distribution.
}
\description{
Computes the variance of the Zipf-PE distribution for given values of \eqn{\alpha} and \eqn{\beta}.
}
\details{
The variance of the distribution only exists for \eqn{\alpha} strictly greater than 3.
}
\examples{
zipfpeVariance(3.5, 1.3)
}
\seealso{
\code{\link{zipfpeMoments}}, \code{\link{zipfpeMean}}.
}
