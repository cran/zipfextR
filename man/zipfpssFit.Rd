% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zipfpssFit.R
\name{zipfpssFit}
\alias{zipfpssFit}
\alias{residuals.zipfpssR}
\alias{fitted.zipfpssR}
\alias{coef.zipfpssR}
\alias{plot.zipfpssR}
\alias{print.zipfpssR}
\alias{summary.zipfpssR}
\alias{logLik.zipfpssR}
\alias{AIC.zipfpssR}
\alias{BIC.zipfpssR}
\title{Zipf-PSS parameters estimation.}
\usage{
zipfpssFit(data, init_alpha = NULL, init_lambda = NULL, level = 0.95,
  isTruncated = FALSE, ...)

\method{residuals}{zipfpssR}(object, isTruncated = FALSE, ...)

\method{fitted}{zipfpssR}(object, isTruncated = FALSE, ...)

\method{coef}{zipfpssR}(object, ...)

\method{plot}{zipfpssR}(x, isTruncated = FALSE, ...)

\method{print}{zipfpssR}(x, ...)

\method{summary}{zipfpssR}(object, isTruncated = FALSE, ...)

\method{logLik}{zipfpssR}(object, ...)

\method{AIC}{zipfpssR}(object, ...)

\method{BIC}{zipfpssR}(object, ...)
}
\arguments{
\item{data}{Matrix of count data in form of table of frequencies.}

\item{init_alpha}{Initial value of \eqn{\alpha} parameter (\eqn{\alpha > 1}).}

\item{init_lambda}{Initial value of \eqn{\lambda} parameter (\eqn{\lambda > 0}).}

\item{level}{Confidence level used to calculate the confidence intervals (default 0.95).}

\item{isTruncated}{Logical; if TRUE, the truncated version of the distribution is returned.(default = FALSE)}

\item{...}{Further arguments to the generic functions. The extra arguments are passing
to the \emph{\link{optim}} function.}

\item{object}{An object from class "zpssR" (output of \emph{zipfpssFit} function).}

\item{x}{An object from class "zpssR" (output of \emph{zipfpssFit} function).}
}
\value{
Returns a \emph{zpssR} object composed by the maximum likelihood parameter estimations jointly
with their standard deviation and confidence intervals and the value of the log-likelihood at the
maximum likelihood estimator.
}
\description{
For a given sample of strictly positive integer numbers,  usually of the type of ranking data or
frequencies of frequencies data, estimates the parameters of the Zipf-PSS distribution by means of
the maximum likelihood method. The input data should be provided as a frequency matrix.
}
\details{
The argument \code{data} is a two column matrix with the first column containing the observations and
the second column containing their frequencies.

The log-likelihood function is equal to:
\deqn{l(\alpha, \lambda, x) = \sum_{i =1} ^{m} f_a(x_i)\, log(P(Y = x_i)),}
where \eqn{m} is the number of different values in the sample, being \eqn{f_{a}(x_i)} is the absolute
frequency of \eqn{x_i}.The probabilities are calculated applying the Panjer recursion.
By default the initial values of the parameters are computed using the function \code{getInitialValues}.
The function \emph{\link{optim}} is used to estimate the parameters.
}
\examples{
data <- rzipfpss(100, 2.5, 1.3)
data <- as.data.frame(table(data))
data[,1] <- as.numeric(as.character(data[,1]))
data[,2] <- as.numeric(as.character(data[,2]))
initValues <- getInitialValues(data, model='zipfpss')
obj <- zipfpssFit(data, init_alpha = initValues$init_alpha, init_lambda = initValues$init_lambda)
}
\references{
{
Panjer, H. H. (1981). Recursive evaluation of a family of compound
distributions. ASTIN Bulletin: The Journal of the IAA, 12(1), 22-26.

Sundt, B., & Jewell, W. S. (1981). Further results on recursive evaluation of
compound distributions. ASTIN Bulletin: The Journal of the IAA, 12(1), 27-39.
 }
}
\seealso{
\code{\link{getInitialValues}}.
}
