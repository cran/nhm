\name{qmatrix.nhm}
\alias{qmatrix.nhm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute transition intensities from a fitted nhm model
}
\description{
Outputs the transition intensities from a non-homogeneous Markov or misclassification type hidden Markov multi-state model fitted using \code{\link{nhm}}.
}
\usage{
qmatrix.nhm(object, time0=0, times=NULL, covvalue=NULL, ci=TRUE, sim=FALSE,
coverage=0.95, B=1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
 Fitted model object produced using \code{\link{nhm}}.
}
  \item{time0}{
 Starting time from which to compute the transition intensities. Defaults to 0.
}
  \item{times}{
 Optional vector of times at which to compute the transition intensities. If omitted, the intensities will be computed at a sequence of times from \code{time0} to the maximum observed time in the data.
}
  \item{covvalue}{
 Optional vector of covariate vectors (should be given in the order specified in the \code{covariate} option in \code{nhm}). If omitted the function will use the mean values of the covariates.
}
  \item{ci}{
 If \code{TRUE} pointwise confidence intervals will be shown in addition to the point estimates.
}
  \item{sim}{
 If \code{TRUE} a simulation Delta method (Mandel, 2013) will be used to calculate the confidence intervals. Otherwise the standard Delta method will be applied.
}
  \item{coverage}{
 Coverage level (should be a value between 0 and 1) for the confidence intervals. Defaults to 0.95.
}
  \item{B}{
 Number of simulations to be performed to compute the simulation Delta method.
}
}
\details{
The \code{intens} function used to fit the model is called to obtain the values of the transition intensities at the supplied times for the supplied covariate value.

If the simulation delta method approach is specified then the function will generate \code{B} parameter vectors from the asymptotic distribution of the MLE and compute the intensities for each of them, before finding pointwise percentile bootstrap confidence intervals from them.}
\value{
Returns a list containing the vector of times at which the intensities are computed, a matrix of probabilities for each state at each of the times. If confidence intervals are requested then the lower and upper limits are also provided.

If transition probability (as opposed to intensity) estimates are required then \code{\link{predict.nhm}} should be used.
}
\references{
Mandel M. Simulation-based confidence intervals for functions with complicated derivatives. 2013. \emph{The American Statistician}, 67. 76-81.}
\author{
Andrew Titman \email{a.titman@lancaster.ac.uk}
}

\seealso{
\code{\link{nhm}}, \code{\link{plot.nhm}}, \code{\link{predict.nhm}}
}
