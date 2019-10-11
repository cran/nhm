\name{predict.nhm}
\alias{predict.nhm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute state occupation or transition probabilities from a fitted nhm model
}
\description{
Outputs the transition probabilites from a non-homogeneous Markov or misclassification type hidden Markov multi-state model fitted using \code{\link{nhm}}.
}
\usage{
\method{predict}{nhm}(object, time0=0, state0=1, times=NULL, covvalue=NULL,
ci=TRUE, sim=FALSE, coverage=0.95, B=1000, rtol=1e-6,
atol=1e-6, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
 Fitted model object produced using \code{\link{nhm}}.
}
  \item{time0}{
 Starting time from which to compute the transition probabilities. Defaults to 0.
}
  \item{state0}{
 Starting state from which to compute the transition probabilities. Defaults to 1.
}
  \item{times}{
 Optional vector of times at which to compute the transition probabilities. If omitted, the probabilities will be computed at a sequence of times from \code{time0} to the maximum observed time in the data.
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
  \item{rtol}{
 Relative tolerance parameter to be used by \code{lsoda} when solving the differential equations
}
  \item{atol}{
 Absolute tolerance parameter to be used by \code{lsoda} when solving the differential equations}
 \item{...}{
 Other items to be passed to the function. Currently not used.
}
}
\details{
The same approach as in the main \code{nhm} function of numerically solving the system of differential equations is used to compute transition probabilities based on the maximum likelihood estimates found in \code{nhm} and assuming a specific vector of covariates.

If the simulation delta method approach is specified then the function will generate \code{B} parameter vectors from the asymptotic distribution of the MLE and solve the system of equations for each of them, before finding pointwise percentile bootstrap confidence intervals from them.}
\value{
Returns a list containing the vector of times at which the probabilities are computed, a matrix of probabilities for each state at each of the times. If confidence intervals are requested then the lower and upper limits are also provided.

If transition intensity (as opposed to probability) estimates are required then \code{\link{qmatrix.nhm}} should be used.
}
\references{
Mandel M. Simulation-based confidence intervals for functions with complicated derivatives. 2013. \emph{The American Statistician}, 67. 76-81.}
\author{
Andrew Titman \email{a.titman@lancaster.ac.uk}
}

\seealso{
\code{\link{nhm}}, \code{\link{plot.nhm}}, \code{\link{qmatrix.nhm}}
}
