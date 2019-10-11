\name{plot.nhm}
\alias{plot.nhm}
\title{
Plot transition probabilities or intensities from a fitted nhm model.
}
\description{
Produces plots of the transition probabilites or intensities from a non-homogeneous Markov or misclassification type hidden Markov multi-state model fitted using \code{\link{nhm}}.
}
\usage{
\method{plot}{nhm}(x, what="probabilities",time0=0, state0=1, times=NULL,
covvalue=NULL, ci=TRUE, sim=FALSE, coverage=0.95, B=1000, rtol=1e-6,
atol=1e-6, main_arg=NULL, xlab="Time", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
 Fitted model object produced using \code{\link{nhm}}.
}
  \item{what}{
 Character string to indicate what should be plotted. Options are \code{probabilities} (the default which produces transition probabilities) or \code{intensities} (to produce a plot of the intensities)
}
  \item{time0}{
 Starting time from which to compute the transition probabilities or intensities. Defaults to 0.
}
  \item{state0}{
 Starting state from which to compute the transition probabilities. Defaults to 1. Not required for transition intensities
}
  \item{times}{
 Optional vector of times at which to compute the transition probabilities or intensities. If omitted, the probabilities/intensities will be computed at a sequence of times of length 100 from \code{time0} to the maximum observed time in the data.
 }
  \item{covvalue}{
 Optional vector of covariate vectors (should be given in the order specified in the \code{covariate} option in \code{nhm}). 
 If omitted the function will use the mean values of the covariates.
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
 Relative tolerance parameter to be used by \code{lsoda} when solving the differential equations for the transition probabilites.
}
  \item{atol}{
 Absolute tolerance parameter to be used by \code{lsoda} when solving the differential equations for the transition probabilites.}
  \item{main_arg}{
 Character string specifying beginning of title to be given to each of the plot panes generated.
}
  \item{xlab}{
 Character string specifying x-axis label to be given to each plot.
}
\item{...}{
 Other items to be passed to the function. Currently not used.
}
}
\details{
Computation is performed by calling \code{\link{predict.nhm}}, for the transition probabilities, or \code{qmatrix.nhm} for the intensities (see for more details).
}
\value{
Generates a multi-pane plot for each state. If values are required they can be obtained using \code{\link{predict.nhm}}.
}
\references{
Mandel M. Simulation-based confidence intervals for functions with complicated derivatives. 2013. \emph{The American Statistician}, 67. 76-81.
}
\author{
Andrew Titman \email{a.titman@lancaster.ac.uk}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{nhm}}, \code{\link{predict.nhm}}, \code{\link{qmatrix.nhm}}
}