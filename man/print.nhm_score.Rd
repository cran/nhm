\name{print.nhm_score}
\alias{print.nhm_score}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Print output from a score test of a nhm object
}
\description{
Print output from a score test based on parameters supplied to \code{\link{nhm}} with \code{score_test=TRUE} specified.
}
\usage{
\method{print}{nhm_score}(x, which_comp = NULL, test_name = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An object of class \code{nhm_code} generated using \code{\link{nhm}}.
}
  \item{which_comp}{
Optional vector to specify which of the parameters are to be tested. If omitted, the function will assume all parameters governing non-homogeneity are to be tested.

Must be supplied if \code{type='bespoke'} was specified when creating the object.
}
  \item{test_name}{
Optional vector to change the heading of the test output.
  }
  \item{\dots}{
Other parameters to be supplied. Currently ignored.
}
}
\details{
The function provides usable output from specifying \code{score_test=TRUE} when using \code{\link{nhm}}. It is most useful to provide a quick(er) test of whether there may be non-homogeneity in a specific model. Note that the model assumes the initial parameters correspond to the constrained maximum likelihood estimate (for instance a model with all the parameters relating to time homogeneity).

The method can be used to compute the local score tests of homogeneity proposed by de Stavola (1988) if \code{type="gompertz"} is specified in \code{\link{nhm}}.

The method can also be used to test for possible informative observation times through the \code{inform} argument in \code{\link{model.nhm}}.

If \code{fisherscore=TRUE} in \code{\link{nhm}} then the expected Fisher information is used. Otherwise, the empirical mean of the squared gradient terms (as used in the BHHH algorithm) is used to estimate the information.

}
\value{
Prints the results of a score test.
}
\references{
de Stavola BL. Testing Departures from Time Homogeneity in Multistate Markov Processes. \emph{Journal of the Royal Statistical Society: Series C (Applied Statistics)} 1988. 37. 242-250.}
\author{Andrew Titman \email{a.titman@lancaster.ac.uk}}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{nhm}}
}

