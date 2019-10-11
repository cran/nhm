\name{initialprob.nhm}
\alias{initialprob.nhm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute the initial probability vector from a fitted nhm model
}
\description{
Outputs the vector of initial state probabilities in a misclassification type hidden Markov multi-state model fitted using \code{\link{nhm}}.
}
\usage{
initialprob.nhm(object, covvalue=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
 Fitted model object produced using \code{\link{nhm}}.
}
  \item{covvalue}{
 Optional vector of covariate vectors (should be given in the order specified in the \code{covariate} option in \code{nhm}). If omitted the function will use the mean values of the covariates.
}
}
\details{
The \code{initp_nhm} function used to fit the model is called to obtain the values of the initial state vector at the supplied times for the supplied covariate value.}
\value{
Returns a list containing a vector of initial state probabilities and a corresponding vector of standard errors computed using the delta method.
}
\author{
Andrew Titman \email{a.titman@lancaster.ac.uk}
}

\seealso{
\code{\link{nhm}}, \code{\link{ematrix.nhm}}
}
