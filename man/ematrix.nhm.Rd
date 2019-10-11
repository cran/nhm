\name{ematrix.nhm}
\alias{ematrix.nhm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute the misclassification probability matrix from a fitted nhm model
}
\description{
Outputs the matrix of misclasification probabilities in a misclassification type hidden Markov multi-state model fitted using \code{\link{nhm}}.
}
\usage{
ematrix.nhm(object, covvalue=NULL)
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
The \code{emat_nhm} function used to fit the model is called to obtain the values of the misclassification probabilities at the supplied times for the supplied covariate value.}
\value{
Returns a list containing a matrix of misclassification probabilities and a matrix of corresponding standard errors computed using the delta method.
}
\author{
Andrew Titman \email{a.titman@lancaster.ac.uk}
}

\seealso{
\code{\link{nhm}}, \code{\link{plot.nhm}}, \code{\link{predict.nhm}}, \code{\link{qmatrix.nhm}}
}
