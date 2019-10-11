\name{example_data2}
\alias{example_data2}
\docType{data}
\title{
Example of data on a progressive 4 state process with state misclassification
}
\description{
The observed states and associated observation times for 1000 patients simulated from a 4 state process non-homogeneous Markov model with misclassification to adjacent transient states.}
\usage{data("example_data1")}
\format{
  A data frame with 3864 rows and 5 variables:
\describe{
   \item{state}{Observed state at the time of observation}
   \item{time}{Time at which the observation occurred}
   \item{id}{Patient identification number}
   \item{cov1}{Binary covariate}
   \item{cov2}{Continuous covariate}
 }
}
\keyword{datasets}
