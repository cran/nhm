\name{state_occupation_probability.nhm}
\alias{state_occupation_probability.nhm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute state occupation probabilities from a fitted nhm model
}
\description{
Outputs the state occupation probabilities from a non-homogeneous Markov or misclassification type hidden Markov multi-state model fitted using \code{\link{nhm}}.}
\usage{
state_occupation_probability.nhm(model, covvalue=NULL, time0 = 0, times = NULL, 
initp = NULL, ltrunc = NULL,  rtol = 1e-06, atol = 1e-06, ci = TRUE, sim = FALSE, 
mode = "main", B = 1000, coverage = 0.95, statemerge = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
 Fitted model object produced using \code{\link{nhm}}.
}
  \item{covvalue}{
 Vector of covariate vectors (should be given in the order specified in the \code{covariate} option in \code{nhm}). If omitted the function will use the mean values of the covariates.
 }
  \item{time0}{
 Starting time from which to compute the state occupation probabilities. Defaults to 0.
}
  \item{times}{
Optional vector of times at which to compute the transition probabilities. If omitted, the probabilities will be computed at a sequence of times from \code{time0} to the maximum observed time in the data.
}
  \item{initp}{
 optional vector of initial state occupation probabilities. If NULL then will use the estimates from the model. If original model was left-truncated will assume probabilities at \code{tstart} correspond to those implied by the left-truncation model. 
 If \code{ltrunc} supplied will similarly calculate based on left-truncation from the value of t0 supplied. If model does not include misclassification, will assume entry in state 1.
}
  \item{ltrunc}{
 Optional list containing \code{ltruncation_time} and \code{ltruncation_states}. If supplied will replace the values in the original model fit object. 
}
  \item{rtol}{
 Relative tolerance parameter to be used by \code{lsoda} when solving the differential equations
}
  \item{atol}{
Absolute tolerance parameter to be used by \code{lsoda} when solving the differential equations
}
  \item{ci}{
 Logical for whether confidence intervals should be calculated for the quantities.
}
  \item{sim}{
 Logical for whether simulation-based (parametric bootstrap) confidence intervals should be used (\code{TRUE}) or delta method-based intervals (\code{FALSE}).
}
  \item{mode}{
 Argument for internal use to faciliate parametric bootstrapping: default \code{"main"} ensures standard errors and calculated, if \code{mode="boot"} then standard errors are not calculated.
}
  \item{B}{
 Integer specifying the number of bootstrap replicates to perform if \code{sim=TRUE}. Defaults to 1000.
}
  \item{coverage}{
 Nominal coverage proportion required for confidence intervals. Defaults to 0.95 implying two-sided 95\% confidence intervals.
 }
  \item{statemerge}{
 Logical to determine whether latent states should be merged for models fitted using the \code{phasemap} option. If \code{TRUE} will return estimates for the number of observable states.
}
}
\details{
The state occupation probabilities are computed by summing over the transition probabilites from \code{time0} with respect to the initial state occupation probabilities.

For models with left-truncation, if \code{time0} is greater than the time origin in the model (\code{ltruncation_time}) then the function will calculate the state occupation probabilites at \code{time0} assuming subjects are sampled conditional on being in the \code{ltruncation_states}

}
\value{
  #out <- list(times= times, ests=ests, ests_cov = ests_cov, ests_low=ests_low, ests_high=ests_high, initp=initp)

Returns a list with the following components:
\item{times}{
 The set of times at which the state occupation probabilities are computed.
}
\item{ests}{
 Matrix of state occupation probability estimates
}
\item{ests_cov}{
 Array of variance-covariance estimates of the state occupation probability estimates at each time.
}
\item{ests_low}{
 Lower limits of component-wise confidence intervals
 }
\item{ests_high}{
Upper limits of component-wise confidence intervals
}
 \item{initp}{
 Value of the initial probability vector at time \code{tstart} calculated from the model.
 }
 }
\author{
Andrew Titman \email{a.titman@lancaster.ac.uk}
}

\seealso{
\code{\link{expected_hitting_time}}, \code{\link{state_life_expectancy}}
}