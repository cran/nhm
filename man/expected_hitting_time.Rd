\name{expected_hitting_time}
\alias{expected_hitting_time}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute the estimated hitting time for a state of a progressive multi-state model.
}
\description{
Uses the estimated state occupation probabilities from a fitted multi-state model using \code{\link{nhm}} to calculate the expected time of entry into a given state, conditional on reaching that state.
}
\usage{
expected_hitting_time(model, state, covvalue=NULL, tstart=0, tmax=NULL,
initp = NULL, npt = 500,  ltrunc = NULL, rtol = 1e-06,  atol = 1e-06, ci = TRUE, 
sim = FALSE, mode = "main", B = 1000, coverage = 0.95)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
 Fitted model object produced using \code{\link{nhm}}.
}
  \item{state}{
The state of interest for the hitting time.
}
  \item{covvalue}{
 Vector of covariate vectors (should be given in the order specified in the \code{covariate} option in \code{nhm}). If omitted the function will use the mean values of the covariates.
}
  \item{tstart}{
 Initial time from which to measure the hitting time and calculate initial state probabilities. Defaults to 0 if omitted.
}
  \item{tmax}{
 Maximum time to integrate over to estimate the mean hitting time. The expected hitting time is effectively conditional on reaching the state before \code{tmax}. Defaults to the maximum observed follow-up time in the data.
}
  \item{initp}{
 optional vector of initial state occupation probabilities. If NULL then will use the estimates from the model. If original model was left-truncated will assume probabilities at \code{tstart} correspond to those implied by the left-truncation model. 
 If \code{ltrunc} supplied will similarly calculate based on left-truncation from the value of t0 supplied. If model does not include misclassification, will assume entry in state 1.
}
  \item{npt}{
 Number of points over which to evaluate state occupation probabilities in order to numerically approximate the integrals via a trapezium rule.
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
}
\details{
Uses the estimated state occupation probabilities and transition intensities from a fitted multi-state model using \code{\link{nhm}} to calculate the expected time of entry into a particular state, conditional on entry.

The method assumes the state of interest can only be entered once and will give incorrect results otherwise. 

The state occupation probabilities are calculated in the same way as in \code{\link{state_occupation_probability.nhm}}. If the estimated \code{initp} has a non-zero probability of being in the hitting state (or states only reachable via that state), then the function computes the expected hitting time from \code{tstart} conditional on not having not yet reached that state.
}
\value{
Returns a list with the following components:
\item{est}{
 Estimated expected hitting time.
}
\item{var_est}{
 Variance of the expected hitting time estimate.
}
\item{var_lest}{
 Variance of the log of the expected hitting time.
}
\item{est_low}{
 Lower limit of symmetric confidence interval
 }
\item{est_high}{
Upper limit of symmetric confidence interval
}
\item{est_low2}{
 Lower limit of confidence interval based on a log-transformation.
 }
\item{est_high2}{
 Upper limit of confidence interval based on a log-transformation.
 } 
 \item{est0}{
 Estimated unconditional probability of reaching the state.
 }
 \item{var_est0}{
 Variance of the estimate of the unconditional probability of reaching the state.
 }
 \item{initp}{
 Value of the initial probability vector at time \code{tstart} calculated from the model.
 }
 }
\author{
Andrew Titman \email{a.titman@lancaster.ac.uk}
}

\seealso{
\code{\link{state_occupation_probability.nhm}}, \code{\link{expected_hitting_time}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
