\name{state_life_expectancy}
\alias{state_life_expectancy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
State-specific life expectancies and quality-adjusted life years
}
\description{
Uses the estimated state occupation probabilities from a fitted multi-state model using \code{\link{nhm}} to calculate estimated life years spent in each state, or weighted sums of the state-wise life expectancies.
}
\usage{
state_life_expectancy(model, covvalue=NULL, tstart=0, tmax=NULL, initp = NULL, 
 npt = 500, discount = NULL, utilities = NULL, ltrunc = NULL, rtol = 1e-06, 
 atol = 1e-06, ci = TRUE, sim = FALSE, mode = "main", B = 1000, coverage = 0.95)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
 Fitted model object produced using \code{\link{nhm}}.
}
  \item{covvalue}{
 Vector of covariate vectors (should be given in the order specified in the \code{covariate} option in \code{nhm}). If omitted the function will use the mean values of the covariates.
}
  \item{tstart}{
 Initial time from which to integrate over state occupancy probabilities. Defaults to 0.
}
  \item{tmax}{
 Maximum time to integrate over state occupancy probabilities. 
}
  \item{initp}{
 optional vector of initial state occupation probabilities. If NULL then will use the estimates from the model. If original model was left-truncated will assume probabilities at \code{tstart} correspond to those implied by the left-truncation model. 
 If \code{ltrunc} supplied will similarly calculate based on left-truncation from the value of \code{ltruncation_time} supplied. If model does not include misclassification, will assume entry in state 1.
}
  \item{npt}{
 Number of points over which to evaluate state occupation probabilities in order to numerically approximate the integral via a trapezium rule.
}
  \item{discount}{
 Optional discounting function (e.g. if discounted LYs are required). See details below.
}
  \item{utilities}{
 Either an \code{nstate} length vector of state-specific utilities in order to produce a single quality-adjusted life years estimate, or an \code{nstate} x \code{m} matrix of utilities. See details below.
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
The function uses the same methods as \code{\link{state_occupation_probability.nhm}} to obtain state occupation probability estimates at a grid of times and then uses the trapezium rule to get estimates of the expected time spent in each state.

The \code{utilties} argument can be used to either create utility weights for quality-adjusted life years estimates, assuming time-constant state-specific health utilities. Alternatively, multiple summaries of the state-specific life expectancies can be obtained by specifying an \code{nstate} x \code{m} matrix.

The \code{discount} option can be used to estimate discounted life years. \code{discount} should be a function that takes a vector of times and converts it to a vector of corresponding discounting factors, or instead returns a matrix with \code{nstate} columns of state-wise discounting factors. If \code{utilities} is also specified then the discounting will be applied to the utility-weighted quantities. Note that, if the state-specific utilities are time-dependent, this can be accommodated by jointly specifying the discounting and utility within the discount function.

If the fitted \code{model} involved left-truncation and a \code{tstart} value is chosen that is greater than the origin time (\code{ltruncation_time}) in the original model then by default the function will calculate state occupation probabilities assuming a random sample of patients who are sampled conditional on being in one of the \code{ltruncation_states} at \code{tstart}. This can be overridden by either providing an initial state vector for \code{initp} or supplying a different \code{ltrunc} object. 
}
\value{
Returns a list containing
\item{est}{
Estimates of state-wise life expectancy for each of the states
}
\item{est_cov}{
Covariance matrix associated with the state-wise life expectancy estimates
}
\item{est_low}{
Lower limits of the individual confidence intervals corresponding to the entries in \code{est}.
}
\item{est_high}{
Upper limits of the individual confidence intervals corresponding to the entries in \code{est}.
}
\item{qaly_est}{
If utilities are specified, provides the weighted summaries.
}
\item{qaly_var}{
Variance or covariance matrix associated with the weighted summaries.
}
\item{qaly_low}{
Lower limits of the individual confidence interval(s) of the weighted summaries.
}
\item{qaly_high}{
Upper limits of the individual confidence interval(s) of the weighted summaries.
}
\item{initp}{
Value of the initial probability vector at time \code{tstart} calculated from the model.
}
\item{ders}{
First derivative of the state-wise life expectancy estimates with respect to each of the parameter values, evaluated at the fitted MLE.
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
