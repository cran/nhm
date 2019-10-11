\name{nhm.control}
\alias{nhm.control}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Ancillary arguments for controlling nhm fits
}
\description{
This is used to set various logical or numeric parameters controlling a non-homogeneous Markov model fit. Usually to be used within a call to \code{nhm}.
}
\usage{
nhm.control(tmax=NULL, coarsen=FALSE, coarsen.vars=NULL, coarsen.lv=NULL,
checks=FALSE,rtol=1e-6, atol=1e-6, fishscore=NULL, linesearch=FALSE, damped=FALSE,
damppar=0,obsinfo=TRUE,splits=NULL,ncores=1,print.level=2, maxLikcontrol=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tmax}{Optional parameter to set the maximum time to which the Kolmogorov Forward equations should be integrated. Defaults to 1+max(time) if left unspecified.}
  \item{coarsen}{If \code{TRUE} the covariate values will be subjected to coarsening using K-means clustering, so there are fewer unique values. This is useful for large datasets with continuous covariates.}
  \item{coarsen.vars}{Vector of the index of covariates which require coarsening. Must be supplied if \code{coarsen=TRUE}.}
  \item{coarsen.lv}{Number of unique covariate values to which the covariates should be coarsened.}
  \item{checks}{If \code{TRUE} some basic checks will be performed to ensure the accuracy of the supplied \code{intens} function. Mainly useful if a user defined \code{type="bespoke"} intensity function is used for which the default is \code{TRUE}, otherwise default is \code{FALSE}}
  \item{rtol}{Relative error tolerance to be passed to lsoda, default is 1e-6}
  \item{atol}{Absolute error tolerance to be passed to lsoda, default is 1e-6}
  \item{fishscore}{If \code{TRUE} then the Fisher scoring algorithm will be used provided the model has no censoring, exact death times or misclassification. This is generally faster, but less robust than the BHHH algorithm.}
  \item{linesearch}{If \code{TRUE} and \code{fishscore=TRUE} then a line search will be performed to find the best step length in the Fisher scoring algorithm.}
  \item{damped}{If \code{TRUE} the Fisher scoring algorithm will be damped (e.g. Levenberg type algorithm). Useful if some parameters are close to being unidentifiable.}
  \item{damppar}{Numerical damping parameter to be applied if \code{damped=TRUE}}
  \item{obsinfo}{If \code{TRUE} the observed Fisher information will be computed in addition to the expected information when the Fisher scoring algorithm is used. For optimization with \code{\link{maxLik}} the observed Fisher information will be used as the Hessian rather than the squared gradient vectors. }
  \item{splits}{Optional vector of intermediate split times for solving the ODEs. Only needed if P(0,t) becomes singular for some t causing the optimization to stop. Should be a set of consecutive values less than tmax.}
  \item{ncores}{Number of cores to use. 1= no parallelization, 2 or more: Uses \code{\link{mclapply}} when solving ODEs with different covariates patterns.}
  \item{print.level}{For \code{\link{maxLik}} optimization; level of detail to print. Integer from 0 to 3. Defaults to 2.}
  \item{maxLikcontrol}{For \code{\link{maxLik}} optimization; optional list of control parameters to be passed to \code{\link{maxLik}}.}
}
\details{
\code{tmax}, \code{rtol} and \code{atol} refer directly to parameters with the \code{lsoda} function in \code{deSolve} and relate to how the Kolmogorov Forward Equations are numerically solved.

\code{coarsen}, \code{coarsen.vars} and \code{coarsen.lv} are useful in situations where it is computationally infeasible (or unattractive) to compute the exact solution for all covariate patterns. Implements an approximate solution in which the covariates are coarsened using K-means clustering (as proposed in Titman (2011)).

\code{linesearch}, \code{damped}, \code{damppar} are specific to the Fisher scoring algorithm.

Setting \code{obsinfo=TRUE} will tend to give more accurate standard error estimates and gives more opportunity to check for non-convergence of the maximum likelihood procedure.

The option \code{splits} modifies the way in which the transition probabilities are computed. By default, \code{\link{nhm}} solves a single system of differential equations starting from 0 to obtain \eqn{P(0,t)} and then uses inversion of the Chapman-Kolmogorov equation \eqn{P(0,t) =P(0,t_0)P(t_0,t)} to find \eqn{P(t_0,t)} for a given \eqn{t_0 > 0}. In some cases \eqn{P(0,t_0)} will be singular or effectively singular. If a split is specified at \eqn{s} then \code{nhm} will find \eqn{P(t_0,t)} for \eqn{t_0 > t*} by solving the system of equations \eqn{P(t*,t)} where \eqn{t*} is the smallest interval start time greater than or equal to \eqn{s} within the data. If \code{\link{nhm}} fails due to the lack of split times, the error message will advise on the interval in which the split should be introduced. Note that the need for splits can also arise if the initial parameters specified are inappropriate. It may often be better to find more appropriate initial parameter estimates,for instance by fitting the analogous homogeneous model in \code{\link[msm]{msm}}, rather than adding multiple split times.

\code{ncores} allows parallel processing to be used, through the \pkg{parallel} package, to simultaneously solve the systems of differential equations for each covariate pattern. If \code{ncores > 1} then \code{ncores} defines the \code{mc.cores} value in \code{mclapply}. Note that the data needs to include multiple covariate patterns for this to successfully increase computation speed.

}
\value{
A list containing the values of each of the above constants
}
\references{
 Titman AC. Flexible Nonhomogeneous Markov Models for Panel Observed Data. \emph{Biometrics}, 2011. 67, 780-787.
}
\author{Andrew Titman \email{a.titman@lancaster.ac.uk}}

\seealso{
\code{\link{nhm}}
}

