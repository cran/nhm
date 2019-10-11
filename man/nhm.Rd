\name{nhm}
\title{Fit a non-homogeneous Markov model using maximum likelihood}
\alias{nhm}
\description{
   Fit a continuous-time Markov or hidden Markov multi-state model
  by maximum likelihood. Observations of the process
  can be made at arbitrary times, or the exact times of
  transition between states can be known.
  Covariates can be fitted to the Markov chain transition intensities or
  to the hidden Markov observation process.
}
\usage{
nhm(model_object, initial=NULL, gen_inits=FALSE,
control, score_test=FALSE, fixedpar=NULL)
}
\arguments{
  \item{model_object}{Model object created using \code{\link{model.nhm}}}

  \item{initial}{Vector of initial parameter values}

  \item{gen_inits}{If \code{TRUE}, then initial values for the
    transition intensities are generated automatically using the method
    in \code{crudeinits.msm} from the \pkg{msm} package.  This is not available for models with misclassified states.
  If \code{FALSE} a BHHH algorithm implemented using \pkg{maxLik} is used.}
  \item{control}{Object of class \code{nhm.control} specifying various settings for the solution of the KFEs and the optimization. See \code{\link{nhm.control}} for default settings.
}
  \item{score_test}{If \code{TRUE} just the gradient and Fisher information at the supplied values will be computed to allow score tests to be performed.}
    \item{fixedpar}{Numerical vector indicating which parameters are taken as fixed at the value specified by \code{initial}.}
  }
\value{

By default returns an object of class \code{nhm} containing model output data such as the estimated parameters, maximized likelihood value, information matrix etc. The object can be used with \code{print}, \code{predict}, \code{plot} and \code{anova}.

If \code{score.test=TRUE} then returns an object of class \code{nhm_score}. See \code{\link{print.nhm_score}} for more details.
}

\details{
  For more details about the methodology behind the \pkg{nhm} package, see Titman (2011) and the package vignette.

}

\examples{
### Example dataset
### For further examples, see the vignette
trans <- rbind(c(0,1,0,0),c(0,0,2,0),c(0,0,0,3),rep(0,4))
nonh <- rbind(c(0,1,0,0),c(0,0,2,0),c(0,0,0,3),rep(0,4))
gomp_model <- model.nhm(state~time, data=example_data1, subject = id,
                        type="gompertz",trans=trans,nonh=nonh)
initial_val <- c(-0.65,-0.45,-0.55,0,0,0)
gomp_fit <- nhm(gomp_model,initial=initial_val,control=nhm.control(obsinfo=FALSE))
gomp_fit
plot(gomp_fit)
plot(gomp_fit,what="intensities")
}

\author{Andrew Titman \email{a.titman@lancaster.ac.uk}}


\references{
 Titman AC. Flexible Nonhomogeneous Markov Models for Panel Observed Data. \emph{Biometrics}, 2011. 67, 780-787.
}


\seealso{
  \code{\link{model.nhm}}, \code{\link{nhm.control}}, \code{\link{plot.nhm}}, \code{\link{predict.nhm}}, \code{\link{print.nhm_score}}
}


\keyword{models}
