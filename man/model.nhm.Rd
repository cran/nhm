\name{model.nhm}
\alias{model.nhm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Model object set up for non-homogeneous Markov models
}
\description{
Sets up a model object in preparation for fitting a non-homogeneous Markov or misclassification type hidden Markov multi-state model.
}
\usage{
model.nhm(formula,  data, subject, covariates=NULL,  type, trans,
nonh=NULL, covm=NULL, centre_time=NULL, emat=NULL, ecovm=NULL,
firstobs=NULL, initp=NULL, initp_value=NULL, initcovm=NULL,
splinelist=NULL,degrees=NULL,censor=NULL,
censor.states=NULL,death=FALSE,death.states=NULL,intens=NULL,
inform = NULL, ltruncation_states=NULL, ltruncation_time=NULL,
ltruncation_entry=NULL, phasemap=NULL)
}

\arguments{
  \item{formula}{
  A formula identifying the state and time variables within \code{data}, for instance \code{state ~ time} would imply the variables are \code{state} and \code{time}, respectively.
  }
  \item{data}{
data frame containing the observed states, observation times, subject identifiers and covariates. Should include initial observation/recruitment times.
}
\item{subject}{
Name of the subject identifier variable within the \code{data} data frame.
}
  \item{covariates}{
  A character vector giving the variable names of the covariates to be used in the model
}
  \item{type}{
  type of intensity model.
  \code{'bespoke'}: user supplied,
  \code{'weibull'}: Model with Weibull transition intensity functions with respect to time.
  \code{'gompertz'}: Gompertz/exponential growth intensity models.
  \code{'bspline'}: b-spline function of time model.
}
  \item{trans}{
  Square matrix of viable transitions with dimension equal to the number of states. Impossible transitions should be 0. Others should be labelled consecutively from 1. Labelling transitions with the same value assumes the parameter is shared. Required even if \code{type='bespoke'}.
}
  \item{nonh}{
  Square matrix to indicate non-homogeneous transitions with dimension equal to the number of states. Impossible transitions or homogeneous transitions should be 0. Otherwise label consecutively from 1. Labelling the same value implies the same non-homogeneity. Not required if \code{type='bespoke'}. If otherwise omitted a time homogeneous model is fitted.
  }
  \item{covm}{
  Either a named list of nstate x nstates indicating the covariate effects with respect to a particular covariate OR
  an nstate x nstate x ncov array to indicate covariate effects, where ncov is the length of the supplied \code{covariates} vector. 0 implies no covariate effect. Otherwise label consecutively from 1. Labelling the same value implies a common covariate effect. Not required if \code{type='bespoke'}.
}
  \item{centre_time}{
 Value by which to centre time for Gompertz models. By default the model is of the form \eqn{h(t) = exp(a + bt)}, centring  by \eqn{c} reparametrizes this to \eqn{h(t) = exp(a + b(t - c))}. Centring can improve the convergence of optimization routines.
}
  \item{emat}{
  Either a square matrix of viable misclassification errors which must be supplied if the model has misclassification. Impossible errors should be 0. Others should be labelled consecutively. Labelling the same implies a common parameter on the logit scale. Alternatively, can supply a function with arguments (z, x) that returns the square matrix of misclassification probabilities. See vignette for more details.
}
  \item{ecovm}{
    Either a named list of nstate x nstates indicating the covariate effects with respect to a particular covariate OR
  an nstate x nstate x ncov array to indicate indicate covariate effects on misclassification, where ncov is the length of the supplied \code{covariates} vector. 0 implies no covariate effect. Otherwise label consecutively from 1. Labelling the same value implies a common covariate effect.
}
  \item{firstobs}{
   For misclassification models: Form of the first observation for each subject in the data.
   \code{'exact'}: Initial state not subject to misclassification (default)
   \code{'absent'}: No initial state. First observation is ignored and state occupied is based on initial probabilities model.
   \code{'misc'}: Initial state is subject to misclassification.
}
 \item{initp}{
   For misclassification models: Numerical vector of length nstate to define the model for the initial probabilities. The first entry should be zero. Should be numbered consecutively. If the same number is repeated implies a shared parameter. If absent then initial probabilities taken from \code{initp_value}.
 }
 \item{initp_value}{
 For misclassification models where \code{firstobs="absent"} or \code{"misc"}: Fixed value of initial probabilities is missing. Should be a numerical vector of length nstate. Ignored if \code{initp} is present. Default if absent is \code{c(1,0,...)}.
 }
 \item{initcovm}{
 For misclassification models; Either a named list of vectors of length nstate, or an nstate x ncovs matrix to specify the covariate effects on misclassification probabilities. 0 implies no covariate effect. Otherwise label consecutively from 1. Labelling the same value implies a common covariate effect.
 }
  \item{splinelist}{
  For bspline models only: list (of length equal to the number of nonhomogeneous transitions) of knot point locations including the boundary knots.
}
  \item{degrees}{
  For bspline models only: optional vector (of length equal to number of nonhomogeneous transitions) of degrees of splines. Defaults to 3 if not specified.
}
  \item{censor}{
  Vector of censor state indicators in the data. Note that censored observations can only occur as the last observation for a subject unless it is a model with misclassification.
}
  \item{censor.states}{
  List of vectors of states in which subject occupy if censored by corresponding censor state indicator. Can be a vector if only one censor state marker is present.
}
  \item{death}{
 Setting \code{TRUE} assumes exact death times are present in the data set
}
  \item{death.states}{
 Vector specifying which states have exact death times. Should only correspond to absorbing states.
}
  \item{intens}{

  Optional supplied intensity function. See below for details.
}
  \item{inform}{
  Square matrix to indicate which transitions to include in a score test for informative observation. Should only be used as part of a score test. 0 implies
  transition not tested. Otherwise label consecutively from 1. Labelling the same value implies assuming common effect. Not required if \code{type='bespoke'}.
  }
  \item{ltruncation_states}{
  If left truncation (delayed entry) exists in the dataset, a vector to specify the subset of states
  possible at recruitment into the dataset. Can be omitted if no left truncation present. Only
  applicable to misclassification models with \code{firstobs='misc'}.
  }
  \item{ltruncation_time}{
  The origin time for left truncation. This is the time (age) at which \code{initp} will be applied.
  If omitted, when \code{ltruncation_states} is given then a value of 0 is assumed.
  }
  \item{ltruncation_entry}{
  The observation number in a subject's sequence to which the left truncation condition applies.
  If omitted defaults to 1. Higher integer values possible if, for instance, data contains only
  subjects with at least \code{m} longitudinal measurements.
  }
  \item{phasemap}{
  Optional vector to specify relationship between observed and latent states to allow latent Markov/
  phase-type semi-Markov models. Should be a vector of length equal to the number of latent states
  containing the integers corresponding to the observed states. Can currently only be used in
  conjunction with \code{type='bespoke'}.
  }
}
\details{
The function allows the model to be specified and creates the metadata needed to use \code{nhm} to fit it. The function automatically generates a function \code{intens} which defines the generator matrix of the model and its first derivatives as a function of time  \code{t}, covariates \code{z} and the underlying parameters \code{x}, provided the model is of Weibull, Gompertz or B-spline type.

Alternatively, \code{type='bespoke'} can be chosen. In which case it is necessary for the user to supply a function \code{intens}. This must have arguments \code{t, z, x} and return a list consisting of a component \code{q} which is the nstate x nstate generator matrix, and \code{dq} which is the nstate x nstate x nparQ first derivatives of the generator matrix with respect to the parameters of the model, where nparQ is the number of parameters in the model for the intensities only (excludes parameters for the emat or initp). Since, by default, unconstrained optimization is used, the parameters should usually take values on \code{-Inf, Inf}. Note that using a hard-coded version via \code{type='bespoke'} can be substantially faster than the analogous automatically generated function, so for large models or datasets it may be advantageous to code directly. See the vignette for further details.

For misclassification type models, the function also automatically creates functions \code{emat_nhm} and \code{initp_nhm}, to allow the misclassification probability matrix and the initial probability vectors and their derivatives to be calculated at given parameter and covariate values. In each case, a multinomial logistic regression is used for the covariate model. User specification of the initial probability vector is not currently possible.




}

\value{

Returns an object of class \code{nhm_model} containing the necessary metadata needed to use \code{nhm} to fit the model.

}


\author{Andrew Titman \email{a.titman@lancaster.ac.uk}}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{nhm}}
}

