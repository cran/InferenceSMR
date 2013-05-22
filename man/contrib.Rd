\name{contrib}
\alias{contrib}
\title{
Weights calculation
}
\description{
\code{contrib} is a function that finds out how much time was contributed for every combination
of the incidence covariates, the survival covariates, and follow up time.
Using these contibutions makes the calculation of the variance of the expected number
of deaths much more efficient than making the calculations on the raw data.
}
\usage{
contrib(start_follow, end_follow, incid_cov, surv_cov, follow_up, increment)
}
\arguments{

  \item{start_follow}{
A vector that contains the amount of follow-up time elapsed when this set of covariate values started.
}
  \item{end_follow}{
A vector that contains the amount of follow-up time that will have elapsed when this set of covariate values changes or when follow-up ends.
}
  \item{incid_cov}{
A vector that contains the values of the covariates for incidence.
}
  \item{surv_cov}{
A vector that contains the values of the covariates for survival.
}
  \item{follow_up}{
A vector that contains the total follow-up time for that individual.
}
  \item{increment}{
The value of the time increment that will be used for Sasieni's estimator and variance (a numeric).
}
}
\details{

If only time independent covariates are used, then all vectors are of dimension \code{n} = sample size.

Otherwise, the function \code{contrib} needs an augmented dataset, in which each 
person-time corresponds to one row. Each row contains the values of the 
time-invariant covariates and the updated values of time-dependent variables.
Therefore, the covariates must have a finite number of values (they must be discrete).

For example, if only the covariate ``age'' is used for incidence and survival, 
an individual followed for 3.4 years that was 54.6 years at the begining 
of the study should be entered as follow:

\code{#start_follow, end_follow, incid_cov, surv_cov, follow_up} \cr
\code{0   0.4 54 54 3.4} \cr
\code{0.4 1.4 55 55 3.4} \cr
\code{1.4 2.4 56 56 3.4} \cr
\code{2.4 3.4 57 57 3.4}

An example of the code required to do such a thing is provided in the \code{screening} dataset help page.

}
\value{

The resulting list is directly usable in functions \code{est.expDeath}, \code{var.expDeath} and \code{inference.SMR}. The list contains the following elements :

\item{contrib}{A matrix giving the total time contributed for every combination of discretized follow-up times, incidence covariates and survival covariates.}

\item{ncov.incid}{The number of covariates for incidence.}

\item{ncov.surv}{The number of covariates for survival.}

\item{increment}{The increment used for discretization of follow-up times.}

}

\author{
Denis Talbot, Thierry Duchesne, Jacques Brisson, Nathalie Vandal.
}

\seealso{
\code{\link{est.expDeath}}, \code{\link{var.expDeath}}, \code{\link{inference.SMR}}, \code{\link{screening}}}

