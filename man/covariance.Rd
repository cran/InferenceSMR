\name{covariance}
\alias{covariance}
\keyword{internal}
\title{
Covariance between two survival probability estimates.
}
\description{
The function \code{covariance} calculates the covariance between two survival probabilty 
obtained with a Cox proportional hazard model. Follow-up times and covariates
may differ for the two survival probabilties. This is mostly an internal function, meant to be used by the function \code{var.expDeath}.
}
\usage{
covariance(t_i, delta_i, beta, time1, z1, time2, z2, S0_i, S1_i, S02_i, Omega_1)
}
\arguments{
  \item{t_i}{
All the follow up times.
}
  \item{delta_i}{
Indicator that the individual was censored: 1 if the individual was censored, 0 otherwise.
}
  \item{beta}{
The coefficients of the Cox model.
}
  \item{time1}{
First survival time.
}
  \item{z1}{
Covariates for first survival.
}
  \item{time2}{
Second survival time.
}
  \item{z2}{
Covariates for second survival.
}
  \item{S0_i}{
Object described in Lin, Flemming and Wei (1994), slightly modified for computationnal efficiency. It is obtained with \code{function1}.
}
  \item{S1_i}{
Object described in Lin et al. (1994), slightly modified for computationnal efficiency. It is obtained with \code{function1}.
}
  \item{S02_i}{
Object described in Lin et al. (1994), slightly modified for computationnal efficiency. It is obtained with \code{function1}.
}
  \item{Omega_1}{
Inverse of the asymptotic variance-covariance matrix of the estimated parameters of the Cox model. It is obtained with \code{function1}.
}
}
\details{
This function was built to be an internal function, only used by the function \code{var.expDeath}.
}
\value{
The covariance between the two estimated survival probability.
}
\references{
Lin D.Y., Fleming T.R., Wei L.J. (1994) \emph{Confidence bands for survival curves under the proportional hazards model}, Biometrika, 81 (1): 73-81.

Talbot, D., Duchesne, T., Brisson, J., Vandal, N. (2011) \emph{Variance estimation and confidence intervals for the standardized mortality ratio with application to the assessment of a cancer screening program}, Statistics in Medicine, 30, 3024-3037. 
}
\author{
Denis Talbot, Thierry Duchesne, Jacques Brisson, Nathalie Vandal.
}

\seealso{
\code{\link{var.expDeath}}, \code{\link{function1}}}

