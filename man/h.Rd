\name{h}
\alias{h}
\keyword{internal}
\title{
h function as defined by Lin et al. (1994).
}
\description{
An internal function meant only to be used by the \code{covariance} function.
}
\usage{
h(t_i, delta_i, beta, time, z, S0_i, S1_i)
}
\arguments{
  \item{t_i}{
The follow up times.
}
  \item{delta_i}{
Indicator that the individual was censored: 1 if the individual was censored, 0 otherwise.
}
  \item{beta}{
The coefficients of the Cox model.}
  \item{time}{
Survival time.
}
  \item{z}{
Covariates for survival time.
}
  \item{S0_i}{
Object described in Lin et al. (1994), slightly modified for computationnal efficiency. It is obtained with \code{function1}.
}
  \item{S1_i}{
Object described in Lin et al. (1994), slightly modified for computationnal efficiency. It is obtained with \code{function1}.
}
}
\details{
This is an internal function, only meant to be used by the \code{covariance} function.
}
\value{
A numeric value described in Lin et al. (1994).
}
\references{
Lin D.Y., Fleming T.R., Wei L.J. (1994) \emph{Confidence bands for survival curves under the proportional hazards model}, Biometrika, 81 (1): 73-81.
}
\author{
Denis Talbot, Thierry Duchense, Jacques Brisson, Nathalie Vandal.
}

\seealso{
\code{\link{covariance}}
}
