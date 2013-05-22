\name{inference.SMR}
\alias{inference.SMR}
\title{
Inference about the standardized mortality ratio (SMR) when evaltuating the effect of a screening program on survival.}
\description{
This function estimates the expected number of deaths, its variance, the SMR and confidence intervals about the SMR.
}
\usage{
inference.SMR(obs.death, normal = "log-smr", alpha = 0.05, contribution,
incid, cox, fuzz = 0.01, Poisson = FALSE, covnames)
}
\arguments{
  \item{obs.death}{
The observed number of deaths for the people participating in the screening program. A numeric value.
}
  \item{normal}{
Indicates at which level should the normality assumption be made, either at the SMR level, at the log-SMR level or at the root-SMR level. A character vector containing one or many of the following elements ``smr'', ``log-smr'' and ``root-smr''.
}
  \item{alpha}{
The nominal error rate of the confidence intervals. A numeric value between 0 and 1, e.g. 0.05 to obtain a 95 \% confidence interval.
}
  \item{contribution}{
An object of contributions produced by the function \code{contrib}.
}
  \item{incid}{
A matrix containing: the incidences, the value of the covariates and the person-years at risk, in that order. It can be obtained with the function \code{incidences}.
}
  \item{cox}{
An oject of class \code{coxph} containing the model that was used to estimate the survival in the cohort of non-participants.
}
  \item{fuzz}{
Numerical precision is problematic when it comes to test equality between objects. The option \code{fuzz} is used to
consider objects not differing by more than \code{fuzz} to be equal. The \code{fuzz} option should be chosen to be a small
positive number, for instance 0.0001. 
}
  \item{Poisson}{
Indicates whether the incidences' variance should be estimated with a Poisson distribution (TRUE) or a binomial distribution (FALSE). The default is FALSE.
}
  \item{covnames}{
An alphanumeric vector containing the names of the covariates used to estimate the survival in the cohort of non-participants, that is, the names of the covariates used to obtain the \code{cox} object.
}
}
\details{
The \code{inference.SMR} function estimates the expected number of deaths as in Sasieni (2003), estimates the variance of the expected number of deaths and builds confidence intervals as in Talbot et al. (2011). As suggested in the latter, the variance of the observed number of deaths is estimated by the observed number of deaths.
}
\value{
\item{expected}{The expected number of deaths}
\item{obs.death}{The observed number of deaths}
\item{variance}{The variance of the expected number of deaths}
\item{smr}{The standardized mortality ratio}
\item{smr.var}{The variance of the SMR. Only returned if ``smr'' was given in the \code{normal} argument.}
\item{smr.ci}{A 1-\code{alpha} confidence interval for the SMR. Only returned if ``smr'' was given in the \code{normal} argument.}
\item{logSMR.var}{The variance of the natural logarithm of the SMR. Only returned if ``log-smr'' was given in the \code{normal} argument.}
\item{logSMR.ci}{A 1-\code{alpha} confidence interval for the log-SMR. Only returned if ``log-smr'' was given in the \code{normal} argument.}
\item{rootSMR.var}{The variance of the square root of the SMR. Only returned if ``root-smr'' was given in the \code{normal} argument.}
\item{rootSMR.ci}{A 1-\code{alpha} confidence interval for the root-SMR. Only returned if ``root-smr'' was given in the \code{normal} argument.}
}
\references{
Sasieni P. (2003) \emph{On the expected number of cancer deaths during follow-up of an initially cancer-free cohort}.
Epidemiology, 14, 108-110.

Talbot, D., Duchesne, T., Brisson, J., Vandal, N. (2011) \emph{Variance estimation and confidence intervals for the standardized mortality ratio with application to the assessment of a cancer screening program}, Statistics in Medicine, 30, 3024-3037. 
}
\author{
Denis Talbot, Thierry Duchesne, Jacques Brisson, Nathalie Vandal.
}
\note{
A complete example of usage is provided in the help page of the \code{screening} dataset.
}

\examples{
#This example uses pre-built objects and shows the simple usage
#of the est.expDeath function when those objects already exist.
#For an example of how to build those objects, refer to the 
#help page of the screening dataset.

#Estimating the variance can be very long even in this small sample example, e.g. a few hours.
#Remove "#" to run example :
#data(req.objects);
#cox.data = req.objects$cox.data;
#results = inference.SMR(obs.death = sum(screening$deathSCN),
# normal = c("smr", "log-smr", "root-smr"),
#	 alpha = 0.05, req.objects$contribution, req.objects$incid,
#	 cox = req.objects$cox, fuzz = 0.01, Poisson = TRUE, req.objects$covnames);


#********  INFERENCE ABOUT THE SMR  ********* 
#
#Observed =  18  Expected =  33.44264 
#Obs.var. =  18  Exp.var. =  39.38153 
#SMR =  0.5382351 
#
# 95 % Confidence intervals with normality assumption at : 
#
#The SMR level : ( 0.2204119 0.8560583 )
#
#The log-SMR level : ( 0.2982118 0.9714471 )
#
#The root-SMR level : ( 0.2673299 0.9029762 )

#results
#
#$expected
#[1] 33.44264
#
#$obs.death
#[1] 18
#
#$variance
#            2
#[1,] 39.38153
#
#$smr
#[1] 0.5400112
#
#$smr.var
#              2
#[1,] 0.02629511
#
#$smr.ci
#[1] 0.2204119 0.8560583
#
#$logSMR.var
#              2
#[1,] 0.09076763
#
#$logSMR.ci
#[1] 0.2982118 0.9714471
#
#$rootSMR.var
#              2
#[1,] 0.01221358
#
#$rootSMR.ci
#[1] 0.2673299 0.9029762



}
\seealso{
\code{\link{est.expDeath}}, \code{\link{var.expDeath}}, \code{\link{screening}}
}

