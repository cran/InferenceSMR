\title{
Variance estimation of the expected number of deaths
}\name{var.expDeath}\alias{var.expDeath}\description{
This function estimates the variance of the expected number of deaths when the latter is estimated using Sasieni's method.
}\usage{var.expDeath(contribution, incid, cox, fuzz, Poisson = FALSE, covnames)}\arguments{
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
consider objects not differing by more then \code{fuzz} to be equal. The \code{fuzz} option should be chosen to be a small
positive number, for instance 0.0001. 
}
  \item{covnames}{
An alphanumeric vector containing the names of the covariates used to estimate the survival in the cohort of non-participants, that is, the names of the covariates used to obtain the \code{cox} object.
}

\item{Poisson}{
Indicates whether the incidences' variance should be estimated with a Poisson distribution (TRUE) or a binomial distribution (FALSE). The default is FALSE.
}
}\value{
The function returns the variance of the expected number of deaths}\note{
A complete example of how to use this function is available in the help page of the \code{screening} dataset.
}\author{
Denis Talbot, Thierry Duchesne, Jacques Brisson, Nathalie Vandal.
}\references{
Talbot, D., Duchesne, T., Brisson, J., Vandal, N. (2011) \emph{Variance estimation and confidence intervals for the standardized mortality ratio with application to the assessment of a cancer screening program}, Statistics in Medicine, 30, 3024-3037. 
}\seealso{
\code{\link{est.expDeath}}, \code{\link{inference.SMR}}, \code{\link{screening}}
}\examples{
#This example uses pre-built objects and shows the simple usage
#of the est.expDeath function when those objects already exists.
#For an example of how to built those object, refer to the 
#help page of the screening dataset.

#Remove "#" to run example. The function can be quite long (a few hours) to run:
#data(req.objects);
#cox.data = req.objects$cox.data;
#var.expDeath(req.objects$contribution,req.objects$incid,req.objects$cox,fuzz = 0.01, 
#req.objects$covnames);


#[1,] 39.31382
}
