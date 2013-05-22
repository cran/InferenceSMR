\name{est.expDeath}
\alias{est.expDeath}
\title{
Estimation of the expected number of deaths.
}
\description{
Estimation of the expected number of deaths in a screening program using the method proposed by Sasieni (2003). 
}
\usage{
est.expDeath(contribution, incid, cox, fuzz, covnames)
}
\arguments{
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
  \item{covnames}{
An alphanumeric vector containing the names of the covariates used to estimate the survival in the cohort of non-participants, that is, the names of the covariates used to obtain the \code{cox} object.
}
}

\value{
Returns the expected number of deaths}
\references{
Sasieni P. (2003) \emph{On the expected number of cancer deaths during follow-up of an initially cancer-free cohort}.
Epidemiology, 14, 108-110.
}
\author{
Denis Talbot, Thierry Duchesne, Jacques Brisson, Nathalie Vandal.
}
\note{
A complete example of usage is provided in the help page of the \code{screening} dataset.
}

\examples{
#This example uses pre-built objects and shows the simple usage
#of the est.expDeath function when those objects already exists.
#For an example of how to built those object, refer to the 
#help page of the screening dataset.

data(req.objects);
cox.data = req.objects$cox.data;
#Remove "#" to run example :
#est.expDeath(req.objects$contribution,req.objects$incid,req.objects$cox,fuzz = 0.01, 
#req.objects$covnames);

#[1] 33.44264

}

\seealso{
\code{\link{var.expDeath}},\code{\link{inference.SMR}}, \code{\link{screening}}
}
