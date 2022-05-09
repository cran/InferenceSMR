\title{
Incidences calculations.
}\name{incidences}\alias{incidences}\description{
A function to calculate incidence rates for every combination of age and calendar years. 
}
\usage{incidences(age_min, age_max, year_min, year_max, follow_up,
start_age, start_year, case)}
\arguments{
  \item{age_min}{
The age at which incidences should begin to be calculated.
}
  \item{age_max}{
The age at which incidences should stop to be calculated.
}
  \item{year_min}{
The calendar year at which indicidences should begin to be calculated.
}
  \item{year_max}{
The calendar year at which incidences should stop to be calculated.
}
  \item{follow_up}{
A vector of dimension \code{n} of follow-up times as non-participants.
}
  \item{start_age}{
A vector of dimension \code{n} of ages at the begining of the follow-up.
}
  \item{start_year}{
A vector of dimension \code{n} of calendar years at the begining of the follow-up.
}
  \item{case}{
A vector of dimension \code{n} where each component equals 1 if the follow-up ended because the individual was infected with the decease and 0 otherwise.
}
}\details{
This function can be used to obtain incidences for the functions \code{est.expDeath}, \code{var.expDeath} and \code{inference.SMR}. A complete example of usage is provided in the help page of the \code{screening} dataset.
}\value{
A matrix whose first column is the number of person-years at risk, the second column is the calendar years, the third column is the ages and the fourth column is the incidence rates. 
}\author{
Denis Talbot, Thierry Duchesne, Jacques Brisson, Nathalie Vandal.
}\seealso{
\code{\link{est.expDeath}}, \code{\link{var.expDeath}}, \code{\link{inference.SMR}}, \code{\link{screening}}
}
