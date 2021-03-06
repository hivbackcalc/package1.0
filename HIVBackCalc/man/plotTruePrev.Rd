% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/results.R
\name{plotTruePrev}
\alias{plotTruePrev}
\title{Function to combine yearly PLWHA prevalence with undiagnosed estimates to provide
true prevalence estimates}
\usage{
plotTruePrev(x, revlevels = TRUE)
}
\arguments{
\item{x}{Results of calcTruePrev}

\item{revlevels}{If the colors of the bars look wrong, try setting to FALSE. The
TRUE setting accommodates an change in factor ordering in newer ggplot versions}
}
\description{
Function to combine yearly PLWHA prevalence with undiagnosed estimates to provide
true prevalence estimates
}
\examples{
# Right now it requires the results to have the Base Case and Upper Bound, 
# named as shown below
undx <- runBackCalc(KCsim, 1, cases=c(`Base Case`='base_case', `Upper Bound`='upper_bound'))
trueprev <- calcTruePrev(undx$results, KCplwh[,c('Year', 'Total')])
plotTruePrev(trueprev)
}

