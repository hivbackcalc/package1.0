% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/results.R
\name{plot.results}
\alias{plot.results}
\title{Plot estimates incidence and undiagnosed cases
 
Overlays the backcalculated incidence on diagnoses in one panel and
undiagnosed counts in another. Cases are indicated by colors.}
\usage{
\method{plot}{results}(x, panel = NULL, valuevar = NULL)
}
\arguments{
\item{x}{List of class "results", the output of combineResults()}

\item{panel}{Character variable name by which to panel the results}

\item{valuevar}{Optional variable name that selects a different group than 
the one represented by the 'value' variable in x$resultsAll}
}
\value{
Paneled ggplot2 object
}
\description{
Plot estimates incidence and undiagnosed cases
 
Overlays the backcalculated incidence on diagnoses in one panel and
undiagnosed counts in another. Cases are indicated by colors.
}
\examples{
undx <- runSubgroups(KCsim, subvar='fakeGroup', intLength=1, 
                         cases=c(`Base Case`='base_case'))
# Right now there's a lot of nesting; you have to access the right element to 
# get the plot
plot(undx[['Group 1']]$results)
# Equivalent: plot Group 2 using the Total-stratified results
plot(undx[['Total-stratified']]$results, valuevar='Group 1')
# Panel
all <- compileSubgroups(undx)
plot(all, panel='Subgroup')
}

