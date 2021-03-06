% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tid.R
\name{plot.TID}
\alias{plot.TID}
\title{Plots data on time from infection to diagnosis
 
Plots the discrete probability and cumulative density functions
from an object of class 'TID' produced by estimateTID(). Overlays
the estimates from all assumptions/cases present in the object}
\usage{
\method{plot}{TID}(x, intLength, cases = NULL, legendpos = NULL)
}
\arguments{
\item{x}{The TID object}

\item{intLength}{Diagnosed interval length (from timeDx in the data)
in years}

\item{cases}{Vector of tidier names for the cases present in the 
TID object, in the same order as names(x). Defaults to
names(x)}

\item{legendpos}{Position for legend, e.g. 'bottom', 'right'}
}
\value{
ggplot object
}
\description{
Plots data on time from infection to diagnosis
 
Plots the discrete probability and cumulative density functions
from an object of class 'TID' produced by estimateTID(). Overlays
the estimates from all assumptions/cases present in the object
}

