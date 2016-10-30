
#' MSM HIV testing histories approximating King County, WA in 2006-2012
#'
#' A dataset containing the date of last negative test as well as other characteristics
#' of HIV diagnoses in a simulated population. The marginal population features approximate
#' HIV diagnoses among MSM in King County, WA between 2006 and 2012.
#'
#' @format A data frame with 1522 rows and 15 columns:
#' \describe{
#'   \item{Population}{ID for the population}
#'   \item{mode}{Mode of HIV transmission}
#'   \item{mode2}{Mode of HIV transmission, intended to have only two values, MSM and non-MSM. In this simulated population, mode and mode2 are MSM for all diagnoses}
#'   \item{race}{Race}
#'   \item{hdx_age}{Age}
#'   \item{agecat5}{5-year age group}
#'   \item{yearDx}{Year of diagnosis}
#'   \item{timeDx}{Quarter-year of diagnosis}
#'   \item{everHadNegTest}{Logical indicating whether the individual ever had a negative test prior to diagnosis}
#'   \item{infPeriod}{Time, in years, between last negative test and HIV diagnosis. For individuals with everHadNegTest as FALSE or never-testers, this is imputed as the minimum of age-16 and 17.98, a max allowed time elapsing between infection and diagnosis based on the age incubation distribution}
#'   \item{infPeriod_imputeNA}{A version of infPeriod in which the missing values are imputed using the same rule stated above for the never-testers}
#'   \item{fakeGroup}{A fake subgroup variable, useful for unit tests}
#'   \item{everHadNegTestM}{A hypothetical everHadNegTest variable with 50% of repeat testers, i.e. everHadNegTest as TRUE, turned into never-testers, i.e. everHadNegTest as FALSE}
#'   \item{infPeriodM}{The corresponding infPeriod to everHadNegTestM, with the 50% of hypothetical never-testers getting the min of age-16 or 17.98 assumption}
#' }
#' @source Simulation based on data from Public Health Seattle King County
"KCsim"


#' MSM living with HIV in KC 2006-2012
#' 
#' A dataset containing rough estimates of MSM PLWH in King County, by race.
#' Estimates are specified by year for unit-tests but reflect a time period
#' average
#' 
#' @format A data frame with 7 rows and 5 columns:
#' \describe{
#'      \item{Year}{Year of estimate}
#'      \item{White}{MSM living with HIV, White race}
#'      \item{Black}{MSM living with HIV, Black race}
#'      \item{Hisp}{MSM living with HIV, Hispanic race}
#'      \item{Total}{MSM living with HIV}
#' }
#' @source Based on estimates from Public Health Seattle King County
"KCplwh"


