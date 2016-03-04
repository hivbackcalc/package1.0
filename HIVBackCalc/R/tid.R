
#' Estimates the continuous CDF of time from infection to diagnosis (TID).
#'  
#' Uses the window of possible infection (infPeriod) and an assumption 
#' about the distribution of the infection hazard in this interval
#' to estimate the the CDF of time from infection to diagnosis (TID).
#'
#' Two infection hazard assumptions are compared in the paper:
#' a "base case" that assumes the probability of infection is uniformly 
#' distributed in the inter-test interval, and 
#' an "upper bound" that assumes infection occurred immediately after
#' the last negative test.
#'  
#' @param infPeriod A vector containing the continuous time from last 
#'          HIV test to diagnosis for each member of the population
#' @param case One of "base_case" or "upper_bound", indicating the distributional
#'          assumption used for the hazard of infection within the infPeriod
#' @param intLength A single number indicating the length in years of discrete 
#'          time intervals by which HIV diagnoses are recorded. The default of 
#'          0.25 represents a quarter-year.
#' @param survivor Logical indicating whether to return S(x)=1-F(x) instead of 
#'          F(x), the cdf
#'  
#' @return A function that takes integer arguments 0 and higher. Will return 0 
#'          for integers greater than floor(max(infPeriod/intLength)) + 1
TID_cdfContinuous <- function(infPeriod,case,intLength,survivor=FALSE) {

    # Remove missing infPeriods and warn if there are zeroes
    infPeriod <- sort(infPeriod[!is.na(infPeriod)]) 
    if(sum(infPeriod==0)!=0) stop('Please resolve infPeriod=0 cases before proceeding, e.g. 
                                  set them to infPeriod=NA (missing) if you cannot identify
                                  a valid infPeriod.')
    n <- length(infPeriod)

    # Define the CDF of time from infection to diagnosis in 
    # the function 'qi'
    switch(case,

         'upper_bound' = {

             # Simply the empirical CDF of the infPeriods
              qi <- function(u) {
                uind <- sum(infPeriod<=u)/n
                if (is.na(uind)) 
                  return(0)
                uind
              }
         },

         'base_case' = {

              # Continuous density of time between infection and diagnosis
              pi <- function(i,eta,infPeriod=infPeriod){
                sapply(i,function(ii){
                  ints <- infPeriod[infPeriod>=ii]
                  sum(1/ints)/length(infPeriod)
                })
              }
              uniqueInf <- unique(infPeriod)
              p<-pi(uniqueInf,,infPeriod) * diff(c(0,uniqueInf))
              cs <- cumsum(p)

              # CDF of density
              qi <- function(u){
                uind <- rev(which(uniqueInf<=u))[1]
                if(is.na(uind))
                  return(0)
                cs[uind]
              }
        })

    if (survivor) {
        Sx <- function(u) { 1-qi(u) }
        return(Sx)
    } else return(qi)
}

#' Evaluates the discrete PDF of time from infection to diagnosis in reporting intervals.
#'  
#' Calls the TID_cdfContinuous function to estimate CDF of the TID, then
#' evaluates the discrete PDF in time intervals that match reporting period
#' units (intLength).
#'  
#' @param infPeriod A vector of continuous times from last HIV test to diagnosis
#'          for a population
#' @param case One of "base_case" or "upper_bound", indicating the 
#'          assumption to apply for when infection occurred within infPeriod
#' @param intLength A single number indicating the length in years of discrete 
#'          time intervals by which HIV diagnoses are recorded. The default of 
#'          0.25 represents a quarter-year.
#'  
#' @return A function that takes integer arguments 0 and higher. Will return 0 
#'          for integers greater than floor(max(infPeriod/intLength)) + 1
TID_pdf <- function(infPeriod,case,intLength) {

    # Remove missing infPeriods 
    infPeriod <- sort(infPeriod[!is.na(infPeriod)]) 
    n <- length(infPeriod)

    # Define the CDF of time from infection to diagnosis in 
    # the function 'qi'
    qi <- TID_cdfContinuous(infPeriod,case,intLength)

    # Use the CDF to define the discrete probability of infection 
    # during interval i to i+1
    pidCalc <- function(i){
    sapply(i,function(ii){
      qi((ii+1)*intLength) - qi(ii*intLength)
    })
    }

    # Calculate the discrete PDF over the m observed intervals, and set
    # probability to zero for longer intervals
    m <- max(infPeriod/intLength) + 1
    pidProbs <- pidCalc(0:m)
    pid <- function(i){
    ifelse(i>m,0,pidProbs[i+1])
    }

    # Return the PDF function
    return(pid)
}

#' Creates a TID (time from infection to diagnosis) object
#'  
#' Calls the TID_pdf() function to estimate the TID from testing history
#' data for each assumption/case. Evaluates the function for the observed
#' time span in the data and returns both the probability and cumulative
#' density distributions.
#'  
#' @param infPeriod The vector containing infection periods, or times 
#'          (in years) between last negative test and diagnosis
#' @param intLength The interval length by which diagnoses are reported
#'          (also in years, 1=1 year)
#' @param cases Cases to estimate; default is c('base_case', 'upper_bound')
#'  
#' @return A nested list. The first tier indicates the assumption used
#'          to estimate the TID. The second tier contains 3 elements:
#'          "pdffxn", the PDF function, and "pdf" and "cdf" which 
#'          indicate the respective distributions
estimateTID <- function(infPeriod, intLength, cases=NULL) {

    # Default cases
    if (is.null(cases)) cases <- c('base_case', 'upper_bound')

    # TID object
    TIDobject <- vector(length=length(cases), mode='list')
    names(TIDobject) <- cases

    # Intervals observed in infPeriod
    maxInt <- max(infPeriod/intLength, na.rm=TRUE)+1

    # Populate with TID functions and actual distributions
    for (c in cases) {
        TIDobject[[c]] <- vector(mode='list', length=3)
        names(TIDobject[[c]]) <- c('pdffxn', 'pdf', 'cdf')
        # PDF function
        TIDobject[[c]]$pdffxn <- TID_pdf(infPeriod, c, intLength)
        # PDF
        TIDobject[[c]]$pdf <- sapply(0:maxInt, TIDobject[[c]]$pdffxn)
        # CDF
        TIDobject[[c]]$cdf <- cumsum(TIDobject[[c]]$pdf)
    }

    class(TIDobject) <- append(class(TIDobject), 'TID')
    return(TIDobject)
}

#' Plots data on time from infection to diagnosis
#'  
#' Plots the discrete probability and cumulative density functions
#' from an object of class 'TID' produced by estimateTID(). Overlays
#' the estimates from all assumptions/cases present in the object
#'  
#' @param x The TID object
#' @param intLength Diagnosed interval length (from timeDx in the data)
#'          in years
#' @param cases Vector of tidier names for the cases present in the 
#'          TID object, in the same order as names(x). Defaults to
#'          names(x)
#' @param legendpos Position for legend, e.g. 'bottom', 'right'
#'  
#' @return ggplot object 
plot.TID <- function(x, intLength, cases=NULL, legendpos=NULL) {

    if (is.null(cases)) cases <- names(x)

    # Construct master data frame
    d <- data.frame(Time=NULL, group=NULL, var=NULL, value=NULL)
    for (c in names(x)) {
        d <- rbind(d, 
                   data.frame(Time=intLength*(0:(length(x[[c]]$pdf)-1)),
                              group='f(x): Probability of Diagnosis',
                              var=cases[which(c==names(x))],
                              value=x[[c]]$pdf),
                   data.frame(Time=intLength*(0:length(x[[c]]$pdf)),
                              group='S(x): Undiagnosed Fraction',
                              var=cases[which(c==names(x))],
                              value=c(1, 1-x[[c]]$cdf)))
    }
    if (is.null(legendpos)) legendpos='bottom'

    # Determine colors needed for scale
    if (length(cases)==2) {
        plotcolors <- c('blue', 'orange2')
    } else if (length(cases)<=12) {
        plotcolors <- c('#a6cee3','#1f78b4','#b2df8a',
                       '#33a02c','#fb9a99','#e31a1c',
                       '#fdbf6f','#ff7f00','#cab2d6',
                       '#6a3d9a','#ffff99','#b15928')[1:length(cases)]
    } else stop('In plot.TID, too many cases to plot')

    p <- ggplot(d) + geom_line(aes(x=Time,y=value,color=var)) + 
      scale_color_manual(values=plotcolors,name="") +
      theme_bw() + 
      ylab("Probability") + 
      xlab("Years Since Infection") +
      scale_x_continuous(expand=c(0,.2)) + 
      theme(legend.position=legendpos) + 
      facet_grid(group~., scales='free_y') 

    return(p)
}


#' Evaluate TID at given time intervals
#'  
#' Returns the discrete probability and cumulative density functions
#' from an object of class 'TID' produced by estimateTID() at specified
#' time points.
#'  
#' @param x The TID object
#' @param intLength Diagnosed interval length (from timeDx in the data)
#'          in years
#' @param times Vector of times to evaluate, in years. Values must be
#'          multiples of intLength
#' @param cases Vector of tidier names for the cases present in the 
#'          TID object, in the same order as names(x). Defaults to
#'          names(x)
#'  
#' @return Data frame
summary.TID <- function(x, intLength, times, cases=NULL) {

    if (sum(times%%intLength)!=0) {
        stop('In summary.TID, times must be multiples of intLength')
    }

    if (is.null(cases)) cases <- names(x)
    ncases <- length(x)

    # Convert intervals to years and identify which times 
    # extract
    years <- intLength*(0:length(x[[1]]$pdf))
    thesetimes <- years%in%times

    # Construct master data frame
    d <- data.frame(Time=times)
    for (i in 1:ncases) {
        d <- cbind(d, x[[i]]$pdf[thesetimes], 1-x[[i]]$cdf[thesetimes])
    }

    colnames(d) <- c('Time', paste(cases, rep(c('f(x)', 'S(x)'), ncases)))
    colnames(d) <- c('Time', paste(as.vector(t(replicate(2,cases))),
                                   c('f(x)', 'S(x)')))

    return(d)
}
