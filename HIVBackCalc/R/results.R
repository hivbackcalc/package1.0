
#' Plots incidence estimates from estimateIncidence()
#'  
#' Plots incidence estimates over time for object of class "backproj" produced 
#' by estimateIncidence(). Optionally can overlay diagnosis counts.
#'  
#' @param x The "backproj" object containing incidence estimates
#' @param time A vector of length two giving the start and end times for the plot.
#'          Allowed to be left unspecified if x$y has names
#' @param showDiagCounts Logical: should diagnosis counts be plotted
#' @param ... Other arguments passed to plot
#'  
#' @return Base plot 
plot.backproj <- function(x,time,showDiagCounts=TRUE, case="", ...){
  obs <- !is.na(x$y)
  if(missing(time)) 
    time <- as.numeric(names(x$y)[obs])
  else
    time <- seq(from=time[1],to=time[2],length.out=sum(obs))
  plot(time,x$lambda[obs],ylim=c(0,100),type="l",
       main=paste("Estimated Incidence", case, sep="\n"),
       ylab="Count",...)
  if(showDiagCounts)
    points(time,x$y[obs],col="red")
}

#' Create an object of class "results" that contains all results
#'  
#' Incidence and undiagnosed estimates for different cases are initially
#' saved in separate objects. This function combines them into one 
#' object of class "results" to facilitate presenting results. It also
#' summarizes results across all time periods and by year.
#'  
#' @param x List with two tiers: the first tier identifies the cases.
#'          Each case is a list of 2: the first element is the "backproj"
#'          object returned by estimateIncidence(), and the second is
#'          the vector of undiagnosed counts returned by estimateUndiagnosed()
#'  
#' @return List object of class "results" 
combineResults <- function(x) {
  # Times with observed diagnoses
  allTimes <- as.numeric(names(x[[1]][[1]]$y))
  obsTimes <- !is.na(allTimes)
  times <- allTimes[obsTimes]
  x$times <- times
  
  # Diagnoses
  x$diagnoses <- x[[1]][[1]]$y[obsTimes]
  
  # Incidence and Undiagnosed organized by case
  cases <- names(x)[!names(x)%in%c('times', 'diagnoses')]
  for (c in cases) {
    incidence <- x[[c]][[1]]$lambda[obsTimes]
    undiagnosed <- x[[c]][[2]][obsTimes]
    x[[c]] <- list(incidence=incidence, undiagnosed=undiagnosed)
  }
  
  # Data frame with all results
  x$resultsAll <- data.frame(time=times, 
                             group='Diagnoses and Incidence', 
                             var='# Diagnosed',
                             value=x$diagnoses)
  for (c in cases) {
    x$resultsAll  <- rbind(x$resultsAll,
                           data.frame(time=times, 
                                      group='Diagnoses and Incidence', 
                                      var=c,
                                      value=x[[c]]$incidence),
                           data.frame(time=times,
                                      group='Undiagnosed Cases',
                                      var=c,
                                      value=x[[c]]$undiagnosed))
  }
  
  # Data frame with summarized results
  x$resultsSummary <- ddply(x$resultsAll, .(var, group), function(x) c(summary(x$value)))
  x$resultsSummary <- within(x$resultsSummary, {
    group <- as.character(group)
    group[group=='Diagnoses and Incidence' &
            var=='# Diagnosed'] <- 'Diagnoses'
    group[group=='Diagnoses and Incidence'] <- 'Incidence'
  })
  colnames(x$resultsSummary)[1:2] <- c('Diagnoses/Case', 'Estimate')
  
  # Data frame with summarized results by year
  x$resultsSummaryYear <- ddply(transform(x$resultsAll, Year=floor(time)), 
                            .(var, group, Year), function(x) c(summary(x$value)))
  x$resultsSummaryYear <- within(x$resultsSummaryYear, {
    group <- as.character(group)
    group[group=='Diagnoses and Incidence' &
            var=='# Diagnosed'] <- 'Diagnoses'
    group[group=='Diagnoses and Incidence'] <- 'Incidence'
  })
  colnames(x$resultsSummaryYear)[1:2] <- c('Diagnoses/Case', 'Estimate')
  
  class(x) <- append(class(x), 'results')
  return(x)
}

#' Plot estimates incidence and undiagnosed cases
#'  
#' Overlays the backcalculated incidence on diagnoses in one panel and
#' undiagnosed counts in another. Cases are indicated by colors.
#'  
#' @param x List of class "results", the output of combineResults()
#'  
#' @return Paneled ggplot2 object 
plot.results <- function(x) {

    d <- x$resultsAll
    p <- ggplot(d,aes(x=time,y=value, linetype=var))  +   
      geom_line(aes(color=var), size=0.5) +
      geom_point(aes(color=var, shape=var), size=2) + 
      scale_alpha_manual(values=c(.5,1,1),name="") + 
      scale_color_manual(name="", values=c("gray3", "blue", "orange2")) + 
      scale_linetype_manual(name="",values=c(6,3,3)) + 
      scale_shape_manual(name="", values=c(3,16,16)) +
      facet_grid(group~.,scales="free_y") +
      xlab("Time") + ylab("Counts") + 
      geom_blank(aes(x=2008,y=0)) +
      scale_y_continuous(expand=c(.15,0)) + 
      theme_bw() + 
      theme(text = element_text(size = 10)) +
      theme(legend.position="bottom",axis.text.x=element_text(angle=90))

    return(p)
}

######################################################################
# runBackCalc
######################################################################

#' Optional wrapper function to run all the backcalculation steps
#' @param testhist Data frame of class 'testinghistories' containing the time of 
#'        diagnosis within "timeDx" and time since last negative test in 
#'        "infPeriod"
#' @param intLength Interval length by which diagnoses are tracked; 1=1 year.
#' @param cases Vector of names of cases to include for the TID assumptions. 
#'        Defaults to all cases computed by estimateTID, which currently 
#'        offers 'base_case' and 'upper_bound'. Names of vector elements will
#'        be used to label results.
#' @param prev  Optional data frame with 1st column 'Year' and a 2nd column with
#'        PLWHA for the population represented in the testhist object
runBackCalc = function(testhist, intLength, cases=NULL, prev=NULL) {
  
  if (is.null(cases)) {
      cases <- c('base_case', 'upper_bound')
      names(cases) <- c('Base Case', 'Upper Bound')
  }

  
  # Estimate TIDs
  TIDs <- estimateTID(testhist$infPeriod, 
                      intLength,
                      cases)
  
  # Diagnoses
  diagCounts = tabulateDiagnoses(testhist, intLength)
  
  # Initialize incidence and undiagnosed count lists
  incidence <- vector(mode='list', length=length(cases))
  names(incidence) <- cases
  undiagnosed <- incidence
  
  # Estimate incidence and undiagnosed
  for (c in cases) {
    cat('\nEstimating case', c, '...\n')
    incidence[[c]] = estimateIncidence(y=diagCounts,
                                      pid=TIDs[[c]]$pdffxn,
                                      gamma=0.1,
                                      verbose=FALSE)
    undiagnosed[[c]] <- estimateUndiagnosed(incidence[[c]])
  }

  # Compile results - this code allows there to be
  # any number of cases
  toCombine <- vector(mode='list',length=length(cases))
  names(toCombine) <- names(cases)
  for (c in cases) {
      toCombine[[which(cases==c)]] <- list(incidence[[c]], undiagnosed[[c]])
  }
  results <- combineResults(toCombine)

  # True prevalence
  if (!is.null(prev)) trueprev <- calcTruePrev(results, prev) else trueprev <- NULL
  
  return(list(TIDs=TIDs, results=results, trueprev=trueprev, N=nrow(testhist)))
}

######################################################################
# runSubgroups
######################################################################

#' Optional wrapper function to run and compile results for subgroups
#' @param testhist Data frame of class 'testinghistories' containing the time of 
#'        diagnosis within "timeDx" and time since last negative test in 
#'        "infPeriod"
#' @param subvar Name of the variable within the testhist data frame that defines
#'        subgroups within which to run the backcalculation
#' @param intLength Interval length by which diagnoses are tracked; 1=1 year.
#' @param cases Vector of names of cases to include for the TID assumptions. 
#'        Defaults to all cases computed by estimateTID, which currently 
#'        offers 'base_case' and 'upper_bound'. Names of vector elements will
#'        be used to label results.
#' @param prev  Optional frame with 1st column 'Year' and a 2nd column with PLWHA
#'        for the population represented in the testhist object
#' @param save  Optional file path to save compiled true prevalence results
runSubgroups = function(testhist, subvar, intLength, cases=NULL, 
                        prev=NULL, save=NULL) {
  
  # Subvar
  if (is.numeric(testhist[,subvar])) {
    warning('Subgroup variable will be coerced to character')
    testhist[,subvar] <- as.character(testhist[,subvar])
  }
  subgroups <- unique(testhist[,subvar])
  numsub <- length(subgroups)
  
  # Cases
  if (is.null(cases)) {
      cases <- c('base_case', 'upper_bound')
      names(cases) <- c('Base Case', 'Upper Bound')
  }

  # Prevalence
  if (!is.null(prev)) {
      # Check that if prevalence is given, it is given for all the 
      # subgroups. If subvar = 'stageGroup', just make fake prev data
      # because we don't need the subgroup results, just the total-weighted
      if (sum(subgroups%in%colnames(prev))!=numsub) stop('In runSubGroups, 
                              prevalence data are insufficient')
  }
  
  # Prepare to store results for each subgroup
  subResults <- vector('list', length=(numsub+1))
  names(subResults) <- c(as.character(subgroups), 'Total-stratified')
  
  # Loop through subgroups
  for (s in subgroups) {
    
    cat('\nSUBGROUP: ', s, '\n')
    # Run the backcalculation for this subgroup, selecting
    # the correct prevalence column if applicable
    if (!is.null(prev)) {
      subPrev <- prev[, c('Year', as.character(s))]
    } else subPrev <- NULL
    subResults[[s]] <- runBackCalc(testhist[testhist[,subvar]==s,], 
                                   intLength,
                                   cases,
                                   prev=subPrev)
    
      # Add a subgroup identifier to the compiled results
      for (r in c('resultsAll', 'resultsSummary', 'resultsSummaryYear')) {
        subResults[[s]]$results[[r]] = data.frame(Subgroup = s,
                                                  subResults[[s]]$results[[r]],
                                                  check.names=FALSE)
      }
  }

  # Extract the results in order to get subgroup-stratified totals
  resultsAllList <- lapply(lapply(subResults, `[[`, 'results'), `[[`, 'resultsAll')

  # Standardize the times common across the groups - some groups may have diagnoses
  # in years or quarters earlier or later than others. Because we're taking averages
  # of incident/undiagnosed cases across time periods rather than sums, let's just
  # remove the extra time periods. Otherwise we would have to impute something reasonable,
  # since we do have to sum across subgroups - can't just put in a zero.

  times <- lapply(subgroups, function(x) subResults[[x]]$results$resultsAll$time)
  mintime <- max(sapply(times, min))
  maxtime <- min(sapply(times, max))
  keeptimes <- seq(mintime, maxtime, by=intLength)

  # Vectors of results for just the keeptimes
  resultsAllList <- lapply(resultsAllList, function(x) x[x$time%in%keeptimes,]$value)

  # Back to extracting results
  resultsAll <- cbind(subResults[[1]]$results$resultsAll[,c('time', 'group', 'var')],
                      do.call(cbind, resultsAllList))
  resultsAll$value <- apply(resultsAll[,as.character(subgroups)],1,sum)
  
  # Summarize subgroup-stratified totals
  # Data frame with summarized results
  resultsSummary <- ddply(resultsAll, .(var, group), function(x) c(summary(x$value)))
  resultsSummary <- within(resultsSummary, {
    group <- as.character(group)
    group[group=='Diagnoses and Incidence' &
            var=='# Diagnosed'] <- 'Diagnoses'
    group[group=='Diagnoses and Incidence'] <- 'Incidence'
  })
  colnames(resultsSummary)[1:2] <- c('Diagnoses/Case', 'Estimate')
  
  # Data frame with summarized results by year
  resultsSummaryYear <- ddply(transform(resultsAll, Year=floor(time)), 
                                .(var, group, Year), function(x) c(summary(x$value)))
  resultsSummaryYear <- within(resultsSummaryYear, {
    group <- as.character(group)
    group[group=='Diagnoses and Incidence' &
            var=='# Diagnosed'] <- 'Diagnoses'
    group[group=='Diagnoses and Incidence'] <- 'Incidence'
  })
  colnames(resultsSummaryYear)[1:2] <- c('Diagnoses/Case', 'Estimate')
  
  # Create the results list of class 'results'
    resultsL=list(resultsAll=resultsAll,
                      resultsSummary=resultsSummary,
                      resultsSummaryYear=resultsSummaryYear)
    class(resultsL) <- append(class(resultsL), 'results')

  # Save in subResults[['Total-stratified']]$results object
  subResults[['Total-stratified']] <- list(results=resultsL)

  if (!is.null(prev)) {
    
    # Calculate total-stratified true prevalence
    subResults[['Total-stratified']]$trueprev <- 
      calcTruePrev(subResults[['Total-stratified']]$results,
                   prev=data.frame(Year=prev$Year,
                                   Total=apply(prev[,as.character(subgroups)],
                                               1,sum)))
  }
  
  if (!is.null(save)) {
    trueprev <- do.call(rbind, lapply(names(subResults), 
                               function(x) {
                                 data.frame(Subgroup=x,
                                            subResults[[x]]$trueprev,
                                            check.names=FALSE)
                                }))
    write.csv(trueprev[, !colnames(trueprev) %in% c('1st Qu.', 'Median', '3rd Qu.')],
              file=save,
              row.names=FALSE)
  }
  
  return(subResults)
}

######################################################################
# calcTruePrev
######################################################################

#' Function to combine yearly PLWHA prevalence with undiagnosed estimates to provide
#' true prevalence estimates
#' @param x Object of class 'results', the output of combineResults
#' @param prev Data frame with 1st column 'Year' and a 2nd column with PLWH
#'        for the population represented in the results
calcTruePrev = function(x, prev) {
  
  # Fix the column name of the prevalence estimate as 'PLWHA'
  colnames(prev)[2] <- 'PLWHA'
  
  # Subset to Undiagnosed estimates and merge on prevalence
  undiag <- subset(x$resultsSummaryYear, Estimate=='Undiagnosed Cases')
  undiag <- merge(undiag, prev, all.y=TRUE, by='Year')
  
  # Now that we've reduced estimates to those years for which we 
  # were given prevalence, lose the PLWHA column and extract the 
  # estimates and PLWHA as a matrix separately
  estMatrix <- undiag[,4:9]
  prevMatrix <- replicate(ncol(estMatrix), undiag$PLWHA)
  undiag$PLWHA <- NULL
  
  # Compute true prevalence and % undiagnosed
  trueprevMatrix <- estMatrix + prevMatrix
  undiagFracMatrix <- 100*(estMatrix/trueprevMatrix)
  
  # Turn those back into data frames
  trueprev <- undiag
  trueprev$Estimate <- 'True Prevalence'
  trueprev[,4:9] <- trueprevMatrix
  undiagFrac <- undiag
  undiagFrac$Estimate <- 'Undiagnosed Fraction (%)'
  undiagFrac[,4:9] <- undiagFracMatrix
  
  # Prepare prevalence to be able to insert it into the trueprev table
  previnsert = transform(prev, diag='PLWHA', est='PLWHA', min=NA, q1=NA, med=NA, q3=NA, max=NA)
  previnsert = previnsert[,c('Year', 'diag', 'est', 'min', 'q1', 'med', 
                             'PLWHA', 'q3', 'max')]
  colnames(previnsert) <- colnames(undiag)
  
  # Combine and sort
  allResults <- rbind(previnsert, undiag, trueprev, undiagFrac)
  allResults[,4:9] <- round(allResults[,4:9], 1)
  allResults <- allResults[order(allResults$Year, allResults[,2]),]
  
  return(allResults)
}


