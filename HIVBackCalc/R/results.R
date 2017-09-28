
#' Check object is of class "results"
#'  
#' @param x Object to check
is.results <- function(x) { inherits(x, 'results') }


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

#' Summarize results by year
#' 
#' Summarize a resultsAll object created by combineResults into 
#' yearly estimates
#
#' @param x a resultsAll data frame
#' @param acrossYears Logical indicating whether to collapse across years
#' @param addclass Logical indicating whether to return x with an 
#'      appended class of "results"
#' @return Original object with resultsSummary or resultsSummaryYear element 
#' @examples
#' 
#' undx <- runBackCalc(KCsim, 1, cases=c(Base='base_case'))
#' # Compare outputs - they're identical
#' undx$results$resultsSummary
#' compare <- sumResults(undx$results, acrossYears=TRUE, addclass=TRUE)
#' compare$resultsSummary
#' # Calculate true prev - need to have year strata
#' compare <- sumResults(undx$results, acrossYears=FALSE, addclass=TRUE)
#' trueprev <- calcTruePrev(compare, KCplwh[,c('Year', 'Total')])
sumResults <- function(x, acrossYears=TRUE, addclass=FALSE) {

    if (!'resultsAll'%in%names(x)) stop('x does not have a resultsAll object')
    # Establish stratum
    byvars <- c('var', 'group')
    if (!acrossYears) {
        byvars <- c(byvars, 'Year')
        x$resultsAll$Year <- floor(x$resultsAll$time)
    }
    # Summarize
    toreturn <- ddply(x$resultsAll, byvars, function(x) c(summary(x$value)))
    toreturn <- within(toreturn,{
        group <- as.character(group)
        group[group=='Diagnoses and Incidence' &
                var=='# Diagnosed'] <- 'Diagnoses'
        group[group=='Diagnoses and Incidence'] <- 'Incidence'
    })
    colnames(toreturn)[1:2] <- c('Diagnoses/Case', 'Estimate')

    if (acrossYears) x$resultsSummary <- toreturn else x$resultsSummaryYear <- toreturn

    if (addclass) class(x) <- append(class(x), 'results')

    return(x)
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
  x <- sumResults(x, acrossYears=TRUE)
  
  # Data frame with summarized results by year
  x <- sumResults(x, acrossYears=FALSE, addclass=TRUE)
  
  return(x)
}

#' Quickly append results from separate subgroups into one results object
#' 
#' Takes runSubgroups output and appends results to facilitate plotting primarily
#' 
#' @param x output from runSubgroups
#' @param includeTotal Logical indicating whether the 'Total-stratified' estimate
#'          should be included
#' @return Data frame of class 'results' 
#' @examples
#' 
#' undx <- runSubgroups(KCsim, subvar='fakeGroup', intLength=1, 
#'                          cases=c(`Base Case`='base_case'))
#' all <- compileSubgroups(undx)
compileSubgroups <- function(x, includeTotal=FALSE) {
    if (!includeTotal) x$`Total-stratified`=NULL
    all <- list(resultsAll=ldply(names(x), function(x) results[[x]]$results$resultsAll))
    class(all)  <- append(class(x), 'results')
    return(all)
}

#' Plot estimates incidence and undiagnosed cases
#'  
#' Overlays the backcalculated incidence on diagnoses in one panel and
#' undiagnosed counts in another. Cases are indicated by colors.
#'  
#' @param x List of class "results", the output of combineResults()
#' @param panel Character variable name by which to panel the results
#' @param valuevar Optional variable name that selects a different group than 
#'          the one represented by the 'value' variable in x$resultsAll
#'  
#' @return Paneled ggplot2 object 
#' @examples
#' undx <- runSubgroups(KCsim, subvar='fakeGroup', intLength=1, 
#'                          cases=c(`Base Case`='base_case'))
#' # Right now there's a lot of nesting; you have to access the right element to 
#' # get the plot
#' plot(undx[['Group 1']]$results)
#' # Equivalent: plot Group 2 using the Total-stratified results
#' plot(undx[['Total-stratified']]$results, valuevar='Group 1')
#' # Panel
#' all <- compileSubgroups(undx)
#' plot(all, panel='Subgroup')
plot.results <- function(x, panel=NULL, valuevar=NULL) {

    # Check that x is of class results
    if (!is.results(x)) stop('x is not of class results; look at x and make sure you have the right list element')
    # Grab the 'resultsAll' object
    d <- x$resultsAll 
    # Define the panel variable
    if (!is.null(panel)) {
        if (!panel%in%colnames(d)) stop('Panel variable not in data')
        d$valuevar  <- d[,panel]
    }
    # Change the value variable if indicated
    # ADD UNIT TEST
    if (!is.null(valuevar)) {
        if (!valuevar%in%names(d)) stop('valuevar variable is not in x$resultsAll')
        d$value <- d[,valuevar]
    }

    p <- ggplot(d,aes(x=time,y=value, linetype=var))  +   
      geom_line(aes(color=var), size=0.5) +
      geom_point(aes(color=var, shape=var), size=2) + 
      scale_alpha_manual(values=c(.5,1,1),name="") + 
      scale_color_manual(name="", values=c("gray3", "blue", "orange2")) + 
      scale_linetype_manual(name="",values=c(6,3,3)) + 
      scale_shape_manual(name="", values=c(3,16,16)) +
      xlab("Time") + ylab("Counts") + 
      geom_blank(aes(x=min(d$time),y=0)) +
      scale_x_continuous(breaks=seq(min(d$time),max(d$time),by=1))+
      scale_y_continuous(expand=c(.15,0)) + 
      theme_bw() + 
      theme(text = element_text(size = 10)) +
      theme(legend.position="bottom",axis.text.x=element_text(angle=90))
   
    if (!is.null(panel)) {
        p <- p+facet_grid(group~Group, scales='free_y')
    } else p <- p+ facet_grid(group~.,scales="free_y") 

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
#' @examples
#' # Example with fake prevalence corresponding to the fake groups
#' # Use check.names=FALSE in the prev data frame construction to ensure
#' # that the column names match the values of the fakeGroup variable
#' undx <- runSubgroups(KCsim, subvar='fakeGroup', intLength=1, 
#'                          cases=c(`Base Case`='base_case'),
#'                          prev=data.frame(Year=KCplwh$Year,
#'                                          `Group 1`=KCplwh$Total/2,
#'                                          `Group 2`=KCplwh$Total/2,
#'                                          Total=KCplwh$Total,
#'                                          check.names=FALSE))
#' # Example with saving results
#' undx <- runSubgroups(KCsim, subvar='fakeGroup', intLength=1, 
#'                          cases=c(`Base Case`='base_case'),
#'                          prev=data.frame(Year=KCplwh$Year,
#'                                          `Group 1`=KCplwh$Total/2,
#'                                          `Group 2`=KCplwh$Total/2,
#'                                          Total=KCplwh$Total,
#'                                          check.names=FALSE),
#'                          save='testfile.csv')
runSubgroups = function(testhist, subvar, intLength, cases=NULL, 
                        prev=NULL, save=NULL) {
  
  # Subvar
  if (is.numeric(testhist[,subvar])) {
    warning('Subgroup variable will be coerced to character')
    testhist[,subvar] <- as.character(testhist[,subvar])
  }
  subgroups <- unique(testhist[,subvar])
  numsub <- length(subgroups)
  
  # Cases (vector must have names)
  if (is.null(cases)) {
      cases <- c('base_case', 'upper_bound')
      names(cases) <- c('Base Case', 'Upper Bound')
  } else if (is.null(names(cases))) stop('Cases vector must have names')

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
  resultsAll$value <- apply(as.matrix(resultsAll[,as.character(subgroups)]),1,sum)
  
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
                                   Total=apply(as.matrix(prev[,as.character(subgroups)]),
                                       1,sum)))
    # Compile and save trueprev
    trueprevAll <- 
        do.call(rbind, lapply(names(subResults), 
                               function(x) {
                                 data.frame(Subgroup=x,
                                            subResults[[x]]$trueprev,
                                            check.names=FALSE)
                                }))
    subResults[['Total-stratified']]$trueprevAll <- 
        trueprevAll[, !colnames(trueprevAll) %in% 
                      c('1st Qu.', 'Median', '3rd Qu.')]
  }
  
  if (!is.null(save)) {
        write.csv(subResults[['Total-stratified']]$results$resultsAll,
                  file=gsub('.csv', '_resultsAll.csv', save),
                  row.names=FALSE)
        if (!is.null(prev)) {
            write.csv(subResults[['Total-stratified']]$trueprevAll,
                      file=gsub('.csv', '_trueprev.csv', save),
                      row.names=FALSE)
        }
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
#' @examples
#' undx <- runBackCalc(KCsim, 1, cases=c(Base='base_case'))
#' trueprev <- calcTruePrev(undx$results, KCplwh[,c('Year', 'Total')])

calcTruePrev = function(x, prev) {
  
    # Check that x is of class results
    if (!is.results(x)) stop('x is not of class results; look at x and make sure you have the right list element')

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

######################################################################
# plotTruePrev
######################################################################

#' Function to combine yearly PLWHA prevalence with undiagnosed estimates to provide
#' true prevalence estimates
#' @param x Results of calcTruePrev
#' @param revlevels If the colors of the bars look wrong, try setting to FALSE. The
#' TRUE setting accommodates an change in factor ordering in newer ggplot versions
#' @examples
#' # Right now it requires the results to have the Base Case and Upper Bound, 
#' # named as shown below
#' undx <- runBackCalc(KCsim, 1, cases=c(`Base Case`='base_case', `Upper Bound`='upper_bound'))
#' trueprev <- calcTruePrev(undx$results, KCplwh[,c('Year', 'Total')])
#' plotTruePrev(trueprev)

plotTruePrev <- function(x, revlevels=TRUE) {

        # Not developed yet
        if (!'Base Case'%in%unique(x[['Diagnoses/Case']]) |
            !'Upper Bound'%in%unique(x[['Diagnoses/Case']])) stop('Must have Base Case and Upper Bound, named that way')

        tp <- subset(x, Estimate=='PLWHA' | Estimate=='Undiagnosed Cases', 
                     select=c('Year', 'Diagnoses/Case', 'Estimate', 'Mean'))
        colnames(tp)[which(colnames(tp)=='Diagnoses/Case')] <- 'Case'
        plwh <- subset(tp, Estimate=='PLWHA')
        tp2 <- rbind(subset(tp, Estimate!='PLWHA'),
                     data.frame(subset(plwh, select=c('Year', 'Mean', 'Estimate')),
                                Case='Base Case'),
                     data.frame(subset(plwh, select=c('Year', 'Mean', 'Estimate')),
                                Case='Upper Bound'))
        tp2$Estimate <- factor(tp2$Estimate, levels=c('PLWHA', 'Undiagnosed Cases'),
                               labels=c('Diagnosed PLWHA', 'Undiagnosed Cases'))
        tp2 = arrange(tp2, Year, Estimate)

        tp3 = ddply(tp2, .(Year, Case), transform, Percent = Mean/sum(Mean) * 100)

        tp3 = ddply(tp3, .(Year, Case), transform, pos = (cumsum(Mean) - 0.5 * Mean))
        tp3$label = paste0(sprintf("%.0f", tp3$Percent), "%")

        if (revlevels) {
            tp3 <- transform(tp3, Estimate = factor(Estimate, 
                                                    levels = rev(levels(Estimate)), 
                                                    labels = rev(levels(Estimate))))
            reverseit <- TRUE
        } else reverseit <- FALSE

        p <- ggplot(tp3, aes(x = factor(Year), y = Mean, fill = Estimate)) +
           geom_bar(stat = "identity", width = .7) +
              geom_text(aes(y = pos, label = label), size = 4) +
                 coord_flip() + 
                 facet_grid(.~Case) + theme_bw() + 
                 scale_x_discrete(name='') + 
                 scale_y_continuous(name='Number of Cases') + 
                 scale_fill_manual(name='', values=rev(c('#a8ddb5', '#43a2ca')),
                                   guide=guide_legend(reverse=reverseit)) + 
                 theme_bw() + 
                 theme(legend.position='bottom') +
                 theme(axis.text = element_text(size = 10))
        
        return(p)
}

######################################################################
# runAnalysis
######################################################################

#' Wrapper function that runs an entire analysis
#' 
#' Ideally this function will successfully produce all the major
#' plots and results in one go
#' 
#' @param testhist Data frame of class 'testinghistories' containing the time of 
#'        diagnosis within "timeDx" and time since last negative test in 
#'        "infPeriod"
#' @param descriptives NAMED vector with descriptive variables to use in
#'          looking at the testing histories. See help page for tabTestHist,
#'          the "variables" argument
#'        subgroups within which to run the backcalculation
#' @param subvar Name of the variable within the testhist data frame that defines
#'        subgroups within which to run the backcalculation
#' @param intLength Interval length by which diagnoses are tracked; 1=1 year.
#' @param cases Vector of names of cases to include for the TID assumptions. 
#'        Defaults to all cases computed by estimateTID, which currently 
#'        offers 'base_case' and 'upper_bound'. Names of vector elements will
#'        be used to label results.
#' @param runEstimation  Logical indicating whether to run a fresh back-calculation
#' @param savedEstimation  Logical indicating whether to load a saved back-calculation. 
#'          Both runEstimation and savedEstimation can be set to FALSE to return
#'          just descriptives and TID but no back-calculation
#' @param prev  Optional frame with 1st column 'Year' and a 2nd column with PLWHA
#'        for the population represented in the testhist object
#' @param save  Optional file path to save compiled results or load saved results
#' @return List object with plot and results list
#' @examples
#' # Basic example with a fresh estimation
#'allRes <- runAnalysis(KCsim, list(c(`Group`='fakeGroup', `Race`='race')),
#'                      'fakeGroup', 1, c(`Base Case`='base_case',
#'                                        `Upper Bound`='upper_bound'),
#'                      runEstimation=TRUE)
#'# Compute true prevalence and save results
#'allRes <- runAnalysis(KCsim, list(c(`Group`='fakeGroup', `Race`='race')),
#'                      'fakeGroup', 1, c(`Base Case`='base_case',
#'                                        `Upper Bound`='upper_bound'),
#'                      runEstimation=TRUE,
#'                      prev=data.frame(Year=2006:2012,
#'                                      `Group 1`=KCplwh$Total/2,
#'                                      `Group 2`=KCplwh$Total/2,
#'                                      Total=KCplwh$Total,
#'                                      check.names=FALSE),
#'                      save='testfile.csv')
#'# Load the saved estimation and plot true prev results from it
#'allRes <- runAnalysis(KCsim, list(c(`Group`='fakeGroup', `Race`='race')),
#'                      'fakeGroup', 1, c(`Base Case`='base_case',
#'                                        `Upper Bound`='upper_bound'),
#'                      runEstimation=FALSE,
#'                      savedEstimation=TRUE,
#'                      save='testfile.csv')
#'plotTruePrev(subset(allRes$results$trueprev, 
#'                    Subgroup=='Total-stratified'))

runAnalysis <- function(testhist, descriptives, 
                        subvar, intLength, cases, 
                        runEstimation=TRUE,
                        savedEstimation=FALSE,
                        prev=NULL, save=NULL) {

    warning('runAnalysis is under development - needs unit testing and special cases for when no subgroup variable is desired, e.g. you want to analyze a single subgroup, like analyzing all KC MSM')

    # Hardcode the subvar as subvar
    testhist$subvar <- testhist[,subvar]
    
    # Diagnoses
    dx <- vector('list')
    dx$plot <- plotDiagnoses(testhist, showlegend=TRUE)

    # Testing histories
    sc <- ifelse(is.list(descriptives), FALSE, TRUE)
    table <- tabTestHist(testhist, descriptives, 
                         supercolumn=sc, fullsample_row=TRUE)
    plot <- plotTestHist(testhist, panel=subvar)
    th <- list(table=table, plot=plot)

    # Inter-test intervals (infPeriod)
    ### Overall summary
    sum <- summary(testhist$infPeriod)
    ### Mean window lengths by year
    means <- c(by(testhist,testhist$yearDx,function(x)mean(x$infPeriod,na.rm=T)))
    ### Test for time trend
    trend <- oneway.test(infPeriod~yearDx,data=testhist)
    infPeriod <- list(sum=sum, means=means, trend=trend)

    # TID
    separate <- sapply(levels(testhist$subvar),
                   function(x) {
                       estimateTID(testhist[testhist$subvar==x,'infPeriod'],
                                   intLength)
                   }, USE.NAMES=TRUE, simplify=FALSE)
    combined <- vector('list')
    for (i in 1:length(separate)) {
        for (j in 1:length(separate[[1]])) {
            index <- 2*i
            if (j==1) index <- index-1
            combined[[index]] <- separate[[i]][[j]]
        }
    }
    names(combined) <- 
        paste(levels(testhist$subvar), 
                             rep(names(separate[[1]]), length(separate)))
    class(combined) <- append(class(combined), 'TID')
    # For some reason, plot.TID doesn't work so I will not include the 
    # combined object in the output
    tid <- list(separate=separate)

    # Results, new or saved
    if (runEstimation) {
        results <- runSubgroups(testhist, subvar, intLength, cases,
                                prev, save)
    } else results <- NULL
    if (runEstimation) {
        resultsCompiled <- list(trueprev=results[['Total-stratified']]$trueprevAll,
                                results=results[['Total-stratified']]$results)
    } else if (savedEstimation) {
        resultsCompiled <- list(trueprev=read.csv(gsub('.csv', '_trueprev.csv', save), 
                                                  header=TRUE, check.names=FALSE), 
                                results=list(resultsAll=read.csv(gsub('.csv', 
                                                                      '_resultsAll.csv', 
                                                                      save), 
                                                                 header=TRUE,
                                                                 check.names=FALSE)))
        class(resultsCompiled$results) <- append(class(results$results), 'results')
        colnames(resultsCompiled$trueprev)[which(colnames(resultsCompiled$trueprev)=='Diagnoses.Case')] <- 
            'Diagnoses/Case'
    } else resultsCompiled <- FALSE

    # Plots don't run for saved estimated, since those are compiled 
    # Could code up mimicing the list structure of fresh back-calc results using the saved 
    # results, and then the plots would run too
    if (runEstimation) {
        resultsPlots <- sapply(results, function(x) {
                                   p <- plot(x$results)
                                   return(p)
                                }, USE.NAMES=TRUE, simplify=FALSE)
    } else if (savedEstimation) {
        subgroups <- colnames(resultsCompiled$results$resultsAll)
        subgroups <- subgroups[!subgroups%in%c('time', 'group', 'var', 'value')]
        resultsPlots <- sapply(c(subgroups,'Total-stratified'),
                               function(x) {
                                   if (x=='Total-stratified') var <- NULL else var <- x
                                   p <- plot(resultsCompiled$results,
                                             valuevar=var)
                                   return(p)

                               }, USE.NAMES=TRUE, simplify=FALSE)
    } else resultsPlots <- NULL
    if (runEstimation & !is.null(prev)) {
        resultsPrevPlots <- sapply(results, function(x) {
                                   p <- plotTruePrev(x$trueprev)
                                   return(p)
                                }, USE.NAMES=TRUE, simplify=FALSE)
    } else if (savedEstimation) {
        resultsPrevPlots <- sapply(levels(resultsCompiled$trueprev$Subgroup),
                                   function(x) {
                                       p <- plotTruePrev(subset(resultsCompiled$trueprev,
                                                                Subgroup==x))
                                       return(p)
                                }, USE.NAMES=TRUE, simplify=FALSE)
    } else resultsPrevPlots <- NULL

    return(list(dx=dx, th=th, infPeriod=infPeriod, tid=tid,
           results=results, 
           resultsCompiled=resultsCompiled,
           resultsPlots=resultsPlots, 
           resultsPrevPlots=resultsPrevPlots))
}

