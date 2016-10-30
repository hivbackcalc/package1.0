
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
