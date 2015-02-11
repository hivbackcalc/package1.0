
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
#' object of class "results" to facilitate presenting results.
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
    p <- ggplot(d,aes(x=time,y=value,group=var)) +
      geom_line(aes(alpha=var, color=var)) +
      geom_point(aes(color=var)) +
      scale_alpha_manual(values=c(.5,1,1,1),name="") + 
      scale_linetype_manual(name="",values=c(3,1,2,2)) + 
      theme_bw()+
      theme(legend.position='bottom',axis.text.x=element_text(angle=90)) + 
      scale_color_hue(name="") + 
      scale_x_continuous(breaks=seq(min(d$time),max(d$time),by=2))+
      xlab("Time") + ylab("Counts") +
      facet_grid(group~., scales='free_y')

    return(p)
}


