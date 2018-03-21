
#' Tabulates HIV diagnoses by interval, as input for estimateIncidence().
#'  
#' HIV diagnosis counts for every interval within the calendar time period
#' of the data are appended to empty intervals that will indicate to 
#' estimateIncidence() how far back to project incidence counts.
#'  
#' @param testhist Object of class 'testinghistories' containing the 
#'          time of diagnosis within "timeDx"
#' @param intLength Interval length by which diagnoses are tracked;
#'          allowed values are 1=1 year, 0.5=half year, 0.25=quarter year,
#'          or 1/12=1 month. Do not specify 1 month as 0.083, but use 1/12.
#' @param nPriorInt Number of incidence intervals to backproject
#'          prior to earliest observed diagnoses. Use NULL to skip
#'          this and simply return observed diagnosis counts. Use the default
#'          to prepare the data for use with estimateIncidence()
#'  
#' @return Vector of NA's plus diagnosis counts
#' @examples
#' # Create KCsimM which has (fake) time steps of months, i.e. 1/12ths in years,
#' # with duplicate rows to get more diagnoses per time step
#' data(KCsim)
#' KCsimM <- rbind(KCsim, KCsim, KCsim, KCsim)
#' KCsimM <- transform(KCsimM, timeDx=yearDx +
#'                    sample((1:10)*(1/12), 
#'                           size=nrow(KCsimM), replace=TRUE))
#' 
#' # Monthly time steps
#' diagCountsM = tabulateDiagnoses(KCsimM, intLength=1/12)
#' # Yearly time steps: the function will aggregate for you
#' diagCounts = tabulateDiagnoses(KCsimM, intLength=1)
#' 
#' # Note that there are too many months with zero counts for 
#' # estimateIncidence() to work. Only if we replace the zero counts
#' # will the estimation run. This demonstrates the value of aggregating
#' # up to larger time steps.
#' diagCountsM[diagCountsM<1] <- 80 # Not good practice, for demonstration only
#' TIDsM <- estimateTID(KCsimM$infPeriod, intLength=1/12)
#' incidenceBaseM = estimateIncidence(y=diagCountsM, 
#'                                   pid=TIDsM[['base_case']]$pdffxn, 
#'                                   gamma=0.1, 
#'                                   verbose=TRUE)

tabulateDiagnoses <- function(testhist, intLength, nPriorInt=100) {

    # Aggregate diagnoses if necessary
    testhist$timeDx <- aggregateDiagnoses(testhist$timeDx, intLength)

    # Vector of all intervals
    allTimes <- seq(min(testhist$timeDx),
                    max(testhist$timeDx), 
                    by=intLength)
    # Tabulation of intervals with observed data
    obsCounts <- table(testhist$timeDx)
    # Account for intervals with 0 cases
    allCounts <- structure(rep(0,length(allTimes)),
                           class='table',
                           names=allTimes)
    allCounts[names(allCounts)%in%names(obsCounts)] <- obsCounts
    # Append empty prior intervals
    if (!is.null(nPriorInt)) allCounts <- c(rep(NA,nPriorInt),allCounts)
    return(allCounts)
}

#' Backcalculate incidence using HIV diagnoses and TID 
#'  
#' Uses the standard convolution equation with a customized EM update step
#' to backcalculate incidence of HIV from diagnoses and the TID aka "pid", time from
#' infection to diagnosis probability function
#'  
#' @param y HIV diagnosis counts vector as formatted by tabulateDiagnoses()
#' @param pid Discrete f(x) giving the probability of diagnosis 
#'        in the current interval given infection x time units prior, like 
#'          that returned by TID_pdf
#' @param gamma Smoothing penalty for EM algorithm
#' @param tol The tolerance for the EM algorithm
#' @param verbose Logical; print intermediate estimates?
#'  
#' @return A list of class "backproj" containing the inputs as well as 
#'          estimated incidence
estimateIncidence <- function(y,pid,gamma=0,tol=10^-5,verbose=FALSE){
  isna <- is.na(y)
  lambda <- y
  impute <- mean(y, na.rm = TRUE)
  for(i in length(y):1){
    if(!isna[i])
      impute <- y[i]
    else
      lambda[i] <- impute
  }
  #lambda <- rep(mean(y,na.rm=TRUE),length(y))
  ll <- lambda
  dev <- Inf
  while(dev>tol){
    lambda <- meanEmUpdate(y,pid,lambda,gamma)
    dev <- sum((ll-lambda)^2/ll)
    if(verbose){
      cat("lambda: ",paste(round(lambda,1),collapse=" "),"\n",
          "parameter change: ",dev,"\n",sep="")
    }
    ll <- lambda
  }
  mod <- list(lambda=lambda,y=y,pid=pid,gamma=gamma,tol=tol)
  class(mod) <- "backproj"
  mod
}

#' EM update for HIV backcalculation
#'  
#' A customized EM update that adjusts for limited HIV surveillance
#' window and imposes a quadratic smoothing penalty
#'  
#' @param y HIV diagnosis counts vector as formatted by tabulateDiagnoses()
#' @param pid Discrete f(x) giving the probability of diagnosis 
#'        in the current interval given infection x time units prior, like 
#'          that returned by TID_pdf
#' @param lambda Vector of starting values, initialized at average # of
#'          observed diagnoses per time interval
#' @param gamma Smoothing penalty for EM algorithm
#'  
#' @return Vector of backcalculated incidences for a given step
meanEmUpdate <- function(y,pid,lambda,gamma=0){
  T <- length(y)
  obs <- !is.na(y)
  a <- b <- c <- lamNew <- rep(NA,T)
  for(k in 1:length(lambda)){
    s <- 0:(T-k)
    b[k] <- sum(pid(s))
    no <- !obs[s+k]
    if(any(no))
      a[k] <- sum(pid(s[no])) / b[k]
    else
      a[k] <- 0
    c[k] <- 0
    for(d in s){
      if(obs[k+d]){
        c[k] <- c[k] + y[k+d]*pid(d) / sum(lambda[1:(k+d)]*pid(k+d-(1:(k+d))))
      }
    }
  }
  if(gamma > 0){
    f <- function(ll){
      (1/ll) * (a*b+c) * lambda - b - 2 * gamma * c(0, ll[2:T] - ll[-T]) - 
        2 * gamma * c(ll[1:(T-1)] - ll[-1] ,0)
    }
    j <- function(ll){
      diag <- (-1/ll^2)* (a*b+c)*lambda - 4 * gamma
      diag[1] <- diag[1] + 2 * gamma
      diag[T] <- diag[T] + 2 * gamma
      off <- rep(0,T)
      off[2:(T-1)] <- 2*gamma
      rbind(-off,diag,off)
    }
    # browser()
    lamNew <- multiroot(f=f,start=lambda,positive=TRUE,jacfunc=j,jactype="bandint")$root
  }else{
    lamNew <- lambda * (a + c / b)
  }
  lamNew
}

#' Print function for "backproj"
#'  
#' Print incidence estimates for class "backproj", the output from
#' estimateIncidence()
#'  
#' @param x Object of class "backproj"
#' @param ... passed to cat
#'  
#' @return Output to console
print.backproj <- function(x,...) {
  cat("Back Projection Incidence Model\n","Estimated incidence: ",x$lambda,...)
}


#' Estimates the number of undiagnosed in the population
#'  
#' Uses backcalculated incidence from estimateIncidence() and the TID to
#' estimate undiagnosed counts per interval
#'  
#' @param mod Object of class "backproj" containing estimated incidence and TID 
#'          function from estimateIncidence()
#' @param nExt the number of time units beyond the end of the observed period to use in calculation
#'
#' @return Vector of undiagnosed counts corresponding to time steps in the vector
#'          created by tabulateDiagnoses() used as input to estimateIncidence()
estimateUndiagnosed <- function(mod,nExt=500){
  pid <- mod$pid
  lambda <- mod$lambda
  y <- mod$y
# nExt <- 500
  obs <- c(!is.na(y),rep(FALSE,nExt))
  T <- length(y)
  P <- T + nExt
  n <- matrix(NA,nrow=T,ncol=P)
  for(i in 1:T){
    for(j in i:P){
      if(!obs[j]){
        n[i,j] <- lambda[i]*pid(j-i)
      }else{
        n[i,j] <- y[j] * lambda[i]*pid(j-i) / sum( lambda[1:j]*pid(j-(1:j)) )
      }
    }
  }
  undiag <- rep(NA,T)
  for(i in 1:T) undiag[i] <- sum(n[1:(i),(i+1):P]) + .5 * sum(n[1:i,i])
  undiag
}



#' Estimates the number of undiagnosed in the population, constant incidence
#'  
#' Uses constant incidence assumption and the TID to estimate undiagnosed 
#' counts per interval
#'  
#' @param infPeriod A vector of continuous times from last HIV test to diagnosis
#'          for a population
#' @param case One of "base_case" or "upper_bound", indicating the 
#'          assumption to apply for when infection occurred within infPeriod
#' @param intLength A single number indicating the length in years of discrete 
#'          time intervals by which HIV diagnoses are recorded. The default of 
#'          0.25 represents a quarter-year.
#' @param incidence Number of new infections per unit intLength. Assumed
#'          to be constant over time
#' @param numSteps Number of discrete time steps used to evaluate the integral
#'
#' @return Scalar representing the number of undiagnosed cases at any
#'          time point. The precise interval length for this estimate 
#'          will be max(infPeriod)/dt
#'          
estimateUndiagnosedConst <- function(infPeriod,case,
                                     intLength,incidence,
                                     numSteps=10000) {

    # Define maximum time and discrete time steps
    maxTime <- max(infPeriod, na.rm=TRUE)
    times <- seq(from=0, to=maxTime, length.out=numSteps)

    # Get continuous survivor function
    Sx <- TID_cdfContinuous(infPeriod,case,intLength,survivor=TRUE)

    # Evaluate the integral
    integral <- sum(sapply(times, function(t) Sx(t) ))

    # Convert time units appropriately
    dt_in_years <- maxTime/numSteps
    incidence_yearly <- incidence/intLength

    # Undiagnosed
    undiag <- incidence_yearly*integral*dt_in_years

    return(undiag)

}
