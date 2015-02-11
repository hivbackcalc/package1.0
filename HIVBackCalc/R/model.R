
#' Tabulates HIV diagnoses by interval, as input for estimateIncidence().
#'  
#' HIV diagnosis counts for every interval within the calendar time period
#' of the data are appended to empty intervals that will indicate to 
#' estimateIncidence() how far back to project incidence counts.
#'  
#' @param testhist Object of class 'testinghistories' containing the 
#'          time of diagnosis within "timeDx"
#' @param intLength Interval length by which diagnoses are tracked;
#'          1=1 year.
#' @param nPriorInt Number of incidence intervals to backproject
#'          prior to earliest observed diagnoses. Use NULL to skip
#'          this and simply return observed diagnosis counts.
#'  
#' @return Vector of NA's plus diagnosis counts
tabulateDiagnoses <- function(testhist, intLength, nPriorInt=100) {

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
  lambda <- rep(mean(y,na.rm=TRUE),length(y))
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



