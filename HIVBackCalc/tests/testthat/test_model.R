
context('Tests for functions in model.R')


test_that('tabulateDiagnoses and aggregateDiagnoses for monthly interval', {

    # Load packages if using this outside of the test framework
    # library(HIVBackCalc)
    # library(testthat)

    # Set up "A", diagInterval = approx one month or 0.083 years, 
    # versus "B", standard diagInterval of 0.25, 
    # and "C", diagInterval = 1/12 exactly one month
    infPeriod = c(0.083, 0.083, 0.083, 1, 1.5, 2, 10, 15, 20, 2, 5, 0.083, 0.5)
    timeDx.A = 2000 + sample((1:10)*.083, size=13, replace=TRUE)
    timeDx.B = 2000 + sample(c(0.25, 0.5, 0.75), size=13, replace=TRUE)
    timeDx.C = 2000 + sample((1:10)*(1/12), size=13, replace=TRUE)
    df.A = data.frame(infPeriod, timeDx=timeDx.A)
    df.B = data.frame(infPeriod, timeDx=timeDx.B)
    df.C = data.frame(infPeriod, timeDx=timeDx.C)


    # For monthly interval, use intLength=1/12 AND specify data as 
    # in 1/12ths, not rounded decimals. I don't quite understand what
    # level of rounding is being employed upon the return of the object,
    # but the test below should pass if the rounding is pretty minimal
    aggCounts <- aggregateDiagnoses(df.C$timeDx, intLength=1/12)
    months <- (aggCounts - floor(aggCounts))/(1/12)
    expect_true(mean(months-round(months))<.0001)

    aggCounts <- aggregateDiagnoses(df.C$timeDx, intLength=0.25)
    months <- (aggCounts - floor(aggCounts))/(1/12)
    expect_true(mean(months-round(months))<.0001)

    # Aggregates to quarterly time steps correctly (mod of 0.25 is 0)
    expect_equal(unique(aggregateDiagnoses(df.C$timeDx, intLength=0.25)%%0.25), 0)


    # Standard intervals
    diagCounts = tabulateDiagnoses(df.B, intLength=0.25)
    diagCounts = tabulateDiagnoses(df.C, intLength=0.25)

    # Monthly interval
    diagCounts = tabulateDiagnoses(df.C, intLength=1/12)

})

test_that('estimateIncidence works for monthly interval', {

    # Create KCsimM, fake data in monthly time steps
    # Repeat the data several times, because otherwise we get too many 
    # months with zero diagnoses and can't estimate
    library(HIVBackCalc)
    data(KCsim)
    KCsimM <- rbind(KCsim, KCsim, KCsim, KCsim)
    KCsimM <- transform(KCsimM, timeDx=yearDx +
                       sample((1:10)*(1/12), 
                              size=nrow(KCsimM), replace=TRUE))


    diagCounts = tabulateDiagnoses(KCsim, intLength=1)
    diagCountsM = tabulateDiagnoses(KCsimM, intLength=1/12)

    # Eh, just remove all the zeroes and replace with 80
    diagCountsM[diagCountsM<1] <- 80

    TIDs <- estimateTID(KCsim$infPeriod, intLength=1)
    TIDsM <- estimateTID(KCsimM$infPeriod, intLength=1/12)

    incidenceBase = estimateIncidence(y=diagCounts, 
                                      pid=TIDs[['base_case']]$pdffxn, 
                                      gamma=0.1, 
                                      verbose=FALSE)
    incidenceBaseM = estimateIncidence(y=diagCountsM, 
                                      pid=TIDsM[['base_case']]$pdffxn, 
                                      gamma=0.1, 
                                      verbose=TRUE)
})
