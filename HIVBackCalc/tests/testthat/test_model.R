
context('Tests for functions in model.R')


test_that('tabulateDiganoses for 0.083 interval', {

    Set up "A", diagInterval = one month or 0.083 years, versus "B", standard diagInterval of 0.25 
    infPeriod = c(0.083, 0.083, 0.083, 1, 1.5, 2, 10, 15, 20, 2, 5, 0.083, 0.5)
    timeDx.A = 2000 + sample((1:10)*.083, size=13, replace=TRUE)
    timeDx.B = 2000 + sample(c(0.25, 0.5, 0.75), size=13, replace=TRUE)
    df.A = data.frame(infPeriod, timeDx=timeDx.A)
    df.B = data.frame(infPeriod, timeDx=timeDx.B)


    # tabulateDiagnoses works for standard intervals
    diagCounts = tabulateDiagnoses(df.B, intLength=0.25)


    # To fix
    # diagCounts = tabulateDiagnoses(df.A, intLength=0.083)
    expect_error(tabulateDiagnoses(df.A, intLength=0.083,
                                   "Error in aggregateDiagnoses(testhist$timeDx, intLength) : object 'timeDxNew' not found"))



    # To fix
    # diagCounts = tabulateDiagnoses(df.A, intLength=0.25)
    expect_error(tabulateDiagnoses(df.A, intLength=0.25, 
                                   "Error in aggregateDiagnoses(testhist$timeDx, intLength) : object 'timeDxNew' not found"))

})
