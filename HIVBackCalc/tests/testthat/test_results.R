
context('Results summaries')


data(KCsim)
data(KCplwh)
plwh <- KCplwh[,c('Year', 'Total')]

# runSubgroups
test_that('Cases vector must have names', {
    expect_error(runSubgroups(KCsim, subvar='fakeGroup', intLength=1, 
                            cases=c('base_case')), 
                 'Cases vector must have names')
})

# plot.results
test_that('x for plot.results must be of class results', {
    undx <- runBackCalc(KCsim, 1, cases=c(Base='base_case'))
    expect_error(plot.results(undx), 'x is not of class results; look at x and make sure you have the right list element')
})
test_that('Panel variable must be in results data', {
    undx <- runSubgroups(KCsim, subvar='fakeGroup', intLength=1, 
                            cases=c(`Base Case`='base_case'))
    toplot <- compileSubgroups(undx)
    expect_error(plot(toplot, panel='sex'), 'Panel variable not in data')
})

# calcTruePrev
test_that('x for calcTruePrev must be of class results', {
    undx <- runBackCalc(KCsim, 1, cases=c(Base='base_case'))
    expect_error(calcTruePrev(undx, plwh), 'x is not of class results; look at x and make sure you have the right list element')
})

# plotTruePrev
test_that('Must have Base Case and Upper Bound, named that way', {
    undx <- runBackCalc(KCsim, 1, cases=c(`Base Case`='base_case'))
    trueprev <- calcTruePrev(undx$results, KCplwh[,c('Year', 'Total')])
    expect_error(plotTruePrev(trueprev), 'Must have Base Case and Upper Bound, named that way')
})

