
# estimateTID
test_that('infPeriod must be a vector', {
    expect_error(estimateTID(KCsim, 1),
                 'infPeriod must be a vector, not a data frame')
})
