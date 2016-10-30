
context('Descriptive Graphics and Tables')

# Test that descriptives work with or without paneling and
# with or without MSM vs non-MSM subgroups
data(KCsim)

# plotDiagnoses
test_that('Plot returns ggplot object, no group or panel', {
    p <- plotDiagnoses(KCsim)
    expect_is(p, 'ggplot')
})
test_that('Plot returns ggplot object, group but no panel', {
    p <- plotDiagnoses(KCsim, showlegend=TRUE) + aes(color=race)
    expect_is(p, 'ggplot')
})
test_that('Plot returns ggplot object, group and panel', {
    p <- plotDiagnoses(KCsim, showlegend=TRUE) + aes(color=race) + 
        facet_grid(.~everHadNegTest)
    expect_is(p, 'ggplot')
})

# tabTestHist 
test_that('Variable vector must have names', {
          expect_error(tabTestHist(KCsim, 'yearDx'))
})

# plotTestHist
test_that('Plot returns ggplot object, no group or panel', {
    p <- plotTestHist(KCsim)
    expect_is(p, 'ggplot')
})
test_that('Plot returns ggplot object, group but no panel', {
    p <- plotTestHist(KCsim, group='race') + aes(lty=Group)
    expect_is(p, 'ggplot')
})
test_that('Plot returns ggplot object, group and panel', {
    p <- plotTestHist(KCsim, panel='race') + facet_grid(.~Panel)
    expect_is(p, 'ggplot')
})


