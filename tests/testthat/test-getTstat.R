context("Faster t-test implementation")

test_that("T-test yields same result as built in version", {
    n = 20
    y1 = rnorm(rpois(1,n)+2)
    y2 = rnorm(rpois(1,n)+2 ,mean = 0.2)
    mm = length(y1)
    nn = length(y2)
    tTestBuiltIn = t.test(y1, y2)
    resultVec = c( tTestBuiltIn[["statistic"]],
                   tTestBuiltIn[["parameter"]])
    names(resultVec) = c("tstat", "df")
    expect_equal(reconsi:::getTstat(y1, y2, mm, nn),
                 resultVec)
})
