context("Faster Wilcoxon rank sum test implementation")

test_that("Wilcoxon rank sum test yields same result as built in version", {
    n = 20; p = 50
    mat = matrix(rnorm(n*p), n, p)
    x = sample(c(0,1), n , replace = TRUE)
    expect_equal(reconsi:::getTestStats(mat, center = FALSE,
                                         x = x, B = 2L, replace = FALSE,
                                        tieBreakRan = FALSE)$statObs,
                 apply(mat, 2, function(y){
                     wilcox.test(y~x)$statistic
                    }))
})
