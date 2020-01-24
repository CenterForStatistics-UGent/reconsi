context("Weights calculation")

test_that("Weights sum to one",{
    n = 20; p = 50; B = 1e2
    mat = matrix(rnorm(n*p), n, p)
    x = sample(c(0,1), n , replace = TRUE)
    expect_equal(sum(reconsi(mat, x)$weights), 1)
})
