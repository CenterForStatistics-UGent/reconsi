context("Weights calculation")

n = 20; p = 50; B = 1e2
mat = matrix(rnorm(n*p), n, p)
x = sample(c(0,1), n , replace = TRUE)

test_that("Weights sum to one",{
    expect_equal(sum(rransi(mat, x)$weights), 1)
})
