context("input to testDAA")

test_that("testDAA deals with all zero rows", {
    n = 20; p = 50; B = 1e2
    mat = matrix(rnorm(n*p), n, p)
    x = sample(c(0,1), n , replace = TRUE)
    mat[1,] = 0
    expect_silent(reconsi(Y = mat, x = x))
})
