context("input to testDAA")

n = 20; p = 50; B = 1e2
mat = matrix(rnorm(n*p), n, p)
x = sample(c(0,1), n , replace = TRUE)
mat[1,] = 0
test_that("testDAA deals with all zero rows", {
    expect_silent(reconsi(Y = mat, x = x))
})
