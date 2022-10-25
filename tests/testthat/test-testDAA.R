context("input to testDAA")

test_that("testDAA deals with all zero rows", {
    n = 20; p = 50; B = 1e2
    mat = matrix(rpois(n*p, lambda = 3), n, p)
    x = sample(c(0,1), n , replace = TRUE)
    FC = rpois(n, 1000)
    mat[1,] = 0
    expect_no_error(testDAA(Y = mat, x = x, FC = FC))
})
