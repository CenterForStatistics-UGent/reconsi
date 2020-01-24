context("rowMultiply")

test_that("Rowmultiply works", {
    n = 50; p = 100
    mat = matrix(rnorm(n*p), n,p)
    vec = rnorm(p)
    expect_equal(rowMultiply(mat, vec), t(t(mat)*vec))
})
