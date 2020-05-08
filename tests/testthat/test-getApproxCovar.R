context("Binned correlation matrix")

test_that("Binned correlation matrix has the right dimensions", {
    n = 20; p = 50; B = 1e2
    mat = matrix(rnorm(n*p), n, p)
    x = sample(c(0,1), n , replace = TRUE)
    nBins = 80L
    reconsiObj = reconsi(Y = mat, x = x)
    expect_equal(dim(reconsi:::getApproxCovar(reconsiObj, nBins = nBins)),
        rep(nBins+2L, 2))
  })
