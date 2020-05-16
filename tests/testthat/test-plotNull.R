context("return of plotNull function")

n = 20; p = 200; B = 1e2
mat = matrix(rnorm(n*p), n, p)
x = sample(c(0,1), n , replace = TRUE)
idDA = seq_len(10) %in% seq_len(p)
mat[x==0, idDA] = mat[x==0, seq_len(10)] + 0.5
reco  = reconsi(Y = mat, x = x)
test_that("testDAA deals with all zero rows", {
    expect_s3_class(plotNull(reco), "ggplot")
    expect_s3_class(plotNull(reco, idDA = idDA), "ggplot")
})
