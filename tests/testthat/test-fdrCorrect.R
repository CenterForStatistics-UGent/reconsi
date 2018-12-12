context("input to fdrCorrect")

n = 20; p = 50; B = 1e2
mat = matrix(rnorm(n*p), n, p)
x = sample(c(0,1), n , replace = TRUE)

test_that("fdrCorrect throws errors to wrong input types", {
  expect_error(fdrCorrect(Y = mat, x = x,
                          zVals = list(zValObs = rnorm(p),
                                         zValsPerm = matrix(rnorm(p*B), B, p)),
                          zValues = FALSE))
    #Zvalues supplied but not used
expect_error(fdrResLm = fdrCorrect(mat, x, B = 5e1,
                                   test = function(x, y){
                                   c(summary(lm(y~x))$coef["x1","t value"],
                                     fit$df.residual)},
                                   quantileFun = function(t){
                                       pt(q = t[1], df = t[2])}))
#Wrong argument name for quantile function (t instead of q)
expect_error(fdrCorrect(Y = mat, x = sample(x, p, replace = TRUE)))
#Array sizes do not correspond
expect_error(fdrCorrect(Y = mat[,1], x = sample(x, p, replace = TRUE)))
#Vector rather than matrix supplied.
})

test_that("fdrCorrect throws warnings where necessary", {
    expect_error(fdrCorrect(Y = mat[,seq_len(25)],
                            x = sample(x, p, replace = TRUE)))
    #Few hypotheses
})
