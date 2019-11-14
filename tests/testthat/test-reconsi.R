context("input to reconsi")

n = 20; p = 50; B = 1e2
mat = matrix(rnorm(n*p), n, p)
x = sample(c(0,1), n , replace = TRUE)

test_that("reconsi throws errors to wrong input types", {
  expect_error(reconsi(Y = mat, x = x,
                          zVals = list(zValObs = rnorm(p),
                                         zValsPerm = matrix(rnorm(p*B), B, p)),
                          zValues = FALSE))
    #Zvalues supplied but not used
expect_error(fdrResLm = reconsi(mat, x, B = 5e1,
                                   test = function(x, y){
                                   c(summary(lm(y~x))$coef["x1","t value"],
                                     fit$df.residual)},
                                   quantileFun = function(t){
                                       pt(q = t[1], df = t[2])}))
#Wrong argument name for quantile function (t instead of q)
expect_error(reconsi(Y = mat, x = sample(x, p, replace = TRUE)))
#Array sizes do not correspond
expect_error(reconsi(Y = mat[,1], x = sample(x, p, replace = TRUE)))
#Vector rather than matrix supplied.
})

test_that("reconsi throws warnings where necessary", {
    expect_error(reconsi(Y = mat[,seq_len(25)],
                            x = sample(x, p, replace = TRUE)))
    #Few hypotheses
})
test_that("reconsi works smoothly when all is fine", {
    expect_silent(reconsi(Y = mat, x = x))
})
test_that("reconsi works when bootstrap is requested", {
              expect_silent(reconsi(Y = mat, test = function(y, x, mu){
                  testRes = t.test(y, mu = mu)
                  c(testRes$statistic, testRes$parameter)}, argList = list(mu = 0), center = FALSE,
                  distFun = function(q){pt(q = q[1],
                                           df = q[2])},
                                  warnConvergence = FALSE))
    expect_silent(reconsi(Y = mat, test = function(y, x){
        wilcox.test(y, mu = 0, exact = FALSE)$statistic
    }, distFun = "psignrank", testPargs = list(n = nrow(mat)),
        warnConvergence = FALSE))
          })
