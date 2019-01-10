
Manual for the use of the pransi functions
==========================================

Introduction
------------

The aim of this *pransi* pacakage is to provide Permutation RAndom Null distributions for simultaneous inference. Wilcoxon rank sum test and two sample t-test are natively implemented, but any other test can be used. The package's aim is multiplicity correction while accounting for correlations among statistical tests.

General use
-----------

We illustrate the general use of the package on a synthetic dataset. The default Wilcoxon rank-sum test is used.

``` r
#Create some synthetic data with 90% true null hypothesis
 p = 200; n = 50
x = rep(c(0,1), each = n/2)
 mat = cbind(
 matrix(rnorm(n*p/10, mean = 5+x),n,p/10), #DA
 matrix(rnorm(n*p*9/10, mean = 5),n,p*9/10) #Non DA
 )
 #Provide just the matrix and grouping factor, and test using the random null
 fdrRes = fdrCorrect(mat, x)
  #The estimated tail-area false discovery rates.
  estFdr = fdrRes$Fdr
```

The method provides an estimate of the proportion of true null hypothesis, which is close to the true 90%.

``` r
fdrRes$p0
```

    ## [1] 0.8562802

The result of the procedure can be represented graphically as follows:

``` r
plotNull(fdrRes)
```

![](README_figs/README-plotNull-1.png)

It is also possible to provide a custom test function, which must accept at least a 'y' response variable and a 'x' grouping factor. Additionally, a quantile function should be supplied for transformation through quantiles to z-values.

``` r
# With another type of test
fdrResLm = fdrCorrect(mat, x, B = 50, test = function(x, y) {
    fit = lm(y ~ x)
    c(summary(fit)$coef["x1", "t value"], fit$df.residual)
}, quantileFun = function(q) {
    pt(q = q[1], df = q[2])
})
```

This framework also accepts more than 2 groups, and additional covariates through the

``` r
 #3 groups
 p = 100; n = 60
x = rep(c(0,1,2), each = n/3)
 mat = cbind(
 matrix(rnorm(n*p/10, mean = 5+x),n,p/10), #DA
 matrix(rnorm(n*p*9/10, mean = 5),n,p*9/10) #Non DA
 )
 #Provide an additional covariate through the 'argList' argument
 z = rpois(n , lambda = 2)
 fdrResLmZ = fdrCorrect(mat, x, B = 5e1,
 test = function(x, y, z){
 fit = lm(y~x+z)
 c(summary(fit)$coef["x1","t value"], fit$df.residual)},
 quantileFun = function(q){pt(q = q[1], df = q[2])},
 argList = list(z = z))
```

If the null distribution of the test statistic is not known, it is also possbile to execute the procedure on the scale of the original test statistics, rather than z-values by setting zValues = FALSE. This may be numerically less stable.

``` r
fdrResKruskal = fdrCorrect(mat, x, B = 50, test = function(x, y) {
    kruskal.test(y ~ x)$statistic
}, zValues = FALSE)
```

Case study
----------

We illustrate the package using an application from microbiology. The species composition of a community of microorganisms can be determined through sequencing. However, this only yields compositional information, and knowledge of the population size can be acquired by cell counting through flow cytometry. Next, the obtained species compositions can multiplied by the total population size to yield approximate absolute cell counts per species. Evidently, this introduces strong correlation between the tests due to the common factor. In other words: random noise in the estimation of the total cell counts will affect all hypotheses. Therefore, we employ permutations to estimate an experiment-specific random null distribution, that will account for this dependence.

The dataset used is taken from Vandeputte *et al.*, 2017 (see [manuscript](https://www.ncbi.nlm.nih.gov/pubmed/29143816)), a study on gut microbiome in healthy and Crohn's disease patients. The test looks for differences in absolute abundance between healthy and diseased patients. It relies on the *phyloseq* package, which is the preferred way to interact with our machinery for microbiome data.

``` r
testVanDePutte = testDAA(Vandeputte, "Health.status", "absCountFrozen")
```

The estimated tail-area false discovery rates can then simply be extracted as

``` r
FdrVDP = testVanDePutte$Fdr
quantile(FdrVDP)
```

    ##           0%          25%          50%          75%         100% 
    ## 5.231074e-16 3.551768e-06 8.838016e-02 6.757317e-01 1.000000e+00
