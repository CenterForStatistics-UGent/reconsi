
``` r
library(corTests)
```

Manual for the use of the corTest functions
===========================================

Introduction
------------

The aim of this pacakage is to provide simultenous inference for correlated hypotheses using random null distributions. These random null distributions are estimated through permutations. Wilcoxon rank sum test and two sample t-test are natively implemented, but any other test can be used.

General use
-----------

We illustrate the general use of the package on a synthetic dataset. The default Wilcoxon rank-sum test is used.

``` r
#Create some synthetic data with 90% true null hypothesis
 p = 100; n = 50
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

    ## [1] 0.8999659

It is also possible to provide a custom test function, and provide a helper function to extract the test statistic, which is all the procedure needs. Note that the helper function must accept the formula interface.

``` r
 #With another type of test
 fdrResLm = fdrCorrect(mat, x, test = "lm", B = 5e1,
 extractFun = function(x){summary(x)$coef["x1","t value"]})
```

Case study
----------

We illustrate the package using an application from microbiology. The species composition of a community of microorganisms can be determined through sequencing. However, this only yields compositional information, and knowledge of the population size can be acquired by cell counting through flow cytometry. Next, the obtained species compositions can multiplied by the total population size to yield approximate absolute cell counts per species. Evidently, this introduces strong correlation between the tests due to the common factor. In other words: random noise in the estimation of the total cell counts will affect all hypotheses. Therefore, we employ permutations to estimate an experiment-specific random null distribution, that will account for this dependence.

The dataset used is taken from Vandeputte *et al.*, 2017 (see ), a study on gut microbiome in healthy and Crohn's disease patients. The test looks for differences in absolute abundance between healthy and diseased patients. It relies on the *phyloseq* package, which is the preferred way to interact with our machinery for microbiome data.

``` r
testVanDePutte = testDAA(Vandeputte, "Health.status", "absCountFrozen")
```
