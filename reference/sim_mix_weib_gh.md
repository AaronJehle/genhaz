# Simulate data from a mixture Weibull GH model

Simulates right-censored survival data from a two-component mixture
Weibull GH model with a binary covariate X ~ Bernoulli(0.5). Censoring
times are drawn from Uniform(0, `tmax`).

## Usage

``` r
sim_mix_weib_gh(
  n,
  p,
  lambda1,
  lambda2,
  gamma1,
  gamma2,
  beta1,
  beta2,
  tmax = 10
)
```

## Arguments

- n:

  Sample size.

- p:

  Mixture proportion.

- lambda1, lambda2:

  Weibull scale parameters.

- gamma1, gamma2:

  Weibull shape parameters.

- beta1:

  AFT coefficient for X.

- beta2:

  PH coefficient for X.

- tmax:

  Administrative censoring time (uniform censoring upper bound).

## Value

Data frame with columns `time`, `X`, `event`, and `T_true`.

## Examples

``` r
set.seed(42)
dat <- sim_mix_weib_gh(n = 200, p = 0.8, lambda1 = 0.1, lambda2 = 0.1,
                       gamma1 = 3, gamma2 = 1.6, beta1 = 0.5, beta2 = 0.5)
head(dat)
#>        time X event   T_true
#> 1 0.2270001 1     0 0.542651
#> 2 1.0010955 1     1 1.001096
#> 3 1.1876568 0     1 1.187657
#> 4 1.0803505 1     1 1.080351
#> 5 1.4679112 1     1 1.467911
#> 6 1.0798707 1     0 1.080861
```
