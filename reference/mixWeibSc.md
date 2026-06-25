# Mixture Weibull model for a named simulation scenario

Convenience wrapper around
[`mixWeib()`](https://aaronjehle.github.io/genhaz/reference/mixWeib.md)
for the three simulation scenarios used in the accompanying study:

- Scenario 1:

  Bathtub-like hazard (p=0.8, lambda1=lambda2=0.1, gamma1=3,
  gamma2=1.6).

- Scenario 2:

  Unimodal hazard (p=0.5, lambda1=lambda2=1, gamma1=1.5, gamma2=0.5).

- Scenario 3:

  Bimodal hazard (p=0.26, lambda1=0.02, lambda2=0.5, gamma1=3,
  gamma2=0.7).

## Usage

``` r
mixWeibSc(scenario, res = "S", t, X = 0, beta1 = 0, beta2 = 0)
```

## Arguments

- scenario:

  Integer (1, 2, or 3) selecting the scenario.

- res:

  Character: `"S"`, `"H"`, or `"h"`.

- t:

  Numeric vector of time points.

- X:

  Covariate value. Default 0.

- beta1:

  AFT coefficient. Default 0.

- beta2:

  PH coefficient. Default 0.

## Value

Numeric vector.

## Examples

``` r
t <- seq(0.01, 5, length.out = 100)
h1 <- mixWeibSc(1, "h", t, X = 1, beta1 = 0.5, beta2 = 0.5)
```
