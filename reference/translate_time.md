# Time-varying acceleration factor between two covariate patterns

For two covariate patterns `covariate0` and `covariate1`, computes the
time-varying acceleration factor phi(t) = h0(t) / h1(tau(t)), where tau
is defined by H0(t) = H1(tau(t)) (equal cumulative-hazard time mapping).
Requires **rstpm2** for the vectorised root-finding.

## Usage

``` r
translate_time(fit, covariate0, covariate1, time)
```

## Arguments

- fit:

  A fitted model object from
  [`fit_genhaz()`](https://aaronjehle.github.io/genhaz/reference/fit_genhaz.md).

- covariate0:

  Numeric vector: covariate pattern for the reference group.

- covariate1:

  Numeric vector: covariate pattern for the comparison group.

- time:

  Numeric vector of evaluation time points.

## Value

Numeric vector of phi(t) values.
