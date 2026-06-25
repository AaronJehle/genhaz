# Time-varying tau mapping between two covariate patterns

Returns tau(t), defined by H0(t) = H1(tau(t)), i.e. the time at which
the comparison group has accumulated the same cumulative hazard as the
reference group at time t. Requires **rstpm2**.

## Usage

``` r
translate_time2(fit, covariate0, covariate1, time)
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

Data frame with columns `time` and `tvta` (tau values).
