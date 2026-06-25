# Predict from a fitted GH model

Evaluates the fitted model at one or more covariate patterns over a time
grid, returning pointwise estimates with delta-method confidence
intervals. Only the requested quantity is computed.

## Usage

``` r
# S3 method for class 'genhaz_fit'
predict(
  object,
  newdata,
  times,
  type = c("hazard", "survival", "cumhaz", "rmst", "surv_diff", "rmst_diff",
    "hazard_ratio", "time_ratio"),
  alpha = 0.05,
  ...
)
```

## Arguments

- object:

  A `"genhaz_fit"` object from
  [`fit_genhaz()`](https://aaronjehle.github.io/genhaz/reference/fit_genhaz.md).

- newdata:

  A `data.frame` of covariate values. Each row defines one covariate
  pattern; row names are used as group labels in the output. Must
  contain all variables referenced in the model formula. For
  `"surv_diff"` and `"rmst_diff"`, exactly two rows are required.

- times:

  Numeric vector of time points. For `"rmst"` and `"rmst_diff"`, these
  are restriction times tau at which RMST(tau) = integral_0^tau S(t) dt
  is evaluated.

- type:

  Quantity to predict: `"hazard"` (default), `"survival"`, `"cumhaz"`
  (cumulative hazard), `"rmst"` (restricted mean survival time),
  `"surv_diff"` (survival difference S1−S2), `"rmst_diff"` (RMST
  difference), `"hazard_ratio"` (h1/h2), or `"time_ratio"` (tau/t where
  S2(tau) = S1(t)). The last four require exactly two rows in `newdata`.

- alpha:

  Significance level for confidence intervals. Default 0.05.

- ...:

  Currently unused.

## Value

A `data.frame` with columns `pattern`, `time`, `estimate`, `lower`,
`upper`. For per-group types rows are grouped by covariate pattern; for
difference types a single block labelled `"group1 - group2"` is
returned.

## Examples

``` r
if (FALSE) { # \dontrun{
t_grid <- seq(0.01, 8, length.out = 200)
nd <- data.frame(X = c(0, 1))
rownames(nd) <- c("Unexposed", "Exposed")

predict(fit, newdata = nd, times = t_grid, type = "survival")
predict(fit, newdata = nd, times = c(1, 2, 5), type = "rmst")
predict(fit, newdata = nd, times = t_grid, type = "surv_diff")
predict(fit, newdata = nd, times = c(1, 2, 5), type = "rmst_diff")
predict(fit, newdata = nd, times = t_grid, type = "hazard_ratio")
predict(fit, newdata = nd, times = t_grid, type = "time_ratio")
} # }
```
