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
    "hazard_ratio", "time_ratio", "acc_factor"),
  interval = c("none", "confidence"),
  level = 0.95,
  tau_max = NULL,
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
  contain all variables referenced in the model formula. For the
  two-group comparison types, exactly two rows are required: **row 1 is
  the exposed group (group 1)** and **row 2 the unexposed/baseline group
  (group 0)**.

- times:

  Numeric vector of time points. For `"rmst"` and `"rmst_diff"`, these
  are restriction times tau at which RMST(tau) = integral_0^tau S(t) dt
  is evaluated.

- type:

  Quantity to predict: `"hazard"` (default), `"survival"`, `"cumhaz"`
  (cumulative hazard), `"rmst"` (restricted mean survival time),
  `"surv_diff"` (survival difference S1−S0), `"rmst_diff"` (RMST
  difference), `"hazard_ratio"` (h1/h0), `"time_ratio"` (tau/t where
  S0(tau) = S1(t): the baseline time tau at which the baseline group
  reaches the comparison group's survival at t), or `"acc_factor"` (the
  time-varying acceleration *effect* on the log scale, f(t) = log
  tau'(t) = log h1(t) - log h0(tau), where tau is defined by H0(tau) =
  H1(t); this satisfies S(t \| x1) = S0(integral_0^t exp(f(u) X) du),
  i.e. tau(t) = integral_0^t exp(f(u) X) du with instantaneous rate
  exp(f(t) X) = tau'(t), so f(t) = beta1 in the constant-AFT limit and
  is directly comparable to beta1). The last five require exactly two
  rows in `newdata`.

- interval:

  Type of interval to return, following
  [`stats::predict.lm()`](https://rdrr.io/r/stats/predict.lm.html).
  `"none"` (default) computes only the point `estimate` (skipping the
  delta-method gradients and any root-solving overhead they require) and
  returns the `lower`/`upper` columns as `NA`; `"confidence"` adds a
  delta-method Wald confidence interval.

- level:

  Confidence level for `interval = "confidence"`. Default 0.95.

- tau_max:

  Upper bound for the equal-cumulative-hazard time tau solved by
  `"time_ratio"` and `"acc_factor"`. Defaults to the model support
  `exp(max(object$knots))` (the largest knot, beyond which the spline is
  pure linear extrapolation). When the baseline group does not reach the
  comparison group's cumulative hazard within `[0, tau_max]`, the warp
  is not identified and the corresponding rows are returned as `NA` with
  a warning; pass a larger value to allow extrapolation. Ignored by all
  other types.

- ...:

  Currently unused.

## Value

A `data.frame` of class `c("genhaz_pred", "data.frame")` with columns
`pattern`, `time`, `estimate`, `lower`, `upper` (the latter two are `NA`
when `interval = "none"`). A `pred_type` attribute records which
quantity was computed;
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) uses this to
set axis labels automatically. For per-group types rows are grouped by
covariate pattern; for two-group types a single block labelled
`"group1 - group0"` (exposed - unexposed) is returned.

## Examples

``` r
if (FALSE) { # \dontrun{
t_grid <- seq(0.01, 8, length.out = 200)
nd <- data.frame(X = c(1, 0))
rownames(nd) <- c("Exposed", "Unexposed")

predict(fit, newdata = nd, times = t_grid, type = "survival",
        interval = "confidence")
predict(fit, newdata = nd, times = c(1, 2, 5), type = "rmst")
predict(fit, newdata = nd, times = t_grid, type = "surv_diff")
predict(fit, newdata = nd, times = c(1, 2, 5), type = "rmst_diff")
predict(fit, newdata = nd, times = t_grid, type = "hazard_ratio")
predict(fit, newdata = nd, times = t_grid, type = "time_ratio")
predict(fit, newdata = nd, times = t_grid, type = "acc_factor")
} # }
```
