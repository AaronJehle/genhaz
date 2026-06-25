# Fit a penalized general hazard model

High-level interface for fitting a penalised cubic regression spline
general hazard (GH) model. Wraps
[`genhaz_work()`](https://aaronjehle.github.io/genhaz/reference/genhaz_work.md)
with automatic knot placement, design matrix construction, and optional
profiling of the smoothing parameter `lambda` via the modified LCV
criterion.

## Usage

``` r
fit_genhaz(
  surv,
  formula,
  data,
  n_knots = 8,
  lambda = 0,
  model_type = "GH",
  control = list(trace = 0, x.tol = 1e-12, rel.tol = 1e-12, max_iter = 100),
  profile = FALSE,
  lcv_method = c("full", "approx", "optimize"),
  lambda_surv = FALSE,
  init = NULL,
  interval = c(0, 10),
  nquad = 25,
  tol_LCV = 0.001,
  knots = NULL,
  reKnot = 1,
  limit_to_95 = TRUE,
  timeIt = FALSE
)
```

## Arguments

- surv:

  A `Surv` object (from the **survival** package). Right- censored
  (`Surv(time, event)`), left-truncated and right-censored
  (`Surv(start, stop, event)`), and interval-censored
  (`Surv(t1, t2, type = "interval2")`) data are all supported.

- formula:

  A one-sided formula specifying the covariate terms, e.g. `~ x + z`. An
  intercept is added internally; do not include one here.

- data:

  A data frame containing all variables referenced in `formula`.

- n_knots:

  Total number of knots for the baseline log-hazard spline, including
  the two boundary knots. Default 8.

- lambda:

  Initial (or fixed) smoothing parameter. Ignored when `profile = TRUE`.
  Default 0.

- model_type:

  Character vector of length 1 or `ncol(X)` selecting the model class
  per covariate: `"GH"` (general hazard, default), `"PH"` (proportional
  hazards), `"AFT"` (accelerated failure time), or `"AH"` (additive
  hazards).

- control:

  List of optimisation controls for
  [`nlminb()`](https://rdrr.io/r/stats/nlminb.html): `trace`, `x.tol`,
  `rel.tol`, `max_iter`.

- profile:

  Logical. If `TRUE`, optimise `lambda` via LCV before returning the
  final fit. Default `FALSE`.

- lcv_method:

  Character string passed to
  [`genhaz_work()`](https://aaronjehle.github.io/genhaz/reference/genhaz_work.md)
  when `profile = TRUE`. Selects the strategy for optimising
  `log(lambda)`: `"full"` (default) uses the full third-derivative LCV
  gradient, `"approx"` uses the faster first-order approximation, and
  `"optimize"` minimises LCV directly without gradient information.
  Ignored when `profile = FALSE`.

- lambda_surv:

  Logical. If `TRUE`, initialise `lambda` from a **survPen**
  proportional hazard fit. Requires **survPen** to be installed. Default
  `FALSE`.

- init:

  Optional numeric vector of starting values. If `NULL`, zeros are used.

- interval:

  Length-2 numeric vector giving the search interval for `log(lambda)`
  during LCV optimisation. Default `c(0, 10)`.

- nquad:

  Number of Gauss-Legendre quadrature points. Default 25.

- tol_LCV:

  Convergence tolerance for LCV optimisation. Default `1e-3`.

- knots:

  Optional numeric vector of knot positions on the log-time scale. If
  `NULL` (default), knots are placed automatically using
  [`knot_pattern()`](https://aaronjehle.github.io/genhaz/reference/knot_pattern.md).

- reKnot:

  Non-negative integer. After the initial fit, re-position knots
  `abs(reKnot)` times using the AFT-adjusted times and refit. Default
  `1` (one re-knot pass). Set to `0` to skip.

- limit_to_95:

  Logical passed to
  [`knot_pattern()`](https://aaronjehle.github.io/genhaz/reference/knot_pattern.md).
  Default `TRUE`.

- timeIt:

  Logical. If `TRUE`, print the wall-clock fitting time. Default
  `FALSE`.

## Value

A list of class `"genhaz_fit"` containing (among others):

- `par`:

  Named vector of parameter estimates.

- `se`:

  Standard errors.

- `z`:

  Wald z-statistics.

- `p_values`:

  Two-sided p-values.

- `var`:

  Estimated covariance matrix of `par`.

- `AIC`:

  Akaike information criterion.

- `edf`:

  Effective degrees of freedom.

- `lambda`:

  Smoothing parameter used.

- `knots`:

  Knot vector on the log-time scale.

- `Z`:

  Spline projection matrix.

## Examples

``` r
if (FALSE) { # \dontrun{
library(survival)
set.seed(42)
dat <- sim_scenario(scenario = 1, beta1 = 0.5, beta2 = 0.5, n = 300)
fit <- fit_genhaz(Surv(dat$time, dat$event), ~ X, data = dat,
                  model_type = "GH", profile = TRUE, n_knots = 6,
                  tol_LCV = 0.01)
fit$par
fit$se
} # }
```
