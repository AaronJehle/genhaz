# Post-estimation from a fitted general hazard model

Evaluates model quantities (hazard, cumulative hazard, survival,
gradients, etc.) at arbitrary covariate patterns and time points, using
the parameters of a fitted model object.

## Usage

``` r
post(fit, X, time, res, event = NULL)
```

## Arguments

- fit:

  A fitted model object returned by
  [`fit_genhaz()`](https://aaronjehle.github.io/genhaz/reference/fit_genhaz.md).

- X:

  Covariate vector or matrix for prediction. If a vector, it is treated
  as a single covariate pattern and replicated across all `time` points.
  If a single-row matrix and `length(time) > 1`, the row is replicated.
  If `NULL`, all covariates are set to zero.

- time:

  Numeric vector of time points at which to evaluate the quantity.

- res:

  Character string; the quantity to return. Supports all non-fitting
  options of
  [`genhaz_work()`](https://aaronjehle.github.io/genhaz/reference/genhaz_work.md):
  `"h"`, `"H"`, `"log_h"`, `"gradient_H"`, `"gradient_h"`, etc.

- event:

  Optional event indicator vector, required for some `res` values (e.g.
  `"scores"`).

## Value

Numeric vector or matrix of the requested quantity evaluated at each
`(time, X)` combination.

## Examples

``` r
if (FALSE) { # \dontrun{
library(survival)
set.seed(42)
dat <- sim_scenario(1, beta1 = 0.5, beta2 = 0.5, n = 300)
fit <- fit_genhaz(Surv(dat$time, dat$event), ~X, data = dat,
                  n_knots = 6, profile = TRUE, tol_LCV = 0.01)

t_grid <- seq(0.01, 8, length.out = 200)
h_est  <- post(fit, X = 0, time = t_grid, res = "h")
plot(t_grid, h_est, type = "l", xlab = "Time", ylab = "Hazard")
} # }
```
