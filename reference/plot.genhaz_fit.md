# Plot estimated curves from a fitted GH model

Calls
[`predict.genhaz_fit()`](https://aaronjehle.github.io/genhaz/reference/predict.genhaz_fit.md)
and passes the result to
[`plot.genhaz_pred()`](https://aaronjehle.github.io/genhaz/reference/plot.genhaz_pred.md).
Supports all eight prediction types; `ylim` is automatically derived
from the full range of CI bands so no line is clipped.

## Usage

``` r
# S3 method for class 'genhaz_fit'
plot(
  x,
  newdata,
  times,
  type = c("hazard", "survival", "cumhaz", "rmst", "surv_diff", "rmst_diff",
    "hazard_ratio", "time_ratio"),
  alpha = 0.05,
  ...
)
```

## Arguments

- x:

  A `"genhaz_fit"` object from
  [`fit_genhaz()`](https://aaronjehle.github.io/genhaz/reference/fit_genhaz.md).

- newdata:

  A `data.frame` of covariate values. One row per group. Row names
  become legend labels.

- times:

  Numeric vector of evaluation time points.

- type:

  Quantity to plot. One of `"hazard"` (default), `"survival"`,
  `"cumhaz"`, `"rmst"`, `"surv_diff"`, `"rmst_diff"`, `"hazard_ratio"`,
  or `"time_ratio"`. The last four require exactly two rows in
  `newdata`.

- alpha:

  Significance level for confidence bands. Default 0.05.

- ...:

  Additional arguments passed to
  [`plot.genhaz_pred()`](https://aaronjehle.github.io/genhaz/reference/plot.genhaz_pred.md)
  and then to
  [`graphics::plot()`](https://rdrr.io/r/graphics/plot.default.html)
  (e.g. `col`, `lty`, `xlab`, `main`, `xlim`).

## Value

Invisibly returns the `"genhaz_pred"` object from
[`predict.genhaz_fit()`](https://aaronjehle.github.io/genhaz/reference/predict.genhaz_fit.md).

## Examples

``` r
if (FALSE) { # \dontrun{
nd <- data.frame(X = c(0, 1))
rownames(nd) <- c("X = 0", "X = 1")
t_grid <- seq(0.01, 8, length.out = 300)

plot(fit, newdata = nd, times = t_grid)
plot(fit, newdata = nd, times = t_grid, type = "survival",
     col = c("steelblue", "firebrick"))
plot(fit, newdata = nd, times = t_grid, type = "hazard_ratio",
     col = "purple")
} # }
```
