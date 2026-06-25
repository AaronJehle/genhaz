# Plot estimated curves from a fitted GH model

Plots hazard, survival, or cumulative hazard curves for one or more
covariate patterns. Each row of `newdata` becomes a separate coloured
line; row names serve as legend labels. Confidence bands are drawn as
dashed lines of the same colour.

## Usage

``` r
# S3 method for class 'genhaz_fit'
plot(
  x,
  newdata,
  times,
  type = c("hazard", "survival", "cumhaz"),
  alpha = 0.05,
  col = NULL,
  lty = 1,
  xlab = "Time",
  ylab = NULL,
  main = NULL,
  legend = TRUE,
  ci = TRUE,
  ...
)
```

## Arguments

- x:

  A `"genhaz_fit"` object from
  [`fit_genhaz()`](https://aaronjehle.github.io/genhaz/reference/fit_genhaz.md).

- newdata:

  A `data.frame` of covariate values. One row per group. Set row names
  for informative legend labels.

- times:

  Numeric vector of evaluation time points.

- type:

  Quantity to plot: `"hazard"` (default), `"survival"`, or `"cumhaz"`.

- alpha:

  Significance level for confidence bands. Default 0.05.

- col:

  Colour vector (recycled over groups). Defaults to `2, 3, 4, ...` (R's
  standard colour sequence, skipping black).

- lty:

  Line type for point-estimate curves. Default `1`.

- xlab:

  X-axis label. Default `"Time"`.

- ylab:

  Y-axis label. Derived from `type` when `NULL`.

- main:

  Plot title. Derived from `type` when `NULL`.

- legend:

  Logical; draw a legend when multiple groups are present? Default
  `TRUE`.

- ci:

  Logical; overlay dashed confidence band lines? Default `TRUE`.

- ...:

  Additional graphical parameters passed to
  [`graphics::plot()`](https://rdrr.io/r/graphics/plot.default.html).

## Value

Invisibly returns the `data.frame` produced by
[`predict.genhaz_fit()`](https://aaronjehle.github.io/genhaz/reference/predict.genhaz_fit.md).

## Examples

``` r
if (FALSE) { # \dontrun{
t_grid <- seq(0.01, 8, length.out = 300)

# Hazard (default) for two groups
nd <- data.frame(X = c(0, 1))
rownames(nd) <- c("X = 0", "X = 1")
plot(fit, newdata = nd, times = t_grid)

# Survival for a single group
plot(fit, newdata = data.frame(X = 0), times = t_grid, type = "survival")
} # }
```
