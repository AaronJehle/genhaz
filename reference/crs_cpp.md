# C++ reimplementation of survPen::crs, including derivatives.

The re-implementation uses the same algorithm but uses different
defaults, and allows for derivatives. For the changes in defaults: the
function here assumes that the intercept is false (as per
[`splines::ns`](https://rdrr.io/r/splines/ns.html)); when the knots are
not specified, the function defaults to using `x` rather than
`unique(x)` for the quantiles; and the `df` argument comes before the
`knots` argument. To use the
[`survPen::crs`](https://rdrr.io/pkg/survPen/man/crs.html) defaults in
R, use `crs_cpp(..., intercept=TRUE, survPen_compatible=TRUE)`.

## Usage

``` r
crs_cpp(
  x,
  df = 10L,
  knots = NULL,
  derivs = 0L,
  intercept = FALSE,
  survPen_compatible = FALSE
)
```

## Arguments

- x:

  a vec of values to provide

- df:

  a non-negative integer for the degrees of freedom (ignored if the
  knots are specified)

- knots:

  a vec of knots (defaults to NULL)

- derivs:

  a non-negative integer for the derivative (default 0)

- intercept:

  a bool for whether to include the intercept

- survPen_compatible:

  a bool for whether to use `unique(x)` for determining the knots if the
  `knots` are not specified

## Value

an `Rcpp::NumericMatrix` with attributes `class="crs_cpp"`, `pen` for
the penalty matrix, `knots` for the knots, `derivs` for the order of
derivatives, and `intercept` for a bool/logical for whether the
intercept is included.

## Examples

``` r
x     <- seq(0, 1, length.out = 11)
knots <- c(1, 3, 5, 7) / 10
crs1  <- crs_cpp(x, knots = knots)
dim(crs1)
#> [1] 11  3
# Comparison with survPen::crs requires the survPen package:
if (FALSE) { # \dontrun{
library(survPen)
crs2 <- crs(x, knots, intercept = FALSE)
range(crs1 - crs2$bs)
range(attr(crs1, "pen") - crs2$pen)
} # }
```
