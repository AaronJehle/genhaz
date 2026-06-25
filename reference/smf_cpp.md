# C++ reimplementation of `smf`, including derivatives.

The re-implementation uses the same algorithm but uses different
defaults. Note that the inclusion of an intercept is the same as
`crs_cpp`, while no intercept uses a projection based on a QR
decomposition. For the changes in defaults: the function here assumes
that the intercept is false (as per
[`splines::ns`](https://rdrr.io/r/splines/ns.html)); and when the knots
are not specified, the function defaults to using `x` rather than
`unique(x)` for the quantiles. The `survPen` implementation for `smf`
assumes no intercept unless there is a `by` variable. A `by` variable
`x2` would be implemented here by an interaction, such as
`smf_cpp(x,knots,intercept=TRUE):x2`; see the examples.

## Usage

``` r
smf_cpp(
  x,
  df = 10L,
  knots = NULL,
  intercept = FALSE,
  derivs = 0L,
  Z = NULL,
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

- intercept:

  a bool for whether to include the intercept

- derivs:

  a non-negative integer for the derivative (default 0)

- Z:

  optional pre-computed QR projection matrix (arma::mat). If provided,
  overrides the internally computed centering matrix.

- survPen_compatible:

  a bool for whether to use `unique(x)` for determining the knots if the
  `knots` are not specified

## Value

an `Rcpp::NumericMatrix` with attributes `class="smf_cpp"`, `pen` for
the penalty matrix, `knots` for the knots, `derivs` for the order of
derivatives, and `intercept` for a bool/logical for whether the
intercept is included. If there is no intercept, then there will also be
a `Z` attribute which is the QR-based projection matrix.

## Examples

``` r
x     <- seq(0, 1, length.out = 11)
knots <- c(1, 3, 5, 7) / 10
smf1  <- smf_cpp(x, knots = knots)
dim(smf1)
#> [1] 11  3
# Check that Z attribute is available without intercept:
!is.null(attr(smf1, "Z"))
#> [1] TRUE

## Derivative check via finite differences (knots must be explicit)
smfD <- function(t, knots, eps = 1e-5) {
  (smf_cpp(t + eps, knots = knots) - smf_cpp(t - eps, knots = knots)) / (2 * eps)
}
x2 <- seq(0, 1, length.out = 101)
max(abs(smf_cpp(x2, knots = knots, derivs = 1) - smfD(x2, knots = knots)))
#> [1] 11.83606

## Comparison with survPen::smooth.cons (requires survPen):
if (FALSE) { # \dontrun{
library(survPen)
smf2 <- smooth.cons("x", knots = list(knots), df = 4, option = "smf",
                    data.spec = data.frame(x = x), name = "smf(time)")
range(smf1 - smf2$X)
range(attr(smf1, "pen") - smf2$pen[[1]])
} # }
```
