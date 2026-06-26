# General hazard model workhorse

Calculates various quantities from a general hazard (GH) model depending
on the `res` argument: negative log-likelihood and its gradient/Hessian,
the cumulative hazard H and its derivatives, effective degrees of
freedom, the LCV criterion, or a fully fitted model object. Should not
be called directly by users; use
[`fit_genhaz()`](https://aaronjehle.github.io/genhaz/reference/fit_genhaz.md)
and [`post()`](https://aaronjehle.github.io/genhaz/reference/post.md)
instead.

## Usage

``` r
genhaz_work(
  theta,
  time,
  X,
  knots,
  Z,
  event = NULL,
  lambda = 1,
  res = c("log_h", "gradient_log_h", "hessian_log_h", "dlogh_dt", "dh_dt", "h",
    "gradient_h", "hessian_h", "H", "gradient_H", "hessian_H", "negll", "negll_upen",
    "scores", "gradient_negll", "hessian_negll", "hessian_negll_upen", "basis", "edf",
    "LCV", "dLCV", "dLCV_full", "fit", "fit_LCV"),
  model_type = "GH",
  time2 = NULL,
  cens_type = "rc",
  gaussified = TRUE,
  nquad = 25,
  control = list(trace = 0, x.tol = 1e-12, rel.tol = 1e-12, max_iter = 100),
  lcv_method = c("full", "approx", "optimize"),
  interval = c(0, 13),
  tol_LCV = 0.001,
  ...
)
```

## Arguments

- theta:

  Numeric vector of parameters: intercept, spline coefficients, then
  beta1 (AFT-style effects) and beta2 (PH-style effects).

- time:

  Numeric vector of observed (or evaluation) times.

- X:

  Numeric matrix of covariates (rows = observations).

- knots:

  Numeric vector of knot positions on the log-time scale, including
  boundary knots.

- Z:

  Projection matrix for the spline basis (from
  [`smf_cpp()`](https://aaronjehle.github.io/genhaz/reference/smf_cpp.md)).

- event:

  Integer event indicator vector (1 = event, 0 = censored). Required for
  fitting.

- lambda:

  Non-negative smoothing parameter for the spline penalty.

- res:

  Character string selecting the quantity to return. One of `"log_h"`,
  `"gradient_log_h"`, `"hessian_log_h"`, `"dlogh_dt"` (derivative of the
  log-hazard w.r.t. time), `"dh_dt"` (derivative of the hazard w.r.t.
  time), `"h"`, `"gradient_h"`, `"hessian_h"`, `"H"`, `"gradient_H"`,
  `"hessian_H"`, `"negll"`, `"negll_upen"`, `"scores"`,
  `"gradient_negll"`, `"hessian_negll"`, `"hessian_negll_upen"`,
  `"basis"`, `"edf"`, `"LCV"`, `"dLCV"`, `"dLCV_full"`, `"fit"`, or
  `"fit_LCV"`. The last two fit a model: `"fit"` at a fixed `lambda`,
  `"fit_LCV"` after profiling `lambda` via LCV. `"dLCV"` returns the
  first-order analytical derivative of LCV w.r.t. `log(lambda)`
  (evaluated at the MLE for the given `lambda`). `"dLCV_full"` returns
  the same derivative including the third-derivative correction term
  from the unpenalised Hessian's dependence on `lambda` (supports `"rc"`
  and `"lt_rc"`; falls back to `"dLCV"` for `"ic"`).

- model_type:

  Character vector specifying the model type for each covariate: `"GH"`
  (general hazard), `"PH"` (proportional hazards, beta1 = 0), `"AFT"`
  (accelerated failure time, beta2 = 0), or `"AH"` (additive hazards,
  beta1 = -beta2). If length 1, recycled to all covariates.

- time2:

  Numeric vector for the second time variable. Used as the truncation
  time for `cens_type = "lt_rc"`, or as the right interval endpoint for
  `cens_type = "ic"`.

- cens_type:

  Censoring type: `"rc"` (right-censoring, default), `"lt_rc"`
  (left-truncation + right-censoring), or `"ic"` (interval censoring).

- gaussified:

  Logical; use Gauss-Legendre quadrature for integrating the cumulative
  hazard? Default `TRUE`. Setting to `FALSE` requires **rstpm2** to be
  installed.

- nquad:

  Number of quadrature points when `gaussified = TRUE`. Default 25.

- control:

  List of optimisation control parameters passed to
  [`nlminb()`](https://rdrr.io/r/stats/nlminb.html): `trace`, `x.tol`,
  `rel.tol`, `max_iter`.

- lcv_method:

  Character string. When `res = "fit_LCV"`, selects the strategy for
  optimising `log(lambda)`. One of:

  `"full"` (default)

  :   Root-find on the full LCV gradient including the third-derivative
      correction (`dLCV_full`). Most accurate.

  `"approx"`

  :   Root-find on the first-order LCV gradient (`dLCV`), which ignores
      the dependence of the Hessian on `lambda`. Faster than `"full"`.

  `"optimize"`

  :   Directly minimise LCV without gradient information via
      [`optimize()`](https://rdrr.io/r/stats/optimize.html). Equivalent
      to the original behaviour before analytic gradients were added.

  Ignored for all other `res` values and falls back to `lcv_gradient`
  for `cens_type = "ic"`.

- interval:

  Length-2 numeric vector giving the search interval for `log(lambda)`
  during LCV optimisation. Default `c(0, 13)`.

- tol_LCV:

  Convergence tolerance for LCV optimisation. Default `1e-3`.

- ...:

  Additional arguments (currently unused).

## Value

Depends on `res`. When `res = "fit"` or `"fit_LCV"`, returns a list with
components `par`, `se`, `z`, `p_values`, `AIC`, `edf`, `lambda`, `var`,
and further fitting details.
