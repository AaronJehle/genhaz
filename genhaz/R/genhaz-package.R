#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @useDynLib genhaz, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom pracma gaussLegendre Trace
#' @importFrom survival Surv
#' @importFrom stats as.formula integrate model.matrix nlminb optimize
#'   pchisq pnorm qnorm quantile rbinom runif uniroot
#' @importFrom utils tail
## usethis namespace: end
NULL

# Package-level mutable state for passing knots into survPen formulas.
.pkg_env <- new.env(parent = emptyenv())
