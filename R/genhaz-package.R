#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @useDynLib genhaz, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom survival Surv
#' @importFrom stats as.formula integrate model.matrix nlminb optimize
#'   pchisq pnorm qnorm quantile rbinom runif uniroot
#' @importFrom utils tail
## usethis namespace: end
NULL

# Package-level mutable state for passing knots into survPen formulas.
.pkg_env <- new.env(parent = emptyenv())

# Pre-seed the GL cache for the default nquad so the first fit_genhaz() call
# pays no eigen() cost.
.onLoad <- function(libname, pkgname) {
  assign("25", gauss_legendre(25L), envir = .gl_cache)
}
