#' Post-estimation from a fitted general hazard model
#'
#' Evaluates model quantities (hazard, cumulative hazard, survival, gradients,
#' etc.) at arbitrary covariate patterns and time points, using the parameters
#' of a fitted model object.
#'
#' @param fit A fitted model object returned by [fit_genhaz()].
#' @param X Covariate vector or matrix for prediction. If a vector, it is
#'   treated as a single covariate pattern and replicated across all `time`
#'   points. If a single-row matrix and `length(time) > 1`, the row is
#'   replicated. If `NULL`, all covariates are set to zero.
#' @param time Numeric vector of time points at which to evaluate the quantity.
#' @param res Character string; the quantity to return. Supports all non-fitting
#'   options of [genhaz_work()]: `"h"`, `"H"`, `"log_h"`, `"gradient_H"`,
#'   `"gradient_h"`, etc.
#' @param event Optional event indicator vector, required for some `res` values
#'   (e.g. `"scores"`).
#'
#' @return Numeric vector or matrix of the requested quantity evaluated at
#'   each `(time, X)` combination.
#' @export
#'
#' @examples
#' \dontrun{
#' library(survival)
#' set.seed(42)
#' dat <- sim_scenario(1, beta1 = 0.5, beta2 = 0.5, n = 300)
#' fit <- fit_genhaz(Surv(dat$time, dat$event), ~X, data = dat,
#'                   n_knots = 6, profile = TRUE, tol_LCV = 0.01)
#'
#' t_grid <- seq(0.01, 8, length.out = 200)
#' h_est  <- post(fit, X = 0, time = t_grid, res = "h")
#' plot(t_grid, h_est, type = "l", xlab = "Time", ylab = "Hazard")
#' }
post <- function(fit, X, time, res, event = NULL) {
  if (is.null(X)) {
    X_mat <- matrix(rep(0, length(time)), ncol = fit$nb1, byrow = TRUE)
    return(genhaz_work(fit$par, time, X_mat, knots = fit$knots, Z = fit$Z,
                       lambda = fit$lambda, event = event,
                       model_type = fit$model_type, res = res,
                       control = list(trace = 0, x.tol = 1e-12,
                                      rel.tol = 1e-12)))
  }
  if (!is.matrix(X)) {
    # single covariate pattern: replicate for each time point
    X_mat <- matrix(rep(X, length(time)), ncol = length(X), byrow = TRUE)
  } else if (is.matrix(X) && length(time) == 1) {
    # single time, multiple covariate patterns
    time <- rep(time, nrow(X))
    X_mat <- X
  } else if (is.matrix(X) && nrow(X) == 1) {
    # single row, replicate for all time points
    X_mat <- matrix(rep(X, each = length(time)), ncol = ncol(X))
  } else {
    X_mat <- X
  }
  genhaz_work(fit$par, time, X_mat, knots = fit$knots, Z = fit$Z,
              lambda = fit$lambda, event = event,
              model_type = fit$model_type, res = res,
              control = list(trace = 0, x.tol = 1e-12, rel.tol = 1e-12))
}
