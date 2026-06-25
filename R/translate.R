#' Time-varying acceleration factor between two covariate patterns
#'
#' For two covariate patterns `covariate0` and `covariate1`, computes the
#' time-varying acceleration factor phi(t) = h0(t) / h1(tau(t)), where tau
#' is defined by H0(t) = H1(tau(t)) (equal cumulative-hazard time mapping).
#' Requires **rstpm2** for the vectorised root-finding.
#'
#' @param fit A fitted model object from [fit_genhaz()].
#' @param covariate0 Numeric vector: covariate pattern for the reference group.
#' @param covariate1 Numeric vector: covariate pattern for the comparison group.
#' @param time Numeric vector of evaluation time points.
#'
#' @return Numeric vector of phi(t) values.
#' @export
translate_time <- function(fit, covariate0, covariate1, time) {
  if (!requireNamespace("rstpm2", quietly = TRUE)) {
    stop("Package 'rstpm2' is required for translate_time().")
  }
  tausimp <- function(t) {
    H0t   <- post(fit, covariate0, t, "H")
    H1    <- function(u) post(fit, covariate1, u, "H")
    Hdiff <- function(u) H1(u) - H0t
    rstpm2::vuniroot(Hdiff, lower = 0, upper = max(time),
                     tol = .Machine$double.eps)$root
  }
  tau  <- Vectorize(tausimp, vectorize.args = "t")
  tau_v <- tau(time)
  h0    <- post(fit, covariate0, time,   "h")
  h1    <- post(fit, covariate1, tau_v,  "h")
  h0 / h1
}

#' Time-varying tau mapping between two covariate patterns
#'
#' Returns tau(t), defined by H0(t) = H1(tau(t)), i.e. the time at which
#' the comparison group has accumulated the same cumulative hazard as the
#' reference group at time t. Requires **rstpm2**.
#'
#' @param fit A fitted model object from [fit_genhaz()].
#' @param covariate0 Numeric vector: covariate pattern for the reference group.
#' @param covariate1 Numeric vector: covariate pattern for the comparison group.
#' @param time Numeric vector of evaluation time points.
#'
#' @return Data frame with columns `time` and `tvta` (tau values).
#' @export
translate_time2 <- function(fit, covariate0, covariate1, time) {
  if (!requireNamespace("rstpm2", quietly = TRUE)) {
    stop("Package 'rstpm2' is required for translate_time2().")
  }
  tausimp <- function(t) {
    H0t   <- post(fit, covariate0, t, "H")
    H1    <- function(u) post(fit, covariate1, u, "H")
    Hdiff <- function(u) H1(u) - H0t
    rstpm2::vuniroot(Hdiff, lower = 0, upper = max(time),
                     tol = .Machine$double.eps)$root
  }
  tau  <- Vectorize(tausimp, vectorize.args = "t")
  tau_v <- tau(time)
  data.frame(time = time, tvta = tau_v)
}
