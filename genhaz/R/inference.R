#' Likelihood ratio test between nested GH models
#'
#' Computes a chi-squared likelihood ratio test statistic comparing a nested
#' (restricted) model against a more general model.
#'
#' @param fit_nested Fitted model object (the restricted model).
#' @param fit_general Fitted model object (the general model).
#'
#' @return Named numeric vector with `LR-statistic` and `p_value`.
#' @export
#'
#' @examples
#' \dontrun{
#' fit_ph <- fit_genhaz(Surv(dat$time, dat$event), ~X,
#'                      data = dat, model_type = "PH", n_knots = 6)
#' fit_gh <- fit_genhaz(Surv(dat$time, dat$event), ~X,
#'                      data = dat, model_type = "GH", n_knots = 6)
#' LR(fit_ph, fit_gh)
#' }
LR <- function(fit_nested, fit_general) {
  lr  <- -2 * (-fit_nested$objective + fit_general$objective)
  pv  <- pchisq(lr, df = fit_general$df - fit_nested$df, lower.tail = FALSE)
  res <- c(lr, pv)
  names(res) <- c("LR-statistic", "p_value")
  res
}

#' Wald confidence interval for a single parameter
#'
#' @param fit A fitted model object from [fit_genhaz()].
#' @param param Character; name of the parameter (must be in `fit$parnames`).
#' @param alpha Significance level. Default 0.05 (95 % CI).
#'
#' @return Named numeric vector with `lower` and `upper` bounds.
#' @export
#'
#' @examples
#' \dontrun{
#' waldCI(fit, "beta2_X")
#' }
waldCI <- function(fit, param, alpha = 0.05) {
  if (!(param %in% fit$parnames)) {
    stop("'", param, "' is not a parameter of the fitted model.")
  }
  z <- qnorm(1 - alpha / 2)
  lower <- fit$par[param] - z * fit$se[param]
  upper <- fit$par[param] + z * fit$se[param]
  c(lower = unname(lower), upper = unname(upper))
}

#' Wald confidence interval for the difference of two parameters
#'
#' Uses the delta method to compute a confidence interval for
#' `param1 - param2` accounting for their covariance.
#'
#' @param fit A fitted model object from [fit_genhaz()].
#' @param param1 Character; name of the first parameter.
#' @param param2 Character; name of the second parameter.
#' @param alpha Significance level. Default 0.05 (95 % CI).
#'
#' @return Named numeric vector with `lower` and `upper` bounds.
#' @export
#'
#' @examples
#' \dontrun{
#' waldCI_minus(fit, "beta1_X", "beta2_X")
#' }
waldCI_minus <- function(fit, param1, param2, alpha = 0.05) {
  cov_m <- fit$var[c(param1, param2), c(param1, param2)]
  var   <- cov_m[1, 1] - 2 * cov_m[1, 2] + cov_m[2, 2]
  z     <- qnorm(1 - alpha / 2)
  diff  <- fit$par[param1] - fit$par[param2]
  c(lower = unname(diff - z * sqrt(var)),
    upper = unname(diff + z * sqrt(var)))
}

#' Pointwise confidence bands for hazard and survival functions
#'
#' Computes pointwise 95 % (or other level) confidence intervals for the
#' cumulative hazard H(t), survival S(t) = exp(-H(t)), and hazard h(t),
#' at a given covariate pattern, using the delta method on the log scale.
#'
#' @param fit A fitted model object from [fit_genhaz()].
#' @param t Numeric vector of evaluation time points.
#' @param covariate Numeric vector of covariate values (one per covariate in
#'   the model).
#' @param alpha Significance level. Default 0.05 (95 % CI).
#'
#' @return A data frame with columns `time`, `H`, `lower_H`, `upper_H`,
#'   `S`, `lower_S`, `upper_S`, `h`, `lower_h`, `upper_h`.
#' @export
#'
#' @examples
#' \dontrun{
#' t_grid <- seq(0.01, 8, length.out = 200)
#' ci     <- CI(fit, t_grid, covariate = 0)
#' plot(ci$time, ci$h, type = "l")
#' lines(ci$time, ci$lower_h, lty = 2)
#' lines(ci$time, ci$upper_h, lty = 2)
#' }
CI <- function(fit, t, covariate, alpha = 0.05) {
  X_mat    <- matrix(rep(covariate, length(t)), ncol = length(covariate),
                     byrow = TRUE)
  var_theta <- fit$var
  z         <- qnorm(1 - alpha / 2)

  # Cumulative hazard
  grad_H  <- post(fit, X_mat, t, "gradient_H")
  H       <- post(fit, X_mat, t, "H")
  logH    <- log(H)
  var_H   <- rowSums((grad_H %*% var_theta) * grad_H)
  SE_logH <- sqrt(var_H) / H

  # Hazard
  grad_h  <- post(fit, X_mat, t, "gradient_h")
  hv      <- post(fit, X_mat, t, "h")
  var_h   <- rowSums((grad_h %*% var_theta) * grad_h)
  SE_logh <- sqrt(var_h) / hv

  data.frame(
    time    = t,
    H       = H,
    lower_H = exp(logH - z * SE_logH),
    upper_H = exp(logH + z * SE_logH),
    S       = exp(-H),
    lower_S = exp(-exp(logH + z * SE_logH)),
    upper_S = exp(-exp(logH - z * SE_logH)),
    h       = hv,
    lower_h = exp(log(hv) - z * SE_logh),
    upper_h = exp(log(hv) + z * SE_logh)
  )
}
