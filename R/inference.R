#' Likelihood ratio tests between nested GH models
#'
#' Compares two or more nested `"genhaz_fit"` models with chi-squared
#' likelihood-ratio tests, in the style of [stats::anova()] for `lm`/`glm`
#' objects. Models are compared sequentially in the order supplied, so list the
#' most restricted model first.
#'
#' @param object A `"genhaz_fit"` object from [fit_genhaz()] (the first, usually
#'   most restricted, model).
#' @param ... One or more further `"genhaz_fit"` objects to compare against
#'   `object`. Non-`"genhaz_fit"` arguments are ignored.
#'
#' @return An analysis-of-deviance table of class `c("anova", "data.frame")`
#'   with one row per model and columns `Df` (model degrees of freedom),
#'   `LogLik`, `Chisq`, `Chi Df`, and `Pr(>Chisq)` (the comparison columns are
#'   `NA` for the first model).
#' @export
#'
#' @examples
#' \dontrun{
#' fit_ph <- fit_genhaz(Surv(dat$time, dat$event), ~X,
#'                      data = dat, model_type = "PH", n_knots = 6)
#' fit_gh <- fit_genhaz(Surv(dat$time, dat$event), ~X,
#'                      data = dat, model_type = "GH", n_knots = 6)
#' anova(fit_ph, fit_gh)
#' }
anova.genhaz_fit <- function(object, ...) {
  dots   <- list(...)
  models <- c(list(object), dots[vapply(dots, inherits, logical(1),
                                        "genhaz_fit")])
  if (length(models) < 2L) {
    stop("anova.genhaz_fit() needs at least two nested 'genhaz_fit' models ",
         "to compare; term-wise ANOVA for a single model is not available.")
  }

  k       <- length(models)
  df_mod  <- vapply(models, function(m) m$df,        numeric(1))
  loglik  <- vapply(models, function(m) -m$objective, numeric(1))
  chisq   <- c(NA_real_, vapply(2:k, function(i)
                 2 * (models[[i - 1L]]$objective - models[[i]]$objective),
                 numeric(1)))
  chi_df  <- c(NA_real_, vapply(2:k, function(i)
                 models[[i]]$df - models[[i - 1L]]$df, numeric(1)))
  pval    <- pchisq(chisq, chi_df, lower.tail = FALSE)

  tab <- data.frame(Df = df_mod, LogLik = loglik, Chisq = chisq,
                    `Chi Df` = chi_df, `Pr(>Chisq)` = pval,
                    check.names = FALSE)
  rownames(tab) <- paste("Model", seq_len(k))

  labels <- vapply(models, function(m) {
    mt <- if (length(unique(m$model_type)) == 1L) unique(m$model_type)
          else paste(m$model_type, collapse = "/")
    sprintf("%s, %s", deparse(m$formula), mt)
  }, character(1))
  heading <- c("Analysis of Deviance Table (likelihood ratio tests)\n",
               paste0(" Model ", seq_len(k), ": ", labels))

  structure(tab, heading = heading, class = c("anova", "data.frame"))
}

#' Wald confidence intervals for model parameters
#'
#' Computes Wald confidence intervals for the parameters of a fitted GH model,
#' following the [stats::confint()] convention. With `diff = TRUE` it instead
#' returns a delta-method interval for the difference of two parameters,
#' accounting for their covariance.
#'
#' @param object A `"genhaz_fit"` object from [fit_genhaz()].
#' @param parm Character vector of parameter names (must be in
#'   `object$parnames`). If missing, all parameters are used. When
#'   `diff = TRUE`, exactly two names must be supplied.
#' @param level Confidence level. Default `0.95`.
#' @param diff Logical. If `TRUE`, return a single interval for
#'   `parm[1] - parm[2]` via the delta method. Default `FALSE`.
#' @param ... Currently unused.
#'
#' @return A matrix with one row per parameter (or a single row for the
#'   difference) and two columns giving the lower and upper bounds.
#' @export
#'
#' @examples
#' \dontrun{
#' confint(fit)                                       # all parameters
#' confint(fit, "beta2_X")                            # one parameter
#' confint(fit, c("beta1_X", "beta2_X"), diff = TRUE) # difference
#' }
confint.genhaz_fit <- function(object, parm, level = 0.95, diff = FALSE, ...) {
  a   <- (1 - level) / 2
  z   <- qnorm(1 - a)
  pct <- paste0(format(100 * c(a, 1 - a), trim = TRUE,
                       scientific = FALSE, digits = 3), " %")

  if (diff) {
    if (length(parm) != 2L) {
      stop("'diff = TRUE' requires exactly two parameter names in 'parm'.")
    }
    unknown <- setdiff(parm, object$parnames)
    if (length(unknown)) {
      stop("Unknown parameter(s): ", paste(unknown, collapse = ", "), ".")
    }
    cov_m <- object$var[parm, parm]
    v     <- cov_m[1, 1] - 2 * cov_m[1, 2] + cov_m[2, 2]
    est   <- object$par[parm[1]] - object$par[parm[2]]
    ci    <- est + c(-1, 1) * z * sqrt(v)
    return(matrix(ci, nrow = 1L,
                  dimnames = list(paste(parm, collapse = " - "), pct)))
  }

  if (missing(parm)) parm <- object$parnames
  unknown <- setdiff(parm, object$parnames)
  if (length(unknown)) {
    stop("Unknown parameter(s): ", paste(unknown, collapse = ", "), ".")
  }
  est <- object$par[parm]
  se  <- object$se[parm]
  out <- cbind(est - z * se, est + z * se)
  dimnames(out) <- list(parm, pct)
  out
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
