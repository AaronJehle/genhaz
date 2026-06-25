# S3 generics for class "genhaz_fit"
#
# Provides print, summary, predict, and plot methods. All user-facing
# functions delegate to post() / CI() so prediction logic stays in one place.


# ---- print ------------------------------------------------------------------

#' Print a fitted GH model
#'
#' Displays a compact summary: model metadata and a covariate coefficient
#' table (spline parameters are suppressed as they are not directly
#' interpretable).
#'
#' @param x A `"genhaz_fit"` object from [fit_genhaz()].
#' @param digits Number of significant digits for numeric output. Default 4.
#' @param alpha Significance level for Wald confidence intervals. Default 0.05.
#' @param ... Currently unused.
#'
#' @return Invisibly returns `x`.
#' @export
#'
#' @examples
#' \dontrun{
#' fit <- fit_genhaz(Surv(dat$time, dat$event), ~X, data = dat,
#'                   model_type = "GH", profile = TRUE, n_knots = 6)
#' print(fit)
#' }
print.genhaz_fit <- function(x, digits = 4, alpha = 0.05, ...) {
  cat("Fitted generalized hazard model (genhaz)\n\n")

  mt <- if (length(unique(x$model_type)) == 1L)
    unique(x$model_type)
  else
    paste(x$model_type, collapse = " / ")

  cat(sprintf("  Model type : %s\n",   mt))
  cat(sprintf("  Knots      : %d (log-time scale)\n", length(x$knots)))
  cat(sprintf("  Lambda     : %s\n",   format(x$lambda, digits = digits)))
  cat(sprintf("  EDF        : %.2f\n", x$edf))
  cat(sprintf("  AIC        : %.2f\n\n", x$AIC))

  cov_idx <- grep("^beta", names(x$par))
  if (length(cov_idx) > 0L) {
    z_crit <- qnorm(1 - alpha / 2)
    est    <- x$par[cov_idx]
    se     <- x$se[cov_idx]
    ci_label <- paste0(round((1 - alpha) * 100), "%")
    tab <- data.frame(
      Estimate              = round(est,                    digits),
      Std.Err               = round(se,                    digits),
      z                     = round(x$z[cov_idx],          digits),
      p.value               = format.pval(x$p_values[cov_idx], digits = 3),
      lower                 = round(est - z_crit * se,     digits),
      upper                 = round(est + z_crit * se,     digits)
    )
    colnames(tab)[5:6] <- paste0(c("lower.", "upper."), ci_label)
    cat(sprintf("Covariate coefficients (Wald %s CI):\n", ci_label))
    print(tab)
  }

  invisible(x)
}


# ---- summary ----------------------------------------------------------------

#' Summarise a fitted GH model
#'
#' Returns a `"summary.genhaz_fit"` object containing model metadata and a
#' covariate coefficient table with Wald confidence intervals and
#' exponentiated estimates.
#'
#' @param object A `"genhaz_fit"` object from [fit_genhaz()].
#' @param alpha Significance level for confidence intervals. Default 0.05.
#' @param ... Currently unused.
#'
#' @return An object of class `"summary.genhaz_fit"` (printed by
#'   `print.summary.genhaz_fit()`).
#' @export
#'
#' @examples
#' \dontrun{
#' summary(fit)
#' }
summary.genhaz_fit <- function(object, alpha = 0.05, ...) {
  cov_idx <- grep("^beta", names(object$par))
  z_crit  <- qnorm(1 - alpha / 2)

  est   <- object$par[cov_idx]
  se    <- object$se[cov_idx]
  lower <- est - z_crit * se
  upper <- est + z_crit * se

  coef_tab <- data.frame(
    Estimate  = est,
    Std.Err   = se,
    z         = object$z[cov_idx],
    p.value   = object$p_values[cov_idx],
    lower     = lower,
    upper     = upper,
    exp.Est   = exp(est),
    exp.lower = exp(lower),
    exp.upper = exp(upper)
  )

  structure(
    list(
      formula    = object$formula,
      model_type = object$model_type,
      n_knots    = length(object$knots),
      lambda     = object$lambda,
      edf        = object$edf,
      AIC        = object$AIC,
      coef_tab   = coef_tab,
      alpha      = alpha
    ),
    class = "summary.genhaz_fit"
  )
}


#' Print a summary of a fitted GH model
#'
#' @param x A `"summary.genhaz_fit"` object from `summary.genhaz_fit()`.
#' @param digits Number of significant digits. Default 4.
#' @param ... Currently unused.
#'
#' @return Invisibly returns `x`.
#' @export
#'
#' @examples
#' \dontrun{
#' print(summary(fit))
#' }
print.summary.genhaz_fit <- function(x, digits = 4, ...) {
  cat("Fitted generalized hazard model (genhaz)\n\n")

  mt <- if (length(unique(x$model_type)) == 1L)
    unique(x$model_type)
  else
    paste(x$model_type, collapse = " / ")

  cat(sprintf("  Formula    : %s\n",   deparse(x$formula)))
  cat(sprintf("  Model type : %s\n",   mt))
  cat(sprintf("  Knots      : %d (log-time scale)\n", x$n_knots))
  cat(sprintf("  Lambda     : %s\n",   format(x$lambda, digits = digits)))
  cat(sprintf("  EDF        : %.2f\n", x$edf))
  cat(sprintf("  AIC        : %.2f\n\n", x$AIC))

  ci_label <- paste0(round((1 - x$alpha) * 100), "% CI")
  cat(sprintf("Covariate coefficients (Wald %s):\n", ci_label))

  tab <- x$coef_tab
  tab$p.value <- format.pval(tab$p.value, digits = 3)
  colnames(tab) <- c(
    "Estimate", "Std.Err", "z", "p.value",
    "lower", "upper",
    "exp(Est.)", "exp(lower)", "exp(upper)"
  )
  print(tab, digits = digits)

  invisible(x)
}


# ---- predict ----------------------------------------------------------------

as_genhaz_pred <- function(df, type) {
  class(df)             <- c("genhaz_pred", "data.frame")
  attr(df, "pred_type") <- type
  df
}

#' Predict from a fitted GH model
#'
#' Evaluates the fitted model at one or more covariate patterns over a time
#' grid, returning pointwise estimates with delta-method confidence intervals.
#' Only the requested quantity is computed.
#'
#' @param object A `"genhaz_fit"` object from [fit_genhaz()].
#' @param newdata A `data.frame` of covariate values. Each row defines one
#'   covariate pattern; row names are used as group labels in the output.
#'   Must contain all variables referenced in the model formula. For
#'   `"surv_diff"` and `"rmst_diff"`, exactly two rows are required.
#' @param times Numeric vector of time points. For `"rmst"` and `"rmst_diff"`,
#'   these are restriction times tau at which RMST(tau) = integral_0^tau S(t) dt
#'   is evaluated.
#' @param type Quantity to predict: `"hazard"` (default), `"survival"`,
#'   `"cumhaz"` (cumulative hazard), `"rmst"` (restricted mean survival time),
#'   `"surv_diff"` (survival difference S1−S2), `"rmst_diff"` (RMST
#'   difference), `"hazard_ratio"` (h1/h2), or `"time_ratio"` (tau/t where
#'   S2(tau) = S1(t)). The last four require exactly two rows in `newdata`.
#' @param alpha Significance level for confidence intervals. Default 0.05.
#' @param ... Currently unused.
#'
#' @return A `data.frame` of class `c("genhaz_pred", "data.frame")` with
#'   columns `pattern`, `time`, `estimate`, `lower`, `upper`. A `pred_type`
#'   attribute records which quantity was computed; `plot()` uses this to set
#'   axis labels automatically. For per-group types rows are grouped by
#'   covariate pattern; for difference types a single block labelled
#'   `"group1 - group2"` is returned.
#' @export
#'
#' @examples
#' \dontrun{
#' t_grid <- seq(0.01, 8, length.out = 200)
#' nd <- data.frame(X = c(0, 1))
#' rownames(nd) <- c("Unexposed", "Exposed")
#'
#' predict(fit, newdata = nd, times = t_grid, type = "survival")
#' predict(fit, newdata = nd, times = c(1, 2, 5), type = "rmst")
#' predict(fit, newdata = nd, times = t_grid, type = "surv_diff")
#' predict(fit, newdata = nd, times = c(1, 2, 5), type = "rmst_diff")
#' predict(fit, newdata = nd, times = t_grid, type = "hazard_ratio")
#' predict(fit, newdata = nd, times = t_grid, type = "time_ratio")
#' }
predict.genhaz_fit <- function(object, newdata, times,
                               type  = c("hazard", "survival", "cumhaz",
                                         "rmst", "surv_diff", "rmst_diff",
                                         "hazard_ratio", "time_ratio"),
                               alpha = 0.05, ...) {
  fit       <- object
  type      <- match.arg(type)
  X_mat     <- model.matrix(fit$formula, newdata)[, -1L, drop = FALSE]
  patterns  <- rownames(newdata)
  var_theta <- fit$var
  z         <- qnorm(1 - alpha / 2)

  # ---- two-group comparison types ------------------------------------------
  if (type %in% c("surv_diff", "rmst_diff", "hazard_ratio", "time_ratio")) {
    if (nrow(X_mat) != 2L)
      stop("'newdata' must have exactly 2 rows for type '", type, "'.")
    pat_label <- paste0(patterns[1L], " - ", patterns[2L])

    if (type == "surv_diff") {
      xi1 <- matrix(rep(X_mat[1L, ], length(times)), ncol = ncol(X_mat), byrow = TRUE)
      xi2 <- matrix(rep(X_mat[2L, ], length(times)), ncol = ncol(X_mat), byrow = TRUE)

      S1        <- exp(-post(fit, xi1, times, "H"))
      S2        <- exp(-post(fit, xi2, times, "H"))
      grad_H1   <- post(fit, xi1, times, "gradient_H")
      grad_H2   <- post(fit, xi2, times, "gradient_H")
      diff_val  <- S1 - S2
      # grad of (S1 - S2) w.r.t. theta: -S1*grad_H1 + S2*grad_H2 (row-wise)
      grad_diff <- -S1 * grad_H1 + S2 * grad_H2
      var_diff  <- rowSums((grad_diff %*% var_theta) * grad_diff)

      return(as_genhaz_pred(data.frame(
                        pattern  = pat_label, time = times,
                        estimate = diff_val,
                        lower    = diff_val - z * sqrt(var_diff),
                        upper    = diff_val + z * sqrt(var_diff),
                        stringsAsFactors = FALSE), type))

    } else if (type == "rmst_diff") {
      gl  <- gauss_legendre(25L)
      out <- t(vapply(times, function(tau) {
        t_q   <- tau * (gl$x + 1) / 2
        w_q   <- tau / 2 * gl$w
        xi1_q <- matrix(rep(X_mat[1L, ], 25L), ncol = ncol(X_mat), byrow = TRUE)
        xi2_q <- matrix(rep(X_mat[2L, ], 25L), ncol = ncol(X_mat), byrow = TRUE)

        S1_q      <- exp(-post(fit, xi1_q, t_q, "H"))
        S2_q      <- exp(-post(fit, xi2_q, t_q, "H"))
        grad_H1_q <- post(fit, xi1_q, t_q, "gradient_H")
        grad_H2_q <- post(fit, xi2_q, t_q, "gradient_H")

        diff_val  <- sum(w_q * S1_q) - sum(w_q * S2_q)
        grad_diff <- colSums(-w_q * S1_q * grad_H1_q) -
                     colSums(-w_q * S2_q * grad_H2_q)
        var_diff  <- as.numeric(crossprod(grad_diff, var_theta %*% grad_diff))
        c(diff_val, diff_val - z * sqrt(var_diff), diff_val + z * sqrt(var_diff))
      }, numeric(3)))

      return(as_genhaz_pred(data.frame(
                        pattern  = pat_label, time = times,
                        estimate = out[, 1L],
                        lower    = out[, 2L],
                        upper    = out[, 3L],
                        stringsAsFactors = FALSE), type))

    } else if (type == "hazard_ratio") {
      xi1 <- matrix(rep(X_mat[1L, ], length(times)), ncol = ncol(X_mat), byrow = TRUE)
      xi2 <- matrix(rep(X_mat[2L, ], length(times)), ncol = ncol(X_mat), byrow = TRUE)

      h1      <- post(fit, xi1, times, "h")
      h2      <- post(fit, xi2, times, "h")
      grad_h1 <- post(fit, xi1, times, "gradient_h")
      grad_h2 <- post(fit, xi2, times, "gradient_h")

      HR         <- h1 / h2
      grad_logHR <- grad_h1 / h1 - grad_h2 / h2   # row-wise: each row scaled by scalar
      var_logHR  <- rowSums((grad_logHR %*% var_theta) * grad_logHR)

      return(as_genhaz_pred(data.frame(
                        pattern  = pat_label, time = times,
                        estimate = HR,
                        lower    = exp(log(HR) - z * sqrt(var_logHR)),
                        upper    = exp(log(HR) + z * sqrt(var_logHR)),
                        stringsAsFactors = FALSE), type))

    } else {
      # time_ratio: TR(t) = tau/t where H2(tau) = H1(t), i.e. S2(tau) = S1(t).
      # Delta method via implicit function theorem:
      #   grad(log TR) = (grad_H1(t) - grad_H2(tau)) / (tau * h2(tau))
      xi1     <- matrix(rep(X_mat[1L, ], length(times)), ncol = ncol(X_mat), byrow = TRUE)
      H1      <- post(fit, xi1, times, "H")
      grad_H1 <- post(fit, xi1, times, "gradient_H")
      t_upper <- max(times) * 20

      out <- t(vapply(seq_along(times), function(j) {
        target <- H1[j]
        tau_j  <- tryCatch(
          uniroot(function(tau_try) {
            post(fit, matrix(X_mat[2L, ], nrow = 1L), tau_try, "H") - target
          }, lower = 1e-8, upper = t_upper, tol = 1e-8)$root,
          error = function(e) NA_real_
        )
        if (is.na(tau_j)) return(c(NA_real_, NA_real_, NA_real_))

        xi2_tau     <- matrix(X_mat[2L, ], nrow = 1L)
        h2_tau      <- post(fit, xi2_tau, tau_j, "h")
        grad_H2_tau <- post(fit, xi2_tau, tau_j, "gradient_H")

        grad_logTR <- (grad_H1[j, ] - grad_H2_tau[1L, ]) / (tau_j * h2_tau)
        var_logTR  <- as.numeric(crossprod(grad_logTR, var_theta %*% grad_logTR))
        TR_j       <- tau_j / times[j]
        c(TR_j, exp(log(TR_j) - z * sqrt(var_logTR)),
                exp(log(TR_j) + z * sqrt(var_logTR)))
      }, numeric(3)))

      return(as_genhaz_pred(data.frame(
                        pattern  = pat_label, time = times,
                        estimate = out[, 1L],
                        lower    = out[, 2L],
                        upper    = out[, 3L],
                        stringsAsFactors = FALSE), type))
    }
  }

  # ---- per-pattern types ---------------------------------------------------
  result_list <- lapply(seq_len(nrow(X_mat)), function(i) {
    xi <- matrix(rep(X_mat[i, ], length(times)),
                 ncol = ncol(X_mat), byrow = TRUE)

    if (type == "hazard") {
      hv      <- post(fit, xi, times, "h")
      grad_h  <- post(fit, xi, times, "gradient_h")
      var_h   <- rowSums((grad_h %*% var_theta) * grad_h)
      SE_logh <- sqrt(var_h) / hv
      data.frame(pattern  = patterns[i], time = times,
                 estimate = hv,
                 lower    = exp(log(hv) - z * SE_logh),
                 upper    = exp(log(hv) + z * SE_logh),
                 stringsAsFactors = FALSE)

    } else if (type %in% c("survival", "cumhaz")) {
      H       <- post(fit, xi, times, "H")
      grad_H  <- post(fit, xi, times, "gradient_H")
      logH    <- log(H)
      var_H   <- rowSums((grad_H %*% var_theta) * grad_H)
      SE_logH <- sqrt(var_H) / H

      if (type == "cumhaz") {
        data.frame(pattern  = patterns[i], time = times,
                   estimate = H,
                   lower    = exp(logH - z * SE_logH),
                   upper    = exp(logH + z * SE_logH),
                   stringsAsFactors = FALSE)
      } else {
        data.frame(pattern  = patterns[i], time = times,
                   estimate = exp(-H),
                   lower    = exp(-exp(logH + z * SE_logH)),
                   upper    = exp(-exp(logH - z * SE_logH)),
                   stringsAsFactors = FALSE)
      }

    } else {
      # RMST(tau) = integral_0^tau S(t) dt via Gauss-Legendre quadrature.
      # grad_RMST = -integral_0^tau S(t) * grad_H(t) dt  (delta method)
      gl <- gauss_legendre(25L)

      out <- t(vapply(times, function(tau) {
        t_q  <- tau * (gl$x + 1) / 2
        w_q  <- tau / 2 * gl$w

        xi_q     <- matrix(rep(X_mat[i, ], 25L), ncol = ncol(X_mat), byrow = TRUE)
        S_q      <- exp(-post(fit, xi_q, t_q, "H"))
        grad_H_q <- post(fit, xi_q, t_q, "gradient_H")

        rmst_val  <- sum(w_q * S_q)
        grad_rmst <- colSums(-w_q * S_q * grad_H_q)
        var_rmst  <- as.numeric(crossprod(grad_rmst, var_theta %*% grad_rmst))

        c(rmst_val, rmst_val - z * sqrt(var_rmst), rmst_val + z * sqrt(var_rmst))
      }, numeric(3)))

      data.frame(pattern  = patterns[i], time = times,
                 estimate = out[, 1L],
                 lower    = out[, 2L],
                 upper    = out[, 3L],
                 stringsAsFactors = FALSE)
    }
  })

  as_genhaz_pred(do.call(rbind, result_list), type)
}


# ---- plot -------------------------------------------------------------------

#' Plot a genhaz prediction
#'
#' Plots the output of [predict.genhaz_fit()] — any of the eight prediction
#' types. Each distinct covariate pattern becomes a coloured line; confidence
#' bands are drawn as dashed lines of the same colour. `ylim` is always
#' derived from the full range of CI bands plus the point estimate, so nothing
#' is clipped. Axis labels and title are set automatically from the prediction
#' type and can be overridden via `ylab` / `main`.
#'
#' @param x A `"genhaz_pred"` object returned by [predict.genhaz_fit()].
#' @param col Colour vector (recycled over groups). Defaults to `2, 3, 4, ...`
#'   (R's standard colour sequence, skipping black).
#' @param lty Line type for point-estimate curves. Default `1`.
#' @param xlab X-axis label. Default `"Time"`.
#' @param ylab Y-axis label. Derived from the prediction type when `NULL`.
#' @param main Plot title. Derived from the prediction type when `NULL`.
#' @param legend Logical; draw a legend when multiple groups are present?
#'   Default `TRUE`.
#' @param ci Logical; overlay dashed confidence band lines? Default `TRUE`.
#' @param ... Additional graphical parameters passed to [graphics::plot()].
#'
#' @return Invisibly returns `x`.
#' @importFrom graphics lines legend
#' @export
#'
#' @examples
#' \dontrun{
#' nd <- data.frame(X = c(0, 1))
#' rownames(nd) <- c("X = 0", "X = 1")
#' t_grid <- seq(0.01, 8, length.out = 200)
#'
#' pred <- predict(fit, newdata = nd, times = t_grid, type = "hazard_ratio")
#' plot(pred, col = "purple")
#' abline(h = 1, lty = 2)
#' }
plot.genhaz_pred <- function(x, col = NULL, lty = 1,
                              xlab = "Time", ylab = NULL, main = NULL,
                              legend = TRUE, ci = TRUE, ...) {
  pred_type <- attr(x, "pred_type")
  patterns  <- unique(x$pattern)
  n_pat     <- length(patterns)

  if (is.null(col)) col <- seq_len(n_pat) + 1L

  if (is.null(ylab)) ylab <- switch(pred_type,
    hazard       = "Hazard h(t)",
    survival     = "Survival S(t)",
    cumhaz       = "Cumulative hazard H(t)",
    rmst         = "RMST",
    surv_diff    = "Survival difference",
    rmst_diff    = "RMST difference",
    hazard_ratio = "Hazard ratio",
    time_ratio   = "Time ratio"
  )
  if (is.null(main)) main <- switch(pred_type,
    hazard       = "Estimated hazard",
    survival     = "Estimated survival",
    cumhaz       = "Estimated cumulative hazard",
    rmst         = "Restricted mean survival time",
    surv_diff    = "Survival difference",
    rmst_diff    = "RMST difference",
    hazard_ratio = "Hazard ratio h1(t) / h2(t)",
    time_ratio   = "Time ratio"
  )

  # Y range covers estimate + both CI bands — nothing is clipped
  ylim <- range(c(x$lower, x$upper, x$estimate), na.rm = TRUE)

  d1 <- x[x$pattern == patterns[1L], ]
  plot(d1$time, d1$estimate, type = "l",
       col = col[1L], lty = lty,
       xlab = xlab, ylab = ylab, main = main, ylim = ylim, ...)
  if (ci) {
    lines(d1$time, d1$lower, col = col[1L], lty = 2L)
    lines(d1$time, d1$upper, col = col[1L], lty = 2L)
  }

  for (i in seq_along(patterns)[-1L]) {
    di <- x[x$pattern == patterns[i], ]
    lines(di$time, di$estimate, col = col[i], lty = lty)
    if (ci) {
      lines(di$time, di$lower, col = col[i], lty = 2L)
      lines(di$time, di$upper, col = col[i], lty = 2L)
    }
  }

  if (legend && n_pat > 1L)
    legend("topright", legend = as.character(patterns),
           col = col[seq_len(n_pat)], lty = lty)

  invisible(x)
}


#' Plot estimated curves from a fitted GH model
#'
#' Calls [predict.genhaz_fit()] and passes the result to `plot.genhaz_pred()`.
#' Supports all eight prediction types; `ylim` is automatically derived from
#' the full range of CI bands so no line is clipped.
#'
#' @param x A `"genhaz_fit"` object from [fit_genhaz()].
#' @param newdata A `data.frame` of covariate values. One row per group.
#'   Row names become legend labels.
#' @param times Numeric vector of evaluation time points.
#' @param type Quantity to plot. One of `"hazard"` (default), `"survival"`,
#'   `"cumhaz"`, `"rmst"`, `"surv_diff"`, `"rmst_diff"`, `"hazard_ratio"`,
#'   or `"time_ratio"`. The last four require exactly two rows in `newdata`.
#' @param alpha Significance level for confidence bands. Default 0.05.
#' @param ... Additional arguments passed to `plot.genhaz_pred()` and then to
#'   [graphics::plot()] (e.g. `col`, `lty`, `xlab`, `main`, `xlim`).
#'
#' @return Invisibly returns the `"genhaz_pred"` object from
#'   [predict.genhaz_fit()].
#' @importFrom graphics lines legend
#' @export
#'
#' @examples
#' \dontrun{
#' nd <- data.frame(X = c(0, 1))
#' rownames(nd) <- c("X = 0", "X = 1")
#' t_grid <- seq(0.01, 8, length.out = 300)
#'
#' plot(fit, newdata = nd, times = t_grid)
#' plot(fit, newdata = nd, times = t_grid, type = "survival",
#'      col = c("steelblue", "firebrick"))
#' plot(fit, newdata = nd, times = t_grid, type = "hazard_ratio",
#'      col = "purple")
#' }
plot.genhaz_fit <- function(x, newdata, times,
                            type  = c("hazard", "survival", "cumhaz",
                                      "rmst", "surv_diff", "rmst_diff",
                                      "hazard_ratio", "time_ratio"),
                            alpha = 0.05, ...) {
  type <- match.arg(type)
  pred <- stats::predict(x, newdata = newdata, times = times,
                         type = type, alpha = alpha)
  plot(pred, ...)
  invisible(pred)
}
