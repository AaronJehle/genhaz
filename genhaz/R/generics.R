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

#' Predict from a fitted GH model
#'
#' Evaluates the fitted model at one or more covariate patterns over a time
#' grid, returning pointwise estimates with delta-method confidence intervals.
#' Delegates to [CI()] for all computations.
#'
#' @param object A `"genhaz_fit"` object from [fit_genhaz()].
#' @param newdata A `data.frame` of covariate values. Each row defines one
#'   covariate pattern; row names are used as group labels in the output.
#'   Must contain all variables referenced in the model formula.
#' @param times Numeric vector of evaluation time points.
#' @param type Quantity to predict: `"hazard"` (default), `"survival"`, or
#'   `"cumhaz"` (cumulative hazard).
#' @param alpha Significance level for confidence intervals. Default 0.05.
#' @param ... Currently unused.
#'
#' @return A `data.frame` with columns `pattern`, `time`, `estimate`,
#'   `lower`, `upper`. Rows are grouped by covariate pattern (one block
#'   per row of `newdata`).
#' @export
#'
#' @examples
#' \dontrun{
#' t_grid <- seq(0.01, 8, length.out = 200)
#'
#' # Single pattern
#' pred0 <- predict(fit, newdata = data.frame(X = 0), times = t_grid)
#'
#' # Two named groups
#' nd <- data.frame(X = c(0, 1))
#' rownames(nd) <- c("Unexposed", "Exposed")
#' pred <- predict(fit, newdata = nd, times = t_grid, type = "survival")
#' }
predict.genhaz_fit <- function(object, newdata, times,
                               type  = c("hazard", "survival", "cumhaz"),
                               alpha = 0.05, ...) {
  type <- match.arg(type)

  # Build covariate matrix using the stored formula (drop intercept column)
  X_mat    <- model.matrix(object$formula, newdata)[, -1L, drop = FALSE]
  patterns <- rownames(newdata)

  col_est <- switch(type, hazard = "h",       survival = "S",       cumhaz = "H")
  col_lo  <- switch(type, hazard = "lower_h", survival = "lower_S", cumhaz = "lower_H")
  col_hi  <- switch(type, hazard = "upper_h", survival = "upper_S", cumhaz = "upper_H")

  result_list <- lapply(seq_len(nrow(X_mat)), function(i) {
    ci_df <- CI(object, times, covariate = X_mat[i, ], alpha = alpha)
    data.frame(
      pattern  = patterns[i],
      time     = times,
      estimate = ci_df[[col_est]],
      lower    = ci_df[[col_lo]],
      upper    = ci_df[[col_hi]],
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, result_list)
}


# ---- plot -------------------------------------------------------------------

#' Plot estimated curves from a fitted GH model
#'
#' Plots hazard, survival, or cumulative hazard curves for one or more
#' covariate patterns. Each row of `newdata` becomes a separate coloured
#' line; row names serve as legend labels. Confidence bands are drawn as
#' dashed lines of the same colour.
#'
#' @param x A `"genhaz_fit"` object from [fit_genhaz()].
#' @param newdata A `data.frame` of covariate values. One row per group.
#'   Set row names for informative legend labels.
#' @param times Numeric vector of evaluation time points.
#' @param type Quantity to plot: `"hazard"` (default), `"survival"`, or
#'   `"cumhaz"`.
#' @param alpha Significance level for confidence bands. Default 0.05.
#' @param col Colour vector (recycled over groups). Defaults to `2, 3, 4, ...`
#'   (R's standard colour sequence, skipping black).
#' @param lty Line type for point-estimate curves. Default `1`.
#' @param xlab X-axis label. Default `"Time"`.
#' @param ylab Y-axis label. Derived from `type` when `NULL`.
#' @param main Plot title. Derived from `type` when `NULL`.
#' @param legend Logical; draw a legend when multiple groups are present?
#'   Default `TRUE`.
#' @param ci Logical; overlay dashed confidence band lines? Default `TRUE`.
#' @param ... Additional graphical parameters passed to [graphics::plot()].
#'
#' @return Invisibly returns the `data.frame` produced by
#'   `predict.genhaz_fit()`.
#' @importFrom graphics lines legend
#' @export
#'
#' @examples
#' \dontrun{
#' t_grid <- seq(0.01, 8, length.out = 300)
#'
#' # Hazard (default) for two groups
#' nd <- data.frame(X = c(0, 1))
#' rownames(nd) <- c("X = 0", "X = 1")
#' plot(fit, newdata = nd, times = t_grid)
#'
#' # Survival for a single group
#' plot(fit, newdata = data.frame(X = 0), times = t_grid, type = "survival")
#' }
plot.genhaz_fit <- function(x, newdata, times,
                            type   = c("hazard", "survival", "cumhaz"),
                            alpha  = 0.05,
                            col    = NULL,
                            lty    = 1,
                            xlab   = "Time",
                            ylab   = NULL,
                            main   = NULL,
                            legend = TRUE,
                            ci     = TRUE,
                            ...) {
  type <- match.arg(type)

  pred     <- predict(x, newdata = newdata, times = times,
                      type = type, alpha = alpha)
  patterns <- unique(pred$pattern)
  n_pat    <- length(patterns)

  if (is.null(col))  col  <- seq_len(n_pat) + 1L
  if (is.null(ylab)) ylab <- switch(type,
    hazard   = "Hazard h(t)",
    survival = "Survival S(t)",
    cumhaz   = "Cumulative hazard H(t)"
  )
  if (is.null(main)) main <- switch(type,
    hazard   = "Estimated hazard",
    survival = "Estimated survival",
    cumhaz   = "Estimated cumulative hazard"
  )

  # Y range spans all CI bands so no line is clipped
  ylim <- range(c(pred$lower, pred$upper), na.rm = TRUE)

  # First group initialises the plot device
  d1 <- pred[pred$pattern == patterns[1L], ]
  plot(d1$time, d1$estimate, type = "l",
       col = col[1L], lty = lty,
       xlab = xlab, ylab = ylab, main = main, ylim = ylim, ...)
  if (ci) {
    lines(d1$time, d1$lower, col = col[1L], lty = 2L)
    lines(d1$time, d1$upper, col = col[1L], lty = 2L)
  }

  # Additional groups added with lines()
  for (i in seq_along(patterns)[-1L]) {
    di <- pred[pred$pattern == patterns[i], ]
    lines(di$time, di$estimate, col = col[i], lty = lty)
    if (ci) {
      lines(di$time, di$lower, col = col[i], lty = 2L)
      lines(di$time, di$upper, col = col[i], lty = 2L)
    }
  }

  if (legend && n_pat > 1L) {
    legend("topright",
           legend = as.character(patterns),
           col    = col[seq_len(n_pat)],
           lty    = lty)
  }

  invisible(pred)
}
