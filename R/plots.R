#' Plot estimated hazard function
#'
#' Quick base-graphics plot of the estimated hazard function from a fitted GH
#' model at a given covariate pattern.
#'
#' @param fit A fitted model object from [fit_genhaz()].
#' @param covariate Numeric vector of covariate values.
#' @param time Numeric vector of time points at which to evaluate the hazard.
#' @param xlim Optional x-axis limits passed to [plot()].
#' @param ylim Optional y-axis limits passed to [plot()].
#' @param ... Additional graphical parameters passed to [plot()].
#'
#' @return Invisibly returns the estimated hazard values.
#' @export
#'
#' @examples
#' \dontrun{
#' t_grid <- seq(0.01, 8, length.out = 300)
#' plot_hazard(fit, covariate = 0, time = t_grid)
#' }
plot_hazard <- function(fit, covariate, time, xlim = NULL, ylim = NULL, ...) {
  X_mat <- matrix(rep(covariate, length(time)), ncol = length(covariate),
                  byrow = TRUE)
  haz   <- post(fit, X_mat, time, "h")
  plot(time, haz, type = "l",
       xlab = "Time", ylab = "Hazard",
       main = "Estimated Hazard Function",
       xlim = xlim, ylim = ylim, ...)
  invisible(haz)
}
