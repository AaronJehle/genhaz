#' Compute knot positions from event times
#'
#' Places `n_knots` equally spaced (in probability) quantiles of the
#' log-event-times. Boundary knots are placed at the 5 % and, optionally,
#' the 95 % or 100 % quantile of the event-time distribution.
#'
#' @param time Numeric vector of observed times.
#' @param event Integer vector of event indicators (1 = event, 0 = censored),
#'   same length as `time`.
#' @param n_knots Total number of knots, including the two boundary knots.
#' @param limit_to_95 Logical. If `TRUE` (default) the upper boundary knot
#'   is placed at the 95 % quantile; if `FALSE` it is placed at the maximum
#'   observed event time.
#'
#' @return Named numeric vector of length `n_knots` with knot positions on the
#'   *log-time* scale, as required by [genhaz_work()] and [fit_genhaz()].
#' @export
#'
#' @examples
#' set.seed(1)
#' t  <- rexp(200)
#' ev <- rbinom(200, 1, 0.8)
#' knot_pattern(t, ev, n_knots = 6)
knot_pattern <- function(time, event, n_knots, limit_to_95 = TRUE) {
  upper_prob <- if (limit_to_95) 0.95 else 1.0
  probs      <- seq(0.05, upper_prob, length.out = n_knots)
  log_event_times <- log(time[event == 1])
  quantile(log_event_times, probs = probs)
}
