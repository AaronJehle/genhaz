#' Fit a penalized general hazard model
#'
#' High-level interface for fitting a penalised cubic regression spline
#' general hazard (GH) model. Wraps [genhaz_work()] with automatic knot
#' placement, design matrix construction, and optional profiling of the
#' smoothing parameter `lambda` via the modified LCV criterion.
#'
#' @param surv A `Surv` object (from the **survival** package). Right-
#'   censored (`Surv(time, event)`), left-truncated and right-censored
#'   (`Surv(start, stop, event)`), and interval-censored
#'   (`Surv(t1, t2, type = "interval2")`) data are all supported.
#' @param formula A one-sided formula specifying the covariate terms, e.g.
#'   `~ x + z`. An intercept is added internally; do not include one here.
#' @param data A data frame containing all variables referenced in `formula`.
#' @param n_knots Total number of knots for the baseline log-hazard spline,
#'   including the two boundary knots. Default 8.
#' @param lambda Initial (or fixed) smoothing parameter. Ignored when
#'   `profile = TRUE`. Default 0.
#' @param model_type Character vector of length 1 or `ncol(X)` selecting the
#'   model class per covariate: `"GH"` (general hazard, default), `"PH"`
#'   (proportional hazards), `"AFT"` (accelerated failure time), or `"AH"`
#'   (additive hazards).
#' @param control List of optimisation controls for [nlminb()]: `trace`,
#'   `x.tol`, `rel.tol`, `max_iter`.
#' @param profile Logical. If `TRUE`, optimise `lambda` via LCV before
#'   returning the final fit. Default `FALSE`.
#' @param lcv_method Character string passed to [genhaz_work()] when
#'   `profile = TRUE`. Selects the strategy for optimising `log(lambda)`:
#'   `"full"` (default) uses the full third-derivative LCV gradient,
#'   `"approx"` uses the faster first-order approximation, and `"optimize"`
#'   minimises LCV directly without gradient information. Ignored when
#'   `profile = FALSE`.
#' @param lambda_surv Logical. If `TRUE`, initialise `lambda` from a
#'   **survPen** proportional hazard fit. Requires **survPen** to be
#'   installed. Default `FALSE`.
#' @param init Optional numeric vector of starting values. If `NULL`, zeros
#'   are used.
#' @param interval Length-2 numeric vector giving the search interval for
#'   `log(lambda)` during LCV optimisation. Default `c(0, 10)`.
#' @param nquad Number of Gauss-Legendre quadrature points. Default 25.
#' @param tol_LCV Convergence tolerance for LCV optimisation. Default `1e-3`.
#' @param knots Optional numeric vector of knot positions on the log-time
#'   scale. If `NULL` (default), knots are placed automatically using
#'   [knot_pattern()].
#' @param reKnot Non-negative integer. After the initial fit, re-position
#'   knots `abs(reKnot)` times using the AFT-adjusted times and refit.
#'   Default `1` (one re-knot pass). Set to `0` to skip.
#' @param limit_to_95 Logical passed to [knot_pattern()]. Default `TRUE`.
#' @param timeIt Logical. If `TRUE`, print the wall-clock fitting time.
#'   Default `FALSE`.
#'
#' @return A list of class `"genhaz_fit"` containing (among others):
#'   \describe{
#'     \item{`par`}{Named vector of parameter estimates.}
#'     \item{`se`}{Standard errors.}
#'     \item{`z`}{Wald z-statistics.}
#'     \item{`p_values`}{Two-sided p-values.}
#'     \item{`var`}{Estimated covariance matrix of `par`.}
#'     \item{`AIC`}{Akaike information criterion.}
#'     \item{`edf`}{Effective degrees of freedom.}
#'     \item{`lambda`}{Smoothing parameter used.}
#'     \item{`knots`}{Knot vector on the log-time scale.}
#'     \item{`Z`}{Spline projection matrix.}
#'   }
#' @export
#'
#' @examples
#' \dontrun{
#' library(survival)
#' set.seed(42)
#' dat <- sim_szenario(szenario = 1, beta1 = 0.5, beta2 = 0.5, n = 300)
#' fit <- fit_genhaz(Surv(dat$time, dat$event), ~ X, data = dat,
#'                   model_type = "GH", profile = TRUE, n_knots = 6,
#'                   tol_LCV = 0.01)
#' fit$par
#' fit$se
#' }
fit_genhaz <- function(surv, formula, data, n_knots = 8, lambda = 0,
                       model_type = "GH",
                       control = list(trace = 0, x.tol = 1e-12,
                                      rel.tol = 1e-12, max_iter = 100),
                       profile = FALSE, lcv_method = c("full", "approx", "optimize"),
                       lambda_surv = FALSE, init = NULL,
                       interval = c(0, 10), nquad = 25, tol_LCV = 1e-3,
                       knots = NULL, reKnot = 1, limit_to_95 = TRUE,
                       timeIt = FALSE) {

  lcv_method <- match.arg(lcv_method)

  if (timeIt) {
    start_time <- Sys.time()
    on.exit({
      end_time <- Sys.time()
      message("Execution time: ",
              round(difftime(end_time, start_time, units = "secs"), 2),
              " seconds")
    }, add = TRUE)
  }

  # --- Parse Surv object ---
  type_surv <- attr(surv, "type")
  if (type_surv == "counting") {
    cens_type <- "lt_rc"
    time  <- surv[, "stop"]    # event / censoring time
    time2 <- surv[, "start"]   # truncation time
    event <- surv[, "status"]
  } else if (type_surv == "interval") {
    cens_type <- "ic"
    time  <- surv[, "time1"]
    time2 <- surv[, "time2"]
    event <- surv[, "status"]
  } else if (type_surv == "right") {
    cens_type <- "rc"
    time  <- surv[, "time"]
    time2 <- time
    event <- surv[, "status"]
  } else {
    stop("Unsupported Surv type: '", type_surv, "'.")
  }

  if (length(time) != nrow(data)) {
    stop("Length of 'surv' and number of rows in 'data' do not match.")
  }
  data$time  <- time
  data$time2 <- time2
  data$event <- event

  # --- Knot placement ---
  if (is.null(knots)) {
    knots <- knot_pattern(time, event, n_knots, limit_to_95 = limit_to_95)
  }
  # Store in package env so survPen formulas can access them
  .pkg_env$knots_smf <- knots

  # --- Design matrix (covariates only, no intercept) ---
  X <- model.matrix(formula, data = data)[, -1, drop = FALSE]

  # --- Fit a survPen PH model to warm-start or borrow lambda (lambda_surv) ---
  ph_mod <- NULL
  if (lambda_surv) {
    if (!requireNamespace("survPen", quietly = TRUE)) {
      stop("Package 'survPen' is required when lambda_surv = TRUE.")
    }
    baseline <- ~ survPen::smf(time, knots = .pkg_env$knots_smf)
    sp_ph    <- as.formula(paste0("~", deparse(formula[[2]]), "+",
                                  deparse(baseline[[2]])))
    ph_mod   <- survPen::survPen(sp_ph, data = data, t1 = time, event = event)
    if (is.null(init)) {
      coef_ph <- ph_mod$coefficients
      n_spl   <- n_knots - 1
      n_cov   <- length(coef_ph) - n_spl - 1
      init    <- c(coef_ph[1], tail(coef_ph, -(n_cov + 1)),
                   rep(0, n_cov), coef_ph[2:(n_cov + 1)]) * 0
    }
  }

  # --- Fall-back zero initialisation ---
  if (is.null(init)) {
    init <- rep(0, 1 + n_knots - 1 + 2 * ncol(X))
  }

  # Spline projection matrix (compute once)
  Z <- attr(smf_cpp(log(time + 1e-10), knots = knots,
                    intercept = FALSE, derivs = 0L), "Z")

  # --- Fit ---
  if (lambda_surv) {
    fit <- genhaz_work(theta = init, time = time, X = X, knots = knots,
                       Z = Z, event = event, lambda = ph_mod$lambda,
                       cens_type = cens_type, model_type = model_type,
                       time2 = time2, res = "fit", control = control,
                       interval = interval, nquad = nquad,
                       tol_LCV = tol_LCV)
    fit$roh_opt <- log(ph_mod$lambda)
  } else {
    fit <- genhaz_work(theta = init, time = time, X = X, knots = knots,
                       Z = Z, event = event, lambda = lambda,
                       cens_type = cens_type, model_type = model_type,
                       time2 = time2,
                       res = ifelse(profile, "fit_LCV", "fit"),
                       lcv_method = lcv_method,
                       control = control, interval = interval,
                       nquad = nquad, tol_LCV = tol_LCV)
  }

  # --- Re-knot iterations ---
  riktning <- sign(reKnot)
  reKnot   <- abs(reKnot)
  knot_hist <- list()
  par_hist  <- list()

  for (i in seq_len(reKnot)) {
    knot_hist[[i]] <- fit$knots
    par_hist[[i]]  <- fit$par
    new_X     <- t(t(X) * fit$par[(fit$nb0 + 1):(fit$nb0 + fit$nb1)])
    new_time  <- time * exp(riktning * rowSums(new_X))
    new_knots <- knot_pattern(new_time, event, n_knots,
                              limit_to_95 = limit_to_95)
    fit <- genhaz_work(theta = fit$par, time = time, X = X,
                       knots = new_knots, Z = Z, event = event,
                       lambda = lambda, cens_type = cens_type,
                       model_type = model_type, time2 = time2,
                       res = ifelse(profile, "fit_LCV", "fit"),
                       lcv_method = lcv_method,
                       control = control, interval = interval,
                       nquad = nquad, tol_LCV = tol_LCV)
  }

  if (reKnot > 0) {
    fit$knot_hist <- knot_hist
    fit$par_hist  <- par_hist
  }

  # --- Name parameters using covariate names ---
  covariate_names <- colnames(X)
  parnames <- c(
    "intercept",
    paste0("s", seq_len(n_knots - 1)),
    paste0("beta1_", covariate_names),
    paste0("beta2_", covariate_names)
  )
  names(fit$par)    <- parnames
  names(fit$se)     <- parnames
  fit$parnames      <- parnames
  if (!is.null(fit$var)) {
    rownames(fit$var) <- colnames(fit$var) <- parnames
  }
  if (!is.null(fit$hessian)) {
    rownames(fit$hessian) <- colnames(fit$hessian) <- parnames
  }

  fit$formula <- formula
  class(fit) <- c("genhaz_fit")
  fit
}
