# Run this script once to produce data/fit_melanoma.rda.
# Avoids long fitting times in the vignette generation process.
# Requires: biostat3, usethis

library(genhaz)
library(survival)
library(biostat3)

set.seed(303567)

mel        <- biostat3::melanoma
mel$X      <- ifelse(mel$stage == "Localised", 0, 1)
mel$event  <- ifelse(mel$status == "Dead: cancer", 1, 0)
mel$time   <- mel$surv_mm
mel$period <- ifelse(mel$year8594 == "Diagnosed 75-84", 0, 1)

fit_melanoma <- fit_genhaz(
  Surv(mel$time, mel$event), ~ X + period + agegrp + sex,
  data       = mel,
  model_type = "GH",
  profile    = TRUE,
  n_knots    = 8,
  tol_LCV    = 0.001,
  lcv_method = "optimize",
  timeIt     = TRUE
)

pkg_root <- normalizePath(file.path(dirname(
  sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))
), ".."))
dir.create(file.path(pkg_root, "data"), showWarnings = FALSE)
save(fit_melanoma, file = file.path(pkg_root, "data", "fit_melanoma.rda"))
message("Saved: ", file.path(pkg_root, "data", "fit_melanoma.rda"))
