pkg <- "C:/Users/Aaron/Documents/02_KI/GH/03_NextGen/genhaz"

# Step 0: Remove compiled artifacts from src/ left by roxygenize debug builds
src_dir <- file.path(pkg, "src")
to_remove <- list.files(src_dir, pattern = "\\.(o|dll|so)$", full.names = TRUE)
if (file.exists(file.path(src_dir, "symbols.rds")))
  to_remove <- c(to_remove, file.path(src_dir, "symbols.rds"))
if (file.exists(file.path(src_dir, "tmp.def")))
  to_remove <- c(to_remove, file.path(src_dir, "tmp.def"))
file.remove(to_remove)

# Step 1: Compile Rcpp attributes (regenerates RcppExports.R / .cpp)
Rcpp::compileAttributes(pkgdir = pkg)
cat("--- compileAttributes done ---\n")

# Step 2: Build
setwd(dirname(pkg))
system(paste0('"', R.home("bin"), '/R" CMD build "', basename(pkg), '"'))
cat("--- build done ---\n")

# Step 3: R CMD check (skip network CRAN check to avoid 10 min timeout)
tgz_file <- rev(sort(list.files(dirname(pkg), pattern = "genhaz_.*.tar.gz",
                                full.names = TRUE)))[1]
cat("Checking:", tgz_file, "\n")
system(paste0('"', R.home("bin"),
              '/R" CMD check --no-manual --no-vignettes --no-multiarch "',
              tgz_file, '"'))

# Step 4: Print full log
log_file <- file.path(dirname(pkg), "genhaz.Rcheck", "00check.log")
cat("\n====== CHECK LOG ======\n")
cat(readLines(log_file), sep = "\n")

install_log <- file.path(dirname(pkg), "genhaz.Rcheck", "00install.out")
if (file.exists(install_log)) {
  cat("\n====== INSTALL LOG ======\n")
  cat(readLines(install_log), sep = "\n")
}
