pkg <- "C:/Users/Aaron/Documents/02_KI/GH/03_NextGen/genhaz"

# Generate man pages with roxygen2 (no devtools needed)
roxygen2::roxygenize(pkg)
cat("--- roxygenize done ---\n")
