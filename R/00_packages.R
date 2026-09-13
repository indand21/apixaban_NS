# 00_packages.R
# Package loading for Apixaban CKD QSP model

required_packages <- c(
  "deSolve",    # ODE integration (lsoda for stiff systems)
  "ggplot2",    # Visualization
  "dplyr",      # Data wrangling
  "tidyr",      # Data reshaping
  "purrr",      # Functional iteration over scenarios

  "pracma",     # Trapezoidal AUC (trapz)
  "patchwork",  # Multi-panel figure composition
  "mgcv"        # GAM fitting for dose nomogram
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Missing dependency: ", pkg, ". Install explicitly before running.")
  }
  library(pkg, character.only = TRUE)
}

rm(required_packages, pkg)
