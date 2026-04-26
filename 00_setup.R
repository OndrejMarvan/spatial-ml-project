# =============================================================================
# 00_setup.R — Package Installation & Project Configuration
# Spatial ML Project: Air Pollution in Poland
# Authors: Ondřej Marvan & Adam Jaworski (WNE UW, 2026)
# =============================================================================

# --- 1. Packages -------------------------------------------------------------

core_pkgs <- base::c("sf","terra","stars","tidyverse","lubridate","glue",
                      "readxl","viridis","corrplot")
analysis_pkgs <- base::c("spdep","spatialreg","ranger","caret","nnet",
                          "gstat","automap","dbscan","cluster")

for (pkg in base::c(core_pkgs, analysis_pkgs)) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste0("Installing ", pkg, "..."))
    tryCatch(utils::install.packages(pkg, dependencies = TRUE),
             warning = function(w) message(paste0("  WARN: ", w$message)),
             error   = function(e) message(paste0("  FAIL: ", e$message)))
  }
}

ok <- sapply(base::c(core_pkgs, analysis_pkgs),
             function(p) requireNamespace(p, quietly = TRUE))
message(paste0("Packages: ", sum(ok), "/", length(ok), " installed"))
if (any(!ok)) message("Missing: ", paste(names(ok[!ok]), collapse = ", "))

# --- 2. Config (bbox as separate values to avoid terra overriding c()) -------

CONFIG <- base::list(
  crs_pl = 2180L, crs_wgs84 = 4326L,
  year = 2024L,
  bbox_xmin = 14.07, bbox_ymin = 49.00, bbox_xmax = 24.15, bbox_ymax = 54.84,
  grid_res = 5000L,
  dir_raw = "data/raw", dir_processed = "data/processed",
  dir_output = "data/output", dir_figures = "figures"
)

for (d in base::c(CONFIG$dir_raw, CONFIG$dir_processed, CONFIG$dir_output,
                   CONFIG$dir_figures,
                   file.path(CONFIG$dir_raw, base::c("tropomi/no2","gios","gus","osm"))))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

# --- 3. Helpers --------------------------------------------------------------

#' Get Poland bbox as named vector — safe from terra's c() override
get_bbox <- function() {
  setNames(c(CONFIG$bbox_xmin, CONFIG$bbox_ymin, CONFIG$bbox_xmax, CONFIG$bbox_ymax),
           c("xmin", "ymin", "xmax", "ymax"))
}

#' Get Poland bbox as sf bbox object
get_bbox_sf <- function(crs = CONFIG$crs_wgs84) {
  sf::st_bbox(get_bbox(), crs = sf::st_crs(crs))
}

safe_read_sf  <- function(path, ...) {
  if (!file.exists(path)) { warning(paste0("Missing: ", path)); return(NULL) }
  sf::st_read(path, quiet = TRUE, ...)
}
safe_read_csv <- function(path, ...) {
  if (!file.exists(path)) { warning(paste0("Missing: ", path)); return(NULL) }
  readr::read_csv(path, show_col_types = FALSE, ...)
}
calc_metrics  <- function(actual, predicted) {
  ok <- !is.na(actual) & !is.na(predicted)
  a <- actual[ok]; p <- predicted[ok]
  if (length(a) < 3) return(NULL)
  data.frame(n=length(a), R2=cor(a,p)^2, RMSE=sqrt(mean((a-p)^2)), MAE=mean(abs(a-p)))
}

message("Setup OK | Year: ", CONFIG$year, " | Grid: ", CONFIG$grid_res/1000, " km")

# --- 4. Data integrity check -------------------------------------------------
# Warns if raw data folders are empty — catches accidental deletions early
check_data <- function() {
  checks <- base::list(
    "TROPOMI (.tif)"  = length(list.files(file.path(CONFIG$dir_raw, "tropomi/no2"),
                                          pattern = "\\.tif$")),
    "GIOŚ (2024/)"    = length(list.files(file.path(CONFIG$dir_raw, "gios/2024"))),
    "GIOŚ (metadata)" = length(list.files(file.path(CONFIG$dir_raw, "gios"),
                                          pattern = "Metadane|stations")),
    "OSM (folders)"   = length(list.dirs(file.path(CONFIG$dir_raw, "osm"),
                                         recursive = FALSE)),
    "GUS (GADM)"      = length(list.files(file.path(CONFIG$dir_raw, "gus"),
                                          pattern = "gadm.*\\.gpkg$")),
    "GUS (population)"= length(list.files(file.path(CONFIG$dir_raw, "gus"),
                                          pattern = "LUDN.*\\.csv$"))
  )

  message("\n--- Data check ---")
  all_ok <- TRUE
  for (nm in names(checks)) {
    n <- checks[[nm]]
    status <- if (n > 0) paste0(n, " files") else "MISSING!"
    flag <- if (n > 0) "OK" else "!!"
    message(sprintf("  [%s] %-20s %s", flag, nm, status))
    if (n == 0) all_ok <- FALSE
  }
  if (!all_ok) message("\n  WARNING: Some data is missing. DO NOT re-download scripts over data folder!")
  return(invisible(all_ok))
}

check_data()
