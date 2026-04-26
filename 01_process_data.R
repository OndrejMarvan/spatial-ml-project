# =============================================================================
# 01_process_data.R — Process Manually Downloaded Raw Data
# =============================================================================
# Reads raw files and produces clean .gpkg / .csv files for 02_build_grid.R
#
# Your raw data layout:
#   data/raw/osm/           16 Geofabrik *-free.shp/ folders
#   data/raw/gios/2024/     2024_NO2_1g.xlsx, 2024_SO2_1g.xlsx, ...
#   data/raw/gios/          Metadane*.xlsx
#   data/raw/gus/           gadm41_POL.gpkg, LUDN*.csv
#   data/raw/tropomi/no2/   *.nc files from S5P-PAL
#
# Run: source("01_process_data.R")
#      process_osm()
#      process_gios()
#      process_tropomi()
#      process_gus()
# =============================================================================

source("00_setup.R")
library(sf)
library(tidyverse)
library(readxl)
library(glue)


# #############################################################################
# A. OSM — Merge Geofabrik voivodeships                                      #
# #############################################################################

process_osm <- function() {
  message("\n========== A. OSM ==========\n")
  t0 <- proc.time()["elapsed"]

  osm_dir <- file.path(CONFIG$dir_raw, "osm")

  # Find Geofabrik .gpkg files (may be inside subfolders from ZIP extraction)
  gpkg_files <- list.files(osm_dir, pattern = "\\.gpkg$",
                           recursive = TRUE, full.names = TRUE)
  # Exclude our own output files
  gpkg_files <- gpkg_files[!basename(gpkg_files) %in%
    c("major_roads.gpkg", "railways.gpkg", "industrial.gpkg", "power_plants.gpkg")]

  message(glue("{length(gpkg_files)} voivodeship .gpkg files found"))
  if (length(gpkg_files) == 0) stop("No Geofabrik .gpkg files in data/raw/osm/")

  merge_layer <- function(layer_name, fclass_vals, out_name) {
    message(glue("  {out_name}: "), appendLF = FALSE)
    parts <- base::list()

    for (i in seq_along(gpkg_files)) {
      tryCatch({
        x <- st_read(gpkg_files[i], layer = layer_name, quiet = TRUE)
        if ("fclass" %in% names(x)) x <- x[x$fclass %in% fclass_vals, ]
        if (nrow(x) > 0) parts[[length(parts) + 1]] <- x
      }, error = function(e) NULL)
    }

    if (length(parts) == 0) { message("no data"); return(NULL) }
    cols <- Reduce(intersect, lapply(parts, names))
    combined <- do.call(rbind, lapply(parts, function(p) p[, cols]))
    combined <- st_transform(combined, CONFIG$crs_pl)
    if ("osm_id" %in% names(combined))
      combined <- combined[!duplicated(combined$osm_id), ]
    st_write(combined, file.path(osm_dir, out_name), delete_dsn = TRUE, quiet = TRUE)
    message(glue("{nrow(combined)} features"))
  }

  merge_layer("gis_osm_roads_free",
              c("motorway", "trunk", "primary", "secondary"), "major_roads.gpkg")
  merge_layer("gis_osm_railways_free",
              c("rail"), "railways.gpkg")
  merge_layer("gis_osm_landuse_a_free",
              c("industrial"), "industrial.gpkg")

  message(glue("\nOSM done ({round(proc.time()['elapsed'] - t0)}s)"))
}


# #############################################################################
# B. GIOŚ — Station metadata + measurements from xlsx                        #
# #############################################################################

process_gios <- function() {
  message("\n========== B. GIOŚ ==========\n")

  gios_dir <- file.path(CONFIG$dir_raw, "gios")

  # --- B1. Station metadata ---
  meta_files <- list.files(gios_dir, pattern = "Metadane.*\\.xlsx$",
                           full.names = TRUE, ignore.case = TRUE)
  if (length(meta_files) == 0) stop("No Metadane xlsx in data/raw/gios/")

  message(glue("Reading: {basename(meta_files[1])}"))
  sheets <- excel_sheets(meta_files[1])
  message(glue("  Sheets: {paste(sheets, collapse = ', ')}"))

  # Read all sheets, find the one with station coordinates
  for (sh in sheets) {
    df <- read_excel(meta_files[1], sheet = sh)
    message(glue("  Sheet '{sh}': {nrow(df)} rows | Cols: {paste(head(names(df),8), collapse=', ')}"))
  }

  # Read the first sheet (usually contains station info)
  meta <- read_excel(meta_files[1], sheet = 1)
  write_csv(meta, file.path(gios_dir, "stations_raw.csv"))
  message(glue("  Saved stations_raw.csv ({nrow(meta)} rows)"))
  message("  >>> Open stations_raw.csv and identify lat/lon + station code columns <<<")

  # --- B2. Measurement xlsx → csv ---
  # GIOŚ xlsx structure:
  #   Row 1: Nr, 1, 2, 3, ...
  #   Row 2: "Kod stacji", station_code1, station_code2, ...
  #   Row 3: Indicator (NO2, SO2, etc.)
  #   Row 4: Averaging time (1g = hourly)
  #   Row 5: Unit (ug/m3)
  #   Row 6+: datetime, value1, value2, ...

  meas_dir <- file.path(gios_dir, "2024")
  target <- base::c("NO2_1g", "SO2_1g", "CO_1g", "PM10_1g", "PM25_1g", "O3_1g")

  for (t in target) {
    f <- file.path(meas_dir, paste0("2024_", t, ".xlsx"))
    if (!file.exists(f)) { message(glue("  Skip: {t} (not found)")); next }

    message(glue("  Reading 2024_{t}.xlsx ..."), appendLF = FALSE)
    tryCatch({
      # Read raw (no headers) to get station codes from row 2
      raw_header <- read_excel(f, col_names = FALSE, n_max = 5, .name_repair = "minimal")
      station_codes <- as.character(raw_header[2, ])  # Row 2 = station codes
      station_codes[1] <- "datetime"  # First column is always datetime

      # Read data (skip 5 header rows)
      df <- read_excel(f, skip = 5, col_names = FALSE, .name_repair = "minimal")
      names(df) <- station_codes[1:ncol(df)]

      poll <- sub("_1g$", "", t)
      out <- file.path(gios_dir, paste0(poll, "_2024.csv"))
      write_csv(df, out)
      message(glue(" {nrow(df)} rows, {ncol(df)-1} stations -> {poll}_2024.csv"))

      if (t == target[1]) {
        message(glue("    Station codes: {paste(head(station_codes[-1], 5), collapse=', ')}..."))
      }
    }, error = function(e) message(glue(" FAILED: {e$message}")))
  }

  message("\nGIOŚ done")
  message(">>> NEXT: inspect stations_raw.csv and NO2_2024.csv column structure <<<")
}


# #############################################################################
# C. TROPOMI — NetCDF → cropped GeoTIFF                                      #
# #############################################################################

process_tropomi <- function() {
  library(terra)
  message("\n========== C. TROPOMI ==========\n")

  nc_dir <- file.path(CONFIG$dir_raw, "tropomi", "no2")

  # Check if TIFs already exist (e.g. from Google Earth Engine)
  tif_files <- list.files(nc_dir, pattern = "\\.tif$")
  if (length(tif_files) > 0) {
    message(glue("{length(tif_files)} GeoTIFFs already present — no processing needed!"))
    message(glue("Files: {paste(head(tif_files, 3), collapse=', ')}..."))
    message("TROPOMI data is ready. Proceed to 02_build_grid.R")
    return(invisible(NULL))
  }

  # Otherwise try converting NetCDF files
  nc_files <- list.files(nc_dir, pattern = "\\.nc$", full.names = TRUE)
  if (length(nc_files) == 0) {
    message("No .nc or .tif files found in data/raw/tropomi/no2/")
    message("Download monthly GeoTIFFs from Google Earth Engine:")
    message("  https://code.earthengine.google.com/")
    return(invisible(NULL))
  }

  # Prefer tropospheric files
  tropo <- nc_files[grepl("tropospheric", nc_files)]
  total <- nc_files[!grepl("tropospheric", nc_files)]
  use <- if (length(tropo) > 0) tropo else total
  message(glue("{length(use)} files to process ({if(length(tropo)>0) 'tropospheric' else 'total column'})"))

  bb <- get_bbox()
  poland <- ext(bb["xmin"], bb["xmax"], bb["ymin"], bb["ymax"])
  index <- base::list()

  for (i in seq_along(use)) {
    f <- use[i]
    fname <- basename(f)

    # Parse start date: fortnight-YYYYMMDD
    m <- regmatches(fname, regexec("fortnight-(\\d{8})", fname))[[1]]
    if (length(m) < 2) { message(glue("  Skip: {fname} (no date)")); next }
    obs_date <- as.Date(m[2], format = "%Y%m%d")

    message(glue("  [{i}/{length(use)}] {obs_date} ..."), appendLF = FALSE)

    tryCatch({
      r <- rast(f)
      r_pl <- crop(r, poland)
      out_name <- paste0("no2_", format(obs_date, "%Y%m%d"), ".tif")
      writeRaster(r_pl, file.path(nc_dir, out_name), overwrite = TRUE)
      index[[i]] <- data.frame(file = out_name, obs_start = obs_date,
                                obs_end = obs_date + 13)
      message(glue(" OK ({out_name})"))
    }, error = function(e) message(glue(" FAIL: {e$message}")))
  }

  if (length(index) > 0) {
    idx <- do.call(rbind, index) |> dplyr::arrange(obs_start)
    write.csv(idx, file.path(nc_dir, "file_index.csv"), row.names = FALSE)
    message(glue("\n{nrow(idx)} TIFs | {min(idx$obs_start)} to {max(idx$obs_end)}"))

    # Warn about gaps
    dates <- sort(idx$obs_start)
    for (g in which(diff(dates) > 16))
      message(glue("  GAP: {dates[g]} → {dates[g+1]}"))
  }

  message("TROPOMI done")
}


# #############################################################################
# D. GUS — Verify GADM + population                                          #
# #############################################################################

process_gus <- function() {
  message("\n========== D. GUS ==========\n")

  gus_dir <- file.path(CONFIG$dir_raw, "gus")

  # GADM
  gadm <- file.path(gus_dir, "gadm41_POL.gpkg")
  if (file.exists(gadm)) {
    layers <- st_layers(gadm)
    message(glue("GADM: {gadm}"))
    message(glue("  Layers: {paste(layers$name, collapse=', ')}"))
    for (ly in layers$name) {
      n <- st_read(gadm, layer = ly, quiet = TRUE) |> nrow()
      message(glue("    {ly}: {n} features"))
    }
  } else {
    message("GADM NOT FOUND — download from gadm.org → Poland → GeoPackage")
  }

  # Population
  pop_files <- list.files(gus_dir, pattern = "LUDN.*\\.csv$", full.names = TRUE)
  if (length(pop_files) > 0) {
    message(glue("\nPopulation: {basename(pop_files[1])}"))
    # Try comma then semicolon delimiter
    tryCatch({
      p <- read.csv(pop_files[1], nrows = 3, fileEncoding = "UTF-8")
      if (ncol(p) <= 2) p <- read.csv2(pop_files[1], nrows = 3, fileEncoding = "UTF-8")
      message(glue("  Columns: {paste(names(p), collapse=', ')}"))
    }, error = function(e) message(glue("  Read error: {e$message}")))
  } else {
    message("Population CSV NOT FOUND")
  }

  message("\nGUS done")
}


# #############################################################################
message("\n=== 01_process_data.R loaded ===")
message("Run in order: process_osm()  process_gios()  process_tropomi()  process_gus()")
