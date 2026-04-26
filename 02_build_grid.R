# =============================================================================
# 02_build_grid.R — Build Analysis Grid & Integrate All Data
# =============================================================================
# Run AFTER 01_process_data.R has produced .gpkg and .csv files.
#
# Creates:
#   data/processed/analysis_grid.gpkg     — grid with OSM + GUS features
#   data/processed/tropomi_grid.csv       — TROPOMI values per cell per period
#   data/processed/stations.gpkg          — GIOŚ stations with coordinates
#   data/processed/no2_station_annual.csv — annual mean NO2 per station
# =============================================================================

source("00_setup.R")
library(sf)
library(terra)
library(tidyverse)
library(readxl)
library(glue)

# =============================================================================
# 1. GRID
# =============================================================================
message("\n========== 1. Building grid ==========\n")

# Poland bbox as sfc
bbox_vec <- get_bbox_sf(CONFIG$crs_wgs84)
bbox_sfc <- st_as_sfc(bbox_vec) |> st_transform(CONFIG$crs_pl)

cs <- as.numeric(rep(CONFIG$grid_res, 2))
grid <- st_make_grid(bbox_sfc, cellsize = cs, what = "polygons") |>
  st_as_sf()

# Rename geometry column
names(grid)[1] <- "geometry"
st_geometry(grid) <- "geometry"
grid$cell_id <- seq_len(nrow(grid))

message(glue("Raw grid: {nrow(grid)} cells"))

# Clip to Poland boundary (GADM)
gadm_path <- file.path(CONFIG$dir_raw, "gus", "gadm41_POL.gpkg")
if (file.exists(gadm_path)) {
  layers <- st_layers(gadm_path)$name
  # Use the coarsest layer (country boundary) for clipping
  poland <- st_read(gadm_path, layer = layers[1], quiet = TRUE) |>
    st_transform(CONFIG$crs_pl) |>
    st_union()

  inside <- st_intersects(grid, poland, sparse = FALSE)[, 1]
  grid <- grid[inside, ]
  grid$cell_id <- seq_len(nrow(grid))
  message(glue("After clipping to Poland: {nrow(grid)} cells"))
} else {
  message("No GADM file — using full bbox grid")
}

# Add centroid coordinates
centroids <- st_coordinates(st_centroid(grid))
grid$x_coord <- centroids[, 1]
grid$y_coord <- centroids[, 2]


# =============================================================================
# 2. OSM FEATURES
# =============================================================================
message("\n========== 2. OSM features ==========\n")

osm_dir <- file.path(CONFIG$dir_raw, "osm")

# --- Roads ---
roads_path <- file.path(osm_dir, "major_roads.gpkg")
if (file.exists(roads_path)) {
  message("Computing road density...")
  roads <- st_read(roads_path, quiet = TRUE)

  tryCatch({
    ri <- st_intersection(roads, grid[, "cell_id"])
    ri$len_m <- as.numeric(st_length(ri))

    road_stats <- ri |>
      st_drop_geometry() |>
      group_by(cell_id) |>
      summarise(total_road_km = sum(len_m) / 1000, .groups = "drop")

    if ("fclass" %in% names(ri)) {
      mw <- ri |> st_drop_geometry() |>
        filter(fclass == "motorway") |>
        group_by(cell_id) |>
        summarise(motorway_km = sum(len_m) / 1000, .groups = "drop")
      road_stats <- left_join(road_stats, mw, by = "cell_id")
    }

    grid <- left_join(grid, road_stats, by = "cell_id")
    grid$total_road_km[is.na(grid$total_road_km)] <- 0
    if ("motorway_km" %in% names(grid)) grid$motorway_km[is.na(grid$motorway_km)] <- 0
    message(glue("  Road density: done ({nrow(road_stats)} cells with roads)"))
  }, error = function(e) warning(paste0("Roads failed: ", e$message)))
} else {
  message("  No major_roads.gpkg — run process_osm() first")
}

# --- Railways ---
rail_path <- file.path(osm_dir, "railways.gpkg")
if (file.exists(rail_path)) {
  message("Computing railway density...")
  tryCatch({
    rail <- st_read(rail_path, quiet = TRUE)
    ri <- st_intersection(rail, grid[, "cell_id"])
    ri$len_m <- as.numeric(st_length(ri))
    rs <- ri |> st_drop_geometry() |>
      group_by(cell_id) |>
      summarise(railway_km = sum(len_m) / 1000, .groups = "drop")
    grid <- left_join(grid, rs, by = "cell_id")
    grid$railway_km[is.na(grid$railway_km)] <- 0
    message(glue("  Railway density: done"))
  }, error = function(e) warning(paste0("Railways failed: ", e$message)))
}

# --- Industrial ---
ind_path <- file.path(osm_dir, "industrial.gpkg")
if (file.exists(ind_path)) {
  message("Computing industrial coverage...")
  tryCatch({
    ind <- st_read(ind_path, quiet = TRUE)
    ii <- st_intersection(ind, grid[, "cell_id"])
    ii$area_m2 <- as.numeric(st_area(ii))
    is_df <- ii |> st_drop_geometry() |>
      group_by(cell_id) |>
      summarise(industrial_m2 = sum(area_m2), .groups = "drop")
    grid <- left_join(grid, is_df, by = "cell_id")
    cell_area <- as.numeric(st_area(grid[1, ]))
    grid$industrial_pct <- ifelse(is.na(grid$industrial_m2), 0,
                                  grid$industrial_m2 / cell_area * 100)
    grid$industrial_m2 <- NULL
    message("  Industrial coverage: done")
  }, error = function(e) warning(paste0("Industrial failed: ", e$message)))
}

# --- Power plants (from earlier Overpass download) ---
pp_path <- file.path(osm_dir, "power_plants.gpkg")
if (file.exists(pp_path)) {
  message("Counting power plants...")
  tryCatch({
    pp <- st_read(pp_path, quiet = TRUE) |> st_transform(CONFIG$crs_pl)
    buffers <- st_buffer(st_centroid(grid), 10000)
    grid$power_plants_10km <- lengths(st_intersects(buffers, pp))
    message("  Power plant count: done")
  }, error = function(e) warning(paste0("Power plants failed: ", e$message)))
}


# =============================================================================
# 3. GUS DEMOGRAPHICS (via GADM spatial join)
# =============================================================================
message("\n========== 3. GUS demographics ==========\n")

if (file.exists(gadm_path)) {
  layers <- st_layers(gadm_path)$name
  # Use the finest available level
  finest_layer <- layers[length(layers)]
  admin <- st_read(gadm_path, layer = finest_layer, quiet = TRUE) |>
    st_transform(CONFIG$crs_pl)
  message(glue("GADM layer '{finest_layer}': {nrow(admin)} units"))

  # Spatial join: each grid centroid → admin unit
  grid_pts <- st_centroid(grid[, "cell_id"])
  joined <- st_join(grid_pts, admin, join = st_within) |> st_drop_geometry()

  # Keep useful columns (NAME_1 = voivodeship, NAME_2 = powiat)
  name_cols <- intersect(names(joined),
                         base::c("cell_id", "NAME_1", "NAME_2", "TYPE_2"))
  if (length(name_cols) > 1) {
    grid <- left_join(grid, joined[, name_cols], by = "cell_id")
    message(glue("  Joined admin names to grid"))
  }

  # Population CSV
  pop_files <- list.files(file.path(CONFIG$dir_raw, "gus"),
                          pattern = "LUDN.*\\.csv$", full.names = TRUE)
  if (length(pop_files) > 0) {
    message(glue("  Population file found: {basename(pop_files[1])}"))
    message("  (Population join requires matching TERYT codes — implement after inspecting CSV)")
  }
} else {
  message("  No GADM — skipping")
}


# =============================================================================
# 4. TROPOMI RASTER EXTRACTION
# =============================================================================
message("\n========== 4. TROPOMI extraction ==========\n")

tif_dir <- file.path(CONFIG$dir_raw, "tropomi", "no2")
tif_files <- list.files(tif_dir, pattern = "\\.tif$", full.names = TRUE)

if (length(tif_files) > 0) {
  message(glue("Found {length(tif_files)} GeoTIFFs"))

  tropomi_list <- base::list()
  for (i in seq_along(tif_files)) {
    f <- tif_files[i]
    fn <- basename(f)
    # Parse date: supports no2_2024_01.tif and no2_20240101.tif
    dm8 <- regmatches(fn, regexec("(\\d{8})", fn))[[1]]
    dm_ym <- regmatches(fn, regexec("(\\d{4})_(\\d{2})", fn))[[1]]

    if (length(dm8) >= 2) {
      obs_date <- as.Date(dm8[2], format = "%Y%m%d")
    } else if (length(dm_ym) >= 3) {
      obs_date <- as.Date(paste0(dm_ym[2], "-", dm_ym[3], "-01"))
    } else {
      message(glue("  Skip: {fn} (no date)"))
      next
    }

    message(glue("  [{i}/{length(tif_files)}] {fn} -> {obs_date}"), appendLF = FALSE)

    tryCatch({
      r <- rast(f)
      if (!identical(crs(r, describe = TRUE)$code, as.character(CONFIG$crs_pl))) {
        r <- project(r, paste0("EPSG:", CONFIG$crs_pl))
      }
      vals <- terra::extract(r, vect(grid), fun = mean, na.rm = TRUE)
      tropomi_list[[i]] <- data.frame(
        cell_id = grid$cell_id, obs_date = obs_date,
        tropomi_no2 = vals[[2]]
      )
      message(" OK")
    }, error = function(e) {
      message(glue(" FAIL: {e$message}"))
    })
  }

  if (length(tropomi_list) > 0) {
    tropomi_df <- bind_rows(tropomi_list)
    write_csv(tropomi_df, file.path(CONFIG$dir_processed, "tropomi_grid.csv"))
    message(glue("  Saved tropomi_grid.csv: {nrow(tropomi_df)} rows"))

    # Annual mean per cell → join to grid
    annual <- tropomi_df |>
      group_by(cell_id) |>
      summarise(tropomi_no2_mean = mean(tropomi_no2, na.rm = TRUE), .groups = "drop")
    grid <- left_join(grid, annual, by = "cell_id")
    message("  Annual NO2 mean joined to grid")
  }
} else {
  message("  No TIFFs — run process_tropomi() first")
}


# =============================================================================
# 5. GIOŚ STATIONS
# =============================================================================
message("\n========== 5. GIOŚ stations ==========\n")

gios_dir <- file.path(CONFIG$dir_raw, "gios")

# Read station metadata
meta_csv <- file.path(gios_dir, "stations_raw.csv")
if (file.exists(meta_csv)) {
  meta <- read_csv(meta_csv, show_col_types = FALSE)
  message(glue("Station metadata: {nrow(meta)} rows"))
  message(glue("  Columns: {paste(names(meta), collapse=', ')}"))

  # Try to find lat/lon columns (common names in GIOŚ metadata)
  lat_col <- names(meta)[grepl("WGS84.*N|lat|szer|geogr.*szer", names(meta), ignore.case = TRUE)]
  lon_col <- names(meta)[grepl("WGS84.*E|lon|dlug|geogr.*dl", names(meta), ignore.case = TRUE)]
  code_col <- names(meta)[grepl("Kod\\.stacji|Kod.stacji|kod.*sta|station.*code", names(meta), ignore.case = TRUE)]

  if (length(lat_col) > 0 && length(lon_col) > 0) {
    message(glue("  Detected lat={lat_col[1]}, lon={lon_col[1]}"))

    stations_sf <- meta |>
      mutate(lat = as.numeric(.data[[lat_col[1]]]),
             lon = as.numeric(.data[[lon_col[1]]])) |>
      filter(!is.na(lat), !is.na(lon)) |>
      st_as_sf(coords = base::c("lon", "lat"), crs = CONFIG$crs_wgs84) |>
      st_transform(CONFIG$crs_pl)

    # Assign to grid cells
    stations_sf <- st_join(stations_sf, grid[, "cell_id"], join = st_within)

    st_write(stations_sf, file.path(CONFIG$dir_processed, "stations.gpkg"),
             delete_dsn = TRUE, quiet = TRUE)
    message(glue("  Saved stations.gpkg: {nrow(stations_sf)} stations"))
  } else {
    message("  Could not auto-detect lat/lon columns.")
    message("  Inspect stations_raw.csv and add manual column mapping.")
  }
} else {
  message("  No stations_raw.csv — run process_gios() first")
}

# Read NO2 measurements and compute annual mean per station
no2_csv <- file.path(gios_dir, "NO2_2024.csv")
if (file.exists(no2_csv)) {
  message("\nProcessing NO2 measurements...")
  no2 <- read_csv(no2_csv, show_col_types = FALSE)
  message(glue("  NO2 data: {nrow(no2)} rows, {ncol(no2)} cols"))
  message(glue("  Columns: {paste(head(names(no2), 8), collapse=', ')}..."))

  # GIOŚ xlsx format is typically: first col = date/time, remaining cols = station codes
  # Each cell is the hourly measurement value
  # We need to pivot to long format

  # Check if first column looks like dates
  col1 <- names(no2)[1]
  message(glue("  First column ('{col1}'): {class(no2[[col1]])[1]}"))

  tryCatch({
    no2_long <- no2 |>
      pivot_longer(cols = -1, names_to = "station_code", values_to = "no2_ugm3") |>
      mutate(no2_ugm3 = as.numeric(no2_ugm3)) |>
      filter(!is.na(no2_ugm3))

    # Annual mean per station
    station_annual <- no2_long |>
      group_by(station_code) |>
      summarise(
        no2_mean = mean(no2_ugm3, na.rm = TRUE),
        no2_median = median(no2_ugm3, na.rm = TRUE),
        no2_max = max(no2_ugm3, na.rm = TRUE),
        n_obs = n(),
        .groups = "drop"
      )

    write_csv(station_annual, file.path(CONFIG$dir_processed, "no2_station_annual.csv"))
    message(glue("  Saved no2_station_annual.csv: {nrow(station_annual)} stations"))
  }, error = function(e) {
    message(glue("  Pivot failed: {e$message}"))
    message("  The xlsx format may differ — inspect NO2_2024.csv manually")
  })
} else {
  message("  No NO2_2024.csv — run process_gios() first")
}


# =============================================================================
# 6. SAVE GRID
# =============================================================================
message("\n========== 6. Saving ==========\n")

st_write(grid, file.path(CONFIG$dir_processed, "analysis_grid.gpkg"),
         delete_dsn = TRUE, quiet = TRUE)

feature_cols <- setdiff(names(grid),
                        base::c("geometry","cell_id","x_coord","y_coord"))
message(glue("Grid: {nrow(grid)} cells"))
message(glue("Features: {paste(feature_cols, collapse=', ')}"))

for (col in feature_cols) {
  if (is.numeric(grid[[col]])) {
    v <- grid[[col]][!is.na(grid[[col]])]
    if (length(v) > 0)
      message(sprintf("  %-20s mean=%.4g  sd=%.4g  n=%d", col, mean(v), sd(v), length(v)))
  }
}

message(glue("\nSaved: {CONFIG$dir_processed}/analysis_grid.gpkg"))
message("02_build_grid.R complete.")
