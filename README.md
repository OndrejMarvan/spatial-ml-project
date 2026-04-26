# Spatial ML Project — Air Pollution Modeling in Poland

**Course:** Spatial Machine Learning in R (WNE UW, 2026)  
**Authors:** Ondřej Marvan & Adam Jaworski  
**Supervisor:** prof. Katarzyna Kopczewska

## Overview

Predicting ground-level NO2 across Poland using Sentinel-5P TROPOMI satellite data, GIOS ground stations, GUS demographics, and OpenStreetMap infrastructure features.

## Files

```
00_setup.R              — Packages, config, helpers, data check
01_process_data.R       — Process downloaded data -> clean .gpkg/.csv
02_build_grid.R         — 5km grid + integrate all data sources
03_analysis.R           — EDA + RF/ANN + DBSCAN/LISA/PAM + Kriging
04_report.Rmd           — Final RPubs report
```

## How to Run (paste into RStudio Console)

```r
# ---- Step 1: Set working directory ----
setwd("~/Documents/GitHub/spatial-ml-project")

# ---- Step 2: Setup + process raw data ----
source("00_setup.R")
source("01_process_data.R")
process_osm()
process_gios()
process_tropomi()
process_gus()

# ---- Step 3: Build analysis grid ----
source("02_build_grid.R")

# ---- Step 4: Run all analysis ----
source("03_analysis.R")
run_eda()
run_supervised()
run_unsupervised()
run_kriging()

# ---- Step 5: Generate report ----
# Open 04_report.Rmd -> click Knit -> publish to RPubs
```

## Data Sources

| Source | Format | Location |
|--------|--------|----------|
| TROPOMI S5P (GEE) | .tif | data/raw/tropomi/no2/ |
| GIOS bulk 2024 | .xlsx -> .csv | data/raw/gios/2024/ |
| Geofabrik OSM | .gpkg | data/raw/osm/ |
| GADM boundaries | .gpkg | data/raw/gus/gadm41_POL.gpkg |
| GUS BDL population | .csv | data/raw/gus/LUDN*.csv |

## Methods

- **Supervised ML:** Random Forest, ANN, (GW-RF if SpatialML available)
- **Unsupervised ML:** DBSCAN hotspots, LISA (Local Moran's I), PAM clustering
- **Kriging:** Ordinary Kriging, Universal Kriging
- **Spatial analysis:** Moran's I, variogram fitting, spatial block CV

## Git — Push Updates to GitHub

```bash
# Open terminal, navigate to project
cd /home/ondrej-marvan/Documents/GitHub/spatial-ml-project

# Check what changed
git status

# Stage all changes
git add .

# Commit with a message
git commit -m "describe what you changed"

# Push to GitHub
git push
```

**Quick one-liner** (stage + commit + push):
```bash
cd ~/Documents/GitHub/spatial-ml-project && git add . && git commit -m "update" && git push
```

**Important:** The `data/` folder is large. Add a `.gitignore` so raw data stays local:
```bash
echo "data/raw/" >> .gitignore
echo "data/processed/" >> .gitignore
echo "data/output/" >> .gitignore
echo "*.RData" >> .gitignore
echo ".Rhistory" >> .gitignore
echo "desktop.ini" >> .gitignore
git add .gitignore
git commit -m "add gitignore"
git push
```

After this, only R scripts, the Rmd, and the README go to GitHub — data stays on your machine. Your colleague Adam can clone with:
```bash
git clone https://github.com/OndrejMarvan/spatial-ml-project.git
```
Then download the raw data separately following the instructions above.
