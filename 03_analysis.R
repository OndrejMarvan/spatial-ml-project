# =============================================================================
# 03_analysis.R — EDA + Supervised ML + Unsupervised ML + Kriging
# =============================================================================
source("00_setup.R")
library(sf); library(terra); library(tidyverse); library(spdep)
library(glue); library(viridis)

grid     <- safe_read_sf(file.path(CONFIG$dir_processed, "analysis_grid.gpkg"))
tropomi  <- safe_read_csv(file.path(CONFIG$dir_processed, "tropomi_grid.csv"))
stations <- safe_read_sf(file.path(CONFIG$dir_processed, "stations.gpkg"))
# st_annual <- safe_read_csv(file.path(CONFIG$dir_processed, "no2_station_annual.csv"))
if (is.null(grid)) stop("Run 02_build_grid.R first")

# ===== A. EDA ================================================================
run_eda <- function() {
  library(corrplot)
  message("\n========== A. EDA ==========\n")
  if ("tropomi_no2_mean" %in% names(grid)) {
    p1 <- ggplot(grid |> filter(!is.na(tropomi_no2_mean))) +
      geom_sf(aes(fill = tropomi_no2_mean), color = NA) +
      scale_fill_viridis(option = "inferno", name = expression(NO[2])) +
      labs(title = "Annual Mean Tropospheric NO2 - Poland 2024") + theme_minimal()
    ggsave(file.path(CONFIG$dir_figures, "no2_map.png"), p1, width=10, height=12, dpi=300)
    message("  Saved: no2_map.png")
  }
  if (!is.null(tropomi)) {
    monthly <- tropomi |> group_by(obs_date) |>
      summarise(mean_no2=mean(tropomi_no2,na.rm=TRUE), sd_no2=sd(tropomi_no2,na.rm=TRUE), .groups="drop")
    p2 <- ggplot(monthly, aes(obs_date, mean_no2)) +
      geom_line(linewidth=1, color="steelblue") +
      geom_ribbon(aes(ymin=mean_no2-sd_no2, ymax=mean_no2+sd_no2), alpha=0.2, fill="steelblue") +
      labs(title="TROPOMI NO2 Over Time", x=NULL, y="NO2 density") + theme_minimal()
    ggsave(file.path(CONFIG$dir_figures, "no2_timeseries.png"), p2, width=10, height=5, dpi=300)
    message("  Saved: no2_timeseries.png")
  }
  num_cols <- grid |> st_drop_geometry() |>
    select(where(is.numeric), -cell_id, -x_coord, -y_coord) |> names()
  if (length(num_cols) >= 3) {
    cor_data <- grid |> st_drop_geometry() |> select(all_of(num_cols)) |> drop_na()
    if (nrow(cor_data) > 30) {
      png(file.path(CONFIG$dir_figures, "correlations.png"), 800, 800, res=100)
      corrplot(cor(cor_data), method="color", type="lower", tl.cex=0.7, addCoef.col="black", number.cex=0.6)
      dev.off()
      message("  Saved: correlations.png")
    }
  }
  if ("tropomi_no2_mean" %in% names(grid)) {
    clean <- grid |> filter(!is.na(tropomi_no2_mean))
    coords <- st_coordinates(st_centroid(clean))
    nb <- knearneigh(coords, k=8) |> knn2nb()
    lw <- nb2listw(nb, style="W")
    mi <- moran.test(clean$tropomi_no2_mean, lw)
    message(sprintf("  Moran's I = %.4f (p = %s)", mi$estimate[1], format.pval(mi$p.value)))
  }
  summ <- grid |> st_drop_geometry() |> select(all_of(num_cols)) |>
    pivot_longer(everything()) |> group_by(name) |>
    summarise(n=sum(!is.na(value)), mean=mean(value,na.rm=TRUE),
              sd=sd(value,na.rm=TRUE), .groups="drop")
  write_csv(summ, file.path(CONFIG$dir_output, "summary_stats.csv"))
  print(summ)
  message("\nEDA complete")
}

# ===== B. SUPERVISED ML ======================================================
run_supervised <- function() {
  library(ranger); library(caret); library(nnet)
  message("\n========== B. Supervised ML ==========\n")
  if (is.null(stations) || is.null(st_annual)) stop("Need stations + measurements")
  st_cols <- names(stations)
  code_col <- st_cols[grepl("kod|code|station", st_cols, ignore.case=TRUE)]
  if (length(code_col) == 0) { message("Cannot find station code column"); return(NULL) }
  code_col <- code_col[1]
  message(glue("Station code column: '{code_col}'"))
  model_df <- stations |> st_drop_geometry() |>
    select(all_of(code_col), cell_id) |>
    inner_join(st_annual, by=setNames("station_code", code_col)) |>
    inner_join(grid |> st_drop_geometry(), by="cell_id") |>
    filter(!is.na(no2_mean))
  message(glue("Model dataset: {nrow(model_df)} stations"))
  if (nrow(model_df) < 20) { message("Too few stations"); return(NULL) }
  exclude <- base::c(code_col, "cell_id", "no2_mean", "no2_median", "no2_max",
                      "n_obs", "NAME_1", "NAME_2", "TYPE_2")
  feat_cols <- setdiff(names(model_df)[sapply(model_df, is.numeric)], exclude)
  message(glue("Features ({length(feat_cols)}): {paste(feat_cols, collapse=', ')}"))
  if (length(feat_cols) < 2) { message("Not enough features"); return(NULL) }
  set.seed(42)
  n_blocks <- max(2, min(5, nrow(model_df) %/% 10))
  km <- kmeans(model_df[, base::c("x_coord","y_coord")], centers=n_blocks)
  model_df$block <- km$cluster
  test_block <- sample(unique(model_df$block), 1)
  train_df <- model_df[model_df$block != test_block, ]
  test_df  <- model_df[model_df$block == test_block, ]
  message(glue("Train: {nrow(train_df)} | Test: {nrow(test_df)}"))
  fml <- as.formula(paste("no2_mean ~", paste(feat_cols, collapse=" + ")))
  message("\nRandom Forest...")
  rf <- ranger(fml, data=train_df, num.trees=500, importance="impurity", seed=42)
  message(sprintf("  OOB R2=%.3f RMSE=%.3f", rf$r.squared, sqrt(rf$prediction.error)))
  rf_pred <- predict(rf, data=test_df)$predictions
  rf_met <- calc_metrics(test_df$no2_mean, rf_pred)
  if (!is.null(rf_met)) message(sprintf("  Test R2=%.3f RMSE=%.3f MAE=%.3f", rf_met$R2, rf_met$RMSE, rf_met$MAE))
  imp <- tibble(var=names(rf$variable.importance), imp=rf$variable.importance) |> arrange(desc(imp))
  p_imp <- ggplot(head(imp,15), aes(reorder(var,imp), imp)) +
    geom_col(fill="steelblue") + coord_flip() +
    labs(title="Variable Importance (RF)", x=NULL) + theme_minimal()
  ggsave(file.path(CONFIG$dir_figures, "rf_importance.png"), p_imp, width=8, height=6, dpi=300)
  message("\nANN...")
  ann_met <- NULL
  tryCatch({
    pp <- preProcess(train_df[feat_cols], method=base::c("center","scale"))
    tr_s <- predict(pp, train_df); te_s <- predict(pp, test_df)
    ctrl <- trainControl(method="cv", number=5)
    ann <- train(fml, data=tr_s, method="nnet", trControl=ctrl,
                 tuneGrid=expand.grid(size=base::c(8,16), decay=base::c(0.01,0.001)),
                 linout=TRUE, maxit=500, trace=FALSE)
    ann_pred <- predict(ann, te_s)
    ann_met <- calc_metrics(test_df$no2_mean, ann_pred)
    if (!is.null(ann_met)) message(sprintf("  ANN R2=%.3f RMSE=%.3f", ann_met$R2, ann_met$RMSE))
  }, error = function(e) message(paste0("  ANN failed: ", e$message)))
  comp <- bind_rows(RF=rf_met, ANN=ann_met, .id="model")
  write_csv(comp, file.path(CONFIG$dir_output, "model_comparison.csv"))
  print(comp)
  gd <- grid |> st_drop_geometry()
  if (all(feat_cols %in% names(gd))) {
    grid$predicted_no2 <- predict(rf, data=gd)$predictions
    st_write(grid, file.path(CONFIG$dir_output, "predicted_no2.gpkg"), delete_dsn=TRUE, quiet=TRUE)
    p <- ggplot(grid |> filter(!is.na(predicted_no2))) +
      geom_sf(aes(fill=predicted_no2), color=NA) +
      scale_fill_viridis(option="inferno", name="NO2") +
      labs(title="Predicted Ground-Level NO2 (RF)") + theme_minimal()
    ggsave(file.path(CONFIG$dir_figures, "predicted_no2.png"), p, width=10, height=12, dpi=300)
  }
  base::list(rf=rf, model_df=model_df, feat_cols=feat_cols)
}

# ===== C. UNSUPERVISED ML ===================================================
run_unsupervised <- function() {
  library(dbscan); library(cluster)
  message("\n========== C. Unsupervised ML ==========\n")
  vc <- if ("tropomi_no2_mean" %in% names(grid)) "tropomi_no2_mean" else NULL
  if (is.null(vc)) { message("No NO2 in grid"); return(NULL) }
  clean <- grid |> filter(!is.na(.data[[vc]]))
  message("DBSCAN...")
  thr <- quantile(clean[[vc]], 0.75, na.rm=TRUE)
  hot <- clean |> filter(.data[[vc]] >= thr)
  co <- st_coordinates(st_centroid(hot))
  db <- dbscan(co, eps=15000, minPts=5)
  hot$cluster <- db$cluster
  nc <- length(unique(db$cluster[db$cluster > 0]))
  message(sprintf("  %d clusters (%d noise)", nc, sum(db$cluster==0)))
  p_hot <- ggplot() + geom_sf(data=clean, fill="grey95", color=NA) +
    geom_sf(data=hot |> filter(cluster>0), aes(fill=factor(cluster)), color=NA, alpha=0.7) +
    scale_fill_brewer(palette="Set1", name="Cluster") +
    labs(title="NO2 Hotspots (DBSCAN)") + theme_minimal()
  ggsave(file.path(CONFIG$dir_figures, "hotspots.png"), p_hot, width=10, height=12, dpi=300)
  message("\nLISA...")
  ca <- st_coordinates(st_centroid(clean))
  nb <- knearneigh(ca, k=8) |> knn2nb()
  lw <- nb2listw(nb, style="W")
  lm <- localmoran(clean[[vc]], lw)
  vz <- scale(clean[[vc]])[,1]; lz <- lag.listw(lw, vz)
  clean$lisa <- case_when(lm[,5]>0.05 ~ "Not significant",
    vz>0 & lz>0 ~ "High-High", vz<0 & lz<0 ~ "Low-Low",
    vz>0 & lz<0 ~ "High-Low",  vz<0 & lz>0 ~ "Low-High")
  lcol <- base::c("High-High"="#d7191c","Low-Low"="#2c7bb6","High-Low"="#fdae61",
                   "Low-High"="#abd9e9","Not significant"="grey90")
  p_l <- ggplot(clean) + geom_sf(aes(fill=lisa), color=NA) +
    scale_fill_manual(values=lcol, name="LISA") +
    labs(title="LISA Cluster Map - NO2") + theme_minimal()
  ggsave(file.path(CONFIG$dir_figures, "lisa.png"), p_l, width=10, height=12, dpi=300)
  message("  Saved: lisa.png")
  clust_cols <- intersect(base::c("tropomi_no2_mean","total_road_km","industrial_pct","railway_km"), names(grid))
  if (length(clust_cols) >= 2) {
    message("\nPAM clustering...")
    cd <- clean |> st_drop_geometry() |> select(cell_id, all_of(clust_cols)) |> drop_na()
    if (nrow(cd) >= 50) {
      sc <- scale(cd |> select(-cell_id))
      sil <- sapply(2:8, function(k) mean(silhouette(pam(sc,k))[,"sil_width"]))
      bk <- (2:8)[which.max(sil)]
      message(sprintf("  Best k=%d (sil=%.3f)", bk, max(sil)))
      cd$cluster <- pam(sc, bk)$clustering
      cj <- clean |> left_join(cd |> select(cell_id, cluster), by="cell_id")
      p_p <- ggplot(cj |> filter(!is.na(cluster))) +
        geom_sf(aes(fill=factor(cluster)), color=NA) +
        scale_fill_brewer(palette="Set2", name="Cluster") +
        labs(title=sprintf("Spatial Typology (PAM k=%d)", bk)) + theme_minimal()
      ggsave(file.path(CONFIG$dir_figures, "pam_clusters.png"), p_p, width=10, height=12, dpi=300)
    }
  }
  message("\nUnsupervised ML complete")
}

# ===== D. KRIGING ============================================================
run_kriging <- function() {
  library(gstat); library(automap)
  message("\n========== D. Kriging ==========\n")
  if (is.null(stations) || is.null(st_annual)) stop("Need stations + measurements")
  st_cols <- names(stations)
  code_col <- st_cols[grepl("kod|code|station", st_cols, ignore.case=TRUE)]
  if (length(code_col)==0) { message("No station code col"); return(NULL) }
  st_pts <- stations |> inner_join(st_annual, by=setNames("station_code", code_col[1])) |>
    filter(!is.na(no2_mean))
  message(glue("Kriging with {nrow(st_pts)} stations"))
  if (nrow(st_pts) < 10) { message("Too few"); return(NULL) }
  message("Variogram...")
  vf <- autofitVariogram(no2_mean ~ 1, input_data=st_pts)
  message(sprintf("  %s | Nug=%.1f Sill=%.1f Range=%.0f", vf$var_model$model[2],
                  vf$var_model$psill[1], sum(vf$var_model$psill), vf$var_model$range[2]))
  png(file.path(CONFIG$dir_figures, "variogram.png"), 600, 400)
  plot(vf); dev.off()
  message("OK kriging...")
  pp <- st_centroid(grid)
  ok <- krige(no2_mean ~ 1, locations=st_pts, newdata=pp, model=vf$var_model)
  grid$ok_pred <- ok$var1.pred; grid$ok_se <- sqrt(ok$var1.var)
  p_ok <- ggplot(grid |> filter(!is.na(ok_pred))) +
    geom_sf(aes(fill=ok_pred), color=NA) +
    scale_fill_viridis(option="inferno", name="NO2") +
    labs(title="Ordinary Kriging - NO2") + theme_minimal()
  ggsave(file.path(CONFIG$dir_figures, "kriging.png"), p_ok, width=10, height=12, dpi=300)
  message("CV...")
  cv <- krige.cv(no2_mean ~ 1, locations=st_pts, model=vf$var_model, nfold=5)
  cv_m <- calc_metrics(cv$observed, cv$observed - cv$residual)
  if (!is.null(cv_m)) message(sprintf("  OK R2=%.3f RMSE=%.3f", cv_m$R2, cv_m$RMSE))
  message("UK...")
  uk_m <- NULL
  tryCatch({
    cr <- st_coordinates(st_pts); st_pts$x <- cr[,1]; st_pts$y <- cr[,2]
    cg <- st_coordinates(pp); pp$x <- cg[,1]; pp$y <- cg[,2]
    vf2 <- autofitVariogram(no2_mean ~ x + y, input_data=st_pts)
    uk <- krige(no2_mean ~ x+y, locations=st_pts, newdata=pp, model=vf2$var_model)
    grid$uk_pred <- uk$var1.pred
    cv2 <- krige.cv(no2_mean ~ x+y, locations=st_pts, model=vf2$var_model, nfold=5)
    uk_m <- calc_metrics(cv2$observed, cv2$observed - cv2$residual)
    if (!is.null(uk_m)) message(sprintf("  UK R2=%.3f RMSE=%.3f", uk_m$R2, uk_m$RMSE))
  }, error = function(e) message(paste0("  UK failed: ", e$message)))
  st_write(grid, file.path(CONFIG$dir_output, "kriging_results.gpkg"), delete_dsn=TRUE, quiet=TRUE)
  comp <- bind_rows(OK=cv_m, UK=uk_m, .id="method")
  write_csv(comp, file.path(CONFIG$dir_output, "kriging_comparison.csv"))
  print(comp)
  message("\nKriging complete")
}

message("\n=== 03_analysis.R loaded ===")
message("Run: run_eda()  run_supervised()  run_unsupervised()  run_kriging()")
