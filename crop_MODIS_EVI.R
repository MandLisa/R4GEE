# =========================================================
# Robust reprojection of MODIS EVI stacks to EPSG:25832
# using gdalwarp with proper quoting + source SRS fallbacks.
# =========================================================

library(terra)

# -------------------- 1) INPUTS --------------------------
mask_gpkg  <- "/mnt/eo/Capra_ibex/_data/CORINE/CORINE_mask_minus_elev.gpkg"
mask_layer <- "corine_minus_elev"

ras_dir <- "/mnt/eo/Capra_ibex/_data/MODIS/EVI"
tif_files <- list.files(ras_dir, pattern = "\\.tif$", full.names = TRUE, recursive = TRUE)
stopifnot(length(tif_files) > 0)
cat("Found", length(tif_files), "rasters.\n")

# Output folders
out_warp_dir <- file.path(ras_dir, "EPSG25832_warped")
dir.create(out_warp_dir, recursive = TRUE, showWarnings = FALSE)

out_mask_dir <- file.path(ras_dir, "EPSG25832_warped_masked")
dir.create(out_mask_dir, recursive = TRUE, showWarnings = FALSE)

# Target settings
target_crs <- "EPSG:25832"
target_res <- 250
noDataVal  <- -9999

# -------------------- 2) CHECK GDAL ----------------------
gdalwarp <- Sys.which("gdalwarp")
if (!nzchar(gdalwarp)) stop("gdalwarp not found in PATH.")
cat("gdalwarp:", gdalwarp, "\n")

# -------------------- 3) LOAD MASK + TARGET EXTENT -------
mask_v <- vect(mask_gpkg, layer = mask_layer)
if (is.null(mask_v)) stop("Could not read mask layer. Check mask_gpkg and mask_layer.")
if (is.na(crs(mask_v))) stop("Mask has NA CRS. Fix CRS metadata first.")

mask_utm <- project(mask_v, target_crs)
e <- ext(mask_utm)
te <- c(e$xmin, e$ymin, e$xmax, e$ymax)
cat("Target extent (EPSG:25832):", paste(round(te, 2), collapse = " "), "\n")

# -------------------- 4) SOURCE SRS CANDIDATES -----------
# Decide if the raster *looks* like degrees; if yes, try EPSG:4326 first.
r0 <- rast(tif_files[1])
e0 <- ext(r0)
is_degrees <- (abs(e0$xmin) <= 180 && abs(e0$xmax) <= 180 && abs(e0$ymin) <= 90 && abs(e0$ymax) <= 90)

# Candidate source SRS strings:
# 1) SR-ORG:6974 is frequently used by GDAL for MODIS Sinusoidal (most robust if available).
# 2) Explicit sinusoidal with spherical radius; use +a/+b rather than +R for maximum PROJ compatibility.
# 3) ESRI world sinusoidal (sometimes works when SR-ORG is unsupported).
# 4) EPSG:4326 fallback (only if the raster is actually lon/lat).
src_candidates <- c(
  "SR-ORG:6974",
  "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs",
  "ESRI:54008",
  "EPSG:4326"
)

# Order candidates based on extent heuristic
if (is_degrees) {
  src_candidates <- c("EPSG:4326", setdiff(src_candidates, "EPSG:4326"))
}

cat("Source SRS candidates (in order):\n")
print(src_candidates)

# -------------------- 5) RUN GDALWARP (QUOTED) -----------
run_gdalwarp <- function(infile, outfile, src_srs) {

  args <- c(
    "-overwrite",
    # Start without -multi to reduce complexity while debugging. You can re-enable later.
    # "-multi",
    "-r", "bilinear",
    "-s_srs", src_srs,
    "-t_srs", target_crs,
    "-tr", target_res, target_res,
    "-tap",
    "-te", te[1], te[2], te[3], te[4],
    "-srcnodata", noDataVal,
    "-dstnodata", noDataVal,
    "-co", "COMPRESS=DEFLATE",
    "-co", "PREDICTOR=2",
    "-co", "ZLEVEL=6",
    "-co", "BIGTIFF=IF_SAFER",
    infile, outfile
  )

  # Build a shell command with proper quoting so the PROJ string stays intact
  cmd <- paste(shQuote(gdalwarp), paste(shQuote(as.character(args)), collapse = " "))
  cmd <- paste(cmd, "2>&1")  # capture stderr

  out <- system(cmd, intern = TRUE)
  status <- attr(out, "status")
  if (is.null(status)) status <- 0

  list(status = status, output = out, cmd = cmd)
}

# -------------------- 6) WARP ALL FILES ------------------
for (f in tif_files) {
  bn <- tools::file_path_sans_ext(basename(f))
  out_f <- file.path(out_warp_dir, paste0(bn, "_EPSG25832.tif"))

  message("\n--- Warping: ", f)

  ok <- FALSE
  for (src_srs in src_candidates) {
    message("Trying -s_srs: ", src_srs)

    res <- run_gdalwarp(f, out_f, src_srs)

    if (res$status == 0 && file.exists(out_f)) {
      message("SUCCESS with -s_srs: ", src_srs)
      message("Wrote: ", out_f)
      ok <- TRUE
      break
    } else {
      message("FAILED (status ", res$status, "). GDAL says:\n",
              paste(res$output, collapse = "\n"))
    }
  }

  if (!ok) {
    stop(
      "All source SRS candidates failed for: ", f,
      "\nThis strongly suggests the input GeoTIFF CRS/geotransform is invalid or inconsistent."
    )
  }
}

cat("\nWarping done.\n")

# -------------------- 7) OPTIONAL: MASK/CROP -------------
warped_files <- list.files(out_warp_dir, pattern = "\\.tif$", full.names = TRUE)

for (f in warped_files) {
  bn <- tools::file_path_sans_ext(basename(f))
  out_f <- file.path(out_mask_dir, paste0(bn, "_CORINEmasked.tif"))

  message("\n--- Masking: ", f)
  r <- rast(f)

  r2 <- crop(r, mask_utm)
  r2 <- mask(r2, mask_utm)

  writeRaster(
    r2, out_f, overwrite = TRUE,
    NAflag = noDataVal,
    gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL=6", "BIGTIFF=IF_SAFER")
  )
  message("Wrote: ", out_f)
}

cat("\nAll done.\n")
