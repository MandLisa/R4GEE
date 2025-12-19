library(terra)
library(data.table)

# ------------------------------------------------------------
# 0) Paths
# ------------------------------------------------------------
ndvi_dir <- "/mnt/eo/Capra_ibex/_data/MODIS/mask/"   # directory containing yearly NDVI stacks
out_csv  <- "/mnt/eo/Capra_ibex/_phenology_output/ndvi_median_2006_2024_full_NP.csv"

# Optional: if you want to restrict to a study area polygon, set this; otherwise leave NULL
# aoi_shp <- "/mnt/eo/Capra_ibex/_data/npbgd_border/npbgd_border.shp"
aoi_shp <- NULL

# Optional: coded missing value (set to NA). If not needed, set to NA_real_ and it is ignored
missing_code <- NA_real_   # e.g., -10000

# ------------------------------------------------------------
# 1) List yearly stack files
# ------------------------------------------------------------
files <- sort(list.files(ndvi_dir, pattern = "\\.tif$", full.names = TRUE))
stopifnot(length(files) > 0)

# If AOI requested, load once
if (!is.null(aoi_shp)) {
  aoi <- vect(aoi_shp)
} else {
  aoi <- NULL
}

# ------------------------------------------------------------
# 2) Helper: parse dates from layer names robustly
#    Expected: contains YYYYMMDD, e.g. NDVI_20060610
# ------------------------------------------------------------
parse_dates_from_names <- function(nm) {
  # extract first 8-digit block that looks like YYYYMMDD
  d8 <- sub(".*?((19|20)\\d{2}(0[1-9]|1[0-2])([0-2]\\d|3[01])).*", "\\1", nm)
  # If no match, d8 will equal original string -> set to NA
  d8[!grepl("^(19|20)\\d{6}$", d8)] <- NA_character_
  as.Date(d8, format = "%Y%m%d")
}

# ------------------------------------------------------------
# 3) Loop over files: for each year in file, select closest-to-15Jun layer and compute global median
# ------------------------------------------------------------
res <- list()
k <- 0L

for (f in files) {
  message("\n--- Processing file: ", f)
  
  r <- rast(f)
  
  # Optional missing code to NA
  if (!is.na(missing_code)) {
    r[r == missing_code] <- NA
  }
  
  # Optional AOI crop/mask
  if (!is.null(aoi)) {
    if (!same.crs(aoi, r)) aoi_use <- project(aoi, crs(r)) else aoi_use <- aoi
    r <- mask(crop(r, aoi_use), aoi_use)
  }
  
  nm <- names(r)
  dates <- parse_dates_from_names(nm)
  
  if (all(is.na(dates))) {
    message("  !! Could not parse ANY dates from layer names in this file. Example names:")
    print(head(nm, 20))
    stop("Fix naming or adapt date parsing.")
  }
  
  yrs <- as.integer(format(dates, "%Y"))
  uyrs <- sort(unique(yrs[!is.na(yrs)]))
  
  # DEBUG: show the years found in this file
  message("  Years found in file: ", paste(uyrs, collapse = ", "))
  
  for (yy in uyrs) {
    target <- as.Date(sprintf("%d-06-15", yy))
    
    idx_y <- which(yrs == yy)
    if (length(idx_y) == 0) next
    
    dists <- abs(as.integer(dates[idx_y] - target))
    min_dist <- min(dists, na.rm = TRUE)
    tied <- idx_y[dists == min_dist]
    
    # tie-break: earlier date
    pick <- tied[which.min(dates[tied])]
    
    pick_date <- dates[pick]
    pick_name <- nm[pick]
    
    # Compute global median over all pixels in the selected layer
    med_df <- terra::global(r[[pick]], fun = median, na.rm = TRUE)
    med_val <- as.numeric(med_df[1, 1])
    
    # DEBUG: show selection and median
    message(sprintf("  Year %d -> selected %s (%s), global median = %.6f",
                    yy, as.character(pick_date), pick_name, med_val))
    
    k <- k + 1L
    res[[k]] <- data.table(
      year = yy,
      selected_date = pick_date,
      layer_name = pick_name,
      NDVI_median = med_val,
      source_file = basename(f)
    )
  }
}

dt_out <- rbindlist(res, use.names = TRUE, fill = TRUE)
setorder(dt_out, year)

# ------------------------------------------------------------
# 4) Sanity checks (important!)
# ------------------------------------------------------------
# (a) Are we really selecting different dates across years?
print(dt_out[, .(n_years = .N,
                 n_unique_dates = uniqueN(selected_date),
                 n_unique_layers = uniqueN(layer_name))])

# (b) Do medians vary at all?
print(dt_out[, .(min = min(NDVI_median, na.rm=TRUE),
                 max = max(NDVI_median, na.rm=TRUE),
                 sd  = sd(NDVI_median,  na.rm=TRUE))])

# If everything is identical, this will show it explicitly:
if (dt_out[, uniqueN(NDVI_median)] == 1) {
  warning("All NDVI_median values are identical across years. Check debug output above: are the selected layers/dates identical?")
}

dt_out[, c(3, 5) := NULL]
# ------------------------------------------------------------
# 5) Save + quick plot
# ------------------------------------------------------------
dir.create(dirname(out_csv), showWarnings = FALSE, recursive = TRUE)
fwrite(dt_out[, .(year, NDVI_median, selected_date)], out_csv)

cat("\nWrote:", out_csv, "\n")
print(dt_out[, .(year, selected_date, NDVI_median)][1:30])

ggplot(dt_out, aes(x = year, y = NDVI_median)) +
  geom_line() +
  geom_point() +
  scale_x_continuous(breaks = sort(unique(dt_out$year))) +
  labs(
    x = "Year",
    y = "Median NDVI (closest to 15 June)",
    title = ""
  ) +
  theme_bw()
