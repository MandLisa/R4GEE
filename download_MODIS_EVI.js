// =======================================================
// MOD13Q1 EVI annual stacks for NPBGD (EPSG:25832)
// - 16-day composites (MOD13Q1)
// - NO QA masking (use all available pixels)
// - EVI scaled and clamped to [-1, 1] (float)
// - One GeoTIFF per year, multi-band (one band per date)
// - Export CRS = EPSG:25832 (matches your shapefile)
// - Explicit NoData = -9999 to prevent GIS min/max artefacts
// =======================================================

// -------------------- 1) AOI ----------------------------
var aoi = ee.FeatureCollection("projects/ee-lisamandl-chm-alps/assets/npbgd_border");
var aoiGeom = aoi.geometry();

Map.centerObject(aoi, 10);
Map.addLayer(aoi, {color: "yellow"}, "AOI outline");

// ---------------- 2) SETTINGS ----------------------------
var outFolder = "MOD13Q1_EVI_STACKS_NPBGD_EPSG25832_NOQA";
var startDate = "2006-01-01";
var endDate   = "2024-12-31";

// Apr–Oct (set to null to disable)
var monthStart = 4;
var monthEnd   = 10;

// CRS must match your shapefile:
var exportCrs = "EPSG:25832";

// Pixel size in meters (MOD13Q1 nominal is 250 m)
var scaleM = 250;

var maxPix = 1e13;
var noDataVal = -9999;

// Resampling during reprojection/export (EVI is continuous)
var resamplingMethod = "bilinear"; // or "nearest"

// Optional sanity mask to avoid fill-value artefacts.
// This is NOT QA masking; it only removes obvious invalid fill values.
// If you truly want EVERYTHING, set to false.
var useEviSanityMask = true;

// ---------------- 3) LOAD MODIS (16-day VI) --------------
var col = ee.ImageCollection("MODIS/061/MOD13Q1")
  .filterBounds(aoiGeom)
  .filterDate(startDate, endDate);

if (monthStart !== null && monthEnd !== null) {
  col = col.filter(ee.Filter.calendarRange(monthStart, monthEnd, "month"));
}

// ---------------- 4) EVI (NO QA) -------------------------
function toEVI(img) {
  // MOD13Q1 EVI is stored as int16 with scale factor 0.0001
  var eviRaw = img.select("EVI");
  var evi = eviRaw.multiply(0.0001).rename("EVI");

  if (useEviSanityMask) {
    // Common MODIS VI fill value is -3000 (before scaling)
    // Keep everything except fill and extreme out-of-range artefacts.
    var m = eviRaw.neq(-3000);
    evi = evi.updateMask(m);
  }

  evi = evi.clamp(-1, 1).toFloat();

  var year = ee.Date(img.get("system:time_start")).get("year");
  return evi.copyProperties(img, ["system:time_start"]).set("year", year);
}

var eviCol = col.map(toEVI);

// ---------------- 5) COUNTS PER YEAR ---------------------
var perYear = ee.Dictionary(eviCol.aggregate_histogram("year"));
print("EVI rasters per year (counts):", perYear);

// ---------------- 6) EXPORT STACKS PER YEAR --------------
var years = ee.List(eviCol.aggregate_array("year")).distinct().sort();
print("Years to export:", years);

years.getInfo().forEach(function(y) {
  y = parseInt(y, 10);

  var yearCol = eviCol.filter(ee.Filter.eq("year", y)).sort("system:time_start");

  // band names EVI_YYYYMMDD (based on system:time_start of the composite)
  var times = ee.List(yearCol.aggregate_array("system:time_start"));
  var bandNames = times.map(function(t) {
    return ee.String("EVI_").cat(ee.Date(t).format("YYYYMMdd"));
  });

  var stack = yearCol.toBands().rename(bandNames);

  // Set resampling behaviour (continuous variable)
  stack = stack.resample(resamplingMethod);

  // Make masked pixels explicit NoData, and tag it in GeoTIFF
  stack = stack.unmask(noDataVal).toFloat();

  Export.image.toDrive({
    image: stack,
    description: "MOD13Q1_EVI_STACK_" + y + "_EPSG25832",
    folder: outFolder,
    fileNamePrefix: "MOD13Q1_EVI_STACK_" + y + "_EPSG25832",
    region: aoiGeom,
    crs: exportCrs,
    scale: scaleM,
    maxPixels: maxPix,
    fileFormat: "GeoTIFF",
    formatOptions: { noData: noDataVal }
  });
});
