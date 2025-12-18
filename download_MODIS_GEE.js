// =======================================================
// MOD09Q1 NDVI annual stacks for NPBGD (EPSG:25832)
// - NO QA masking (use all available pixels)
// - NDVI clamped to [-1, 1] (float)
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
var outFolder = "MOD09Q1_NDVI_STACKS_NPBGD_EPSG25832_NOQA";
var startDate = "2006-01-01";
var endDate   = "2024-12-31";

// Apr–Oct (set to null to disable)
var monthStart = 4;
var monthEnd   = 10;

// CRS must match your shapefile:
var exportCrs = "EPSG:25832";

// Pixel size in meters (in EPSG:25832)
var scaleM = 250;

var maxPix = 1e13;
var noDataVal = -9999;

// Resampling during reprojection/export (NDVI is continuous)
var resamplingMethod = "bilinear"; // or "nearest"

// Optional sanity mask to avoid NDVI artefacts from invalid reflectance.
// If you truly want EVERYTHING, set to false.
var useReflectanceSanityMask = true;

// ---------------- 3) LOAD MODIS -------------------------
var col = ee.ImageCollection("MODIS/061/MOD09Q1")
  .filterBounds(aoiGeom)
  .filterDate(startDate, endDate);

if (monthStart !== null && monthEnd !== null) {
  col = col.filter(ee.Filter.calendarRange(monthStart, monthEnd, "month"));
}

// ---------------- 4) NDVI (NO QA) ------------------------
function toNDVI(img) {
  var red = img.select("sur_refl_b01").multiply(0.0001);
  var nir = img.select("sur_refl_b02").multiply(0.0001);
  var denom = nir.add(red);

  var ndvi = nir.subtract(red).divide(denom).rename("NDVI");

  if (useReflectanceSanityMask) {
    var m =
      red.gte(0).and(nir.gte(0))
        .and(red.lte(1)).and(nir.lte(1))
        .and(denom.neq(0));
    ndvi = ndvi.updateMask(m);
  }

  ndvi = ndvi.clamp(-1, 1).toFloat();

  var year = ee.Date(img.get("system:time_start")).get("year");
  return ndvi.copyProperties(img, ["system:time_start"]).set("year", year);
}

var ndviCol = col.map(toNDVI);

// ---------------- 5) COUNTS PER YEAR ---------------------
var perYear = ee.Dictionary(ndviCol.aggregate_histogram("year"));
print("NDVI rasters per year (counts):", perYear);

// ---------------- 6) EXPORT STACKS PER YEAR --------------
var years = ee.List(ndviCol.aggregate_array("year")).distinct().sort();
print("Years to export:", years);

years.getInfo().forEach(function(y) {
  y = parseInt(y, 10);

  var yearCol = ndviCol.filter(ee.Filter.eq("year", y)).sort("system:time_start");

  // band names NDVI_YYYYMMDD
  var times = ee.List(yearCol.aggregate_array("system:time_start"));
  var bandNames = times.map(function(t) {
    return ee.String("NDVI_").cat(ee.Date(t).format("YYYYMMdd"));
  });

  var stack = yearCol.toBands().rename(bandNames);

  // Crop to AOI at export (region), and set reprojection behaviour
  stack = stack.resample(resamplingMethod);

  // Make masked pixels explicit NoData, and tag it in GeoTIFF
  stack = stack.unmask(noDataVal).toFloat();

  Export.image.toDrive({
    image: stack,
    description: "MOD09Q1_NDVI_STACK_" + y + "_EPSG25832",
    folder: outFolder,
    fileNamePrefix: "MOD09Q1_NDVI_STACK_" + y + "_EPSG25832",
    region: aoiGeom,
    crs: exportCrs,
    scale: scaleM,
    maxPixels: maxPix,
    fileFormat: "GeoTIFF",
    formatOptions: { noData: noDataVal }
  });
});

