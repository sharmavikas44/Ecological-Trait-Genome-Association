**#########################################################**

**#######// Canopy cover/ forest cover #################**

**#########################################################**

// Load your points (coordinates with population IDs)

var points = ee.FeatureCollection("projects/betula123/assets/Individual\_cordinates");



// Load the latest Hansen Global Forest Change dataset

var gfc = ee.Image("UMD/hansen/global\_forest\_change\_2024\_v1\_12");



// Extract the 2000 tree canopy cover percentage (0–100)

var treeCover2000 = gfc.select('treecover2000');



// Calculate current forest cover, accounting for loss (2001–2024)

var lossyear = gfc.select('lossyear');

var gain = gfc.select('gain');



// Mask areas where forest has been lost

var lossMask = lossyear.gt(0);

var cover2024 = treeCover2000.where(lossMask, 0).add(gain);



// Extract % forest cover for each point

var coverPerPoint = cover2024.reduceRegions({

  collection: points,

  reducer: ee.Reducer.mean(),

  scale: 30

});



// Display results and export

print(coverPerPoint);

Export.table.toDrive({

  collection: coverPerPoint,

  description: 'ForestCover\_PerPoint\_2024',

  fileFormat: 'CSV'

});







**###################################################**

**########// Canopy height #######################**

**###################################################**

// Preferred Lang et al. canopy height asset (10 m)

var langAsset = "users/nlang/ETH\_GlobalCanopyHeight\_2020\_10m\_v1";

// Fallback (example) - replace with any other accessible canopy-height asset you have

var fallbackAsset = "NASA/GSFC/Global\_Forest\_Canopy\_Height\_2019";

// output filenames

var exportPointsName = "CanopyHeight\_PerPoint";

var exportPopName = "CanopyHeight\_PerPopulation";

// scale to sample at (use dataset native resolution: 10 for Lang, 30 for some globals)

var sampleScale = 10;

// ----------------------------------------------------------------



// load preferred asset; if not available, switch to fallback

var gch;

try {

  gch = ee.Image(langAsset);

  // test existence by reading band names (will throw if asset not available)

  var bands = gch.bandNames();

  print('Using Lang et al. canopy height asset:', langAsset);

} catch (e) {

  print('Preferred asset not available, attempting fallback:', fallbackAsset);

  gch = ee.Image(fallbackAsset);

  print('Using fallback canopy height asset:', fallbackAsset);

}



// Optionally inspect band names and choose the appropriate band if needed:

print('Canopy height bands:', gch.bandNames());



// If the dataset band name is not "canopy\_height" or similar, set band accordingly.

// Example adjustments (uncomment and set the correct band name if needed):

// gch = gch.select('canopy\_height');        // for NASA global product

// gch = gch.select('height');               // if the Lang asset uses 'height'

// gch = gch.select('band\_name\_here');



// Sample canopy height at each point (mean of pixel(s) touching point)

var heightPerPoint = gch.reduceRegions({

  collection: points,

  reducer: ee.Reducer.mean(),

  scale: sampleScale,

  tileScale: 4

});



// Add lon/lat coordinates and keep population id (assumes 'population' property exists)

var heightWithCoords = heightPerPoint.map(function(f) {

  var geom = f.geometry();

  var coords = geom.coordinates();

  return f.set({

    longitude: coords.get(0),

    latitude: coords.get(1),

    // If your point property for population is named differently (e.g. 'pop'), copy it:

    population: f.get('population')

  });

});



// Replace missing values (null) with a sentinel, e.g., -9999

var heightClean = heightWithCoords.map(function(f) {

  var h = ee.Number(f.get('mean'));

  // If mean is null, set to -9999 (or keep null)

  var h2 = ee.Algorithms.IsEqual(h, null).evaluate ? ee.Algorithms.If(ee.Algorithms.IsEqual(h, null), -9999, h) : h;

  // safer: set only if null using despatch on server side

  h2 = ee.Algorithms.If(ee.Algorithms.IsEqual(f.get('mean'), null), -9999, f.get('mean'));

  return f.set('canopy\_height\_m', h2);

});



// Print sample rows

print('Point-level canopy height (first 20):', heightClean.limit(20));



// Aggregate to population-wise mean (ignores sentinel -9999 values)

var populations = heightClean.aggregate\_array('population').distinct();

var popMeanFC = ee.FeatureCollection(

  populations.map(function(pop) {

    pop = ee.String(pop);

    var sub = heightClean.filter(ee.Filter.eq('population', pop));

    // convert canopy\_height\_m to numeric array and filter out sentinel -9999

    var valid = sub.filter(ee.Filter.neq('canopy\_height\_m', -9999));

    var meanCover = valid.aggregate\_mean('canopy\_height\_m');

    var count = valid.size();

    return ee.Feature(null, {

      population: pop,

      mean\_canopy\_height\_m: meanCover,

      n\_points: count

    });

  })

);

print('Population-wise canopy height:', popMeanFC);



// Export point-level table to Drive (CSV)

Export.table.toDrive({

  collection: heightClean.select(\['population','longitude','latitude','canopy\_height\_m']),

  description: exportPointsName,

  fileFormat: 'CSV'

});



// Export population-wise table to Drive (CSV)

Export.table.toDrive({

  collection: popMeanFC,

  description: exportPopName,

  fileFormat: 'CSV'

});



// Map for visualisation (optional)

Map.centerObject(points, 7);

Map.addLayer(gch, {min:0, max:40, palette:\['lightblue','yellow','green','darkgreen']}, 'Canopy height (raster)');

Map.addLayer(points, {color: 'red'}, 'Points');









**##################################################################**

**########### VV:VH ratio based Ecological traits #################**

**##################################################################**

// ----------------------

// Load your population points (must have property "popID")

// ----------------------

var points = ee.FeatureCollection("projects/betula123/assets/Individual\_cordinates");



// ----------------------

// Sentinel-1 collection

// ----------------------

var s1 = ee.ImageCollection("COPERNICUS/S1\_GRD")

  .filter(ee.Filter.eq('instrumentMode','IW'))

  .filter(ee.Filter.listContains('transmitterReceiverPolarisation','VV'))

  .filter(ee.Filter.listContains('transmitterReceiverPolarisation','VH'))

  .filterDate('2015-01-01','2024-12-31')

  .filter(ee.Filter.eq('orbitProperties\_pass','ASCENDING')); // optional filter



// ----------------------

// Function to compute VV/VH ratio^2

// ----------------------

var addRatio = function(img) {

  var vv = img.select('VV');

  var vh = img.select('VH');

  var ratio = vv.divide(vh).pow(2).rename('vv\_vh\_ratio2');

  return img.addBands(ratio).select(\['vv\_vh\_ratio2'])

            .copyProperties(img, \['system:time\_start']);

};



var s1\_ratios = s1.map(addRatio);



// ----------------------

// Extract time series for each point

// ----------------------

var bufferRadius = 200; // meters



var extract\_ts = function(feature) {

  var geom = feature.geometry().buffer(bufferRadius);

 

  var perImage = s1\_ratios.map(function(img) {

    var mean\_val = img.reduceRegion({

      reducer: ee.Reducer.mean(),

      geometry: geom,

      scale: 10,

      maxPixels: 1e8

    }).get('vv\_vh\_ratio2');

 

    return ee.Feature(null, {

      'date': ee.Date(img.get('system:time\_start')).format('YYYY-MM-dd'),

      'value': mean\_val,

      'popID': feature.get('popID')

    });

  });

 

  return ee.FeatureCollection(perImage).filter(ee.Filter.notNull(\['value']));

};



// Map over points and flatten results

var ts = points.map(extract\_ts).flatten();



print('Extracted time series (sample):', ts.limit(10));



// ----------------------

// Export results to Drive

// ----------------------

Export.table.toDrive({

  collection: ts,

  description: 'S1\_VV\_VH\_Ratio2\_TimeSeries',

  fileFormat: 'CSV'

});









//-----------------------------------

// New script for Ecological traits

// ----------------------------------

// ----------------------------------------------------------

// INPUT: Your point feature collection (must contain popID)

// ----------------------------------------------------------

var points = ee.FeatureCollection("projects/betula123/assets/Pop\_cordinates");



// Optional: apply buffer around each point (in meters)

var bufferRadius = 200; // set to 0 if you want EXACT point only

var bufferedPoints = points.map(function(f){

  return f.buffer(bufferRadius);

});



// ----------------------------------------------------------

// SENTINEL-1 COLLECTION \& VV/VH^2 RATIO

// ----------------------------------------------------------

var s1 = ee.ImageCollection("COPERNICUS/S1\_GRD")

  .filter(ee.Filter.eq('instrumentMode','IW'))

  .filter(ee.Filter.listContains('transmitterReceiverPolarisation','VV'))

  .filter(ee.Filter.listContains('transmitterReceiverPolarisation','VH'))

  .filterDate('2015-01-01','2024-12-31')

  .filter(ee.Filter.eq('orbitProperties\_pass','ASCENDING')); // optional



var addRatio = function(img) {

  var ratio = img.select('VV').divide(img.select('VH')).pow(2).rename('vv\_vh\_ratio2');

  return img.addBands(ratio).select('vv\_vh\_ratio2');

};



var s1\_ratio = s1.map(addRatio);



// ----------------------------------------------------------

// COMPUTE TRAIT METRICS ACROSS ENTIRE TIME SERIES

// ----------------------------------------------------------

var mean\_img      = s1\_ratio.mean().rename('mean\_ratio');

var stddev\_img    = s1\_ratio.reduce(ee.Reducer.stdDev()).rename('stddev\_ratio');

var min\_img       = s1\_ratio.min().rename('min\_ratio');

var max\_img       = s1\_ratio.max().rename('max\_ratio');

var amplitude\_img = max\_img.subtract(min\_img).rename('amplitude\_ratio');



// Combine into a single multi-band trait image

var traits = mean\_img.addBands(\[stddev\_img, min\_img, max\_img, amplitude\_img]);



// ----------------------------------------------------------

// EXTRACT TRAITS PER POINT BUFFER (CORRECT FUNCTION)

// ----------------------------------------------------------

var traitTable = traits.reduceRegions({

  collection: bufferedPoints,

  reducer: ee.Reducer.mean(),   // average inside buffer

  scale: 10

});



// View in Console

print("Trait Table:", traitTable.limit(10));



// ----------------------------------------------------------

// EXPORT TO GOOGLE DRIVE (CSV)

// ----------------------------------------------------------

Export.table.toDrive({

  collection: traitTable,

  description: "S1\_VV\_VH2\_Structural\_Canopy\_Traits",

  fileFormat: "CSV"

});







**############################################################**

**################## Leaf Area Index  ########################**

**############################################################**

// 1. Load your points

var pts = ee.FeatureCollection('projects/betula123/assets/Individual\_cordinates');

print('First feature:', pts.first());

print('Total points:', pts.size());



Map.centerObject(pts, 7);

Map.addLayer(pts, {}, 'Points');



// 2. Load MODIS LAI (MOD15A2H v6.1) and select correct band

var laiAll = ee.ImageCollection('MODIS/061/MOD15A2H')

&nbsp; .select('Lai\_500m');  // IMPORTANT: Lai\_500m



print('LAI collection size (all):', laiAll.size());

print('Bands of first image:', ee.Image(laiAll.first()).bandNames());



// 3. Filter by date

var laiCol = laiAll.filterDate('2003-01-01', '2023-12-31');

print('Filtered LAI size:', laiCol.size());



// 4. Mean LAI over time

var laiMean = laiCol.mean()

&nbsp; .multiply(0.1)           // scale factor

&nbsp; .rename('LAI\_mean');



Map.addLayer(laiMean, {min: 0, max: 7}, 'Mean LAI');



// 5. Extract LAI at each individual point

var ptsWithLai = laiMean.sampleRegions({

&nbsp; collection: pts,

&nbsp; properties: \['id', 'popID'],  // must match your table

&nbsp; scale: 500,

&nbsp; tileScale: 2

});



print('Points with LAI:', ptsWithLai.limit(5));

print('Number of sampled points:', ptsWithLai.size());



// 6. Pop-wise mean LAI

var grouped = ptsWithLai.reduceColumns({

&nbsp; selectors: \['popID', 'LAI\_mean'],

&nbsp; reducer: ee.Reducer.mean().group({

&nbsp;   groupField: 0,

&nbsp;   groupName: 'popID'

&nbsp; })

});



print('Grouped dictionary:', grouped);



var groupList = ee.List(grouped.get('groups'));



var laiByPop = ee.FeatureCollection(groupList.map(function(g) {

&nbsp; g = ee.Dictionary(g);

&nbsp; return ee.Feature(null, {

&nbsp;   popID: g.get('popID'),

&nbsp;   LAI\_mean: g.get('mean')

&nbsp; });

}));



print('Pop-wise LAI:', laiByPop);



// 7a. Export POP-wise LAI

Export.table.toDrive({

&nbsp; collection: laiByPop,

&nbsp; description: 'BU\_popwise\_LAI\_MOD15A2H\_2003\_2023',

&nbsp; fileFormat: 'CSV'

});



// 7b. Export POINT-wise LAI (no coords; join by id later)

Export.table.toDrive({

&nbsp; collection: ptsWithLai,

&nbsp; description: 'BU\_pointwise\_LAI\_MOD15A2H\_2003\_2023',

&nbsp; fileFormat: 'CSV'

});





**##################################################################**

**##################   Chlorophyll content  ########################**

**##################################################################**

// ===============================

// 1. Load your population points

// ===============================

var pts = ee.FeatureCollection('projects/betula123/assets/Individual\_cordinates');

print('First feature:', pts.first());

print('Total points:', pts.size());



Map.centerObject(pts, 7);

Map.addLayer(pts, {}, 'Points');





// =====================================

// 2. Cloud masking for Sentinel-2 SR

// =====================================

// Using COPERNICUS/S2\_SR\_HARMONIZED

// Mask clouds \& cloud shadows using SCL band

function maskS2clouds(image) {

  var scl = image.select('SCL');

  // Keep only "good" classes (vegetation, bare, water, etc.)

  var good = scl.eq(4)   // vegetation

    .or(scl.eq(5))       // not vegetated

    .or(scl.eq(6))       // water

    .or(scl.eq(7))       // unclassified

    .or(scl.eq(11));     // snow/ice (optional; remove if not needed)

 

  return image.updateMask(good);

}





// =====================================

// 3. Load Sentinel-2 and compute indices

// =====================================

// Adjust dates to your field sampling period / growing season

var startDate = '2017-06-01';

var endDate   = '2023-09-30';



// Area of interest around your points

var aoi = pts.geometry().buffer(5000);



var s2 = ee.ImageCollection('COPERNICUS/S2\_SR\_HARMONIZED')

  .filterDate(startDate, endDate)

  .filterBounds(aoi)

  .filter(ee.Filter.lte('CLOUDY\_PIXEL\_PERCENTAGE', 40))

  .map(maskS2clouds);



print('S2 images after filtering:', s2.size());



// Function to add chlorophyll indices

var addChlIndices = function(img) {

  // Sentinel-2 bands:

  // B4 = red (665 nm)

  // B5 = red edge 1 (705 nm)

  // B6 = red edge 2 (740 nm)

  // B8A = narrow NIR (865 nm)

 

  var red   = img.select('B4');

  var re1   = img.select('B5');

  var re2   = img.select('B6');

  var nirN  = img.select('B8A');

 

  // 1) Chlorophyll Index red-edge: CIre = (NIR / RE1) - 1

  var CIre = nirN.divide(re1).subtract(1).rename('CIre');

 

  // 2) MERIS Terrestrial Chlorophyll Index (MTCI-like):

  //    MTCI ≈ (RE2 - RE1) / (RE1 - Red)

  var mtci = re2.subtract(re1)

                .divide(re1.subtract(red))

                .rename('MTCI');

 

  return img.addBands(\[CIre, mtci]);

};



var s2Chl = s2.map(addChlIndices);



// Create a median composite over the period

var chlComposite = s2Chl.median().select(\['CIre', 'MTCI']);

print('Chlorophyll composite:', chlComposite);



Map.addLayer(chlComposite.select('CIre'), {min:0, max:3}, 'CIre (median)');

Map.addLayer(chlComposite.select('MTCI'), {min:-1, max:5}, 'MTCI (median)');





// =====================================

// 4. Extract chlorophyll indices at points

// =====================================

var ptsWithChl = chlComposite.sampleRegions({

  collection: pts,

  properties: \['id', 'popID'],

  scale: 20,       // Sentinel-2 bands at 20 m (B5, B6, B8A)

  tileScale: 2

});



print('Points with chlorophyll indices:', ptsWithChl.limit(5));

print('Total sampled points:', ptsWithChl.size());





// =====================================

// 5. Pop-wise mean chlorophyll indices

// =====================================

var grouped = ptsWithChl.reduceColumns({

  selectors: \['popID', 'CIre', 'MTCI'],

  reducer: ee.Reducer.mean().repeat(2).group({

    groupField: 0,    // index of 'popID'

    groupName: 'popID'

  })

});



print('Grouped dictionary:', grouped);



var groupsList = ee.List(grouped.get('groups'));



var chlByPop = ee.FeatureCollection(groupsList.map(function(g) {

  g = ee.Dictionary(g);

  var popID = g.get('popID');

  var means = ee.List(g.get('mean'));  // \[mean\_CIre, mean\_MTCI]

 

  return ee.Feature(null, {

    'popID': popID,

    'CIre\_mean': means.get(0),

    'MTCI\_mean': means.get(1)

  });

}));



print('Pop-wise chlorophyll indices:', chlByPop);





// =====================================

// 6. Export pop-wise chlorophyll to CSV

// =====================================

Export.table.toDrive({

  collection: chlByPop,

  description: 'BU\_popwise\_Chlorophyll\_S2\_2017\_2023',

  fileFormat: 'CSV'

});



// (Optional) Export point-wise chlorophyll too:

Export.table.toDrive({

  collection: ptsWithChl,

  description: 'BU\_pointwise\_Chlorophyll\_S2\_2017\_2023',

  fileFormat: 'CSV'

});





**################################################################**

**############### NDRE-based Nitrogen content ####################**

**################################################################**

// 1. Load your population points

var pts = ee.FeatureCollection('projects/betula123/assets/Individual\_cordinates');

print('First feature:', pts.first());

print('Total points:', pts.size());



Map.centerObject(pts, 7);

Map.addLayer(pts, {}, 'Points');





// 2. Cloud mask for Sentinel-2 SR

function maskS2clouds(image) {

  var scl = image.select('SCL');

  var good = scl.eq(4)   // vegetation

    .or(scl.eq(5))       // not vegetated

    .or(scl.eq(6))       // water

    .or(scl.eq(7))       // unclassified

    .or(scl.eq(11));     // snow/ice (optional)

  return image.updateMask(good);

}





// 3. Load Sentinel-2 and compute NDRE (nitrogen proxy)



// Adjust to your preferred period (e.g. sampling years, growing season)

var startDate = '2017-06-01';

var endDate   = '2023-09-30';



var aoi = pts.geometry().buffer(5000);



var s2 = ee.ImageCollection('COPERNICUS/S2\_SR\_HARMONIZED')

  .filterDate(startDate, endDate)

  .filterBounds(aoi)

  .filter(ee.Filter.lte('CLOUDY\_PIXEL\_PERCENTAGE', 40))

  .map(maskS2clouds);



print('S2 images after filtering:', s2.size());



// Add NDRE band: NDRE = (B8A - B5) / (B8A + B5)

var addNDRE = function(img) {

  var nirN = img.select('B8A'); // narrow NIR

  var re1  = img.select('B5');  // red-edge 1

  var ndre = nirN.subtract(re1).divide(nirN.add(re1)).rename('NDRE');

  return img.addBands(ndre);

};



var s2Ndre = s2.map(addNDRE);



// Median NDRE composite over period

var ndreComposite = s2Ndre.median().select('NDRE');

print('NDRE composite:', ndreComposite);



Map.addLayer(ndreComposite, {min: -0.2, max: 0.6}, 'NDRE (median)');





// 4. Extract NDRE (nitrogen proxy) at each point

var ptsWithN = ndreComposite.sampleRegions({

  collection: pts,

  properties: \['id', 'popID'],

  scale: 20,       // B5, B8A are 20 m

  tileScale: 2

});



print('Points with NDRE (N proxy):', ptsWithN.limit(5));

print('Total sampled points:', ptsWithN.size());





// 5. Pop-wise mean nitrogen content (NDRE\_mean)

var grouped = ptsWithN.reduceColumns({

  selectors: \['popID', 'NDRE'],

  reducer: ee.Reducer.mean().group({

    groupField: 0,   // index of 'popID'

    groupName: 'popID'

  })

});



print('Grouped dictionary:', grouped);



var groupsList = ee.List(grouped.get('groups'));



var N\_byPop = ee.FeatureCollection(groupsList.map(function(g) {

  g = ee.Dictionary(g);

  return ee.Feature(null, {

    'popID':      g.get('popID'),

    'NDRE\_mean':  g.get('mean')   // your nitrogen content proxy

  });

}));



print('Pop-wise nitrogen content (proxy):', N\_byPop);





// 6. Export pop-wise nitrogen proxy to CSV

Export.table.toDrive({

  collection: N\_byPop,

  description: 'BU\_popwise\_Nitrogen\_NDRE\_2017\_2023',

  fileFormat: 'CSV'

});



// (Optional) export point-wise NDRE too

Export.table.toDrive({

  collection: ptsWithN,

  description: 'BU\_pointwise\_Nitrogen\_NDRE\_2017\_2023',

  fileFormat: 'CSV'

});

