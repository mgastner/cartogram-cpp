# Sources

## bahamas_by_district_since_1999.geojson
Global Administrative Areas (2012). GADM database of Global Administrative Areas, version 4.0. URL: www.gadm.org.
Downloaded from: https://geodata.ucdavis.edu/gadm/gadm4.1/shp/gadm41_BHS_shp.zip on 28 July 2022.
The shapefile is converted into a .geojson file using the code below.

### Code for conversion
```
library(rmapshaper)
library(geojsonio)
library(sf)
library(mapview)
library(tidyverse)
temp <- tempfile()
download.file(
  "https://geodata.ucdavis.edu/gadm/gadm4.1/shp/gadm41_BHS_shp.zip",
  temp
)
bhs_unzip <- unzip(temp)
bhs <- read_sf("gadm41_BHS_1.shp")
geojson_write(
  bhs,
  lat = NULL,
  lon = NULL,
  geometry = "point",
  group = NULL,
  file = "bahamas_by_district_since_1999.geojson",
  overwrite = TRUE,
  precision = NULL,
  convert_wgs84 = FALSE,
  crs = NULL
)
unlink(bhs_unzip)
```

## bahamas_population_2010.csv
Source: City Population. Downloaded 24 July 2022 from https://www.citypopulation.de/Bahamas.html. The data is spread across the various pages under the ARCHIPELAGOS: Districts, Settlements & Populated Islands segement of the page.



