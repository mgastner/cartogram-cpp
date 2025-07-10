# Sources

## NUTS_by_region_until_2020.geojson
Eurostat. Accessed via GISCO Data Distribution.
Downloaded from: https://gisco-services.ec.europa.eu/distribution/v2/nuts/download/ref-nuts-2021-03m.geojson.zip on 29 July 2022.

### Code for simplification
```
library(rmapshaper)
library(geojsonio)
library(sf)
library(mapview)
library(tidyverse)
temp <- tempfile()
download.file(
  "https://gisco-services.ec.europa.eu/distribution/v2/nuts/download/ref-nuts-2021-03m.geojson.zip",
  temp
)
nuts_unzip <- unzip(temp)
nuts <- geojson_sf("NUTS_RG_03M_2021_4326_LEVL_2.geojson")
target_n_pts_in_output <- 47000
nuts_simp <- ms_simplify(nuts, keep = target_n_pts_in_output / npts(nuts),keep_shapes = T)
npts(nuts_simp)
geojson_write(
  nuts_simp,
  lat = NULL,
  lon = NULL,
  geometry = "point",
  group = NULL,
  file = "NUTS_by_region_until_2020.geojson",
  overwrite = TRUE,
  precision = NULL,
  convert_wgs84 = FALSE,
  crs = NULL
)
unlink(nuts_unzip)
  ```
  
## NUTS_population_2021.csv
Source: Eurostat. Population on 1 January by age, sex and NUTS 2 region.
Downloaded on 29 July 2022 from https://appsso.eurostat.ec.europa.eu/nui/show.do?dataset=demo_pjan&lang=en.


