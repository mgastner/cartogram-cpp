# Sources

## vietnam_by_province_since_1996.geojson
Global Administrative Areas (2012). GADM database of Global Administrative Areas, version 2.0. [online] URL: www.gadm.org.
Downloaded from https://geodata.ucdavis.edu/gadm/gadm4.0/shp/gadm40_VNM_shp.zip on 4 July 2022.

### Code for simplification
```
library(rmapshaper)
library(geojsonio)
library(sf)
library(mapview)
vnm <- geojson_sf("~/vietnam_by_province_since_1996/gadm40_VNM_1.geojson")
target_n_pts_in_output <- 48500
npts(vnm)
vnm_simp <- ms_simplify(vnm, keep = target_n_pts_in_output/npts(vnm))
geojson_write(
  vnm_simp,
  lat = NULL,
  lon = NULL,
  geometry = "point",
  group = NULL,
  file = "vietnam_by_province_since_1996.geojson",
  overwrite = TRUE,
  precision = NULL,
  convert_wgs84 = FALSE,
  crs = NULL)
  ```
  
## vietnam_population_2019.csv
Source: General Statistics Office of Vietnam. Downloaded 30 June 2022 from https://www.gso.gov.vn/en/data-and-statistics/2020/11/completed-results-of-the-2019-viet-nam-population-and-housing-census/.


