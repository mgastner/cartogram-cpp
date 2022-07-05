# Sources

## vietnam_by_province_since_2021.geojson
Global Administrative Areas ( 2012). GADM database of Global Administrative Areas, version 2.0. [online] URL: www.gadm.org.


## Code for simplification
library(rmapshaper)
library(geojsonio)
library(sf)
library(mapview)
vnm <- geojson_sf("~/vietnam_by_province_since_2022_FAULTY/vietnam_by_province_since_2021.geojson")
target_n_pts_in_output <- 50000
npts(vnm)
vnm_simp <- ms_simplify(vnm, keep = target_n_pts_in_output/npts(vnm))
geojson_write(
  vnm_simp,
  lat = NULL,
  lon = NULL,
  geometry = "point",
  group = NULL,
  file = "vietnam_by_province_since_2021.geojson",
  overwrite = TRUE,
  precision = NULL,
  convert_wgs84 = FALSE,
  crs = NULL)

### vietnam_population_2019.csv
Source: General Statistics Office of Vietnam. Downloaded 30 June 2022 from https://www.gso.gov.vn/en/data-and-statistics/.


