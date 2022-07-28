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
library(tidyverse)
temp <- tempfile()
download.file(
  "https://geodata.ucdavis.edu/gadm/gadm4.0/shp/gadm40_VNM_shp.zip",
  temp
)
vnm_unzip <- unzip(temp)
vnm <- read_sf("gadm40_VNM_1.shp")
target_n_pts_in_output <- 45000
vnm_simp <- ms_simplify(vnm, keep = target_n_pts_in_output / npts(vnm))
npts(vnm_simp)
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
  crs = NULL
)
unlink(vnm_unzip)
  ```
  
## vietnam_population_2019.csv
Source: General Statistics Office of Vietnam. Downloaded 30 June 2022 from https://www.gso.gov.vn/wp-content/uploads/2019/12/Ket-qua-toan-bo-Tong-dieu-tra-dan-so-va-nha-o-2019.pdf, and accessed via https://data.humdata.org/dataset/bf6858bf-a5a8-403b-8380-f78cd72c9fee/resource/90642c40-f854-49a2-bdd2-fdea7db4130d/download/vnm_admpop_2019.xlsx.


