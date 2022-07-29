# Sources

## contiguous_usa_by_county_since_2001.geojson
United States Census Bureau. Cartographic Boundary Files.
Downloaded from: https://www2.census.gov/geo/tiger/GENZ2021/shp/cb_2021_us_county_500k.zip on 18 July 2022.
Lines within the .json file containing shapes for Alaska, Puerto Rico, US Virgin Islands, American Samoa, Northern Mariana Islands, Guam, and US Virgin Islands were manually deleted.

### Code for simplification
```
library(rmapshaper)
library(geojsonio)
library(sf)
library(mapview)
library(tidyverse)
temp <- tempfile()
download.file(
  "https://www2.census.gov/geo/tiger/GENZ2021/shp/cb_2021_us_county_500k.zip",
  temp
)
usa_unzip <- unzip(temp)
usa <- read_sf("cb_2021_us_county_500k.shp") |>
  filter(!(STATE_NAME %in% c(
    "Alaska",
    "Hawaii",
    "Commonwealth of the Northern Mariana Islands",
    "Puerto Rico",
    "Guam",
    "American Samoa",
    "United States Virgin Islands"
  )))
target_n_pts_in_output <- 48500
usa_simp <- ms_simplify(usa, keep = target_n_pts_in_output / npts(usa), keep_shapes = T)
geojson_write(
  usa_simp,
  lat = NULL,
  lon = NULL,
  geometry = "point",
  group = NULL,
  file = "contiguous_usa_by_county_since_2001.geojson",
  overwrite = TRUE,
  precision = NULL,
  convert_wgs84 = FALSE,
  crs = NULL
)
unlink(usa_unzip)

  ```
  
## usa_population_2021.csv
Source: United States Census Bureau. County Population Totals: 2020-2021.
Downloaded on 18 July 2022 from https://www2.census.gov/programs-surveys/popest/tables/2020-2021/counties/totals/co-est2021-pop.xlsx.


