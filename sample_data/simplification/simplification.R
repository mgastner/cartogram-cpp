library(rmapshaper)
library(geojsonio)
library(sf)
library(mapview)
bra <- geojson_sf("C:/Users/edwar/Documents/GitHub/cartogram-cpp/sample_data/brazil_by_state_2022_NEW/brazil_by_state_2022.geojson")
target_n_pts_in_output <- 50000
npts(bra)
bra_simp <- ms_simplify(bra, keep = target_n_pts_in_output/npts(bra))
geojson_write(
  bra_simp,
  lat = NULL,
  lon = NULL,
  geometry = "point",
  group = NULL,
  file = "brazil_by_state_2022.geojson",
  overwrite = TRUE,
  precision = NULL,
  convert_wgs84 = FALSE,
  crs = NULL)
