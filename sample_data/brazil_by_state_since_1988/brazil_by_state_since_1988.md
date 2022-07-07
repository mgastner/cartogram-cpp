# Sources

## brazil_by_state_since_1988.geojson
Runfola D, Anderson A, Baier H, Crittenden M, Dowker E, Fuhrig S, et al. (2020) 
geoBoundaries: A global database of political administrative boundaries. 
PLoS ONE 15(4): e0231866. https://doi.org/10.1371/journal.pone.0231866. 

## Simplification
library(rmapshaper)
library(geojsonio)
library(sf)
library(mapview)
bra <- geojson_sf("~/brazil_by_state_since_1988/brazil_by_state_since_1988.geojson")
target_n_pts_in_output <- 50000
npts(bra)
bra_simp <- ms_simplify(bra, keep = target_n_pts_in_output/npts(bra))
geojson_write(
  bra_simp,
  lat = NULL,
  lon = NULL,
  geometry = "point",
  group = NULL,
  file = "brazil_by_state_since_1988.geojson",
  overwrite = TRUE,
  precision = NULL,
  convert_wgs84 = FALSE,
  crs = NULL)

### brazil_population_2019.csv
Source: Population estimates, Brazilian Institute of Geography and Statistics, 2019. Downloaded 30 June 2022 from 
https://www.ibge.gov.br/en/statistics/social/population/18448-estimates-of-resident-population-for-municipalities-and-federation-units.html?=&t=resultados

