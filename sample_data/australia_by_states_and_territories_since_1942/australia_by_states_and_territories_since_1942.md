# Sources

## australia_by_states_and_territories_since_1942.geojson
Runfola D, Anderson A, Baier H, Crittenden M, Dowker E, Fuhrig S, et al. (2020) 
geoBoundaries: A global database of political administrative boundaries. 
PLoS ONE 15(4): e0231866. https://doi.org/10.1371/journal.pone.0231866. 
Downloaded from: https://github.com/wmgeolab/geoBoundaries/tree/main/releaseData/gbOpen/AUS/ADM1 on 4 July 2022.

### Code for simplification
``` 
library(rmapshaper)
library(geojsonio)
library(sf)
library(mapview)
aus <- geojson_sf("~/australia_by_states_and_territories_since_1942/geoBoundaries-AUS-ADM1.geojson")
target_n_pts_in_output <- 48500
aus_simp <- ms_simplify(aus, keep = target_n_pts_in_output/npts(aus))
geojson_write(
  aus_simp,
  lat = NULL,
  lon = NULL,
  geometry = "point",
  group = NULL,
  file = "australia_by_states_and_territories_since_2021.geojson",
  overwrite = TRUE,
  precision = NULL,
  convert_wgs84 = FALSE,
  crs = NULL) 
```

## australia_population_2021.csv
Source: Australian Bureau of Statistics, National, state and territory population December 2021. Downloaded 30 June 2022 from https://www.abs.gov.au/census. 


