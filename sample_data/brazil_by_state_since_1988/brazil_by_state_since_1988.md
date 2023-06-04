# Sources

## brazil_by_state_since_1988.geojson
Runfola D, Anderson A, Baier H, Crittenden M, Dowker E, Fuhrig S, et al. (2020) 
geoBoundaries: A global database of political administrative boundaries. 
PLoS ONE 15(4): e0231866. https://doi.org/10.1371/journal.pone.0231866. 
Downloaded from: https://github.com/wmgeolab/geoBoundaries/blob/main/releaseData/gbOpen/BRA/ADM1/geoBoundaries-BRA-ADM1_simplified.geojson on 4 July 2022.

### Simplification
```
library(rmapshaper)
library(geojsonio)
library(sf)
library(mapview)
bra <- geojson_sf("https://raw.githubusercontent.com/wmgeolab/geoBoundaries/main/releaseData/gbOpen/BRA/ADM1/geoBoundaries-BRA-ADM1_simplified.geojson")
target_n_pts_in_output <- 48500
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
```

## brazil_population_2021.csv
Source: Population estimates, Brazilian Institute of Geography and Statistics, 2019. Downloaded 30 June 2022 from 
https://ftp.ibge.gov.br/Estimativas_de_Populacao/Estimativas_2021/POP2021_20220711.pdf

