# Sources

## japan_by_prefecture_since_1888.geojson
Runfola D, Anderson A, Baier H, Crittenden M, Dowker E, Fuhrig S, et al. (2020) 
geoBoundaries: A global database of political administrative boundaries. 
PLoS ONE 15(4): e0231866. https://doi.org/10.1371/journal.pone.0231866. 
Downloaded from: https://github.com/wmgeolab/geoBoundaries/blob/main/releaseData/gbOpen/JPN/ADM1/geoBoundaries-JPN-ADM1_simplified.geojson on 4 July 2022.

### Code for simplification
```
library(rmapshaper)
library(geojsonio)
library(sf)
library(mapview)
jpn <- geojson_sf("https://raw.githubusercontent.com/wmgeolab/geoBoundaries/main/releaseData/gbOpen/JPN/ADM1/geoBoundaries-JPN-ADM1_simplified.geojson")
target_n_pts_in_output <- 48500
jpn_simp <- ms_simplify(jpn, keep = target_n_pts_in_output/npts(jpn))
geojson_write(
  jpn_simp,
  lat = NULL,
  lon = NULL,
  geometry = "point",
  group = NULL,
  file = "japan_by_prefecture_since_1888.geojson",
  overwrite = TRUE,
  precision = NULL,
  convert_wgs84 = FALSE,
  crs = NULL)
  ```

## japan_population_2020.csv
Source: Portal Site of Official Statistics of Japan website (https://www.e-stat.go.jp/)
Source: 2020 Population Census - Basic Complete Tabulation on Population and Households. Downloaded 30 June 2022 from https://www.e-stat.go.jp/en/stat-search/file-download?statInfId=000032142402&fileKind=0.


