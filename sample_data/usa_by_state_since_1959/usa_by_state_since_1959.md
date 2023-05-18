# Sources

## usa_by_state_since_1959.geojson
Runfola D, Anderson A, Baier H, Crittenden M, Dowker E, Fuhrig S, et al. (2020) 
geoBoundaries: A global database of political administrative boundaries. 
PLoS ONE 15(4): e0231866. https://doi.org/10.1371/journal.pone.0231866.
Downloaded from: https://github.com/wmgeolab/geoBoundaries/blob/main/releaseData/gbOpen/USA/ADM1/geoBoundaries-USA-ADM1_simplified.geojson on 4 July 2022.

### Code for simplification
```
usa <- geojson_sf("https://raw.githubusercontent.com/wmgeolab/geoBoundaries/main/releaseData/gbOpen/USA/ADM1/geoBoundaries-USA-ADM1_simplified.geojson")
target_n_pts_in_output <- 48500
npts(usa)
usa_simp <- ms_simplify(usa, keep = target_n_pts_in_output/npts(usa))
geojson_write(
  usa_simp,
  lat = NULL,
  lon = NULL,
  geometry = "point",
  group = NULL,
  file = "usa_by_state_since_1959.geojson",
  overwrite = TRUE,
  precision = NULL,
  convert_wgs84 = FALSE,
  crs = NULL)
```
## usa_population_2020.csv
Source: The United States Census Bureau. Downloaded 30 June 2022 from https://www2.census.gov/programs-surveys/popest/tables/2020-2021/counties/totals/co-est2021-pop.xlsx.






