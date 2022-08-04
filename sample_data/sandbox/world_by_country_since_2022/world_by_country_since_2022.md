# Sources

## world_by_country_since_2022.geojson
Runfola D, Anderson A, Baier H, Crittenden M, Dowker E, Fuhrig S, et al. (2020) 
geoBoundaries: A global database of political administrative boundaries. 
PLoS ONE 15(4): e0231866. https://doi.org/10.1371/journal.pone.0231866.
Downloaded from: https://github.com/wmgeolab/geoBoundaries/blob/main/releaseData/CGAZ/geoBoundariesCGAZ_ADM0.geojson on 19 July 2022.

### Code for simplification
```
wld <- geojson_sf("https://github.com/wmgeolab/geoBoundaries/raw/main/releaseData/CGAZ/geoBoundariesCGAZ_ADM0.geojson")
target_n_pts_in_output <- 48500
npts(wld)
wld_simp <- ms_simplify(wld, keep = target_n_pts_in_output/npts(wld), keep_shapes=T)
geojson_write(
  wld_simp,
  lat = NULL,
  lon = NULL,
  geometry = "point",
  group = NULL,
  file = "world_by_country_since_2022.geojson",
  overwrite = TRUE,
  precision = NULL,
  convert_wgs84 = FALSE,
  crs = NULL)
```
## world_population_2010.csv
Source: The United Nations Department of Economic and Social Affairs Population Division. World Population Prospects. Downloaded 19 July 2022 from https://population.un.org/wpp/Download/Standard/CSV/.






