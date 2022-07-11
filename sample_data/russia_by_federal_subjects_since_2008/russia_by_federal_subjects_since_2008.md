# Sources

## russia_by_federal_subjects_since_2008.geojson
Runfola D, Anderson A, Baier H, Crittenden M, Dowker E, Fuhrig S, et al. (2020) 
geoBoundaries: A global database of political administrative boundaries. 
PLoS ONE 15(4): e0231866. https://doi.org/10.1371/journal.pone.0231866. 
Downloaded from: https://github.com/wmgeolab/geoBoundaries/tree/main/releaseData/gbOpen/RUS/ADM1 on 4 July 2022.

### Code for simplification
```
rus <- geojson_sf("~/russia_by_federal_subjects_since_2008/geoBoundaries-RUS-ADM1.geojson")
target_n_pts_in_output <- 48500
npts(rus)
rus_simp <- ms_simplify(rus, keep = target_n_pts_in_output/npts(rus))
geojson_write(
  rus_simp,
  lat = NULL,
  lon = NULL,
  geometry = "point",
  group = NULL,
  file = "russia_by_federal_subjects_since_2008.geojson",
  overwrite = TRUE,
  precision = NULL,
  convert_wgs84 = FALSE,
  crs = NULL)
```

## russia_population_2010.csv
Source: Russian Federal State Statistics Service. 2010 All-Russian Population Census. Downloaded 30 June 2022 from https://en.wikipedia.org/wiki/List_of_federal_subjects_of_Russia_by_population#cite_note-2010Census-3.


