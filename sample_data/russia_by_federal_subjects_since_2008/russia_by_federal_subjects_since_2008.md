# Sources

## russia_by_federal_subjects_since_2008.geojson
Runfola D, Anderson A, Baier H, Crittenden M, Dowker E, Fuhrig S, et al. (2020) 
geoBoundaries: A global database of political administrative boundaries. 
PLoS ONE 15(4): e0231866. https://doi.org/10.1371/journal.pone.0231866. 

## Code for simplification
rus <- geojson_sf("~/russia_by_federal_subjects_since_2008/russia_by_federal_subjects_since_2008.geojson")
target_n_pts_in_output <- 50000
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

### russia_population_2010.csv
Source: Russian Federal State Statistics Service. 2010 All-Russian Population Census. Downloaded 30 June 2022 from http://en.rian.ru/infographics/20111222/170405728.html.


