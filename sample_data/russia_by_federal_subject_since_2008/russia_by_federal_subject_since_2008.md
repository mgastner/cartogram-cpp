# Sources

## russia_by_federal_subject_since_2008.geojson
Global Administrative Areas (2012). GADM database of Global Administrative Areas, version 2.0. [online] URL: www.gadm.org. 
Downloaded from: https://geodata.ucdavis.edu/gadm/gadm4.0/shp/gadm40_RUS_shp.zip on 4 July 2022.

### Code for simplification
```
rus <- geojson_sf("~/russia_by_federal_subject_since_2008/gadm40_RUS_1.geojson")
target_n_pts_in_output <- 48500
npts(rus)
rus_simp <- ms_simplify(rus, keep = target_n_pts_in_output/npts(rus))
geojson_write(
  rus_simp,
  lat = NULL,
  lon = NULL,
  geometry = "point",
  group = NULL,
  file = "russia_by_federal_subject_since_2008.geojson",
  overwrite = TRUE,
  precision = NULL,
  convert_wgs84 = FALSE,
  crs = NULL)
```

## russia_population_2010.csv
Source: Russian Federal State Statistics Service. 2010 All-Russian Population Census. Downloaded 30 June 2022 from https://en.wikipedia.org/wiki/List_of_federal_subjects_of_Russia_by_population#cite_note-2010Census-3.


