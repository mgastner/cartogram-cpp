#Adding and simplifying data for Russia map (GADM)
The source for boundaries is GADM. 
Downloaded 27.12.2021

Here is the code for simplification and conversion of the data (.sf) to GeoJSON.
```
library(mapview)  # Find number of points in sf geometry
library(rmapshaper)  # For ms_simplify()
library(sf) 
russia <- readRDS("~/Downloads/gadm36_RUS_1_sf.rds")
npts(russia)
n_pts_in_output <- 50000
russia_smp <- ms_simplify(russia, keep = n_pts_in_output / npts(russia))
npts(russia_smp)
plot(russia_smp$geometry)
st_area(russia_smp)
st_write(russia_smp, "russia.geojson")
```

## Bibliography for the CSV file
“Russia Cities Database | Simplemaps.Com.” Interactive HTML5 and JavaScript Maps for Websites | Simplemaps.Com, simplemaps, https://simplemaps.com/data/ru-cities. Accessed 2 Jan. 2022.