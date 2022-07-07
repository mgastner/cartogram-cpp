#Adding and simplifying data for Germany map (GADM)
The source for boundaries is GADM. 
Downloaded 27.12.2021

Here is the code for simplification and conversion of the data (.sf) to GeoJSON.
```
library(mapview)  # Find number of points in sf geometry
library(rmapshaper)  # For ms_simplify()
library(sf) 
germany <- readRDS("~/Downloads/gadm36_DEU_1_sf.rds")
npts(germany)
n_pts_in_output <- 50000
germany_smp <- ms_simplify(germany, keep = n_pts_in_output / npts(germany))
npts(germany_smp)
plot(germany_smp$geometry)
st_area(germany_smp)
st_write(germany_smp, "germany.geojson")
#recreating the object
st_crs(germany$geometry) <- st_crs(germany$geometry)
```

## Bibliography for the CSV file
“Germany Cities Database | Simplemaps.Com.” Interactive HTML5 and JavaScript Maps for Websites | Simplemaps.Com, simplemaps, https://simplemaps.com/data/de-cities. Accessed 2 Jan. 2022.