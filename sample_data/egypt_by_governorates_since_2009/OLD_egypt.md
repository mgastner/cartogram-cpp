#Adding and simplifying data for Egypt map (GADM)
The source for boundaries is GADM. 
Downloaded 27.12.2021

Here is the code for simplification and conversion of the data (.sf) to GeoJSON.
```
library(mapview)  # Find number of points in sf geometry
library(rmapshaper)  # For ms_simplify()
library(sf) 
egypt <- readRDS("~/Downloads/gadm36_EGY_1_sf.rds")
npts(egypt)
n_pts_in_output <- 50000
egypt_smp <- ms_simplify(egypt, keep = n_pts_in_output / npts(egypt))
npts(egypt_smp)
plot(egypt_smp$geometry)
st_area(egypt_smp)
st_write(egypt_smp, "egypt.geojson")
#recreating the object
st_crs(egypt$geometry) <- st_crs(egypt$geometry)
```

## Bibliography for the CSV file
“Egypt Cities Database | Simplemaps.Com.” Interactive HTML5 and JavaScript Maps for Websites | Simplemaps.Com, simplemaps, https://simplemaps.com/data/eg-cities. Accessed 2 Jan. 2022.