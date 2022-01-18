#Adding and simplifying data for Bahamas map (GADM)
The source for boundaries is GADM. 
Downloaded 27.12.2021

Here is the code for simplification and conversion of the data (.sf) to GeoJSON.
```
library(mapview)  # Find number of points in sf geometry
library(rmapshaper)  # For ms_simplify()
library(sf) 
bahamas <- readRDS("~/Downloads/gadm36_BHS_1_sf.rds")
npts(bahamas)
n_pts_in_output <- 50000
bahamas_smp <- ms_simplify(bahamas, keep = n_pts_in_output / npts(bahamas))
npts(bahamas_smp)
plot(bahamas_smp$geometry)
st_area(bahamas_smp)
st_write(bahamas_smp, "bahamas.geojson")
#recreating the object
st_crs(bahamas$geometry) <- st_crs(bahamas$geometry)
```

## Bibliography for the CSV file
“Bahamas Cities Database | Simplemaps.Com.” Interactive HTML5 and JavaScript Maps for Websites | Simplemaps.Com, simplemaps, https://simplemaps.com/data/bs-cities. Accessed 2 Jan. 2022.