#Adding and simplifying data for Vietnam map (GADM)
The source for boundaries is GADM. 
Downloaded 27.12.2021

Here is the code for simplification and conversion of the data (.sf) to GeoJSON.
```
library(mapview)  # Find number of points in sf geometry
library(rmapshaper)  # For ms_simplify()
library(sf) 
vietnam <- readRDS("~/Downloads/gadm36_VNM_1_sf.rds")
npts(vietnam)
n_pts_in_output <- 50000
vietnam_smp <- ms_simplify(vietnam, keep = n_pts_in_output / npts(vietnam))
npts(vietnam_smp)
plot(vietnam_smp$geometry)
st_area(vietnam_smp)
st_write(vietnam_smp, "vietnam.geojson")
#recreating the object
st_crs(vietnam$geometry) <- st_crs(vietnam$geometry)
```

## Bibliography for the CSV file
“Vietnam Cities Database | Simplemaps.Com.” Interactive HTML5 and JavaScript Maps for Websites | Simplemaps.Com, simplemaps, https://simplemaps.com/data/vn-cities. Accessed 2 Jan. 2022.