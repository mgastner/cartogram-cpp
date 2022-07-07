#Adding and simplifying data for Japan map (GADM)
The source for boundaries is GADM. 
Downloaded 27.12.2021

Here is the code for simplification and conversion of the data (.sf) to GeoJSON.
```
library(mapview)  # Find number of points in sf geometry
library(rmapshaper)  # For ms_simplify()
library(sf) 
japan <- readRDS("~/Downloads/gadm36_JPN_1_sf.rds")
npts(japan)
n_pts_in_output <- 50000
japan_smp <- ms_simplify(japan, keep = n_pts_in_output / npts(japan))
npts(japan_smp)
plot(japan_smp$geometry)
st_area(japan_smp)
st_write(japan_smp, "japan.geojson")
#recreating the object
st_crs(japan$geometry) <- st_crs(japan$geometry)
```

## Bibliography for the CSV file
“Japan Cities Database | Simplemaps.Com.” Interactive HTML5 and JavaScript Maps for Websites | Simplemaps.Com, simplemaps, https://simplemaps.com/data/jp-cities. Accessed 2 Jan. 2022.