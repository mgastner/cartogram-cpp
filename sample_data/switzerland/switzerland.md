#Adding and simplifying data for Switzerland map (GADM)
The source for boundaries is GADM. 
Downloaded 27.12.2021

Here is the code for simplification and conversion of the data (.sf) to GeoJSON.
```
library(mapview)  # Find number of points in sf geometry
library(rmapshaper)  # For ms_simplify()
library(sf) 
switzerland <- readRDS("~/Downloads/gadm36_CHE_1_sf.rds")
npts(switzerland)
plot(switzerland$geometry)
st_area(switzerland)
st_write(switzerland, "switzerland.geojson")
```

## Bibliography for the CSV file
“Switzerland Cities Database | Simplemaps.Com.” Interactive HTML5 and JavaScript Maps for Websites | Simplemaps.Com, simplemaps, https://simplemaps.com/data/ch-cities. Accessed 2 Jan. 2022.