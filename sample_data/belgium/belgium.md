#Adding and simplifying data for Belgium map (GADM)
The source for boundaries is GADM. 
Downloaded 27.12.2021

Here is the code for simplification and conversion of the data (.sf) to GeoJSON.
```
library(mapview)  # Find number of points in sf geometry
library(rmapshaper)  # For ms_simplify()
library(sf) 
belgium <- readRDS("~/Downloads/gadm36_BEL_1_sf.rds")
npts(belgium)
plot(belgium$geometry)
st_area(belgium)
st_write(belgium, "belgium.geojson")
```

## Bibliography for the CSV file
“Belgium Cities Database | Simplemaps.Com.” Interactive HTML5 and JavaScript Maps for Websites | Simplemaps.Com, simplemaps, https://simplemaps.com/data/be-cities. Accessed 2 Jan. 2022.