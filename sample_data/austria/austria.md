#Adding and simplifying data for Austria map (GADM)
The source for boundaries is GADM. 
Downloaded 27.12.2021

Here is the code for simplification and conversion of the data (.sf) to GeoJSON.
```
library(mapview)  # Find number of points in sf geometry
library(rmapshaper)  # For ms_simplify()
library(sf) 
austria <- readRDS("~/Downloads/gadm36_AUT_1_sf.rds")
npts(austria)
plot(austria$geometry)
st_area(austria)
st_write(austria, "austria.geojson")
```

## Bibliography for the CSV file
“Austria Cities Database | Simplemaps.Com.” Interactive HTML5 and JavaScript Maps for Websites | Simplemaps.Com, simplemaps, https://simplemaps.com/data/at-cities. Accessed 2 Jan. 2022.