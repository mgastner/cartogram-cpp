#Adding and simplifying data for Algeria map (GADM)
The source for boundaries is GADM. 
Downloaded 27.12.2021

Here is the code for simplification and conversion of the data (.shp) to GeoJSON.
```
install.packages("mapview")
install.packages("rmapshaper")
install.packages("sf")
library(mapview) 
library(rmapshaper)  
library(sf) 
algeria <- readRDS("~/Downloads/gadm36_DZA_1_sf.rds") 
npts(algeria)
n_pts_in_output <- 50000
algeria_smp <- ms_simplify(algeria, keep = n_pts_in_output / npts(algeria))
npts(algeria_smp)
plot(algeria_smp$geometry)
st_area(algeria_smp)
st_write(algeria_smp, "algeria.geojson")
```

## Bibliography for the CSV file
“Algeria Cities Database | Simplemaps.Com.” Interactive HTML5 and JavaScript Maps for Websites | Simplemaps.Com, simplemaps, https://simplemaps.com/data/dz-cities. Accessed 2 Jan. 2022.
