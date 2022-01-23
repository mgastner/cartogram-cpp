#Adding and simplifying data for USA map (GADM)
The source for boundaries is GADM. 
Downloaded 27.12.2021

Here is the code for simplification and conversion of the data (.sf) to GeoJSON.
```
library(mapview)  # Find number of points in sf geometry
library(rmapshaper)  # For ms_simplify()
library(sf) 
usa <- readRDS("~/Downloads/gadm36_USA_1_sf.rds")
npts(usa)
n_pts_in_output <- 50000
usa_smp <- ms_simplify(usa, keep = n_pts_in_output / npts(usa))
npts(usa_smp)
plot(usa_smp$geometry)
st_area(usa_smp)
st_write(usa_smp, "usa.geojson")
#recreating the object
st_crs(usa$geometry) <- st_crs(usa$geometry)
```

## Bibliography for the CSV file
“US States - Ranked by Population 2021.” 2021 World Population by Country, https://worldpopulationreview.com/states. Accessed 7 Jan. 2022.