library(sf)
library(tmap)

phi <- seq(0, 2 * pi, by = 0.01 * pi)
phi[length(phi)] <- 0
outer <- cbind(2 * cos(phi), 2 * sin(phi))
hole <- round(cbind(cos(phi), -sin(phi)), digits = 6)
inner <- round(cbind(cos(phi), sin(phi)), digits = 6)
pol1 <- st_polygon(list(outer, hole))
pol2 <- st_polygon(list(inner))
sf_obj <-
  st_sf(name = c("outer", "inner"),
        geometry = st_sfc(pol1, pol2))
tm_shape(sf_obj) +
  tm_polygons(col = "name")
st_write(sf_obj, "concentric_circles.geojson", delete_dsn = TRUE)
