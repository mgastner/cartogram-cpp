library(sf)
library(tmap)

phi <- seq(0, 2*pi, by = 0.01 * pi)
phi[length(phi)] <- 0
outer <- cbind(2 * cos(phi), 2 * sin(phi))
hole <- cbind(sin(phi), cos(phi))
inner <- cbind(cos(phi), sin(phi))
pol1 <- st_polygon(list(outer, hole))
pol2 <- st_polygon(list(inner))
sf_obj <- st_sf(name = c("outer", "inner"), geometry = st_sfc(pol1, pol2))
tm_shape(sf_obj) +
  tm_polygons(col = "name")
st_write(sf_obj, "concentric_circles.geojson")

# outer <- matrix(c(0, 0, 10, 0, 10, 10, 0, 10, 0, 0), ncol = 2, byrow = TRUE)
# hole1 <- matrix(c(1, 1, 1, 2, 2, 2, 2, 1, 1, 1), ncol = 2, byrow = TRUE)
# hole2 <- matrix(c(5, 5, 5, 6, 6, 6, 6, 5, 5, 5), ncol = 2, byrow = TRUE)

# pol1 <- list(outer, hole1, hole2)
# pol2 <- list(outer + 12, hole1 + 12)
# pol3 <- list(outer + 24)
# mp <- list(pol1, pol2, pol3)
# (mp1 <- st_multipolygon(mp))
# mp_sf <- st_sf(name = "bla", geometry = st_sfc(mp1))
# 
# tm_shape(mp_sf) +
#   tm_polygons()
