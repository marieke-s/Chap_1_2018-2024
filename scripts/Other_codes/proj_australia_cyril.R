library(sf)
library(ggplot2)

m <-st_read("~/Downloads/shp_cyril/meow_ecos.shp")
st_crs(m)
plot(m[1])

# m <- st_buffer(m, 0)  ## fixes some issues
# plot(m[1])
# m


#> Warning in st_buffer.sfc(st_geometry(x), dist, nQuadSegs, endCapStyle =
#> endCapStyle, : st_buffer does not correctly buffer longitude/latitude data
#> dist is assumed to be in decimal degrees (arc_degrees).
m = st_transform(st_crop(m, st_bbox(c(xmin = -30, xmax = -30, ymin = -90, ymax = 90))),  crs = '+proj=robin +ellps=WGS84 +datum=WGS84 +no_defs')

#> although coordinates are longitude/latitude, st_intersection assumes that they are planar
#> Warning: attribute variables are assumed to be spatially constant throughout all
#> geometries

plot(m[1])




#-----------

library(sf)

polys <- m


# Assume coordinates are lon/lat, so build a big rectangle west of -30
# adjust -180/180, -90/90 if your data extent is different
clip_poly <- st_sfc(
  st_polygon(list(rbind(
    c(-180, -90),
    c(-30,  -90),
    c(-30,   90),
    c(-180,  90),
    c(-180, -90)
  ))),
  crs = st_crs(polys)
)

# Intersect polygons with this clip polygon
polys_west <- st_intersection(polys, clip_poly)

plot(polys_west[1])

clip_poly_east <- st_sfc(
  st_polygon(list(rbind(
    c(-30, -90),
    c(180, -90),
    c(180, 90),
    c(-30, 90),
    c(-30, -90)
  ))),
  crs = st_crs(polys)
)

polys_east <- st_intersection(polys, clip_poly_east)
plot(polys_east[1])





























