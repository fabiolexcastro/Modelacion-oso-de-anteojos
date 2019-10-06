
# Load libraries
require(pacman)
pacman::p_load(raster, rgdal, rgeos, stringr, velox, dismo, usdm, gtools, tidyverse)

g <- gc(reset = TRUE)
rm(list = ls())
options(scipen = 999, 
        stringsAsFactors = FALSE)

# Functions to use


# Load climate data
ecu <- list.files('../datos/ecu/ascii', full.names = TRUE, pattern = '.asc$') %>% 
  mixedsort() %>% 
  stack()
per <- list.files('../datos/per/ascii', full.names = TRUE, pattern = '.asc$') %>% 
  mixedsort() %>%
  stack()
occ1 <- read_csv('../datos/ecu/gbif/OSO_EXCEL.csv')
occ2 <- read_csv('../datos/per/gbif/Osos_Peru_UNIGIS.csv')
shp1 <- shapefile('../datos/ecu/shp/ECU_CONT.shp')
shp2 <- shapefile('../datos/per/shp/PERU.shp')

# Mosaicking the climate raster, Ecuador and Peru into only one raster
stk <- raster::mosaic(ecu, per, fun = 'mean')
plot(stk[[1]])

# Joining the occurrences tables into only one
occ <- rbind(occ1, occ2)

# Initial plot of the ocurrences data
png(filename = '../png/maps/points_v1.png', width = 9, height = 7, units = 'in', res = 300)
plot(stk[[1]])
plot(shp1, add = TRUE, border = 'black')
plot(shp2, add = TRUE, border = 'black')
points(occ$Longitude, occ$Latitude, pch = 16, col = 'red', cex = .65)
dev.off()

# Removing the duplicated by cell











