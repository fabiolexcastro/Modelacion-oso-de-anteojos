
# Load libraries
require(pacman)
pacman::p_load(raster, rgdal, rgeos, stringr, velox, dismo, usdm, gtools, tidyverse, rJava)

g <- gc(reset = TRUE)
rm(list = ls())
options(scipen = 999, 
        stringsAsFactors = FALSE)

# Functions to use
dup_cell <- function(mask, df){
  cellNum <- raster::extract(mask, df[,c('Lon', 'Lat')], cellnumbers = T) 
  cells <- xyFromCell(mask, cellNum[,'cells'])
  dupvec <- duplicated(cells[,c('x', 'y')])
  occ_rmDupCell <- tbl_df(df[!dupvec,])
  occ_DupCell <- tbl_df(df[dupvec,])
  return(occ_rmDupCell)
}

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

occ <- occ %>%
  setNames(c('Specie', 'Lat', 'Lon')) %>% 
  dplyr::select(Specie, Lon, Lat)

# Removing the duplicated by cell
occ_rmv <- dup_cell(mask = stk[[1]], df = occ)
nrow(occ) - nrow(occ_rmv)

# Extracting the values for the ocurrences 
occ_swd <- raster::extract(stk, occ_rmv[,2:3]) %>% 
  as.data.frame() %>% 
  setNames(c(paste0('bio_', 1:19))) %>% 
  mutate(Lon = pull(occ_rmv, 2),
         Lat = pull(occ_rmv, 3)) %>% 
  dplyr::select(Lon, Lat, everything())

# VIF Analysis
vif.res <- vif(x = occ_swd[,3:ncol(occ_swd)])
vif.step <- vifstep(x = occ_swd[,3:ncol(occ_swd)], th=10)
vrs <- vif.step@results$Variables

names(stk) <- paste0('bio_', 1:19)
stk_sub <- raster::subset(stk, grep(paste0(vrs, collapse = '|'), names(stk), value = FALSE))

occ_swd <- occ_swd %>% 
  dplyr::select(Lon, Lat, vrs) %>% 
  as_tibble()

# Generating the background (pseudo-absences)
back_raster <- stk_sub[[1]]
speciescell <- raster::extract(stk_sub[[1]], occ_swd[,1:2], cellnumber = TRUE)
back_raster[speciescell[,1]]  <- NA #remove the cell with presences
back <- randomPoints(back_raster, nrow(occ_swd)) %>%
  as_data_frame()
coordinates(back) <- ~ x + y
back_swd  <- raster::extract(stk_sub, back) %>% 
  cbind(coordinates(back), .)
write.csv(back_swd, '../model/occ/bck_run1.csv', row.names = FALSE)


















