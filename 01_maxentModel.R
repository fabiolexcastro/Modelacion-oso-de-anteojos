
# Load libraries
require(pacman)
pacman::p_load(raster, rgdal, rgeos, stringr, velox, dismo, usdm, gtools, tidyverse, rJava)

g <- gc(reset = TRUE)
rm(list = ls())
options(scipen = 999, 
        stringsAsFactors = FALSE)

options(java.parameters = '-Xmx4g')

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
  cbind(coordinates(back), .) %>% 
  as.data.frame %>% 
  rename(Lon = x,
         Lat = y)
write.csv(back_swd, '../model/occ/bck_run4.csv', row.names = FALSE)

# MaxEnt Model
# Evaluations using K-fold partioning
pres.covs <- rbind(occ_swd, back_swd)
fold <- kfold(pres.covs, k = 25)
occtest <- pres.covs[fold == 1,]
occtrain <- pres.covs[fold != 1,]

dir.create('../model/maxent/run4')
spp <- occ_swd
coordinates(spp) <- ~ Lon + Lat
bck <- back_swd
coordinates(bck) <- ~ Lon + Lat
mxn <- maxent(x = stk_sub, p = spp, a = bck, args = c('addsamplestobackground=true'), path = '../model/maxent/run4') 
prd <- predict(mxn, stk_sub)



y <- c(rep(1, nrow(occtrain)), rep(0, nrow(back_swd)))
env.values <- data.frame(rbind(occtrain, back_swd))
dir.create('../model/maxent/run1', recursive = TRUE)

# Now, for all folds
auc <- rep(NA, 25)
max.tss <- rep(NA,5)
maps <- list()
e <- list()
me <- list()
dir_out <- paste0('../model/maxent/run1/', 1:25)

for (i in 1:25){
  occtest <- pres.covs[fold == i, ]
  occtrain <- pres.covs[fold != i, ]
  env.values <- data.frame(rbind(occtrain, back_swd))
  y  <- c(rep(1, nrow(occtrain)), rep(0, nrow(back_swd)))
  me[[i]] <- maxent(env.values[,3:ncol(env.values)], y, args = c('addsamplestobackground=true'), path = dir_out[[i]])
  maps[[i]] <- predict(me[[i]], stk_sub)
  e[[i]] <- evaluate(me[[i]], p = data.frame(occtest[,3:ncol(occtest)]), a = data.frame(back_swd[,3:ncol(back_swd)]))
  auc[i] <- e[[i]]@auc
  lines((1 - e[[i]]@TNR), e[[i]]@TPR)
  tss <- e[[i]]@TPR + e[[i]]@TNR-1
  max.tss[i] <- e[[i]]@t[which.max(tss)]
}

map_avg <- mean(stack(maps))
Map('writeRaster', x = maps, filename = paste0('../model/maxent/avg/map_run_', 1:25, '.asc'))
writeRaster(map_avg, '../model/maxent/run1/avg/map_avg.asc')

dir.create('../model/rds')
th_tss <- mean(max.tss)
























