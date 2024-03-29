##all necessary packages
library(lattice)
library(ggplot2)
library(tidyr)
library(dplyr)
library(readr)
library(mobr)
library(vegan)
library(lidR)
library(remotes)
library(lidR)
library(sp)
library(terra)
library(rgdal)
library(sf)
library(geometry)
library(gridExtra)
library(lme4)
library(GGally)

## LiDAR Data Manipulations
## This LiDAR analysis is broken up by site (STONO, HALIDON) due to the size of the .las
## files. The code is exactly the same for both, until we merge the two
## outputs on lines 253-255.

# lidar -----

##STONO

las <- readLAS("Stono_Preserve_all_returns.las", select='xyzcr' )
wetland_coords <- read.csv("raw_data - wetland_coords.csv")

## coordinate change to utm for Stono
coords_ll <- las@data[ , c('X', 'Y')]
coords_ll <- SpatialPoints(coords_ll, 
                           proj4string =  CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"))

# plot(coords_ll)
coords_utm <- sp::spTransform(coords_ll,
                              CRSobj = CRS("+proj=utm +zone=17 +datum=WGS84"))
# plot(coords_utm)

las@data$X <- coordinates(coords_utm)[ , 1]
las@data$Y <- coordinates(coords_utm)[ , 2]

#https://epsg.io/32617#:~:text=WGS%2084%20%2F%20UTM%20zone%2017N%20%2D%20EPSG%3A32617
st_crs(las) <- 32617

#checking to see if the coordinates changed to UTM properly
head(las@data)
table(las@data$Classification)

plot(las)
# str(las)


# normalizing height by creating a dtm, 
dtm <- rasterize_terrain(las, 1, tin())
# plot(dtm, col = gray(1:50/50))

nlas <- las - dtm
# plot(nlas, size = 4, bg = "white")

# checking to see if we actually normalized the height
hist(filter_ground(nlas)$Z, breaks = seq(-0.6, 0.6, 0.01), main = "", xlab = "Elevation")

#now normalized need to rasterize the canopy

chm_stono <- rasterize_canopy(nlas, 1, p2r(0.15))

ttops <- locate_trees(nlas, lmf(4), uniqueness = "incremental")

algo <- dalponte2016(chm_stono, ttops)
head(algo)
ctg_segmented <- segment_trees(nlas, algo)

# try to run a for loop that derives canopy variables 
stono_indices <- c(21:31, 36:37)
circle_LAS_by_wet <- clip_circle(ctg_segmented,
                                 x = wetland_coords$utm_easting[stono_indices],
                                 y = wetland_coords$utm_northing[stono_indices],
                                 radius = 50)

names(circle_LAS_by_wet) <- wetland_coords$wetland_id[stono_indices]

stono_crown_metrics <- lapply(circle_LAS_by_wet, crown_metrics, .stdtreemetrics , geom = "convex")
stono_tree_height_mean <- unlist(lapply(stono_crown_metrics, function(x) mean(x$Z)))
stono_tree_height_mean <- as.data.frame(stono_tree_height_mean)

stono_tree_height_sd <- unlist(lapply(stono_crown_metrics, function(x) sd(x$Z)))
stono_tree_height_sd <- as.data.frame(stono_tree_height_sd)

stono_crown_hull_area_mean <- unlist(lapply(stono_crown_metrics, function(x) mean(x$convhull_area)))
stono_crown_hull_area_mean <- as.data.frame(stono_crown_hull_area_mean)

stono_tree_metrics <- cbind(stono_tree_height_mean, stono_tree_height_sd, stono_crown_hull_area_mean)
as.data.frame(stono_tree_metrics)
colnames(stono_tree_metrics) <- c("tree_height_mean", "tree_height_sd", "crown_hull_mean")


# lidar_out <- data.frame(wetland_coords[stono_indices, ], mean_tree, sd_tree)
# lidar_out_SP01 <- lidar_out[ 1,]

##MAKE THIS INTO A FUNCTION

# Estimating forest cover using simulation model and crown hulls
n_pts <- 1e4

all_canopy_cover <- data.frame()

for (i in seq_along(stono_crown_metrics)){
  tst <- as(stono_crown_metrics[[i]], 'Spatial')
  box <- bbox(tst)
  rnd_pts <- data.frame(x = runif(n_pts, box[1,1], box[1,2]),
                        y = runif(n_pts, box[2,1], box[2,2]))
  
  # need to either only simulate points within the circle or 
  # drop those points greater than radius of circle from center
  center <- as.data.frame(t(rowMeans(box)))
  radius <- 50
  pt_dist <- raster::pointDistance(center, rnd_pts, lonlat = FALSE)
  rnd_pts <- rnd_pts[pt_dist < radius, ]
  # update number of random points
  n_pts <- nrow(rnd_pts)
  # set projection which is needed for over function
  rnd_pts <- SpatialPoints(rnd_pts)
  proj4string(rnd_pts) <- crs(tst)
  
  out <- over(rnd_pts, tst)
  
  canopy_cover <- data.frame(names(stono_crown_metrics)[i],
                             (sum(!is.na(out$treeID)) / n_pts) * 100)
  
  all_canopy_cover <- rbind(all_canopy_cover, canopy_cover)
}

colnames(all_canopy_cover) <- c('wetland_id', 'canopy_cover_sim')


##LiDAR Data Manipulations
## Halidon 

las_hh <- readLAS('Hallidon_Hill_all_returns.las', select='xyzc')


## coordinate change for Halidon
coords_hh <- las_hh@data[ , c('X', 'Y')]
coords_hh <- SpatialPoints(coords_hh, 
                           proj4string =  CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"))

#plot(coords_hh)
coords_utm_hh <- sp::spTransform(coords_hh,
                                 CRSobj = CRS("+proj=utm +zone=17 +datum=WGS84"))
#plot(coords_utm_hh)

las_hh@data$X <- coordinates(coords_utm_hh)[ , 1]
las_hh@data$Y <- coordinates(coords_utm_hh)[ , 2]

# Dan added the following because this was used with the stono las
st_crs(las_hh) <- 32617

# class(las_hh@data$Y)
# plot(las_hh)

# normalize height by creating a dtm, seems like maybe we don't have to do this
dtm_hh <- rasterize_terrain(las_hh, 1, tin())

nlas_hh <- las_hh - dtm_hh

chm_hh <- rasterize_canopy(nlas_hh, 1, p2r(0.15))

ttops_hh <- locate_trees(nlas_hh, lmf(4), uniqueness = "incremental")

algo_hh <- dalponte2016(chm_hh, ttops_hh)
head(algo)
ctg_segmented_hh <- segment_trees(nlas_hh, algo_hh)

# plot(nlas_hh, size = 4, bg = "white")

# segmenting trees
# las_segmented_trees_hh <- segment_trees(nlas_hh, li2012(R = 3, speed_up = 5))

# plot(las_segmented_trees_hh, color = "treeID")


#trying to clip and subset based on plot locations for halidon


hh_indices <- c(1:20, 32:35)
circle_LAS_by_wet_hh <- clip_circle(ctg_segmented_hh,
                                    x = wetland_coords$utm_easting[hh_indices],
                                    y = wetland_coords$utm_northing[hh_indices],
                                    radius = 25)
names(circle_LAS_by_wet_hh) <- wetland_coords$wetland_id[hh_indices]

# mean_tree_hh <- NULL
# sd_tree_hh <- NULL
# for (i in seq_along(circle_LAS_by_wet_hh)) {
# extract canopy variables
# mean_tree_hh[i] <- cloud_metrics(circle_LAS_by_wet_hh[[i]], func = ~mean(Z))
# sd_tree_hh[i] <- cloud_metrics(circle_LAS_by_wet_hh[[i]], func = ~sd(Z))
#}

# lidar_out_hh <- data.frame(wetland_coords[hh_indices, ], mean_tree_hh, sd_tree_hh)

hh_crown_metrics <- lapply(circle_LAS_by_wet_hh, crown_metrics, .stdtreemetrics , geom = "convex")
hh_tree_height_mean <- unlist(lapply(hh_crown_metrics, function(x) mean(x$Z)))
hh_tree_height_mean <- as.data.frame(hh_tree_height_mean)

hh_tree_height_sd <- unlist(lapply(hh_crown_metrics, function(x) sd(x$Z)))
hh_tree_height_sd <- as.data.frame(hh_tree_height_sd)

hh_crown_hull_area_mean <- unlist(lapply(hh_crown_metrics, function(x) mean(x$convhull_area)))
hh_crown_hull_area_mean <- as.data.frame(hh_crown_hull_area_mean)

hh_tree_metrics <- cbind(hh_tree_height_mean, hh_tree_height_sd, hh_crown_hull_area_mean)

colnames(hh_tree_metrics) <- c("tree_height_mean", "tree_height_sd", "crown_hull_mean")

lidar_tree_metrics <- rbind(stono_tree_metrics, hh_tree_metrics)

# Estimating forest cover using simulation model and crown hulls
n_pts <- 1e4

all_canopy_cover_hh <- data.frame()

for (i in seq_along(hh_crown_metrics)){
  tst <- as(hh_crown_metrics[[i]], 'Spatial')
  box <- bbox(tst)
  rnd_pts <- data.frame(x = runif(n_pts, box[1,1], box[1,2]),
                        y = runif(n_pts, box[2,1], box[2,2]))
  
  # need to either only simulate points within the circle or 
  # drop those points greater than radius of circle from center
  center <- as.data.frame(t(rowMeans(box)))
  radius <- 50
  pt_dist <- raster::pointDistance(center, rnd_pts, lonlat = FALSE)
  rnd_pts <- rnd_pts[pt_dist < radius, ]
  # update number of random points
  n_pts <- nrow(rnd_pts)
  # set projection which is needed for over function
  rnd_pts <- SpatialPoints(rnd_pts)
  proj4string(rnd_pts) <- crs(tst)
  
  out <- over(rnd_pts, tst)
  
  canopy_cover_hh <- data.frame(names(hh_crown_metrics)[i],
                                (sum(!is.na(out$treeID)) / n_pts) * 100)
  
  all_canopy_cover_hh <- rbind(all_canopy_cover_hh, canopy_cover_hh)
}

colnames(all_canopy_cover_hh) <- c('wetland_id', 'canopy_cover_sim')

## Combine both canopy cover simulations from Halidon and Stono

all_canopy_cover <- rbind(all_canopy_cover, all_canopy_cover_hh)

lidar_tree_metrics <- cbind(lidar_tree_metrics, all_canopy_cover)
lidar_tree_metrics <- lidar_tree_metrics[ , c(-6)]

write.csv(lidar_tree_metrics, file = 'lidar_tree_metrics.csv', row.names = FALSE)