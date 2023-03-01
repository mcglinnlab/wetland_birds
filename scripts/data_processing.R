library(dplyr)
library(tidyr)

# read in data --------------------------------------------------------------
bird_lng <- read.csv("./data/raw_data - all_bird_data.csv")

# create unique wetland_id date field so that each unique
# sample (i.e., point count) has a separate identifier
bird_lng$uni_id_date <- with(bird_lng, paste(wetland_id, date, sep='_'))

# check that all sites visited 3 times
tapply(bird_lng$uni_id_date, list(bird_lng$wetland_id, bird_lng$year), function(x) length(unique(x)))
# SP03 was visited 7 times
# 3 times in 2021 and 4 times in 2022
bird_lng[bird_lng$wetland_id == 'SP03' & bird_lng$year == '2022', ]
# just keep the first 3 visits in 2022
bird_lng <- subset(bird_lng, subset = uni_id_date != 'SP03_6/8/2022')

#cols_to_drop <- c('method', 'X25m', 'X50m', 'breeding.codes')
#bird_wid <- subset(bird_lng, select = !(names(bird_lng) %in% cols_to_drop)) %>% 
#  pivot_wider(names_from = species, values_from = total, values_fill = 0,
#              values_fn = sum)
cols_to_drop <- c('method', 'X50m', 'total', 'breeding.codes')
bird_lng$X25m <- ifelse(is.na(bird_lng$X25m), 0, bird_lng$X25m)
bird_wid <- subset(bird_lng, select = !(names(bird_lng) %in% cols_to_drop)) %>% 
  pivot_wider(names_from = species, values_from = X25m, values_fill = 0,
              values_fn = sum)
dim(bird_wid)
length(unique(bird_lng$uni_id_date))

# env matrix
veg <- read.csv("./data/raw_data - veg_2022.csv")

# fix mistake in veg for site SP01 where it is coded as an upland
veg$site_type[veg$wetland_id == 'SP01']
veg$site_type[veg$wetland_id == 'SP01'] <- rep('wetland', 4)

# wetland assessment env matrix

wet_env <- read.csv("./data/raw_data - just_wet_attribute.csv")
summary(wet_env$block)
# note in this data file a lot of the blocks are incomplete listed as NAs 
# drop this block column as it is misleading
wet_env <- subset(wet_env, select = -block)

lidar <- read.csv("./data/lidar_tree_metrics.csv")

## aggregate vegetation data across 4 tree samples ----------------------------
table(veg$wetland_id) # ok so each row of this matrix has 4 closest tree species
# this dataset needs to be aggregated
# note site HH17 has 8 trees recorded - maybe from year 1 and year 2 - not clear
# also note that the column spp.number isn't really informative
veg_agg <- veg %>% 
  group_by(wetland_id) %>%
  summarize(site = site[1], site_type = site_type[1], block = block[1],
            tree_sr = length(unique(spp)), tree_dist = mean(dist),
            tree_dbh = sum(dbh), tree_ba = sum(pi * (dbh / 2)^2),
            tupelo_dbh = sum(dbh[spp == 'tupelo']),
            pine_dbh = sum(dbh[spp %in% c('loblolly', 'loblolly/slash', 'slash/loblolly', 'longleaf')]),
            other_dbh = sum(dbh[!(spp %in% c('tupelo','loblolly', 'loblolly/slash', 'slash/loblolly', 'longleaf'))]))
# quick checks for mistakes
sum(veg$dbh)
sum(veg_agg$tree_dbh)
sum(veg_agg$tupelo_dbh + veg_agg$pine_dbh + veg_agg$other_dbh)


all(bird_lng$wetland_id %in% veg_agg$wetland_id)
all(veg_agg$wetland_id %in% bird_lng$wetland_id)

head(wet_env)
with(wet_env, tapply(date, list(wetland_id), length))

# merge bird counts with field measured vegetation attributes 
dim(bird_wid)
dim(veg_agg)
dat <- merge(bird_wid, veg_agg, by = c('wetland_id', 'site', 'site_type', 'block'))
dim(dat)
head(dat)

# merge that matrix with field measured environmental attributes
dat <- merge(dat, wet_env, by = 'wetland_id')
dim(dat)
# date.x is the date of bird survey
# date.y is the date of the env wetland survey

# merge lidar data

dat <- merge(dat, lidar, by = 'wetland_id')
dim(dat)

# pull out just the species matrix
sp_col_indices <- match(unique(bird_lng$species), names(dat))
comm <- dat[ , sp_col_indices]
names(comm)
rownames(comm) <- dat$uni_id_date
head(comm)

# export data files
dir.create('./data/filtered_data')
write.csv(dat, file = './data/filtered_data/clean_bird_dat.csv', row.names = FALSE)
write.csv(comm, file = './data/filtered_data/clean_bird_comm.csv', row.names = TRUE)
