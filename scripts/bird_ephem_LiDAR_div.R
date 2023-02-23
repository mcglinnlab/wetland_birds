library(dplyr)
library(vegan)
library(mobr)
library(pwr)
library(ggplot2)
library(remotes)
install.packages("sf")
library(rgdal)
library(sf)
#LiDAR data manipulation
library(treetop)
library(lidR)
launchApp()
??treetop
treetop::launchApp()
#reading the file in
stono_las <- readLAS("merged Stono all returns_clipped.las")
print(stono_las)
las_check(stono_las)
plot(hill_las)
hill_las <- readLAS("Halidon_Hill_clipped.las")

SC_peedee_lid <- readLAS("./USGS_LPC_SC_SavannahPeeDee_2019_B19_10493971.laz")

#need to create shapefile to calculate means within each polygon
map <-readOGR(dsn=".", layer="ephem_bird_LiDAR_poly_shape")
plots <- sf::st_read("bird_shape", quiet = TRUE)
class(plots)
center <- read.csv("~/raw_data - Sheet5.csv")
class(center)
metrics <- plot_metrics(stono_las, .stdmetrics_z, center, radius = 50)

las <- classify_ground(stono_las, algorithm = pmf(ws = 5, th = 3))

#normalizing the surface model
stono_norm <- normalize_height(stono_las, knnidw())

ttops <- locate_trees(SC_peedee_lid, lmf(ws = 5))
x <- plot(SC_peedee_lid)
add_treetops3d(x, ttops)
plot(x)
#trying to segment trees
stono_seg_tree <- segment_trees(stono_norm, li2012())
#trying to actually get the mean metrics for tree heights
tree_metrics(stono_seg_tree, func = ~max(Z))

#trying to make a sf data frame to use as the centroid
center$geometry <- as.list(center$geometry)

mysf <- st_as_sf(x = center$wetland_id, geometry = center$geometry)

?st_sfc
cloud_metrics(las, func = ~mean(Z)) # calculate mean height
cloud_metrics(hill_las, func = ~mean(Z))

polygon_metrics(stono_las, func = ~mean(Z))
#trying to do tree detection and plotting


chm <- rasterize_canopy(stono_las, 0.5, pitfree(subcircle = 0.2))
plot(stono_las, bg = "white", size = 4)


?ggmap
remotes::install_github("gavinsimpson/ggvegan")
library(ggvegan)

# need to fix block variable in excel, 
# block 3 replicated across block 4

dat <- read.csv("~/raw_bird_data.csv")

lid_stono <- read.csv("~/clipped_stono.csv", head = TRUE, sep = " ")



head(dat)
wet_attribute_only <- read.csv("~/raw_data - just_wet_attribute.csv")

summary(aov(bird_sr ~ as.factor(canopy_cover) + canopy_dom + 
             as.factor(midstory_cover) + herb_layer
            , data = wet_attribute_only))

summary(aov(bird_sr ~ canopy_dom, data = wet_attribute_only))

# sr ~ Midstory cover
boxplot(bird_sr ~ property, data = wet_attribute_only)

canopy_bird_sr <- aov(bird_sr ~ canopy_cover, data = wet_attribute_only)

plot(canopy_bird_sr$residuals, canopy_bird_sr$fitted.values)
shapiro.test(wet_attribute_only$shan_sr)
# sr ~ 
ggplot(data = wet_attribute_only, 
       mapping = aes(x = canopy_cover, y = bird_sr, colour = canopy_cover)) + 
  geom_point() + geom_boxplot() + theme_bw()

summary(aov(bird_sr ~ herb_layer, data = wet_attribute_only))

wet_attribute_only$block <- factor(wet_attribute_only$block)

#shannon div and canopy cover
summary(aov(shan_sr ~ canopy_cover + herb_layer + midstory_cover, data = wet_attribute_only))


ggplot(data = wet_attribute_only, 
       mapping = aes(x = canopy_cover, y = shan_sr, colour = canopy_cover)) + 
  geom_point() + geom_boxplot() + theme_bw() + 
  xlab("Canopy Cover (%), broken into 4 categories") +
  ylab("Shannon-Weiner Index of Biodiversity at wetland sites")

plot(canopy_bird_sr$residuals, canopy_bird_sr$fitted)
sr_pivot <- read.csv("~/raw_data - pivot.csv")
head(sr_pivot)
head(com_sum)


# Shannon-Weiner diversity index
shan_comm <- diversity(com_sum, index = "shannon")
shan_comm
shan_comm <- as.data.frame(shan_comm)

site_type <- c("wetland", "wetland", "wetland", "wetland", "wetland",
              "wetland", "wetland", "wetland", "wetland", "wetland",
              "wetland", "wetland", "wetland", "wetland", "wetland",
              "wetland", "wetland", "wetland", "wetland", "upland",
              "upland", "upland", "upland", "upland", "upland")
shan_comm2 <- cbind(site_type, shan_comm)
head(shan_comm2)
summary(aov(shan_comm ~ site_type, data = shan_comm2))

#boxplot of shannon-weiner biodiv
ggplot(data = wet_attribute_only, 
       mapping = aes(x = canopy_cover, y = bird_sr, colour = canopy_cover)) + 
  geom_point() + geom_boxplot() + theme_bw()

dim(wet_attribute)

bird_sr <- dat_sr$sr
print(dat_sr)
print(bird_sr)

comm_2 <- as.data.frame(comm)
comm_2
shan_comm <- diversity(comm, index = "shannon")
shan_comm
row.names(comm_2)
total_shan_div <- colSums(matrix(shan_comm, nrow=3))
comm_total <- rowSums(comm_2,rep(1:3,each=43))

?geom_boxplot
head(comm)
print(comm_total)
total_shan_div <- as.data.frame(total_shan_div)
total_shan_div <- total_shan_div/3
head(total_shan_div)
print(total_shan_div)
rownames(total_shan_div) <- c('HH02', 'HH04', 'HH07', 'HH13', 'HH14',
                         'HH17', 'HH22', 'HH33', 'HH34', 'HH48', 
                         'HH49', 'HH50', 'SP01', 'SP02', 'SP03', 
                         'SP05', 'SP07', 'SP08', 'SP09', 'UP01',
                         'UP02', 'UP03', 'UP04', 'UP05', 'UP06')
length(shan_comm)

n1 <- 3
head(comm)

?reshape
dim(comm)
dim(dat_sr)
head(dat_sr)
head(dat)

# assumming each species only recorded once per point count then
# species richness (sr) is just number of rows which function n() computes
dat_sr <- dat %>%
    group_by(date, site, site_type, block, wetland_id) %>%
    summarise(sr = n())

# richness is normalish in distribution
plot(density(dat_sr$sr))
hist(dat_sr$sr)

sd(dat_sr$sr)
?pwr.f2.test
# 3 treatments + control - 3 df
# wet vs dry - 1 df
# stono vs halidon - 1 df
# yr 1 vs yr 2 = 1 df
# total of 6 df
# f2 = R2 / (1-R2)
# 10 % effect
pwr.f2.test(u = 6, 246, f2=0.1/(1-0.1),sig.level=0.05)
# 5 % effect
pwr.f2.test(u = 6, 246, f2=0.05/(1-0.05), sig.level=0.05)

# just looking at Halidon hill for just year 2 and 3 where we have control wetlands
# power to detect treatment effects at point count scale
pwr.f2.test(u = 3, 120-4, f2=0.1/(1-0.1), sig.level=0.05) # power is 87% which is great
# power to detect treatment effects at wetland
pwr.f2.test(u = 3, 20, f2=0.1/(1-0.1), sig.level=0.05) # power is 20% which is poor
pwr.f2.test(u = 3, 20, f2=0.35/(1-0.35), sig.level=0.05) # power is 79% 
# ok so if psudeo-reps are summarized to wetland scale then need 35% effect to be detected



summary(lm(sr ~ block + site_type + site, data = dat_sr))
summary(glm(sr ~ block + site_type + site, data = dat_sr, family='poisson'))

boxplot(sr ~ site_type, data = dat_sr)

id <- with(dat, paste(wetland_id, date))
comm <- with(dat, tapply(total, list(id, species), function(x) sum(x)))
comm <- ifelse(is.na(comm), 0, comm)


# trying to just get wetland

id_2 <- with(dat, paste(wetland_id))
comm_new <- with(dat, tapply(total, list(id_2, species), function(x) sum(x)))
comm_new <- ifelse(is.na(comm), 0, comm)


env <- aggregate(dat[ , c('block', 'site_type', 'site')], list(id), function(x) x[1])
head(id)
head(comm)
head(env)
env$Group.1 == row.names(comm)

bird_rda <- rda(comm ~  as.factor(site_type) + as.factor(site), data = env)
bird_rda
plot(bird_rda, display = c('cn', 'sp'))

RsquareAdj(bird_rda)

autoplot(bird_rda, display = c('cn', 'sp'))
autoplot(bird_rda)

bird_wet_rda <- rda(comm ~ as.factor(canopy_cover) + as.factor(midstory_cover) +
                      as.factor(herb_layer), data = wet_attribute)

bird_in <- make_mob_in(comm, env)

bird_wet_mob <- make_mob_in(comm, wet_attribute)

mob_test <- get_delta_stats(bird_in, env_var = 'site_type', ref_level = 'upland',
                            type = 'discrete', tests = c('SAD', 'N'), n_perm = 199)
                                
plot(mob_test)

par(mfrow=c(2,2))
plot_rarefaction(bird_in, group_var = 'site_type', method = 'IBR')
plot_rarefaction(bird_in, group_var = 'site_type', method = 'SBR')
plot_rarefaction(bird_in, group_var = 'site', method = 'IBR')
plot_rarefaction(bird_in, group_var = 'site', method = 'SBR')
plot_rarefaction(subset(bird_in, subset = site == 'halidon'), group_var = 'block', method = 'SBR')
plot_rarefaction(subset(bird_in, subset = site == 'stono'), group_var = 'block', method = 'SBR')



