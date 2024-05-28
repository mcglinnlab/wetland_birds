#' ---
#' title: "Within Wetland Bird Analysis"
#' author: "Dan McGlinn and Jackson Barratt Heitmann"
#' date: '`r paste("First created on 2022-02-24. Updated on", Sys.Date())`'
#' output: 
#'  pdf_document
#'---
library(vegan)
library(mobr)
library(dplyr)
library(nlme)
library(lme4)
library(MuMIn)
library(scales)
library(ggplot2)
library(gridExtra)
library(tidyr)
library(performance)

dat <- read.csv('./data/filtered_data/clean_bird_dat.csv')
comm <- read.csv('./data/filtered_data/clean_bird_comm.csv')
coords <- read.csv('./data/raw_data - wetland_coords.csv')

row.names(comm) <- comm[ , 1]
comm <- comm[ , -1]

## First let's try an RDA to see which env variables are driving community comp
comm_wet <- comm[1:150, ]
dat_wet <- dat[1:150, ]


# creating the RDA
rda_tree <- rda(comm_wet ~ canopy_cover_sim + crown_hull_mean + pine_dbh + tupelo_dbh + arc_area + 
                  tree_height_sd + tree_height_mean + groundcover, data=dat_wet)

rda_tree

RsquareAdj(rda_tree)

plot(rda_tree, type='n', scaling=1)
orditorp(rda_tree, display='sp', cex=0.5, scaling=1, col='blue')
text(rda_tree, display='cn', col='red')


anova(rda_tree, by='margin', permutations=1000)


pseudo_r2 = function(glm_mod) {
  1 -  glm_mod$deviance / glm_mod$null.deviance
}

# calculate diversity metrics
dat_wet$N <- rowSums(comm_wet)
dat_wet$S <- rowSums(comm_wet > 0 )
dat_wet$S_n <- apply(comm_wet, 1, rarefaction, 'IBR', effort = 5, extrapolate = F,
                 quiet_mode = TRUE)
# singletons will result in rarefaction that results in 1 species which isn't useful
dat_wet$S_n <- ifelse(dat_wet$S_n == 1, NA, dat_wet$S_n)
dat_wet$S_PIE <- calc_PIE(comm_wet)
dat_wet$S_asymp <- apply(comm_wet, 1, calc_chao1)

# Problem with PIE here, the variable has a ton of values at 1 = perfectly even
# creating the glm models 
div_mods_env <- list()
div_mods_env$N <- glm(N ~ canopy_cover_sim + crown_hull_mean + pine_dbh + tupelo_dbh + arc_area + 
                        tree_height_sd + tree_height_mean + groundcover, data = dat_wet, family = 'poisson')
div_mods_env$S <- glm(S ~ canopy_cover_sim + crown_hull_mean + pine_dbh + tupelo_dbh + arc_area + 
                        tree_height_sd + tree_height_mean + groundcover, data = dat_wet, family = 'poisson')
div_mods_env$S_n <- glm(S_n ~ canopy_cover_sim + crown_hull_mean + pine_dbh + tupelo_dbh + arc_area + 
                          tree_height_sd + tree_height_mean + groundcover, data = dat_wet, family = 'poisson')
div_mods_env$S_PIE <- glm(S_PIE ~ canopy_cover_sim + crown_hull_mean + pine_dbh + tupelo_dbh + arc_area + 
                            tree_height_sd + tree_height_mean + groundcover, data = dat_wet, family = 'quasipoisson')
div_mods_env$S_asymp <- glm(S_asymp ~ canopy_cover_sim + crown_hull_mean + pine_dbh + tupelo_dbh + arc_area + 
                              tree_height_sd + tree_height_mean + groundcover, data = dat_wet, family = 'poisson')

lapply(div_mods_env, summary)

#check the models for assumptions of normality
lapply(div_mods_env, check_model)

lapply(div_mods_env, anova)
lapply(div_mods_env, pseudo_r2)

# mixed effect models to predict env variables after controlling for wetland_id
# year and block
indices_env <- c('N', 'S', 'S_n', 'S_PIE', 'S_asymp')
div_mods_env <- null_mods_env <- vector('list', length(indices))
names(div_mods_env) <- indices
for(i in seq_along(indices)) {
  div_mods_env[[i]] <- lme(as.formula(paste(indices_env[i], "~ canopy_cover_sim + crown_hull_mean + tupelo_dbh + wet_area + pine_dbh")),
                           random = ~1 |  year / block / wetland_id, data = dat_wet,
                           na.action = na.omit)
  null_mods_env[[i]] <- lme(as.formula(paste(indices_env[i], "~ 1")),
                            random = ~1 |  year / block / wetland_id, data = dat_wet,
                            na.action = na.omit)
}

lapply(div_mods_env, summary)
lapply(div_mods_env, r.squaredGLMM)
lapply(div_mods_env, anova)
mapply(anova, div_mods_env, null_mods_env)

# individual species abundance models
# subsetting to get individual species
#NOPA ----

comm_NOPA <- comm_wet["NOPA"]
as.data.frame(comm_NOPA)
class(comm_NOPA$NOPA)

div_mods_NOPA <- (lme(comm_NOPA$NOPA ~ canopy_cover_sim + 
                        crown_hull_mean + pine_dbh + tupelo_dbh + wet_area,
                      random = ~1 | year / block / wetland_id, data = dat_wet,
                      na.action = na.omit))

alpha_div_NOPA <- (glm(comm_NOPA$NOPA ~ canopy_cover_sim + 
                         crown_hull_mean + pine_dbh + tupelo_dbh + wet_area,
                       data = dat_wet, family = binomial))
summary(alpha_div_NOPA)
?glm
histogram(comm_NOPA$NOPA)

# WEVI ----

comm_WEVI <- comm_wet["WEVI"]
as.data.frame(comm_WEVI)
class(comm_WEVI$WEVI)

div_mods_WEVI <- (lme(comm_WEVI$WEVI ~ canopy_cover_sim + 
                        crown_hull_mean + pine_dbh + tupelo_dbh + wet_area,
                      random = ~1 | year / block / wetland_id, data = dat_wet,
                      na.action = na.omit))

alpha_div_WEVI <- (glm(comm_WEVI$WEVI ~ canopy_cover_sim + 
                         crown_hull_mean + pine_dbh + tupelo_dbh + wet_area, 
                       data = dat_wet, family = "poisson"))

summary(alpha_div_WEVI)
?glm
histogram(comm_WEVI$WEVI)

# COYE ----

comm_COYE <- comm_wet["COYE"]
as.data.frame(comm_COYE)
class(comm_COYE$COYE)

alpha_div_COYE <- (glm(comm_COYE$COYE ~ canopy_cover_sim + 
                         crown_hull_mean + pine_dbh + tupelo_dbh + wet_area, 
                       data = dat_wet, family = poisson))

summary(alpha_div_COYE)
?glm
histogram(comm_COYE$COYE)