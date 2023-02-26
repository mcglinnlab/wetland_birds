#' ---
#' title: "Wetland Bird Analysis"
#' author: "Dan McGlinn"
#' date: '`r paste("First created on 2015-01-29. Updated on", Sys.Date())`'
#' output: 
#'  html_document
#'---

library(mobr)
library(dplyr)
library(nlme)
library(lme4)
library(MuMIn)


#+ echo=FALSE
knitr::opts_knit$set(root.dir='../', tidy = TRUE)

#' Read in data
dat <- read.csv('./data/filtered_data/clean_bird_dat.csv')
comm <- read.csv('./data/filtered_data/clean_bird_comm.csv')
row.names(comm) <- comm[ , 1]
comm <- comm[ , -1]

#head(dat)
#head(comm)

dim(dat)
dim(comm)

#' #Q1: is bird diversity higher in wetlands and uplands 
#div <- calc_biodiv(comm, dat$uni_id_date, effort = 5, extrapolate = TRUE)

dat$N <- rowSums(comm)
dat$S <- rowSums(comm > 0 )
dat$S_n <- apply(comm, 1, rarefaction, 'IBR', effort = 5, extrapolate = F,
                 quiet_mode = TRUE)
# singletons will result in rarefaction that results in 1 species which isn't useful
dat$S_n <- ifelse(dat$S_n == 1, NA, dat$S_n)
dat$S_PIE <- calc_SPIE(comm)
dat$S_asymp <- apply(comm, 1, calc_chao1)

par(mfrow=c(2,3)) 
plot(density(dat$N))
plot(density(dat$S))
plot(density(dat$S_n, na.rm = TRUE))
plot(density(dat$S_PIE, na.rm = TRUE))
plot(density(dat$S_asymp, na.rm = TRUE))


par(mfrow=c(2,3)) 
boxplot(N ~ site_type, data = dat)
boxplot(S ~ site_type, data = dat)
boxplot(S_n ~ site_type, data = dat)
boxplot(S_PIE ~ site_type, data = dat)
boxplot(S_asymp ~ site_type, data = dat)

div_mods <- list()
div_mods$N <- glm(N ~ site_type + site + block, data = dat, family = 'poisson')
div_mods$S <- glm(S ~ site_type + site + block, data = dat, family = 'poisson')
div_mods$S_n <- glm(S_n ~ site_type + site + block, data = dat, family = 'quasipoisson')
div_mods$S_PIE <- glm(S_PIE ~ site_type + site + block, data = dat, family = 'quasipoisson')
div_mods$S_asymp <- glm(S_asymp ~ site_type + site + block, data = dat, family = 'quasipoisson')

lapply(div_mods, summary)
lapply(div_mods, anova)

#' use mixed effect model to account for pseudo-replicates

indices <- c('N', 'S', 'S_n', 'S_PIE', 'S_asymp')
div_mods_me <- vector('list', length(indices))
names(div_mods_me) <- indices
for(i in seq_along(indices)) {
  div_mods_me[[i]] <- lme(as.formula(paste(indices[i], "~ site_type + site")),
                        random = ~1 | year / block / wetland_id, data = dat,
                        na.action = na.omit)
}

lapply(div_mods_me, summary)
lapply(div_mods_me, r.squaredGLMM)
lapply(div_mods_me, anova)

#' # Model Diagnostics
# fitted values vs residuals
par(mfrow=c(3,2))
for(i in seq_along(div_mods_me)) { 
  plot(predict(div_mods_me[[i]], newdata = dat), 
       residuals(div_mods_me[[1]]))
  abline(h= 0, col='red', lwd = 2)
}

# obs vs predicted plots
par(mfrow=c(3,2))
for(i in seq_along(div_mods_me)) { 
  plot(dat[ , indices[i]], 
       predict(div_mods_me[[i]], newdata = dat))
  abline(a = 0, b= 1, col='red', lwd = 2)
}


#' run again but after aggregating across replicates within a year
#' for all strings just take first value because these do not change during
#' repeat visits
dat_chr_agg <- dat %>%
  group_by(wetland_id, year) %>% 
  summarise_if(is.character, first)

#' for numeric values this mostly don't change visit to visit but a few do
#' like temp and wind speed so take an average
#' for the species we want to sum their counts though. 
names(dat)
sp_cols <- 12:66
dat_num_agg <- dat[ , -(12:66)] %>%
  group_by(wetland_id, year) %>% 
  summarise_if(is.numeric, mean, na.rm =TRUE)
dat_sp_agg <- dat[ , c('wetland_id', 'year', names(dat)[sp_cols])] %>% 
  group_by(wetland_id, year) %>% 
  summarize_all(sum)
# note on next line some of the names get changed with column binding
dat_agg <- cbind(dat_chr_agg, dat_num_agg, dat_sp_agg)
dim(dat_agg)
comm_agg <- dat_sp_agg[ , -(1:2)]
comm_agg

dat_agg$S <- rowSums(comm_agg > 0 )
dat_agg$N <- rowSums(comm_agg)
summary(dat_agg$N)
dat_agg$S_n <- apply(comm_agg, 1, rarefaction, 'IBR', effort = 10, extrapolate = F,
                 quiet_mode = TRUE)
dat_agg$S_n <- ifelse(dat_agg$S_n == 1, NA, dat_agg$S_n)
dat_agg$S_PIE <- calc_SPIE(comm_agg)
dat_agg$S_asymp <- apply(comm_agg, 1, calc_chao1)
dat_agg$pct_rare <- apply(comm_agg, 1, calc_div, 'pct_rare', rare_thres = 0.2)

par(mfrow=c(2,3)) 
plot(density(dat_agg$N))
plot(density(dat_agg$S))
plot(density(dat_agg$S_n, na.rm = TRUE))
plot(density(dat_agg$S_PIE, na.rm = TRUE))
plot(density(dat_agg$S_asymp, na.rm = TRUE))


par(mfrow=c(2,3)) 
boxplot(N ~ site_type, data = dat_agg)
boxplot(S ~ site_type, data = dat_agg)
boxplot(S_n ~ site_type, data = dat_agg)
boxplot(S_PIE ~ site_type, data = dat_agg)
boxplot(S_asymp ~ site_type, data = dat_agg)

div_mods <- list()
div_mods$N <- glm(N ~ site_type + site + block, data = dat_agg, family = 'poisson')
div_mods$S <- glm(S ~ site_type + site + block, data = dat_agg, family = 'poisson')
div_mods$S_n <- glm(S_n ~ site_type + site + block, data = dat_agg, family = 'quasipoisson')
div_mods$S_PIE <- glm(S_PIE ~ site_type + site + block, data = dat_agg, family = 'quasipoisson')
div_mods$S_asymp <- glm(S_asymp ~ site_type + site + block, data = dat_agg, family = 'quasipoisson')

lapply(div_mods, summary)
lapply(div_mods, anova)

#' use mixed effect model to account for pseudo-replicates

indices <- c('N', 'S', 'S_n', 'S_PIE', 'S_asymp')
div_mods_me <- vector('list', length(indices))
names(div_mods_me) <- indices
for(i in seq_along(indices)) {
  div_mods_me[[i]] <- lme(as.formula(paste(indices[i], "~ site_type + site")),
                          random = ~1 |  wetland_id...1, data = dat_agg,
                          na.action = na.omit)
}

lapply(div_mods_me, summary)
lapply(div_mods_me, r.squaredGLMM)
lapply(div_mods_me, anova)



