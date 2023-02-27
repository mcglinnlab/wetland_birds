#' ---
#' title: "Wetland Bird Analysis"
#' author: "Dan McGlinn"
#' date: '`r paste("First created on 2022-02-24. Updated on", Sys.Date())`'
#' output: 
#'  pdf_document
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
lapply(div_mods_me, function(x) intervals(x, which = 'fixed'))

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
# fix wonky names
names(dat_agg) <- c('wetland_id', 'year', names(dat_agg)[-(1:2)])
dim(dat_agg)
comm_agg <- dat_sp_agg[ , -(1:2)]
comm_agg
sum(comm_agg)
sum(comm)

rowSums(comm_agg)
rowSums(comm_agg[dat_agg$wetland_id == "UP03", ])


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
lapply(div_mods_me, function(x) intervals(x, which = 'fixed'))

#' use mixed effect model to account for pseudo-replicates

indices <- c('N', 'S', 'S_n', 'S_PIE', 'S_asymp')
div_mods_me <- vector('list', length(indices))
names(div_mods_me) <- indices
for(i in seq_along(indices)) {
  div_mods_me[[i]] <- lme(as.formula(paste(indices[i], "~ site_type + site")),
                          random = ~1 |  wetland_id, data = dat_agg,
                          na.action = na.omit)
}

lapply(div_mods_me, summary)
lapply(div_mods_me, r.squaredGLMM)
lapply(div_mods_me, anova)

# beta diversity analysis ------------
# loop through each block and year and pull 1 upland and 1 wetland
# Note:  two blocks do not have an upland and a wetlands sites, blocks 7 & 8
# Note: block 3 the upland site UP03 only has 1 individual across the 3 visits
# in the 0-25m range so drop this site and block from analysis

nboot <- 500
uni_yrs <- unique(dat_agg$year)
uni_blocks <- 1:6 # not including blocks 7 & 8 here b/c these only had wetland sites
#uni_blocks <- 1:4 # to only do halidon sites

betas <- data.frame()
curves <- data.frame()

Navg_up <- sum(comm_agg[dat_agg$site_type == 'upland', ]) / 
           nrow(comm_agg[dat_agg$site_type == 'upland', ])
Navg_we <- sum(comm_agg[dat_agg$site_type == 'wetland', ]) /
           nrow(comm_agg[dat_agg$site_type == 'wetland', ])

#+ eval = FALSE
for (i in 1:nboot) {
  for (j in seq_along(uni_yrs)) {
    uplands <- data.frame()
    wetlands <- data.frame()
    for (k in seq_along(uni_blocks)) { 
       good_rows <- dat_agg$year == uni_yrs[j] & 
                    dat_agg$block == uni_blocks[k] 
       sample_ids <- unique(dat_agg$wetland_id[good_rows])
       # from this list draw a single wetland and a single upland
       upland_id <- sample_ids[grep('UP', sample_ids)]
       wetland_id <- sample(sample_ids[!(sample_ids %in% upland_id)], 1)
       upland_samples <- comm_agg[good_rows & dat_agg$wetland_id == upland_id, ]
       wetland_samples <- comm_agg[good_rows & dat_agg$wetland_id == wetland_id, ]
       # keep just one of the sites in the wetland samples
       #random_wetland_id <- sample(unique(wetland_samples$wetland_id), 1)
       #wetland_samples <- subset(wetland_samples, wetland_id == random_wetland_id)
       # ok now we have 3 samples from the upland and 3 samples from a wetland
       # in a specific year
       uplands <- rbind(uplands, upland_samples)
       wetlands <- rbind(wetlands, wetland_samples)
    }
    # now capture rarefaction curve and compute beta div at two scales point count to block & block to study
    #up_N_min <- min(rowSums(uplands))
    #we_N_min <- min(rowSums(wetlands))
    betas <- rbind(betas, 
                   data.frame(boot = i , site_type = 'upland',
                     calc_comm_div(uplands, index = c('S', 'S_n', 'S_PIE'),
                       effort = 10, scale = 'beta')))
    betas <- rbind(betas, 
                   data.frame(boot = i , site_type = 'wetland',
                     calc_comm_div(wetlands, index = c('S', 'S_n', 'S_PIE'),
                       effort = 10, scale = 'beta')))
    nmin <- max(min(rowSums(uplands)), 5)
    S_up <- apply(uplands, 1, rarefaction, 'IBR', effort  = 1:nmin, extrapolate = TRUE,
                  quiet_mode = TRUE)
    if (is.matrix(S_up)) S_up <- rowMeans(S_up)
    curves <- rbind(curves,
                    data.frame(boot = i, site_type = 'upland', scale = 'alpha',
                               effort = 1:nmin, S = S_up))
    curves <- rbind(curves,
                    data.frame(boot = i, site_type = 'upland', scale = 'gamma',
                               effort = 1:sum(uplands), S = rarefaction(uplands, 'IBR')))
    
    nmin <- max(min(rowSums(wetlands)), 5)
        S_we <- apply(wetlands, 1, rarefaction, 'IBR', effort  = 1:nmin, extrapolate = TRUE,
                  quiet_mode = TRUE)
    if (is.matrix(S_we)) S_we <- rowMeans(S_we)
    

    curves <- rbind(curves,
                    data.frame(boot = i, site_type = 'wetland', scale = 'alpha',
                               effort = 1:nmin, S = S_we))
    curves <- rbind(curves,
                    data.frame(boot = i, site_type = 'wetland', scale = 'gamma',
                               effort = 1:sum(wetlands), S = rarefaction(wetlands, 'IBR')))
  }
}  

#+ eval = FALSE
save(betas, curves, file = './results/div_boostrap_results.Rdata')


#+ eval = TRUE
load(file = './results/div_boostrap_results.Rdata')

#' aggregated across boostraps
#+ eval = TRUE
beta_sum <- betas %>% group_by(index, site_type) %>%
  summarize(beta_avg = mean(value), beta_lo = quantile(value, 0.025),
            beta_hi = quantile(value, 0.975))
curves_sum <- curves %>% group_by(site_type, scale, effort) %>%
  summarize(S_avg = mean(S), S_lo = quantile(S, 0.025), S_hi = quantile(S, 0.975))


#' Here are the computed beta-diversity metrics with 95% CI 
beta_sum
# make parplot
par(mfrow=c(1,1))
tmp <- barplot(height = beta_sum$beta_avg, col = c('pink', 'dodgerblue'), ylim = c(0, 5),
               names = beta_sum$index, cex.names = 0.75)
arrows(x0 = tmp, 
       y0 = beta_sum$beta_lo,
       y1 = beta_sum$beta_hi, 
       angle = 90,
       code = 3,
       length = 0.1)
abline(h = 1, col='grey', lty =2, lwd =2)


nmax <- curves %>% group_by(site_type, scale, boot) %>%
  summarize(N_max = max(effort)) %>%
  summarize(N_min = min(N_max))

plot(S_hi ~ effort, curves_sum, subset = site_type == 'wetland' & scale == 'gamma',
     type = 'l', col= 'dodgerblue', log = 'xy')
lines(S_lo ~ effort, curves_sum, subset = site_type == 'wetland' & scale == 'gamma',
      col = 'dodgerblue')
lines(S_hi ~ effort, curves_sum, subset = site_type == 'wetland' & scale == 'alpha',
      col ='dodgerblue')
lines(S_lo ~ effort, curves_sum, subset = site_type == 'wetland' & scale == 'alpha',
      col = 'dodgerblue')

lines(S_hi ~ effort, curves_sum, subset = site_type == 'upland' & scale == 'gamma',
      col = 'pink')
lines(S_lo ~ effort, curves_sum, subset = site_type == 'upland' & scale == 'gamma',
      col = 'pink')
lines(S_hi ~ effort, curves_sum, subset = site_type == 'upland' & scale == 'alpha', 
      col = 'pink')
lines(S_lo ~ effort, curves_sum, subset = site_type == 'upland' & scale == 'alpha',
      col = 'pink')


plot(S_avg ~ effort, curves_sum, subset = site_type == 'wetland' & scale == 'gamma',
     type = 'l', col= 'dodgerblue', log = 'xy')
lines(S_avg ~ effort, curves_sum, subset = site_type == 'wetland' & scale == 'alpha',
      col ='dodgerblue')

lines(S_avg ~ effort, curves_sum, subset = site_type == 'upland' & scale == 'gamma',
      col = 'pink')
lines(S_avg ~ effort, curves_sum, subset = site_type == 'upland' & scale == 'alpha', 
      col = 'pink')

#' classic un-balanced rarefaction comparison
dat_mob_in <- make_mob_in(dat[ , 12:66], dat)
par(mfrow=c(1,3))
plot_rarefaction(dat_mob_in, group_var = 'site_type', method = 'IBR', avg= TRUE,
                 log='xy')


tst <- aggregate(dat_agg[ , 69:123], list(dat_agg$site_type), sum)
rowSums(tst[ , -1])
calc_comm_div(tst[, -1], index = c('S', 'S_n', "S_PIE"), effort = 90, scales = 'beta')




