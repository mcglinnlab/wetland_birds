#' ---
#' title: "Wetland-Upland Bird Biodiversity Analysis"
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
#+ echo=FALSE
knitr::opts_knit$set(root.dir='../', tidy = TRUE)

#' Read in data
dat <- read.csv('./data/filtered_data/clean_bird_dat.csv')
comm <- read.csv('./data/filtered_data/clean_bird_comm.csv')
coords <- read.csv('./data/raw_data - wetland_coords.csv')

dat <- merge(dat, coords[ , c('wetland_id', 'utm_easting', 'utm_northing')],
             all.x = TRUE, all.y = FALSE)

# filtering out blocks 7 + 8
dat <- dat %>%
  filter(block < 7)


# filtering out 5 sites that don't have corresponding upland counts in the block
dat <- dat %>%
  filter(wetland_id != "SP05", wetland_id != "SP13", wetland_id != "SP12", wetland_id != "SP16", 
         wetland_id != "SP05", wetland_id != "SP11", wetland_id != "HH34", wetland_id != "HH05",
         wetland_id != "HH07")

dat$block <- as.factor(dat$block)

# getting rid of those same rows in the community matrix
comm <- comm %>%
  filter(X %in% dat$uni_id_date)

row.names(comm) <- comm[ , 1]
comm <- comm[ , -1]

## REDUNDANCY ANALYSIS to GET COMMUNITY COMP between SITE_TYPE ----------
## after controlling for SITE

# first aggregate comm4 at the wetland id year scale so that we reduce
# psuedo replication

bird_rda <- rda(comm ~ site_type + site, 
                data = dat)
bird_rda
RsquareAdj(bird_rda)

par(mfrow=c(1, 1))

plot(bird_rda, display = c('cn', 'sp'))
orditorp(bird_rda, display= 'sp', col = "red", air = 0.65)

pdf('fig4.pdf')
plot(bird_rda, type='n', scaling=1)
orditorp(bird_rda, display='sp', scaling = 1, cex=1, col='black', air = 1.5)
text(bird_rda, display='bp', col='red', labels = c('Upland', 'Wetland','Halidon', 'Stono'))
dev.off()

# parsing through both variables
anova(bird_rda, by='margin', permutations=1000)

#overall model fit
anova(bird_rda)

#head(dat)
#head(comm)

dim(dat)
dim(comm)

test <- dat %>%
  group_by(site_type) %>%
  summarize_if(is.numeric, sum)

test <- test[, c(1, 6:60)]

rowSums(test >0)


#' #Q1: is bird diversity higher in wetlands and uplands 
#div <- calc_biodiv(comm, dat$uni_id_date, effort = 5, extrapolate = TRUE)

dat$N <- rowSums(comm)
dat$S <- rowSums(comm > 0 )
dat$S_n <- apply(comm, 1, rarefaction, 'IBR', effort = 5, extrapolate = F,
                 quiet_mode = TRUE)
# singletons will result in rarefaction that results in 1 species which isn't useful
dat$S_n <- ifelse(dat$S_n == 1, NA, dat$S_n)
dat$S_PIE <- calc_PIE(comm)
dat$S_asymp <- apply(comm, 1, calc_chao1)

# glm models
div_mods <- list()
div_mods$N <- glm(N ~ site_type + site + block, data = dat, family = 'poisson')
div_mods$S <- glm(S ~ site_type + site + block, data = dat, family = 'poisson')
div_mods$S_n <- glm(S_n ~ site_type + site + block, data = dat, family = 'quasipoisson')
div_mods$S_PIE <- glm(S_PIE ~ site_type + site + block, data = dat, family = 'quasipoisson')
div_mods$S_asymp <- glm(S_asymp ~ site_type + site + block, data = dat, family = 'quasipoisson')

lapply(div_mods, summary)
lapply(div_mods, anova)


######### These are the model outputs reported in manuscript
#' use mixed effect model to account for pseudo-replicates

indices <- c('N', 'S', 'S_n', 'S_PIE', 'S_asymp')
div_mods_me <- vector('list', length(indices))
names(div_mods_me) <- indices
for(i in seq_along(indices)) {
  div_mods_me[[i]] <- lme(as.formula(paste(indices[i], "~ site_type")),
                        random = ~1 | site / year / block / wetland_id, data = dat,
                        na.action = na.omit)
}

lapply(div_mods_me, summary)
lapply(div_mods_me, r.squaredGLMM)
lapply(div_mods_me, anova)
lapply(div_mods_me, function(x) intervals(x, which = 'fixed'))


# those confidence intervals for the parameters are almost
# what we want use Bolker's function to get predicted means
# and standard errors or 95% CI's if desired. 
bolker_ci <- function(model, newdat, pred_int = FALSE, conf_level = 0.95) {
  if(class(model) != "lme") {
    stop("works for lme-models only")
  }
  z <- round(qnorm((1-conf_level)/2, lower.tail = FALSE), 2)
  newdat$pred <- predict(model, newdat, level = 0)
  Designmat <- model.matrix(formula(model)[-2], newdat)
  predvar <- diag(Designmat %*% vcov(model) %*% t(Designmat))
  newdat$se <- sqrt(predvar)
  newdat$ci_l <- newdat$pred - z*newdat$se
  newdat$ci_h <- newdat$pred + z*newdat$se
  if(pred_int == TRUE) {
    newdat$se2 <- sqrt(predvar+model$sigma^2)
    newdat$predint_l <- newdat$pred - z*newdat$se2
    newdat$predint_h <- newdat$pred + z*newdat$se2
  }
  newdat
}

newdata <- data.frame(site_type = as.factor(c('upland','wetland')))

#coefficients for ggplot figure construction
table_output <- list()
for(i in seq_along(div_mods_me)) {
  table_output[[i]] <- (bolker_ci(div_mods_me[[i]], newdata))
}
table_output2 <- as.data.frame(table_output)
write.csv(table_output, file = "./results/div_mods_me_ci.csv")

# creating the final PLOTS for upland v wetland-----

plt <- list()
ylabs <- c('Mean abundance (N)', 'Mean species richness (S)',
           expression(Mean~rarefied~richness~(S[n])),
           expression(Mean~species~evenness~(S[PIE])),
           expression(Mean~asymptotic~richness~(S[asym])))
xlabs <- c('', '', 'Habitat type', 'Habitat type', 'Habitat type')
for(i in seq_along(div_mods_me)) {
  plt[[i]] <- ggplot(bolker_ci(div_mods_me[[i]], newdata)) +
    geom_bar(aes(x=site_type, y = pred, fill = site_type), stat = "identity",
             width = 0.7) +
    geom_errorbar( aes(x=site_type, ymin=pred + se, ymax=pred - se), width=0.15, 
                   colour="black", alpha=0.7, size=0.5) + theme_bw() + 
    xlab(xlabs[i]) + ylab(ylabs[i]) +
    theme(legend.position = 'none')
}

## arranging plots
pdf('results/div_bar_plots.pdf')
grid.arrange(arrangeGrob(plt[[1]], plt[[2]], plt[[3]], plt[[4]]), ncol = 1, nrow = 1) 
dev.off()

# sup plot
pdf('results/asymS_barplot.pdf')
plt[[5]]
dev.off()



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
  group_by(wetland_id, year, block) %>% 
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
dat_agg$S_PIE <- calc_PIE(comm_agg)
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

# beta diversity and rarefaction analysis ------------
# the beta diversity analysis uses the aggregated community data which is the
# accumulation of three visits to the same site in a given year. This was done
# so that a reasonable minimum number of individuals is sampled at the site to compute
# a beta-diversity statistic for. 
# loop through each block and year and pull 1 upland and 1 wetland


# the rarefaction analysis follows a similar algo but it uses the 
# unaggregated community data so the grain is a single point count

nboot <- 500
uni_yrs <- unique(dat_agg$year)
uni_blocks <- 1:6 # not including blocks 7 & 8 here b/c these only had wetland sites
#uni_blocks <- 1:4 # to only do halidon sites

betas <- data.frame()

Navg_up <- sum(comm_agg[dat_agg$site_type == 'upland' & dat_agg$block %in% uni_blocks, ]) / 
           nrow(comm_agg[dat_agg$site_type == 'upland' & dat_agg$block %in% uni_blocks, ])
Navg_we <- sum(comm_agg[dat_agg$site_type == 'wetland' & dat_agg$block %in% uni_blocks, ]) /
           nrow(comm_agg[dat_agg$site_type == 'wetland' & dat_agg$block %in% uni_blocks, ])

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
    }
}  

# bootstrapped sample based rarefactions that go across spatial and temporal variation
# the grain is a single point count
curves <- data.frame()
#+ eval = FALSE
for (i in 1:nboot) {
  uplands <- data.frame()
  wetlands <- data.frame()
  for (j in seq_along(uni_yrs)) {
    for (k in seq_along(uni_blocks)) { 
      good_rows <- dat$year == uni_yrs[j] & 
        dat$block == uni_blocks[k] 
      sample_ids <- unique(dat$wetland_id[good_rows])
      # from this list draw a single wetland and a single upland
      upland_id <- sample_ids[grep('UP', sample_ids)]
      wetland_id <- sample(sample_ids[!(sample_ids %in% upland_id)], 1)
      upland_samples <- comm[good_rows & dat$wetland_id == upland_id, ]
      wetland_samples <- comm[good_rows & dat$wetland_id == wetland_id, ]
      # keep just one of the sites in the wetland samples
      #random_wetland_id <- sample(unique(wetland_samples$wetland_id), 1)
      #wetland_samples <- subset(wetland_samples, wetland_id == random_wetland_id)
      # ok now we have 3 samples from the upland and 3 samples from a wetland
      # in a specific year
      uplands <- rbind(uplands, upland_samples)
      wetlands <- rbind(wetlands, wetland_samples)
    }
  }
  # compute upland sample-based rarefaction curves
  S_up <- vegan::specaccum(uplands, conditioned = FALSE)
  S_up_pool <- specpool(uplands)
  curves <- rbind(curves,
                  data.frame(boot = i, site_type = 'upland', scale = 'gamma',
                             effort = c(S_up$sites, 100),
                             S = c(S_up$richness, S_up_pool$chao),
                             S_hi = c(S_up$richness + S_up$sd, 
                                      S_up_pool$chao + S_up_pool$chao.se),
                             S_lo = c(S_up$richness - S_up$sd, 
                                      S_up_pool$chao - S_up_pool$chao.se)))
  # compute wetland sample-based rarefaction curves
  S_we <- vegan::specaccum(wetlands, conditioned = FALSE)
  S_we_pool <- specpool(wetlands)
  curves <- rbind(curves,
                  data.frame(boot = i, site_type = 'wetland', scale = 'gamma',
                             effort = c(S_we$sites, 100),
                             S = c(S_we$richness, S_we_pool$chao),
                             S_hi = c(S_we$richness + S_we$sd, 
                                      S_we_pool$chao + S_we_pool$chao.se),
                             S_lo = c(S_we$richness - S_we$sd, 
                                      S_we_pool$chao - S_we_pool$chao.se)))
  
  S_study <- vegan::specaccum(rbind(wetlands, uplands), conditioned = FALSE)
  S_study_pool <- specpool(rbind(wetlands, uplands))
  curves <- rbind(curves,
                  data.frame(boot = i, site_type = 'both', scale = 'study',
                             effort = c(S_study$sites, 100),
                             S = c(S_study$richness, S_study_pool$chao),
                             S_hi = c(S_study$richness + S_study$sd,
                                      S_study_pool$chao + S_study_pool$chao.se),
                             S_lo = c(S_study$richness - S_study$sd,
                                      S_study_pool$chao - S_study_pool$chao.se)))
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

#' Here are the computed beta-diversity metrics with 95% CI 
beta_sum
beta_sum$site_type <- as.factor(beta_sum$site_type)
as.data.frame(beta_sum)

#re-ordering
beta_order <- rbind(beta_sum[3:8, ], beta_sum[1:2 , ])

# make barplot
par(mfrow=c(1,1))
tmp <- barplot(height = beta_order$beta_avg, col = c( "#F8766D", "#00BFC4"),
               border = NA, ylim = c(0, 5), 
               ylab = expression(italic(beta) * '-diversity'),
               names = beta_order$index, cex.names = 0.75, legend.text = TRUE)
arrows(x0 = tmp, 
       y0 = beta_order$beta_lo,
       y1 = beta_order$beta_hi, 
       angle = 90,
       code = 3,
       length = 0.1)
abline(h = 1, col='grey', lty =2, lwd =2)

# plot rarefaction results
# first aggregate across the bootstrap samples
curves$effort <- as.integer(curves$effort)
curves_sum <- curves %>% group_by(site_type, scale, effort) %>%
  summarize(S = mean(S), S_lo = mean(S_lo), S_hi = mean(S_hi))


#pdf('./results/sample_based_bootstrapped_rarefaction.pdf')
par(mfrow=c(1,1))
plot(S ~ effort , data = curves_sum, subset = scale == 'study' & effort < 100,
     xlim = c(1, 100), type = 'l', lwd = 2, log = 'xy', ylim = c(1,150), 
     xlab = 'Number of point counts', ylab = 'Species richness', 
     bty = 'n', axes=F)
axis(side=1, at = c(1, 2, 5, 10, 20, 50))
axis(side=2)
with(subset(curves_sum, site_type == 'wetland' & effort < 100), 
     polygon(c(1:36, 36:1), c(S_hi, rev(S_lo)), 
             col =  rgb(0, 191, 196, maxColorValue = 255, alpha = 255/2),
             border = NA))
with(subset(curves_sum, site_type == 'upland' & effort < 100), 
     polygon(c(1:36, 36:1), c(S_hi, rev(S_lo)), 
             col = rgb(248, 118, 109, maxColorValue = 255, alpha = 255/2), 
             border = NA))
with(subset(curves_sum, effort == 100),
     points(c(80, 93, 110), S, pch = 19, cex = 1.1, 
     col = c(1, "#00BFC4", "#F8766D")))
with(subset(curves_sum, scale == 'study' & effort == 100), 
     arrows(80, S_hi, 80, S_lo, angle = 90, code = 3, length = 0.075,
            lwd = 2))
with(subset(curves_sum, site_type == 'wetland' & effort == 100), 
     arrows(93, S_hi, 93, S_lo, angle = 90, code = 3, length = 0.075,
            col =  "#00BFC4", lwd = 2))
with(subset(curves_sum, site_type == 'upland' & effort == 100), 
     arrows(110, S_hi, 110, S_lo, angle = 90, code = 3, length = 0.075,
            col = "#F8766D", lwd = 2))
legend('bottomright', c('Wetlands and Uplands', 'Wetlands', 'Uplands'),
       lty = 1, lwd = c(2, 10, 10), 
       col = c(1, "#00BFC4", "#F8766D"), bty = 'n', inset = 0.2)
#dev.off()






########## Rarefaction for all point counts in both habitats
#' classic un-balanced rarefaction comparison
dat_mob_in <- make_mob_in(dat[ , 12:66], dat,
                          coord_names = c('utm_easting', 'utm_northing'))

par(mfrow=c(1,3))
plot_rarefaction(dat_mob_in, group_var = 'site_type',
                 method = 'IBR', avg= TRUE,
                 log='xy')
# sample based rarefaction unbalanced
Sstudy <- rarefaction(dat_mob_in$comm, 'SBR')
Sup <- rarefaction(dat_mob_in$comm[dat_mob_in$env$site_type == 'upland', ],
                   method = 'SBR')
Swet <- rarefaction(dat_mob_in$comm[dat_mob_in$env$site_type == 'wetland', ],
                    method = 'SBR')
par(mfrow=c(1,1))
plot(Sstudy, xlab = 'Number of Point Counts', ylab = 'Number of Species',
     type ='l', lwd = 2, frame.plot = F, xlim = c(0, 200), ylim = c(0, 50))
lines(Sup, col = '#F8766D', lwd = 2)
lines(Swet, col = '#00BFC4', lwd = 2)
legend('bottomright', c('Wetlands & Uplands', 'Wetlands', 'Uplands'),
       col=c('black','#00BFC4','#F8766D'), lty = 1, lwd = 2, bty='n', inset = 4)
?legend




####### SUPPLEMENTAL FIGURES
## Halidon v Stono






