################################################################################
######### 07: Modeling saliva cortisol ########
################################################################################

######################################################################
##### 1.0 Loading data #####
######################################################################

########## 1.1 Set working directory & download packages ##########

rm(list = ls())
setwd("~/Documents/R/Saliva")
options(stringsAsFactors = F)
library(GGally)
library(lme4)
library(lmerTest)
library(MuMIn)
library(multcomp)
library(DHARMa)
library(sjPlot)
library(performance)
library(viridis)
library(patchwork)
library(rptR)
library(tidyverse)

#Set plot colors using viridis
viridis_2 <- viridis(7)[-c(1,3,5,6,7)]
viridis_3 <- viridis(7)[-c(2,4,6,7)]

#Load data
load("07.saliva_behavior.Rdata")

#Reformat time
saliva.cortisol$start_time <- as.POSIXct(format(saliva.cortisol$start_time, format = "%H:%M"), 
                                         format = "%H:%M")
saliva.cortisol$stop_time <- as.POSIXct(format(saliva.cortisol$stop_time, format = "%H:%M"), 
                                        format = "%H:%M")


######################################################################
##### 2.0 Exploring data #####
######################################################################

########## 2.1 Summarize data ##########

nrow(saliva.cortisol)     #262 samples
length(unique(saliva.cortisol$hyena_id))     #69 hyenas

#Age
hist(saliva.cortisol$age)
round(mean(saliva.cortisol$age), dig = 1)     #7.3
round(median(saliva.cortisol$age), dig = 1)     #6.9
round(min(saliva.cortisol$age), dig = 1)     #2.2
round(max(saliva.cortisol$age), dig = 1)     #23.4

#Sex
table(saliva.cortisol$sex)
# f   m 
# 112 149 

#Maternal social rank
hist(saliva.cortisol$rank)
round(mean(saliva.cortisol$rank), dig = 2)     #0.06
round(median(saliva.cortisol$rank), dig = 2)     #0.11
round(min(saliva.cortisol$rank), dig = 2)     #-1
round(max(saliva.cortisol$rank), dig = 2)     #1

#Litter status
table(saliva.cortisol$litter)
# singleton   dominant    subordinate 
# 85          92          85 

##Number of samples per individual
cort.by.id <- saliva.cortisol %>% group_by(hyena_id) %>% dplyr::summarize(num_samples = length(saliva_sample_id))
nrow(cort.by.id)     # 69 individuals
min(cort.by.id$num_samples)     # 1 sample
max(cort.by.id$num_samples)     # 20 samples
round(mean(cort.by.id$num_samples), 2)     # 3.8 samples
median(cort.by.id$num_samples)     # 2 samples


########## 2.2 Look at data ##########

#Check ICC for hyena_id - keep as random effect
icc(lme4::lmer(data = saliva.cortisol, log_cortisol_ug_dl ~ (1 | hyena_id)))
# Adjusted ICC: 0.479
# Conditional ICC: 0.479

#Check ICC for clan - not important, don't include
icc(lme4::lmer(data = saliva.cortisol, log_cortisol_ug_dl ~ (1 | clan)))
# Adjusted ICC: 0.012
# Conditional ICC: 0.012

#Look at distribution of data
hist(saliva.cortisol$log_cortisol_ug_dl)
hist(saliva.cortisol$date, breaks = 36)
hist(saliva.cortisol$start_time, breaks = 20)
hist(saliva.cortisol$time_lag)
hist(saliva.cortisol$chew_time)
hist(saliva.cortisol$ln2_diff)
hist(saliva.cortisol$cortisol_assay_date, breaks = 36)   #3 rounds of assays
hist(saliva.cortisol$cortisol_assay_diff)
hist(saliva.cortisol$freeze_thaw)
hist(saliva.cortisol$temp_min)
hist(saliva.cortisol$temp_max)
hist(saliva.cortisol$precip)
hist(saliva.cortisol$prey_density)
hist(saliva.cortisol$age)
hist(saliva.cortisol$rank)
hist(saliva.cortisol$number_littermates)
hist(saliva.cortisol$litrank)
table(saliva.cortisol$clan)
table(saliva.cortisol$ampm)
table(saliva.cortisol$sex)
table(saliva.cortisol$litter)

#Look at potential outliers
dotchart(saliva.cortisol$log_cortisol_ug_dl)
dotchart(saliva.cortisol$time_lag)
dotchart(saliva.cortisol$chew_time)
dotchart(saliva.cortisol$ln2_diff)
dotchart(saliva.cortisol$cortisol_assay_diff)
dotchart(saliva.cortisol$freeze_thaw)
dotchart(saliva.cortisol$temp_min)
dotchart(saliva.cortisol$temp_max)
dotchart(saliva.cortisol$precip)
dotchart(saliva.cortisol$prey_density)
dotchart(saliva.cortisol$age)
dotchart(saliva.cortisol$rank)
dotchart(saliva.cortisol$number_littermates)
dotchart(saliva.cortisol$litrank)


######################################################################
##### 3.0 Modeling data #####
######################################################################

########## 3.1 Prepare data ##########

#Standardize variables
saliva.cortisol$time_lag_z <- as.numeric(scale(saliva.cortisol$time_lag, center = TRUE, scale = TRUE))
saliva.cortisol$chew_time_z <- as.numeric(scale(saliva.cortisol$chew_time, center = TRUE, scale = TRUE))
saliva.cortisol$ln2_diff_z <- as.numeric(scale(saliva.cortisol$ln2_diff, center = TRUE, scale = TRUE))
saliva.cortisol$cortisol_assay_diff_z <- as.numeric(scale(saliva.cortisol$cortisol_assay_diff, 
                                                          center = TRUE, scale = TRUE))
saliva.cortisol$freeze_thaw_z <- as.numeric(scale(saliva.cortisol$freeze_thaw, center = TRUE, scale = TRUE))
saliva.cortisol$temp_min_z <- as.numeric(scale(saliva.cortisol$temp_min, center = TRUE, scale = TRUE))
saliva.cortisol$temp_max_z <- as.numeric(scale(saliva.cortisol$temp_max, center = TRUE, scale = TRUE))
saliva.cortisol$precip_z <- as.numeric(scale(saliva.cortisol$precip, center = TRUE, scale = TRUE))
saliva.cortisol$prey_density_z <- as.numeric(scale(saliva.cortisol$prey_density, center = TRUE, scale = TRUE))
saliva.cortisol$age_z <- as.numeric(scale(saliva.cortisol$age, center = TRUE, scale = TRUE))
saliva.cortisol$rank_z <- as.numeric(scale(saliva.cortisol$rank, center = TRUE, scale = TRUE))
saliva.cortisol$number_littermates_z <- as.numeric(scale(saliva.cortisol$number_littermates, 
                                                         center = TRUE, scale = TRUE))
saliva.cortisol$litrank_z <- as.numeric(scale(saliva.cortisol$litrank, center = TRUE, scale = TRUE))

#Look at potential collinearity of independent variables
head(saliva.cortisol[,c(8,9,10,12,15:25,27)])
ggcorr(saliva.cortisol[,c(8,9,10,12,15:25,27)], label = T)
# number_littermates and litrank are highly correlated -> use litter


########## 3.2 Bivariate models ##########

summary(lm(data = saliva.cortisol, log_cortisol_ug_dl ~ ampm))     #significant
summary(lm(data = saliva.cortisol, log_cortisol_ug_dl ~ time_lag_z))     #significant

summary(lm(data = saliva.cortisol, log_cortisol_ug_dl ~ chew_time_z))
round(mean(saliva.cortisol$chew_time, na.rm = T), dig = 1)    #3.5 mins
round(min(saliva.cortisol$chew_time, na.rm = T), dig = 1)     #1
round(max(saliva.cortisol$chew_time, na.rm = T), dig = 1)     #8

summary(lm(data = saliva.cortisol, log_cortisol_ug_dl ~ ln2_diff_z))     #significant
round(mean(saliva.cortisol$ln2_diff, na.rm = T)/60, dig = 1)    #2.1 hours
round(min(saliva.cortisol$ln2_diff, na.rm = T)/60, dig = 1)     #0.6
round(max(saliva.cortisol$ln2_diff, na.rm = T)/60, dig = 1)     #5

summary(lm(data = saliva.cortisol, log_cortisol_ug_dl ~ cortisol_assay_diff_z))     #significant
round(mean(saliva.cortisol$cortisol_assay_diff, na.rm = T), dig = 1)    #8.9 months
round(min(saliva.cortisol$cortisol_assay_diff, na.rm = T), dig = 1)     #2.3
round(max(saliva.cortisol$cortisol_assay_diff, na.rm = T), dig = 1)     #31.6

summary(lm(data = saliva.cortisol, log_cortisol_ug_dl ~ freeze_thaw_z))
round(mean(saliva.cortisol$freeze_thaw, na.rm = T), dig = 1)    #2.2
round(min(saliva.cortisol$freeze_thaw, na.rm = T), dig = 1)     #2
round(max(saliva.cortisol$freeze_thaw, na.rm = T), dig = 1)     #4

summary(lm(data = saliva.cortisol, log_cortisol_ug_dl ~ temp_min_z))
summary(lm(data = saliva.cortisol, log_cortisol_ug_dl ~ temp_max_z))     #significant
summary(lm(data = saliva.cortisol, log_cortisol_ug_dl ~ precip_z)) 
summary(lm(data = saliva.cortisol, log_cortisol_ug_dl ~ prey_density_z))

summary(lm(data = saliva.cortisol, log_cortisol_ug_dl ~ age_z))
summary(lm(data = saliva.cortisol, log_cortisol_ug_dl ~ sex))
summary(lm(data = saliva.cortisol, log_cortisol_ug_dl ~ rank_z))
summary(lm(data = saliva.cortisol, log_cortisol_ug_dl ~ litter))     #significant


########## 3.3 Additive model ##########

#Full model
mod.cort <- lmer(data = saliva.cortisol, log_cortisol_ug_dl ~ ampm + time_lag_z + 
                   chew_time_z + ln2_diff_z + cortisol_assay_diff_z + freeze_thaw_z + 
                   temp_min_z + temp_max_z + precip_z + prey_density_z + 
                   age_z + sex + rank_z + litter + (1|hyena_id))
check_collinearity(mod.cort)     #all below 3
check_model(mod.cort)

#Use AIC criterion to determine best model
options(na.action = "na.fail")
head(saliva.cortisol[,c(8,9,10,12,15:22,25,27)])
tmp <- filter(saliva.cortisol, complete.cases(saliva.cortisol[,c(8,9,10,12,15:22,25,27)]))     #217
mod.cort <- lmer(data = tmp, log_cortisol_ug_dl ~ ampm + time_lag_z + 
                   chew_time_z + ln2_diff_z + cortisol_assay_diff_z + freeze_thaw_z + 
                   temp_min_z + temp_max_z + precip_z + prey_density_z + 
                   age_z + sex + rank_z + litter + (1|hyena_id))
results.dredge <- dredge(mod.cort)
get.models(subset(results.dredge, delta == 0), subset = T)
# log_cortisol_ug_dl ~ ampm + cortisol_assay_diff_z + litter + temp_max_z + time_lag_z + (1 | hyena_id)
importance(subset(results.dredge, delta <= 5 & !nested(.)))    
#                      ampm cortisol_assay_diff_z temp_max_z time_lag_z litter
# Sum of weights:      1.00 1.00                  1.00       1.00       0.63  
# N containing models:    2    2                     2          2          1  
options(na.action = "na.omit")


########## 3.4 Interactive model ##########

#Full model
mod.cort <- lmer(data = saliva.cortisol, log_cortisol_ug_dl ~ ampm + time_lag_z + 
                   chew_time_z + ln2_diff_z + cortisol_assay_diff_z + freeze_thaw_z + 
                   temp_min_z + temp_max_z + precip_z + prey_density_z + 
                   age_z + sex + rank_z + litter + 
                   ampm * time_lag_z + #if slopes different in AM and PM
                   temp_min_z * temp_max_z + #if temperature change matters more than max temp
                   rank_z * prey_density_z + #food availability based on maternal rank
                   litter * prey_density_z + #food availability based on intra-litter rank
                   age_z * sex + #different developmental trajectories based on sex
                   age_z * rank_z + #different developmental trajectories based on maternal rank
                   age_z * litter + #different developmental trajectories based on litter status
                   sex * litter + #following Benhaiem et al. 2013
                   (1|hyena_id))
check_collinearity(mod.cort)
mod.cort <- lmer(data = saliva.cortisol, log_cortisol_ug_dl ~ ampm + time_lag_z + 
                   chew_time_z + ln2_diff_z + cortisol_assay_diff_z + freeze_thaw_z + 
                   temp_min_z + temp_max_z + precip_z + prey_density_z + 
                   age_z + sex + rank_z + litter + 
                   ampm * time_lag_z + #if slopes different in AM and PM
                   temp_min_z * temp_max_z + #if temperature change matters more than max temp
                   rank_z * prey_density_z + #food availability based on maternal rank
                   age_z * sex + #different developmental trajectories based on sex
                   age_z * rank_z + #different developmental trajectories based on maternal rank
                   (1|hyena_id))
check_collinearity(mod.cort)     #all below/at 3
check_model(mod.cort)

#Use AIC criterion to determine best model
options(na.action = "na.fail")
mod.cort <- lmer(data = tmp, log_cortisol_ug_dl ~ ampm + time_lag_z + 
                   chew_time_z + ln2_diff_z + cortisol_assay_diff_z + freeze_thaw_z + 
                   temp_min_z + temp_max_z + precip_z + prey_density_z + 
                   age_z + sex + rank_z + litter + 
                   ampm * time_lag_z + #if slopes different in AM and PM
                   temp_min_z * temp_max_z + #if temperature change matters more than max temp
                   rank_z * prey_density_z + #food availability based on maternal rank
                   age_z * sex + #different developmental trajectories based on sex
                   age_z * rank_z + #different developmental trajectories based on maternal rank
                   (1|hyena_id))
results.dredge <- dredge(mod.cort)
get.models(subset(results.dredge, delta == 0), subset = T)
# log_cortisol_ug_dl ~ ampm + cortisol_assay_diff_z + litter + temp_max_z + time_lag_z + (1 | hyena_id)
importance(subset(results.dredge, delta <= 5 & !nested(.)))    
#                      ampm cortisol_assay_diff_z temp_max_z time_lag_z litter
# Sum of weights:      1.00 1.00                  1.00       1.00       0.63  
# N containing models:    2    2                     2          2          1  
options(na.action = "na.omit")


######################################################################
##### 4.0 Final model #####
######################################################################

########## 4.1 Create dataset ##########

head(saliva.cortisol[,c(8,9,15,17,25)])
saliva.cort.final <- filter(saliva.cortisol, complete.cases(saliva.cortisol[,c(8,9,15,17,25)]))     #256
saliva.cort.final$time_lag_z <- as.numeric(scale(saliva.cort.final$time_lag, center = TRUE, scale = TRUE))
saliva.cort.final$cortisol_assay_diff_z <- as.numeric(scale(saliva.cort.final$cortisol_assay_diff, 
                                                            center = TRUE, scale = TRUE))
saliva.cort.final$temp_max_z <- as.numeric(scale(saliva.cort.final$temp_max, center = TRUE, scale = TRUE))


########## 4.2 Describe dataset ##########

cort <- saliva.cort.final
# 256 samples
length(unique(cort$hyena_id))  
# 69 individuals


########## 4.3 Final model ##########

cortisol <- lme4::lmer(data = saliva.cort.final, log_cortisol_ug_dl ~ ampm + 
                         cortisol_assay_diff_z + litter + temp_max_z + time_lag_z + (1 | hyena_id))
summary(cortisol)
#                        Estimate Std. Error        df t value Pr(>|t|)    
# (Intercept)            -1.17765    0.18372  65.46062  -6.410 1.84e-08 ***
# ampmPM                  0.49461    0.12523 235.62786   3.950 0.000103 ***
# cortisol_assay_diff_z  -0.43975    0.06556 183.04556  -6.708 2.38e-10 ***
# litterdominant         -0.66594    0.22087  53.99281  -3.015 0.003909 ** 
# littersubordinate      -0.19853    0.24337  49.52068  -0.816 0.418552    
# temp_max_z              0.26382    0.06132 248.95383   4.303 2.43e-05 ***
# time_lag_z             -0.20051    0.06124 236.15995  -3.274 0.001219 ** 

summary(glht(cortisol, linfct = mcp(litter = "Tukey")))
#                              Estimate Std. Error z value Pr(>|z|)   
# dominant - singleton == 0     -0.6659     0.2209  -3.015  0.00719 **
# subordinate - singleton == 0  -0.1985     0.2434  -0.816  0.69284   
# subordinate - dominant == 0    0.4674     0.2290   2.041  0.10228   

#Check model
check_collinearity(cortisol)    #all below 3
check_model(cortisol)
simulationOutput <- simulateResiduals(fittedModel = cortisol, n = 250)
plot(simulationOutput)   #KS test: p = 0.00183, Dispersion test: p = 0.608, Outlier test: p = 0.46314
model_performance(cortisol)
# AIC     |     BIC | R2 (cond.) | R2 (marg.) |   ICC |  RMSE | Sigma
# 691.113 | 723.020 |      0.586 |      0.352 | 0.361 | 0.704 | 0.773

# ## Remove outliers - changes significance of litter but not direction
# check_outliers(cortisol)
# tmp <- saliva.cort.final[-c(66),]    #remove outliers
# cortisol <- lmer(data = tmp, log_cortisol_ug_dl ~ ampm + cortisol_assay_diff_z +
#                    litter + temp_max_z + time_lag_z + (1 | hyena_id))
# check_outliers(cortisol)
# summary(cortisol)
# summary(glht(cortisol, linfct = mcp(litter = "Tukey")))
# #                              Estimate Std. Error z value Pr(>|z|)
# # dominant - singleton == 0    -0.55933    0.20504  -2.728   0.0175 *
# # subordinate - singleton == 0 -0.09174    0.22484  -0.408   0.9122
# # subordinate - dominant == 0   0.46759    0.21140   2.212   0.0690 .

#Calculate residuals
saliva.cort.final$residuals <- resid(cortisol)


########## 4.4 Plots ##########

#Litter
myaxis <- saliva.cort.final %>% group_by(litter) %>% 
  dplyr::summarize(count = length(unique(saliva_sample_id))) %>% 
  mutate(myaxis = paste0(litter, "\n", "n=", count))
dat.litter <- ggeffects::ggpredict(cortisol, type = "fe", terms = c("litter"), full.dat.rsa = FALSE)
dat.litter <- left_join(dat.litter, myaxis, by = c("x" = "litter"))
plot.litter <- ggplot(dat.litter, aes(x = myaxis, y = predicted)) +
  geom_point(size = 5) +
  geom_line(data = data.frame(x = c(dat.litter$myaxis, dat.litter$myaxis),
                              predicted = c(dat.litter$predicted, dat.litter$predicted),
                              ci = c(dat.litter$conf.low, dat.litter$conf.high)),
            aes(x = x, y = ci, group= x), size = 1.5)+
  xlab('Litter') +
  ylab("Predicted salivary cortisol (log)")+
  scale_x_discrete(limits=c("singleton\nn=84", "subordinate\nn=85", "dominant\nn=87"))+
  theme_classic() +
  theme(legend.title = element_blank(), legend.text = element_text(size = 12),
        axis.title = element_text(size = 20), axis.text = element_text(size = 16)) + 
  ylim(c(-3.7,0.3))
pdf('09.plot.litter.pdf', width = 7, height = 5)
plot.litter
dev.off()

#Maximum temperature
summary(saliva.cort.final$temp_max)
dat.temp <- ggeffects::ggpredict(cortisol, terms = c("temp_max_z [all]"), full.data = FALSE)
plot.temp <- ggplot(dat.temp, aes(x = dat.temp$x*sd(saliva.cort.final$temp_max, na.rm = T) + 
                                    mean(saliva.cort.final$temp_max, na.rm = T), y = predicted)) + 
  geom_line() + 
  geom_ribbon(aes(x = dat.temp$x*sd(saliva.cort.final$temp_max, na.rm = T) + mean(saliva.cort.final$temp_max, na.rm = T), 
                  ymin = conf.low, ymax = conf.high), alpha = 0.2, inherit.aes = FALSE) +
  xlab(expression(paste("Max temperature (",degree,"C)"))) +
  ylab("Predicted salivary cortisol (log)")+
  theme_classic() + 
  theme(legend.title = element_blank(), legend.text = element_text(size = 12), 
        axis.title = element_text(size = 20), axis.text = element_text(size = 16)) + 
  ylim(c(-3.7,0.3))
pdf('09.plot.temp.pdf', width = 7, height = 5)
plot.temp
dev.off()

#Circadian rhythms
summary(saliva.cort.final$time_lag)
dat.time_lag <- ggeffects::ggpredict(cortisol, terms = c("time_lag_z [all]", "ampm"), full.data = FALSE)
dat.time_lag.am <- filter(dat.time_lag, group == "AM")
dat.time_lag.pm <- filter(dat.time_lag, group == "PM")
plot.time_lag.am <- ggplot(dat.time_lag.am, aes(x = dat.time_lag.am$x*sd(saliva.cort.final$time_lag, na.rm = T) + 
                                                  mean(saliva.cort.final$time_lag, na.rm = T), y = predicted)) + 
  geom_line() + 
  geom_ribbon(aes(x = dat.time_lag.am$x*sd(saliva.cort.final$time_lag, na.rm = T) + 
                    mean(saliva.cort.final$time_lag, na.rm = T), 
                  ymin = conf.low, ymax = conf.high), alpha = 0.2, inherit.aes = FALSE) +
  xlab('Time to sunrise (mins)') +
  ylab("Predicted salivary cortisol (log)")+
  theme_classic() + 
  theme(legend.title = element_blank(), legend.text = element_text(size = 12), 
        axis.title = element_text(size = 20), axis.text = element_text(size = 16)) + 
  ylim(c(-3.7,0.3))
plot.time_lag.pm <- ggplot(dat.time_lag.pm, aes(x = dat.time_lag.pm$x*sd(saliva.cort.final$time_lag, na.rm = T) + 
                                                  mean(saliva.cort.final$time_lag, na.rm = T), y = predicted)) + 
  geom_line() + 
  geom_ribbon(aes(x = dat.time_lag.pm$x*sd(saliva.cort.final$time_lag, na.rm = T) + 
                    mean(saliva.cort.final$time_lag, na.rm = T), 
                  ymin = conf.low, ymax = conf.high), alpha = 0.2, inherit.aes = FALSE) +
  xlab('Time to sunset (mins)') +
  ylab("Predicted salivary cortisol (log)")+
  theme_classic() + 
  theme(legend.title = element_blank(), legend.text = element_text(size = 12), 
        axis.title = element_text(size = 20), axis.text = element_text(size = 16), 
        axis.title.y = element_blank()) + 
  ylim(c(-3.7,0.3))
plot.time_lag <- plot.time_lag.am + plot.time_lag.pm
pdf('09.plot.time_lag.pdf', width = 10, height = 5)
plot.time_lag
dev.off()

#Time collection to assay
summary(saliva.cort.final$cortisol_assay_diff)
dat.assaydiff <- ggeffects::ggpredict(cortisol, terms = c("cortisol_assay_diff_z [all]"), full.data = FALSE)
plot.assaydiff <- ggplot(dat.assaydiff, aes(x = dat.assaydiff$x*sd(saliva.cort.final$cortisol_assay_diff, na.rm = T) + 
                                              mean(saliva.cort.final$cortisol_assay_diff, na.rm = T), y = predicted)) + 
  geom_line() + 
  geom_ribbon(aes(x = dat.assaydiff$x*sd(saliva.cort.final$cortisol_assay_diff, na.rm = T) + 
                    mean(saliva.cort.final$cortisol_assay_diff, na.rm = T), 
                  ymin = conf.low, ymax = conf.high), alpha = 0.2, inherit.aes = FALSE) +
  xlab('Collection to assay (mos)') +
  ylab("Predicted salivary cortisol (log)")+
  theme_classic() + 
  theme(legend.title = element_blank(), legend.text = element_text(size = 12), 
        axis.title = element_text(size = 20), axis.text = element_text(size = 16), 
        axis.title.y = element_blank()) + 
  ylim(c(-3.7,0.3))
pdf('09.plot.assaydiff.pdf', width = 7, height = 5)
plot.assaydiff
dev.off()

layout <- 
  "AB
   CD"
plots.all <- plot.temp + plot.assaydiff + plot.time_lag.am + plot.time_lag.pm + plot_layout(design = layout)
pdf('09.plots.all.pdf', width = 7, height = 5)
plots.all
dev.off()

#Re-order model terms for plotting
cortisol <- lmer(data = saliva.cort.final, log_cortisol_ug_dl ~ litter + temp_max_z + 
                   ampm + time_lag_z + cortisol_assay_diff_z + (1 | hyena_id))
summary(cortisol)

#Table
tab_model(cortisol, show.se = T, show.ci = F, show.re.var = F, show.intercept = F,
          pred.labels = c("Litter [dominant]",
                          "Litter [subordinate]",
                          "Daily maximum temperature",
                          "Time of day of collection [PM]",
                          "Time to sunrise/set",
                          "Time between collection and assay"),
          dv.labels = c("Log salivary cortisol (ug/dL)"),
          string.se = "SE", file = "09.table.cort.model.doc")

#Model plot
set_theme(base = theme_classic(), axis.textcolor = "black", axis.title.color = "black", 
          axis.textsize.y = 1.5, axis.textsize.x = 1.2, axis.title.size = 1.7)
plot.model <- plot_model(cortisol, type = "est", transform = NULL,
                         axis.labels = c("Time between collection and assay",
                                         "Time to sunrise/set",
                                         "Time of day of collection [PM]",
                                         "Daily maximum temperature",
                                         "Litter [subordinate]",
                                         "Litter [dominant]"),
                         vline.color = "black", title = "", dot.size = 4.5, line.size = 1.5,
                         show.values = TRUE, show.p = TRUE, digits = 2, value.offset = 0.3, 
                         value.size = 6, colors = viridis_2, axis.lim = c(-2,2))
pdf('09.plot.cort.model.pdf', width = 7, height = 5)
plot.model
dev.off()


########## 4.5 Calculate repeatability ##########

#Only with random effect
rep_c <- rpt(log_cortisol_ug_dl ~ (1 | hyena_id), grname = "hyena_id", data = saliva.cort.final)
print(rep_c)
#Repeatability for hyena_id: R = 0.48, SE = 0.077, CI = [0.306, 0.614], P = 1.54e-17

#Including fixed effects
rep_c2 <- rpt(log_cortisol_ug_dl ~ ampm + cortisol_assay_diff_z + litter + 
                temp_max_z + time_lag_z + (1 | hyena_id),
              grname = "hyena_id", data = saliva.cort.final)
print(rep_c2)
#Repeatability for hyena_id: R = 0.361, SE = 0.08, CI = [0.202, 0.512], P = 9.45e-09


######################################################################
##### 5.0 Save final data #####
######################################################################

saliva.cortisol <- saliva.cortisol[,c(1:27)]
saliva.cortisol <- left_join(saliva.cortisol, saliva.cort.final[,c("saliva_sample_id", "residuals")])

save(file = "10.saliva_cortisol.Rdata", list = c("saliva.cortisol"))





