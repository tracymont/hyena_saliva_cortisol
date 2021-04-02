################################################################################
######### 11: Cortisol and behavior ########
################################################################################

######################################################################
##### 1.0 Loading data #####
######################################################################

########## 1.1 Set working directory & download packages ##########

rm(list = ls())
setwd("~/Documents/R/Saliva")
options(stringsAsFactors = F)
library(lme4)
library(lmerTest)
library(sjPlot)
library(viridis)
library(tidyverse)

load("07.saliva_behavior.Rdata")
load("10.saliva_cortisol.Rdata")

#Set plot colors using viridis
viridis_2 <- viridis(7)[-c(1,3,5,6,7)]
viridis_3 <- viridis(7)[-c(2,4,6,7)]

#Fix times
saliva.behavior$time <- as.POSIXct(format(saliva.behavior$time, format = "%H:%M"), format = "%H:%M")
saliva.cortisol$start_time <- as.POSIXct(format(saliva.cortisol$start_time, format = "%H:%M"), format = "%H:%M")


######################################################################
##### 2.0 Analysing data #####
######################################################################

########## 2.1 Assign behavior to hormone samples ##########

#Filter to hyenas in saliva dataset
saliva.behavior <- filter(saliva.behavior, date %in% saliva.cortisol$date & id %in% saliva.cortisol$hyena_id)
summary(as.factor(saliva.behavior$behav))

#Assign behavior
saliva.cortisol.behav <- saliva.cortisol
saliva.cortisol.behav$behav.all <- ""
for(i in 1:nrow(saliva.cortisol.behav)){
  clan.i <- saliva.cortisol.behav$clan[i]
  date.i <- saliva.cortisol.behav$date[i]
  time.i <- saliva.cortisol.behav$start_time[i]
  id.i <- saliva.cortisol.behav$hyena_id[i]
  #Include behaviors that occurred 15-45 min before saliva sample
  behav.i <- filter(saliva.behavior, clan == clan.i & date == date.i & id == id.i & 
                      time <= (time.i - 15*60) & time >= (time.i - 45*60))$behav
  if(length(behav.i) == 1){
    saliva.cortisol.behav$behav.all[i] <- behav.i
  }
  if(length(behav.i) > 1){
    saliva.cortisol.behav$behav.all[i] <- paste0(behav.i, collapse = ",")
  }
}


########## 2.2 Assign behavior category ##########

summary(as.factor(saliva.behavior$behav))

saliva.cortisol.behav$behav.cat <- NA

saliva.cortisol.behav[grepl("aggressor", saliva.cortisol.behav$behav.all) &
                        !grepl("play", saliva.cortisol.behav$behav.all) &
                        !grepl("recipient", saliva.cortisol.behav$behav.all),]$behav.cat <- "2aggressor"

saliva.cortisol.behav[grepl("recipient", saliva.cortisol.behav$behav.all) &
                        !grepl("play", saliva.cortisol.behav$behav.all) &
                        !grepl("aggressor", saliva.cortisol.behav$behav.all),]$behav.cat <- "3recipient"

saliva.cortisol.behav[grepl("play", saliva.cortisol.behav$behav.all) &
                        !grepl("aggressor", saliva.cortisol.behav$behav.all) &
                        !grepl("recipient", saliva.cortisol.behav$behav.all),]$behav.cat <- "4play"

saliva.cortisol.behav[grepl("active", saliva.cortisol.behav$behav.all) & 
                        !grepl("aggressor", saliva.cortisol.behav$behav.all) &
                        !grepl("play", saliva.cortisol.behav$behav.all) &
                        !grepl("recipient", saliva.cortisol.behav$behav.all),]$behav.cat <- "1active"

summary(as.factor(saliva.cortisol.behav$behav.cat))


########## 2.3 Model residuals ##########

#Model residuals from cortisol model in script 09.
# Model: log_cortisol_ug_dl ~ ampm + cortisol_assay_diff_z + litter + temp_max_z + time_lag_z + (1 | hyena_id)

behavior.dat <- filter(saliva.cortisol.behav, !is.na(behav.cat) & !is.na(residuals))
behavior.dat$behav.cat <- as.factor(behavior.dat$behav.cat)
behavior.dat$behav.cat <- relevel(behavior.dat$behav.cat, ref = "1active")
mod.cort <- lm(data = behavior.dat, residuals ~ behav.cat)
summary(mod.cort)
#                      Estimate Std. Error t value Pr(>|t|)  
# (Intercept)         -0.10771    0.09524  -1.131   0.2613  
# behav.cat2aggressor -0.45230    0.21822  -2.073   0.0413 *
# behav.cat3recipient  0.21711    0.19977   1.087   0.2803  
# behav.cat4play       0.22681    0.24590   0.922   0.3590  

aov.cort <- aov(mod.cort)
TukeyHSD(aov.cort)
# Tukey multiple comparisons of means
# 95% family-wise confidence level
#                               diff         lwr       upr     p adj
# 2aggressor-1active    -0.452302120 -1.02443750 0.1198333 0.1706123
# 3recipient-1active     0.217105510 -0.30667037 0.7408814 0.6985636
# 4play-1active          0.226814845 -0.41790869 0.8715384 0.7929675
# 3recipient-2aggressor  0.669407629 -0.02122942 1.3600447 0.0609346
# 4play-2aggressor       0.679116965 -0.10720791 1.4654418 0.1148699
# 4play-3recipient       0.009709336 -0.74216105 0.7615797 0.9999859


########## 2.4 Plots ##########

#Plot
behavior.dat %>% group_by(behav.cat) %>% summarize(count = length(unique(saliva_sample_id)))
boxplot.behav <- ggplot(behavior.dat, aes(x = behav.cat, y = residuals)) +
  geom_boxplot() + 
  xlab('Behavior 15-45 min prior to sample') +
  ylab("Cortisol residuals") +
  theme_classic() +
  theme(legend.title = element_blank(), legend.text = element_text(size = 12),
        axis.title = element_text(size = 20), axis.text = element_text(size = 16)) + 
  scale_x_discrete(limits = c("1active", "2aggressor", "3recipient", "4play"),
                   labels = c("other\nn=51", "aggressor\nn=12", "recipient\nn=15", "play\nn=9")) +
  ylim(c(-2,2))
pdf('11.plot.cort.behav.pdf', width = 7, height = 5)
boxplot.behav
dev.off()

#Table
tab_model(mod.cort, show.se = T, show.ci = F, show.intercept = F, show.reflvl = TRUE,
          pred.labels = c("Aggressor", "Recipient", "Play"),
          dv.labels = c("Cortisol residuals"),
          string.se = "SE")#, file = "11.table.cort.behav.doc")

#Model plot
set_theme(base = theme_classic(), axis.textcolor = "black", axis.title.color = "black", 
          axis.textsize.y = 1.5, axis.textsize.x = 1.2, axis.title.size = 1.7)
plot.model <- plot_model(mod.cort, type = "est", transform = NULL,
                         axis.labels = c("Play", "Recipient", "Aggressor"),
                         vline.color = "black", title = "", dot.size = 4.5, line.size = 1.5,
                         show.values = TRUE, show.p = TRUE, digits = 2, value.offset = 0.3, 
                         value.size = 6, colors = viridis_2, axis.lim = c(-2,2))
pdf('11.plot.behav.model.pdf', width = 7, height = 5)
plot.model
dev.off()


