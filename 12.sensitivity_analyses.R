################################################################################
######### 12: Sensitivity analyses requested by reviewers ########
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
##### 2.0 Littermates: is the pattern the same within vs. between litters #####
######################################################################

########## 2.1 Create dataset ##########

head(saliva.cortisol[,c(8,9,15,17,25,26)])
twins <- filter(saliva.cortisol, complete.cases(saliva.cortisol[,c(8,9,15,17,25,26)]))     #114

#Remove twin litters with only one littermate in dataset for comparison
twins %>% group_by(litter_id) %>% summarize(num_twins = length(unique(hyena_id))) %>% ungroup(.) %>% filter(num_twins != 2) 
twins <- filter(twins, litter_id != "app-momo" & litter_id != "goli-orb")

#Standardize variables
twins$litter_status <- as.factor(as.character(twins$litter_status))
twins$litter_status <- relevel(twins$litter_status, ref = "subordinate")
twins$time_lag_z <- as.numeric(scale(twins$time_lag, center = TRUE, scale = TRUE))
twins$cortisol_assay_diff_z <- as.numeric(scale(twins$cortisol_assay_diff, 
                                                center = TRUE, scale = TRUE))
twins$temp_max_z <- as.numeric(scale(twins$temp_max, center = TRUE, scale = TRUE))

#Describe dataset
cort <- twins
# 110 samples
length(unique(cort$hyena_id))  
# 28 individuals
table(twins$litter_status)
# subordinate    dominant 
# 75          35 


########## 2.2 Final model with litter status ##########

cortisol.twins <- lmer(data = twins, log_cortisol_ug_dl ~ ampm + cortisol_assay_diff_z + 
                         litter_status + temp_max_z + time_lag_z + (1 | hyena_id))
summary(cortisol.twins)
#                        Estimate Std. Error        df t value Pr(>|t|)    
# ampmPM                  0.78631    0.19418 100.57031   4.049 0.000101 ***
# cortisol_assay_diff_z  -0.28013    0.10058  54.57705  -2.785 0.007340 ** 
# litter_statusdominant  -0.38706    0.28181  24.25793  -1.374 0.182158    
# temp_max_z              0.26099    0.08401 103.35855   3.107 0.002443 ** 
# time_lag_z             -0.09833    0.09972 101.72051  -0.986 0.326437    

#Check model
check_collinearity(cortisol.twins)    #all below 3
check_model(cortisol.twins)
simulationOutput <- simulateResiduals(fittedModel = cortisol.twins, n = 250)
plot(simulationOutput)   #KS test: p = 0.02564, Dispersion test: p = 0.744, Outlier test: p = 0.58522
model_performance(cortisol.twins)
# AIC     |     BIC | R2 (cond.) | R2 (marg.) |   ICC |  RMSE | Sigma
# 298.180 | 319.784 |      0.549 |      0.320 | 0.336 | 0.686 | 0.754

#Re-order model terms for plotting
cortisol.twins <- lmer(data = twins, log_cortisol_ug_dl ~ litter_status + temp_max_z +
                           ampm + time_lag_z + cortisol_assay_diff_z + (1 | hyena_id))
summary(cortisol.twins)

#Table
tab_model(cortisol.twins, show.se = T, show.ci = F, show.re.var = F, show.intercept = F,
          pred.labels = c("Litter status [dominant]",
                          "Daily maximum temperature",
                          "Time of day of collection [PM]",
                          "Time to sunrise/set",
                          "Time between collection and assay"),
          dv.labels = c("Log salivary cortisol (ug/dL)"),
          string.se = "SE", file = "12.table.twins.model.doc")

#Model plot
set_theme(base = theme_classic(), axis.textcolor = "black", axis.title.color = "black", 
          axis.textsize.y = 1.5, axis.textsize.x = 1.2, axis.title.size = 1.7)
plot.model <- plot_model(cortisol.twins, type = "est", transform = NULL,
                         axis.labels = c("Time between collection and assay",
                                         "Time to sunrise/set",
                                         "Time of day of collection [PM]",
                                         "Daily maximum temperature",
                                         "Litter status [dominant]"),
                         vline.color = "black", title = "", dot.size = 4.5, line.size = 1.5,
                         show.values = TRUE, show.p = TRUE, digits = 2, value.offset = 0.3, 
                         value.size = 6, colors = viridis_2, axis.lim = c(-2,2))
pdf('12.plot.twins.model.pdf', width = 7, height = 5)
plot.model
dev.off()


########## 2.3 Final model without litter status ##########

cortisol.twins2 <- lmer(data = twins, log_cortisol_ug_dl ~ ampm + 
                         cortisol_assay_diff_z + temp_max_z + time_lag_z + (1 | hyena_id))
summary(cortisol.twins2)
#                        Estimate Std. Error        df t value Pr(>|t|)    
# ampmPM                  0.76948    0.19371 101.29582   3.972 0.000133 ***
# cortisol_assay_diff_z  -0.29429    0.10238  57.95063  -2.875 0.005650 ** 
# temp_max_z              0.24986    0.08407 103.71751   2.972 0.003676 ** 
# time_lag_z             -0.10751    0.09947 102.20813  -1.081 0.282298    

#Calculate random effects
twins.ranef <- as.data.frame(ranef(cortisol.twins2, condVar = TRUE, whichel = "hyena_id"))
twins.ranef <- twins.ranef[,3:5]
colnames(twins.ranef)[1] <- "hyena_id"
colnames(twins.ranef)[2] <- "ranef"
colnames(twins.ranef)[3] <- "ranef.sd"

#Random effects
twins.ranef <- left_join(twins.ranef, twins[,c("hyena_id", "litter_id", "litter_status")])
twins.ranef <- unique(twins.ranef)
twins.ranef <- pivot_wider(twins.ranef, id_cols = "litter_id", names_from = "litter_status", values_from = "ranef")
twins.ranef$direction <- ifelse(twins.ranef$dominant < twins.ranef$subordinate,
                                "subordinate has higher cortisol", "dominant has higher cortisol")

plot.twins.ranef <- ggparcoord(twins.ranef, columns = 3:2, groupColumn = 4, showPoints = T) + 
  xlab('Litter status') +
  ylab("Individual random effect estimate")+
  theme_classic() +
  theme(legend.title = element_blank(), legend.text = element_text(size = 12),
        axis.title = element_text(size = 20), axis.text = element_text(size = 16)) + 
  scale_color_manual(values = viridis_2)
pdf('12.plot.twins.ranef.pdf', width = 7, height = 5)
plot.twins.ranef
dev.off()


######################################################################
##### 3.0 Weaning #####
######################################################################

########## 3.1 Create dataset ##########

head(saliva.cortisol[,c(8,9,15,17,25,27)])
weaning.data <- filter(saliva.cortisol, complete.cases(saliva.cortisol[,c(8,9,15,17,25,27)]))     #248
weaning.data$time_lag_z <- as.numeric(scale(weaning.data$time_lag, center = TRUE, scale = TRUE))
weaning.data$cortisol_assay_diff_z <- as.numeric(scale(weaning.data$cortisol_assay_diff, 
                                                       center = TRUE, scale = TRUE))
weaning.data$temp_max_z <- as.numeric(scale(weaning.data$temp_max, center = TRUE, scale = TRUE))
weaning.data$weaning_status <- as.factor(weaning.data$weaning_status)
weaning.data$weaning_status <- relevel(weaning.data$weaning_status, ref = "nursing")


########## 3.2 Describe dataset ##########

cort <- weaning.data
# 248 samples
length(unique(cort$hyena_id))  
# 66 individuals
table(weaning.data$weaning_status)
# nursing  weaned 
# 199      49 


########## 3.3 Final model ##########

cortisol.weaning <- lmer(data = weaning.data, log_cortisol_ug_dl ~ ampm + 
                   cortisol_assay_diff_z + litter_status + temp_max_z + time_lag_z + 
                   weaning_status + (1 | hyena_id))
summary(cortisol.weaning)
#                           Estimate Std. Error        df t value Pr(>|t|)    
# ampmPM                     0.49061    0.12653 223.24049   3.877 0.000139 ***
# cortisol_assay_diff_z     -0.44379    0.06919 150.01216  -6.414 1.75e-09 ***
# litter_statusdominant     -0.73142    0.23664  49.93198  -3.091 0.003260 ** 
# litter_statussubordinate  -0.27776    0.25875  46.88565  -1.073 0.288560    
# temp_max_z                 0.25847    0.06209 239.83326   4.163 4.38e-05 ***
# time_lag_z                -0.18309    0.06187 225.06372  -2.959 0.003415 ** 
# weaning_statusweaned       0.35500    0.16424 227.05087   2.161 0.031704 *  
  
summary(glht(cortisol.weaning, linfct = mcp(litter_status = "Tukey")))
#                              Estimate Std. Error z value Pr(>|z|)   
# dominant - singleton == 0     -0.7314     0.2366  -3.091  0.00555 **
# subordinate - singleton == 0  -0.2778     0.2588  -1.073  0.52996   
# subordinate - dominant == 0    0.4537     0.2362   1.921  0.13239   

#Check model
check_collinearity(cortisol.weaning)    #all below 3
check_model(cortisol.weaning)
simulationOutput <- simulateResiduals(fittedModel = cortisol.weaning, n = 250)
plot(simulationOutput)   #KS test: p = 0.02483, Dispersion test: p = 0.648, Outlier test: p = 0.04981
model_performance(cortisol.weaning)
# AIC     |     BIC | R2 (cond.) | R2 (marg.) |   ICC |  RMSE | Sigma
# 669.378 | 704.512 |      0.616 |      0.368 | 0.391 | 0.689 | 0.761


########## 3.4 Plots ##########

#Re-order model terms for plotting
cortisol.weaning <- lmer(data = weaning.data, log_cortisol_ug_dl ~ litter_status + temp_max_z + 
                   ampm + time_lag_z + cortisol_assay_diff_z + weaning_status + (1 | hyena_id))
summary(cortisol.weaning)

#Table
tab_model(cortisol.weaning, show.se = T, show.ci = F, show.re.var = F, show.intercept = F,
          pred.labels = c("Litter status [dominant]",
                          "Litter status [subordinate]",
                          "Daily maximum temperature",
                          "Time of day of collection [PM]",
                          "Time to sunrise/set",
                          "Time between collection and assay", 
                          "Weaning status [weaned]"),
          dv.labels = c("Log salivary cortisol (ug/dL)"),
          string.se = "SE", file = "12.table.weaning.model.doc")

#Model plot
set_theme(base = theme_classic(), axis.textcolor = "black", axis.title.color = "black", 
          axis.textsize.y = 1.5, axis.textsize.x = 1.2, axis.title.size = 1.7)
plot.model <- plot_model(cortisol.weaning, type = "est", transform = NULL,
                         axis.labels = c("Weaning status [weaned]", 
                                         "Time between collection and assay",
                                         "Time to sunrise/set",
                                         "Time of day of collection [PM]",
                                         "Daily maximum temperature",
                                         "Litter status [subordinate]",
                                         "Litter status [dominant]"),
                         vline.color = "black", title = "", dot.size = 4.5, line.size = 1.5,
                         show.values = TRUE, show.p = TRUE, digits = 2, value.offset = 0.3, 
                         value.size = 6, colors = viridis_2, axis.lim = c(-2,2))
pdf('12.plot.weaning.model.pdf', width = 7, height = 5)
plot.model
dev.off()


