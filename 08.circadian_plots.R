################################################################################
######### 08: Plotting raw saliva cortisol vs. time ########
################################################################################

######################################################################
##### 1.0 Loading data #####
######################################################################

########## 1.1 Set working directory & download packages ##########

rm(list = ls())
setwd("~/Documents/R/Saliva")
options(stringsAsFactors = F)
library(patchwork)
library(tidyverse)

#Set plot colors using viridis
viridis_2 <- viridis(7)[-c(1,3,5,6,7)]
viridis_3 <- viridis(7)[-c(2,4,6,7)]

#Load data
load("07.saliva_behavior.Rdata")

#Reformat time
saliva.cortisol$start_time <- as.POSIXct(format(saliva.cortisol$start_time, format = "%H:%M"), format = "%H:%M")
saliva.cortisol$stop_time <- as.POSIXct(format(saliva.cortisol$stop_time, format = "%H:%M"), format = "%H:%M")


######################################################################
##### 2.0 Circadian boxplots #####
######################################################################

########## 2.1 AM and PM samples (raw data in facet plot) ##########

#AMPM (facet)
cort.time.facet <- ggplot(saliva.cortisol, aes(x = start_time, y = log_cortisol_ug_dl)) +
  geom_smooth(method = "lm", se = FALSE, size = 2, color = "darkgray") + 
  geom_point() + 
  facet_wrap(vars(ampm), scales = "free_x", strip.position = "top") + 
  xlab('Time') +
  ylab("Salivary cortisol (log)") +
  theme_classic() +
  theme(legend.title = element_blank(), legend.text = element_text(size = 12),
        axis.title = element_text(size = 20), axis.text = element_text(size = 16), 
        strip.text.x = element_text(size = 12))
pdf('08.plot.cort.time.facet.pdf', width = 7, height = 5)
cort.time.facet
dev.off()


########## 2.2 Boxplot by hour of cortisol concentrations and maternal den activity ##########

#Sunrise - sunset:   6:36 - 18:43 (see script 04 for these calculations)
# round(36/60, dig = 2) = 0.6  -> sunrise at hour 6.6
# round(43/60, dig = 2) = 0.72 -> sunset at hour 18.72

#Boxplot by hour
saliva.cortisol$hour <- format(saliva.cortisol$start_time, "%H")
cort.hour <- ggplot(saliva.cortisol, aes(x = hour, y = log_cortisol_ug_dl)) +
  geom_rect(xmin = 0, xmax = (6.6+1), ymin = -5, ymax = 3, fill = "gray80") + 
  geom_rect(xmin = (18.72+1), xmax = 25, ymin = -5, ymax = 3, fill = "gray80") + 
  geom_boxplot(position = position_nudge(0.5)) + 
  xlab('Hour') +
  ylab("Juvenile salivary cortisol (log)") +
  theme_classic() +
  theme(legend.title = element_blank(), legend.text = element_text(size = 12),
        axis.title = element_text(size = 20), axis.text = element_text(size = 16)) + 
  scale_x_discrete(limits = c("00","01","02","03","04","05","06","07","08","09","10","11","12",
                              "13","14","15","16","17","18","19","20","21","22","23"))


########## 2.3 Boxplot by hour of maternal den activity (GPS collar data 2012-14) ##########

##### CODE AND DATA FROM JULIA GREENBERG (Greenberg 2017 dissertation) #####

propperhour.all <- read.csv("00.raw_data/proportionperhourwithinfo_cleaned.csv")

# exclude moms with cubs older than 1 year
propperhour <- subset(propperhour.all, cubage < 365 & park == "C")
propperhour$LMT_Date <- as.Date(propperhour$LMT_Date, format = "%m/%d/%y")
summary(propperhour$LMT_Date)
length(unique(propperhour$hyena))     #10

# graph comparing Serena by hour. Make summary table; graph mean.
error <- plyr::ddply(propperhour, c("hour"), summarize,
                     N    = length(presentbinary),
                     mean = mean(presentbinary),
                     sd   = sd(presentbinary),
                     se   = sd / sqrt(N) )
error$hour <- as.character(error$hour)
error[error$hour == "0",]$hour <- "00"
error[error$hour == "1",]$hour <- "01"
error[error$hour == "2",]$hour <- "02"
error[error$hour == "3",]$hour <- "03"
error[error$hour == "4",]$hour <- "04"
error[error$hour == "5",]$hour <- "05"
error[error$hour == "6",]$hour <- "06"
error[error$hour == "7",]$hour <- "07"
error[error$hour == "8",]$hour <- "08"
error[error$hour == "9",]$hour <- "09"

#Plot modified by TMM to add daylight/nighttime shading
mom.den <- ggplot(error, aes(x = hour, y = mean)) + 
  geom_rect(xmin = 0, xmax = (6.6+1), ymin = -5, ymax = 3, fill = "gray80") + 
  geom_rect(xmin = (18.72+1), xmax = 25, ymin = -5, ymax = 3, fill = "gray80") + 
  geom_point(stat = "identity", position = position_nudge(0.5), size = 3) + 
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se), width = 0, position = position_nudge(0.5), size = 1.25) +
  xlab('Hour') +
  ylab("Proportion of GPS fixes at den (mothers)") +
  theme_classic() + 
  theme(legend.title = element_blank(), legend.text = element_text(size = 12), 
        axis.title = element_text(size = 20), axis.text = element_text(size = 16)) + 
  ylim(0.2,0.5) + 
  scale_x_discrete(limits = c("00","01","02","03","04","05","06","07","08","09","10","11","12",
                              "13","14","15","16","17","18","19","20","21","22","23"))

##### END CODE FROM JULIA GREENBERG #####

rm(error)
rm(propperhour)
rm(propperhour.all)


########## 2.4 Combine boxplots by hour ##########

layout <- "A
           B"
cort.by.hour <- cort.hour + mom.den + plot_layout(design = layout)
pdf('08.plot.cort.by.hour.pdf', width = 10, height = 10)
cort.by.hour
dev.off()


