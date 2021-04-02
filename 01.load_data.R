################################################################################
##########        01: Loading saliva hormones & behavior data          #########
################################################################################

########## 1.1 Set working directory & download packages ##########

rm(list = ls())
setwd("~/Documents/R/Saliva")
options(stringsAsFactors = F)
library(tidyverse)


########## 1.2 Saliva behavior data ##########

saliva.behav <- read.csv("00.raw_data/SalivaBehaviorEntry_all.csv")
saliva.behav$sampleID <- as.numeric(saliva.behav$sampleID)
saliva.behav$clan <- gsub("HZ", "happy.zebra", saliva.behav$clan)
saliva.behav$clan <- gsub("S", "serena.s", saliva.behav$clan)
saliva.behav$clan <- gsub("N", "serena.n", saliva.behav$clan)
saliva.behav$clan <- as.factor(gsub(" ", "", saliva.behav$clan))
saliva.behav$date <- as.Date(saliva.behav$date, format = "%d-%b-%y")
saliva.behav$hyenaID <- gsub(" ", "", saliva.behav$hyenaID)
saliva.behav$time <- as.POSIXct(saliva.behav$time, format = "%H:%M")
saliva.behav$meridian <- format(saliva.behav$time, '%p')


########## 1.3 Saliva hormone data ##########

#Repository
saliva.repository <- read.csv("00.raw_data/tblSalivaRepository.csv")
saliva.repository$saliva_sample_id <- as.numeric(saliva.repository$saliva_sample_id)
saliva.repository$clan <- as.factor(saliva.repository$clan)
saliva.repository$hyena_id <- gsub(" ", "", saliva.repository$hyena_id)
saliva.repository$date <- as.Date(saliva.repository$date, format = "%d-%b-%y")
saliva.repository$start_time <- as.POSIXct(saliva.repository$start_time, format = "%H:%M")
saliva.repository$stop_time <- as.POSIXct(saliva.repository$stop_time, format = "%H:%M")
saliva.repository$ln2_time <- as.POSIXct(saliva.repository$ln2_time, format = "%H:%M")
saliva.repository$repeated <- as.logical(saliva.repository$repeated)

#All samples
nrow(saliva.repository)     #302
length(unique(saliva.repository$hyena_id))     #81

#Hormone data
##If updated hormone data, must change column names in Excel prior to loading
saliva.hormones <- read.csv("00.raw_data/tblSalivaHormones.csv")
saliva.hormones$saliva_sample_id <- as.numeric(saliva.hormones$saliva_sample_id)
saliva.hormones$cortisol_assay_date <- as.Date(saliva.hormones$cortisol_assay_date, 
                                               format = "%Y-%m-%d")


########## 1.4 Final dataset ##########

#Join tables
saliva.final <- left_join(saliva.repository, saliva.hormones, by = "saliva_sample_id")
saliva.final <- saliva.final[,1:11]

rm(saliva.repository)
rm(saliva.hormones)

save(file = "02.cleaned_data.Rdata", list = c("saliva.final", "saliva.behav"))





