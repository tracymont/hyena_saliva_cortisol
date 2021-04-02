################################################################################
######### 04: Creating saliva.cortisol datasets ########
################################################################################

######################################################################
##### 1.0 Loading data #####
######################################################################

########## 1.1 Set working directory & download packages ##########

rm(list = ls())
setwd("~/Documents/R/Saliva")
options(stringsAsFactors = F)
library(tidyverse)
hyenadata::update_tables("1.2.88")
library(hyenadata)

load("02.cleaned_data.Rdata")

#Remove Talek samples (only 3 - different ecological context than rest of samples)
saliva.final <- filter(saliva.final, clan != "talek.w")

#Add freeze-thaw data (estimated)
saliva.freezethaw <- read.csv("00.raw_data/saliva_sample_freezethaw.csv")
saliva.freezethaw <- unique(saliva.freezethaw[,c(2:3)])
saliva.freezethaw$assay_date <- as.Date(saliva.freezethaw$assay_date, format = "%d-%b-%y")
saliva.freezethaw <- saliva.freezethaw %>% group_by(sample_id) %>% 
  mutate(freeze_thaw = as.numeric(min_rank(assay_date) + 1))  
#add 1 due to transfer from cryotubes to microcentrifuge tubes [removal of Kimbo]


########## 1.2 hyenadata tables ##########

#tblLifeHistory
data(tblLifeHistory)
tblLifeHistory$id <- gsub(" ", "", tblLifeHistory$id)
tblLifeHistory$event_code <- as.factor(tblLifeHistory$event_code)
tblLifeHistory$error <- as.numeric(tblLifeHistory$error)    
tblLifeHistory$event_status <- as.factor(tblLifeHistory$event_status)
tblLifeHistory <- tblLifeHistory[,1:7]
tblLifeHistory <- filter(tblLifeHistory, !is.na(id) & !is.na(event_code))   #remove 0

#tblHyenas
data("tblHyenas")
tblHyenas$id <- gsub(" ", "", tblHyenas$id)
tblHyenas$sex <- as.factor(tblHyenas$sex)
tblHyenas$status <- as.factor(tblHyenas$status)
tblHyenas$mom <- gsub(" ", "", tblHyenas$mom)
tblHyenas$dad <- gsub(" ", "", tblHyenas$dad)
tblHyenas$number_littermates <- as.numeric(tblHyenas$number_littermates)
tblHyenas$litrank <- as.numeric(tblHyenas$litrank)
tblHyenas <- tblHyenas[,c(1,4,7:13,15)]
tblHyenas <- filter(tblHyenas, !is.na(id))   #remove 0

#Adding litter variable
tblHyenas$litter <- NA
tblHyenas[tblHyenas$number_littermates != '0' & !is.na(tblHyenas$number_littermates) & 
            !is.na(tblHyenas$litrank) & tblHyenas$litrank == '1',]$litter <- 'dominant'
tblHyenas[tblHyenas$number_littermates != '0' & !is.na(tblHyenas$number_littermates) & 
            !is.na(tblHyenas$litrank) & tblHyenas$litrank == '2',]$litter <- 'subordinate'
tblHyenas[is.na(tblHyenas$number_littermates) | tblHyenas$number_littermates == '0',]$litter <- 'singleton'
tblHyenas$litter <- as.factor(tblHyenas$litter)

#Adding clan variable - because all natal animals, use natal clan
tblHyenas$clan <- NA
for(i in 1:nrow(tblHyenas)){
  id.i <- tblHyenas$id[i]
  if(id.i %in% tblLifeHistory$id){
    tblHyenas$clan[i] <- filter(tblLifeHistory, id == id.i & event_code == "DOB")$event_data
  }
}
rm(i)
tblHyenas$clan <- as.factor(tblHyenas$clan)

#tblRanks
data("tblFemaleRanks")
tblFemaleRanks <- filter(tblFemaleRanks, clan == "serena.n" | clan == "serena.s" | 
                           clan == "happy.zebra")

#tblWeather
data("tblWeather")
tblWeather <- filter(tblWeather, park == "Conservancy")
tblWeather[,3:5] <- sapply(tblWeather[,3:5], as.numeric)

#tblPreyCensus
data("tblPreyCensus")
tblPreyCensus <- filter(tblPreyCensus, region == "Conservancy" & !is.na(clan))
tblPreyCensus$year <- as.numeric(format(tblPreyCensus$date, "%Y"))
tblPreyCensus$month <- as.numeric(format(tblPreyCensus$date, "%m"))
tblPreyCensus[,c(4,9:39)] <- sapply(tblPreyCensus[,c(4,9:39)], as.numeric)


########## 1.3 Calculate prey density ##########

#Data: each clan has 2-4 transects
# Transects are 1.45-5.40 km in length
# Prey is surveyed for 100m on either side of transect
# Transects are surveyed 2x per month (once in first half of month, once in second half)

#Set up column categories
prey_all <- c("thomsons", "impala", "zebra", "wildebeest", "topi", "warthog", "hartebeest", "grants", 
              "buffalo", "hippo", "giraffe", "ostrich", "eland", "elephant", "oribi", "reedbuck", 
              "waterbuck", "baboon", "bushbuck")
to.sum <- c("distance", prey_all)

#Calculate prey summary data
prey.summary <- tblPreyCensus %>% group_by(region, clan, year, month) %>%
  summarise_at(vars(all_of(to.sum)), sum)     #prey density calculated per month 
prey.summary$total_prey_count <- rowSums(prey.summary[,prey_all])
prey.summary$area <- prey.summary$distance*0.2    #prey censused for 100m on either side of road (e.g., census of distance = 1km has an area of 0.2 km2)
prey.summary$prey_density <- as.numeric(prey.summary$total_prey_count/prey.summary$area)     #prey density calculated as # animals per km2
tblPreyDensity <- prey.summary[,c(1:5,25:27)]

#Sanity check
boxplot(tblPreyDensity$prey_density ~ tblPreyDensity$month, ylim = c(0,1000))
abline(h = median(tblPreyDensity$prey_density))


######################################################################
##### 2.0 Combining data #####
######################################################################

########## 2.1 Combine tables ##########

#Add weather
saliva.horm <- left_join(saliva.final, tblWeather[,c(2:5)], by = "date")

#Add prey density
saliva.horm$year <- as.numeric(format(saliva.horm$date, "%Y"))
saliva.horm$month <- as.numeric(format(saliva.horm$date, "%m"))
saliva.horm <- left_join(saliva.horm, tblPreyDensity[,c(2:4,8)], by = c("clan", "year", "month"))

#Add age
saliva.horm <- left_join(saliva.horm, tblHyenas[,c(1,4)], by = c("hyena_id" = "id"))
saliva.horm$age <- (as.numeric(saliva.horm$date - saliva.horm$birthdate)/365)*12    #age in months

#Add sex, litter
saliva.horm <- left_join(saliva.horm, tblHyenas[,c(1,3,6,8,9,11)], by = c("hyena_id" = "id"))

#Add maternal rank
saliva.horm$rank <- NA
for(i in 1:nrow(saliva.horm)){
  mom.i <- saliva.horm$mom[i]
  year.i <- as.numeric(format(saliva.horm$date[i], "%Y"))
  if(!is.na(mom.i) & !is.na(year.i)){
    rank.i <- filter(tblFemaleRanks, id == mom.i & year == year.i)$stan_rank
    #Add mom's rank from current year
    if(length(rank.i) == 1){
      saliva.horm$rank[i] <- rank.i
    }
    #Add mom's rank from previous year if current year rank is not available
    if(length(rank.i) == 0){
      rank.i <- filter(tblFemaleRanks, id == mom.i & year == (year.i-1))$stan_rank
      if(length(rank.i) == 1){
        saliva.horm$rank[i] <- rank.i
      }
      if(length(rank.i) == 0){
        saliva.horm$rank[i] <- NA
      }
    }
  }
}

#Reformat clan
saliva.horm$clan <- as.factor(as.character(saliva.horm$clan))

#Fix sex
saliva.horm$sex <- as.character(saliva.horm$sex)
saliva.horm[saliva.horm$sex == "u" & !is.na(saliva.horm$sex),]$sex <- NA    #if unknown sex - 3 samples
saliva.horm$sex <- as.factor(saliva.horm$sex)

#Fix litter
saliva.horm$litter <- relevel(saliva.horm$litter, ref = "singleton")
saliva.horm[saliva.horm$number_littermates == 0,]$litrank <- 0   #if singleton

#Add AM/PM column for time
saliva.horm$ampm <- format(saliva.horm$start_time, '%p')
saliva.horm$ampm <- as.factor(saliva.horm$ampm)

#Reformat times
saliva.horm$start_time <- as.POSIXct(format(saliva.horm$start_time, format = "%H:%M"), 
                                     format = "%H:%M")
saliva.horm$stop_time <- as.POSIXct(format(saliva.horm$stop_time, format = "%H:%M"), 
                                    format = "%H:%M")
saliva.horm$ln2_time <- as.POSIXct(format(saliva.horm$ln2_time, format = "%H:%M"), 
                                   format = "%H:%M")

#Calculate time lag from sunrise/sunset
# Data on sunrise/sunset from esrl.noaa.gov
#   June 21 2020: 6:40    18:42  @ Talek Gate
#   Dec  21 2020: 6:31    18:43
#   June 21 2020: 6:41    18:43  @ South Mara Bridge
#   Dec  21 2020: 6:32    18:45
#Take average time of sunrise/sunset
#   Final time:   6:36    18:43

#Add time lag data
saliva.horm$time_lag <- NA
for(i in 1:nrow(saliva.horm)){
  if(saliva.horm$ampm[i] == "AM"){
    saliva.horm$time_lag[i] <- as.numeric(difftime(saliva.horm$start_time[i], 
                                                   as.POSIXct("06:36", format = "%H:%M"), 
                                                   units = "mins"))
  }
  if(saliva.horm$ampm[i] == "PM"){
    saliva.horm$time_lag[i] <- as.numeric(difftime(saliva.horm$start_time[i], 
                                                   as.POSIXct("18:43", format = "%H:%M"), 
                                                   units = "mins"))
  }
}

#Calculate time differences
saliva.horm$chew_time <- as.numeric(difftime(saliva.horm$stop_time, saliva.horm$start_time, 
                                             units = "mins"))   #time chewed on rope in mins
saliva.horm$chew_time <- ifelse(saliva.horm$chew_time == 0, 1, saliva.horm$chew_time)
saliva.horm$cortisol_assay_diff <- (as.numeric(saliva.horm$cortisol_assay_date - 
                                                 saliva.horm$date)/365)*12   #time collection to assay in months
saliva.horm$ln2_diff <- as.numeric(difftime(saliva.horm$ln2_time, saliva.horm$start_time, 
                                            units = "mins"))    #time collection to freezing in mins
saliva.horm$ln2_diff <- ifelse(saliva.horm$ln2_diff < 0, NA, saliva.horm$ln2_diff)


########## 2.2 Select columns ##########

#Cortisol
saliva.cortisol <- saliva.horm[,c("saliva_sample_id", "repeated", "clan", "hyena_id", 
                                  "date", "start_time", "stop_time", "ampm", "time_lag", 
                                  "chew_time", "ln2_time", "ln2_diff", "cortisol_ug_dl", 
                                  "cortisol_assay_date", "cortisol_assay_diff", "temp_min", 
                                  "temp_max", "precip", "prey_density", "age", "sex", 
                                  "rank", "number_littermates", "litrank", "litter")]
saliva.cortisol <- filter(saliva.cortisol, !is.na(cortisol_ug_dl))

#Log-transform to achieve normality
qqnorm(saliva.cortisol$cortisol_ug_dl)
shapiro.test(saliva.cortisol$cortisol_ug_dl)
# Shapiro-Wilk normality test
# W = 0.41647, p-value < 2.2e-16

saliva.cortisol$log_cortisol_ug_dl <- log(saliva.cortisol$cortisol_ug_dl)
qqnorm(saliva.cortisol$log_cortisol_ug_dl)
shapiro.test(saliva.cortisol$log_cortisol_ug_dl)
# Shapiro-Wilk normality test
# W = 0.97428, p-value = 8.969e-05
# approximates normal 

#Add freeze-thaw cycles to data
saliva.cortisol <- left_join(saliva.cortisol, saliva.freezethaw, 
                             by = c("saliva_sample_id" = "sample_id", 
                                    "cortisol_assay_date" = "assay_date"))
summary(saliva.cortisol)


########## 2.3 Final clean of dataset ##########

#Filter to only juveniles
saliva.cortisol <- filter(saliva.cortisol, age < 24)     #remove 1


######################################################################
##### 3.0 Save data #####
######################################################################

save(file = "05.saliva_hormones.Rdata", list = c("saliva.cortisol"))




