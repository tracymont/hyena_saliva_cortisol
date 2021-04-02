################################################################################
######### 06: Clean behavior data ########
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
load("05.saliva_hormones.Rdata")


########## 1.2 hyenadata tables ##########

#tblSessions
data("tblSessions")
tblSessions$clan <- tolower(gsub(" ", "", tblSessions$clan))
tblSessions$location <- as.factor(tblSessions$location)
tblSessions$start <- as.POSIXct(tblSessions$start, format = "%H:%M:%S")
tblSessions$stop <- as.POSIXct(tblSessions$stop, format = "%H:%M:%S")
tblSessions$hyenas <- tolower(gsub(" ", "", tblSessions$hyenas))
tblSessions$unidhyenas <- tolower(tblSessions$unidhyenas)
tblSessions <- tblSessions[,1:9]
tblSessions <- filter(tblSessions, !is.na(session))   #remove 0
tblSessions <- filter(tblSessions, clan == "happy.zebra" | clan == "serena.s" | clan == "serena.n")

#tblHyenasPerSession
data("tblHyenasPerSession")
tblHyenasPerSession <- tblHyenasPerSession[,1:2]
tblHyenasPerSession <- filter(tblHyenasPerSession, !is.na(session) & !is.na(id))   #remove 0
tblHyenasPerSession <- filter(tblHyenasPerSession, session %in% tblSessions$session)

#tblLifeHistory
data("tblLifeHistory")
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

#Adding clan variable - because all natal animals, use natal clan
tblHyenas$clan <- NA
for(i in 1:nrow(tblHyenas)){
  id.i <- tblHyenas$id[i]
  if(id.i %in% tblLifeHistory$id){
    tblHyenas$clan[i] <- filter(tblLifeHistory, id == id.i & event_code == "DOB")$event_data
  }
}
rm(i)
rm(id.i)
tblHyenas$clan <- as.factor(tblHyenas$clan)


######################################################################
##### 2.0 tblAggression #####
######################################################################

########## 2.1 Load data ########## 

#Load data
data("tblAggression")
tblAggression$clan <- tolower(gsub(" ", "", tblAggression$clan))
tblAggression$aggid <- as.numeric(tblAggression$aggid)
tblAggression$session <- tolower(gsub(" ", "", tblAggression$session))
tblAggression$time <- as.POSIXct(tblAggression$time, format = "%H:%M:%S")
tblAggression$aggressor <- tolower(gsub(" ", "", tblAggression$aggressor))
tblAggression$recip <- tolower(gsub(" ", "", tblAggression$recip))
tblAggression <- unique(tblAggression)   #remove 1

#Combine tblAggression and tblSessions
tblAggression <- left_join(tblAggression[,c(1:17)], tblSessions[,c(1:6)], by = "session")
tblAggression <- left_join(tblAggression, tblHyenas[,c(1,11)], by = c("aggressor" = "id"))
tblAggression <- left_join(tblAggression, tblHyenas[,c(1,11)], by = c("recip" = "id"))
colnames(tblAggression)[1] <- "clan.tblAgg"
colnames(tblAggression)[4] <- "date.tblAgg"
colnames(tblAggression)[18] <- "clan.tblSes"
colnames(tblAggression)[20] <- "date.tblSes"
colnames(tblAggression)[23] <- "clan.aggressor"
colnames(tblAggression)[24] <- "clan.recip"

#Filter to Serena sessions May 2015 to June 2018
tblAggression <- filter(tblAggression, (clan.tblAgg == "happy.zebra" | clan.tblAgg == "serena.n" | 
                                          clan.tblAgg == "serena.s") | 
                          (clan.tblSes == "happy.zebra" | clan.tblSes == "serena.n" | clan.tblSes == "serena.s"))
nrow(filter(tblAggression, is.na(date.tblAgg) & is.na(date.tblSes)))    #0
tblAggression <- filter(tblAggression, (date.tblAgg > "2015-04-30") | (date.tblSes > "2015-04-30"))
tblAggression <- filter(tblAggression, (date.tblAgg < "2018-07-01") | (date.tblSes < "2018-07-01"))


########## 2.2 Data checking and cleaning ########## 

#Fix clan 
# View(filter(tblAggression, clan.tblAgg != clan.tblSes))     #tblAgg is most reliable
tblAggression$clan.tblAgg <- ifelse(is.na(tblAggression$clan.tblAgg), as.character(tblAggression$clan.tblSes), 
                                    as.character(tblAggression$clan.tblAgg))
tblAggression$clan.tblSes <- NULL
tblAggression$clan.aggressor <- NULL
tblAggression$clan.recip <- NULL
colnames(tblAggression)[1] <- "clan"

#Fix date
# View(filter(tblAggression, date.tblAgg != date.tblSes))     #tblSes is most reliable
tblAggression$date.tblSes <- ifelse(is.na(tblAggression$date.tblSes), as.character(tblAggression$date.tblAgg), 
                                    as.character(tblAggression$date.tblSes))
tblAggression$date.tblAgg <- tblAggression$date.tblSes
tblAggression$date.tblSes <- NULL
colnames(tblAggression)[4] <- "date"
tblAggression$date <- as.Date(tblAggression$date, format = "%Y-%m-%d")

#Filter to Serena sessions May 2015 to June 2018
tblAggression <- filter(tblAggression, (clan == "happy.zebra" | clan == "serena.n" | clan == "serena.s"))
tblAggression <- filter(tblAggression, (date > "2015-04-30") & (date < "2018-07-01"))

#Clean IDs - aggressors
tblAggression[tblAggression$aggressor == "alvn" & !is.na(tblAggression$aggressor),]$aggressor <- "avln"
tblAggression[tblAggression$aggressor == "bazr" & !is.na(tblAggression$aggressor),]$aggressor <- "gazr"
tblAggression[tblAggression$aggressor == "berd" & !is.na(tblAggression$aggressor),]$aggressor <- "nerd"
tblAggression[tblAggression$aggressor == "bmn" & !is.na(tblAggression$aggressor),]$aggressor <- "tbmn"
tblAggression[tblAggression$aggressor == "c;ay" & !is.na(tblAggression$aggressor),]$aggressor <- "clay"
tblAggression[tblAggression$aggressor == "drast" & !is.na(tblAggression$aggressor),]$aggressor <- "rast"
tblAggression[tblAggression$aggressor == "fava" & !is.na(tblAggression$aggressor),]$aggressor <- "java"
tblAggression[tblAggression$aggressor == "fses" & !is.na(tblAggression$aggressor),]$aggressor <- "dawa"
tblAggression[tblAggression$aggressor == "hthnd" & !is.na(tblAggression$aggressor),]$aggressor <- "thnd"
tblAggression[tblAggression$aggressor == "inid" & !is.na(tblAggression$aggressor),]$aggressor <- "unid"
tblAggression[tblAggression$aggressor == "inidblackcub1" & !is.na(tblAggression$aggressor),]$aggressor <- "unidblackcub1"
tblAggression[tblAggression$aggressor == "inidcub" & !is.na(tblAggression$aggressor),]$aggressor <- "unidcub"
tblAggression[tblAggression$aggressor == "jat" & !is.na(tblAggression$aggressor),]$aggressor <- "taj"
tblAggression[tblAggression$aggressor == "jlry" & !is.na(tblAggression$aggressor),]$aggressor <- "jlyr"
tblAggression[tblAggression$aggressor == "jylr" & !is.na(tblAggression$aggressor),]$aggressor <- "jlyr"
tblAggression[tblAggression$aggressor == "kava" & !is.na(tblAggression$aggressor),]$aggressor <- "java"
tblAggression[tblAggression$aggressor == "kiti" & !is.na(tblAggression$aggressor),]$aggressor <- "kita"
tblAggression[tblAggression$aggressor == "lamc" & !is.na(tblAggression$aggressor),]$aggressor <- "lanc"
tblAggression[tblAggression$aggressor == "lunc" & !is.na(tblAggression$aggressor),]$aggressor <- "lanc"
tblAggression[tblAggression$aggressor == "msnt" & !is.na(tblAggression$aggressor),]$aggressor <- "mnst"
tblAggression[tblAggression$aggressor == "ngon" & !is.na(tblAggression$aggressor),]$aggressor <- "nogn"
tblAggression[tblAggression$aggressor == "pesky" & !is.na(tblAggression$aggressor),]$aggressor <- "whiz"
tblAggression[tblAggression$aggressor == "pici" & !is.na(tblAggression$aggressor),]$aggressor <- "pixi"
tblAggression[tblAggression$aggressor == "sapr" & !is.na(tblAggression$aggressor),]$aggressor <- "spar"
tblAggression[tblAggression$aggressor == "sbi" & !is.na(tblAggression$aggressor),]$aggressor <- "wsbi"
tblAggression[tblAggression$aggressor == "tili" & !is.na(tblAggression$aggressor),]$aggressor <- "gili"
tblAggression[tblAggression$aggressor == "wcr" & !is.na(tblAggression$aggressor),]$aggressor <- "wncr"
check <- anti_join(filter(tblAggression, !is.na(aggressor)), tblHyenas, by = c("aggressor" = "id"))
check <- filter(check, !grepl("/", check$aggressor) & !grepl("\\?", check$aggressor) & !grepl("unid", check$aggressor))
check <- left_join(check, tblSessions[,c(2,7)], by = "session")
sort(unique(check$aggressor))     #none
rm(check)

#Clean IDs - recip
tblAggression[tblAggression$recip == "al1057" & !is.na(tblAggression$recip),]$recip <- "a1057"
tblAggression[tblAggression$recip == "al1060" & !is.na(tblAggression$recip),]$recip <- "a1060"
tblAggression[tblAggression$recip == "alvn" & !is.na(tblAggression$recip),]$recip <- "avln"
tblAggression[tblAggression$recip == "bhcr" & !is.na(tblAggression$recip),]$recip <- "bchr"
tblAggression[tblAggression$recip == "csb" & !is.na(tblAggression$recip),]$recip <- "csby"
tblAggression[tblAggression$recip == "foot" & !is.na(tblAggression$recip),]$recip <- "boot"
tblAggression[tblAggression$recip == "gazrn" & !is.na(tblAggression$recip),]$recip <- "gazr"
tblAggression[tblAggression$recip == "jylr" & !is.na(tblAggression$recip),]$recip <- "jlyr"
tblAggression[tblAggression$recip == "mnsy" & !is.na(tblAggression$recip),]$recip <- "mnst"
tblAggression[tblAggression$recip == "mrdk" & !is.na(tblAggression$recip),]$recip <- "mdrk"
tblAggression[tblAggression$recip == "palal" & !is.na(tblAggression$recip),]$recip <- "pala"
tblAggression[tblAggression$recip == "pix" & !is.na(tblAggression$recip),]$recip <- "pixi"
tblAggression[tblAggression$recip == "qhiz" & !is.na(tblAggression$recip),]$recip <- "whiz"
tblAggression[tblAggression$recip == "rru" & !is.na(tblAggression$recip),]$recip <- "tru"
tblAggression[tblAggression$recip == "sapr" & !is.na(tblAggression$recip),]$recip <- "spar"
tblAggression[tblAggression$recip == "sino" & !is.na(tblAggression$recip),]$recip <- "wino"
tblAggression[tblAggression$recip == "stkl" & !is.na(tblAggression$recip),]$recip <- "sktl"
tblAggression[tblAggression$recip == "urckt" & !is.na(tblAggression$recip),]$recip <- "rckt"
tblAggression[tblAggression$recip == "vemo" & !is.na(tblAggression$recip),]$recip <- "veni"
tblAggression[tblAggression$recip == "wchr" & !is.na(tblAggression$recip),]$recip <- "wncr"
tblAggression[tblAggression$recip == "wtfg" & !is.na(tblAggression$recip),]$recip <- "wtg"
check <- anti_join(filter(tblAggression, !is.na(recip)), tblHyenas, by = c("recip" = "id"))
check <- filter(check, !grepl("/", check$recip) & !grepl("\\?", check$recip) & !grepl("unid", check$recip))
check <- left_join(check, tblSessions[,c(2,7)], by = "session")
sort(unique(check$recip))    #all alien hyenas
rm(check)

#Final dataset
tblAggression <- unique(tblAggression)
summary(tblAggression)


######################################################################
##### 3.0 tblPlayScans #####
######################################################################

########## 3.1 Load data ########## 

#Load data
HZ_play <- read.csv("00.raw_data/HZPlayBehav_15April19_playscans.csv")
N_play <- read.csv("00.raw_data/NPlayBehav_15Apr19_playscans.csv")
S_play <- read.csv("00.raw_data/SPlayBehav_23Mar19_playscans.csv")
play.scans <- rbind(HZ_play, N_play, S_play)
rm(HZ_play)
rm(N_play)
rm(S_play)

#Clean data
play.scans$SessionID <- NULL
play.scans$Clan <- gsub(" ", "", play.scans$Clan)
play.scans$Clan <- gsub("HZ", "happy.zebra", play.scans$Clan)
play.scans$Clan <- gsub("S", "serena.s", play.scans$Clan)
play.scans$Clan <- gsub("N", "serena.n", play.scans$Clan)
play.scans$Date <- gsub(" ", "", play.scans$Date)
play.scans$Date <- as.Date(play.scans$Date, format= '%d-%b-%y')
play.scans$Time <- gsub(" ", "", play.scans$Time)
play.scans$Time <- as.POSIXct(play.scans$Time, format = "%H:%M")
play.scans$PlayGroup <- gsub(" ", "", play.scans$PlayGroup)
play.scans <- unique(play.scans[,1:4])    #remove 2

#Rename columns
colnames(play.scans)[1] <- "clan"
colnames(play.scans)[2] <- "date"
colnames(play.scans)[3] <- "time"
colnames(play.scans)[4] <- "play.group"


########## 3.2 Data checking and cleaning ########## 

#Fix mistyped hyena_ids
play.scans$play.group <- gsub("melh", "meln", play.scans$play.group)   #typo
play.scans$play.group <- gsub("fegr", "ferg", play.scans$play.group)   #typo
play.scans$play.group <- gsub("btsi", "bsti", play.scans$play.group)   #typo
play.scans$play.group <- gsub("dela", "bela", play.scans$play.group)   #typo
play.scans$play.group <- gsub("shrt", "shtr", play.scans$play.group)   #typo
play.scans$play.group <- gsub("pdnt", "pdtn", play.scans$play.group)   #typo
play.scans$play.group <- gsub("dnag", "dang", play.scans$play.group)   #typo
play.scans$play.group <- gsub("gig", "gili", play.scans$play.group)   #typo
play.scans$play.group <- gsub(",pix,", ",pixi,", play.scans$play.group)   #typo
play.scans$play.group <- gsub("ferri", "frri", play.scans$play.group)   #typo
play.scans$play.group <- gsub("king", "kng", play.scans$play.group)   #name changed later
play.scans$play.group <- gsub("ragn", "rgnk", play.scans$play.group)   #typo

#Fix dates/times
play.scans$time <- format(play.scans$time, "%H:%M")
play.scans[play.scans$date == "2017-08-07",]$date <- "2016-08-07"
play.scans[play.scans$date == "2016-10-14" & play.scans$play.group == ",bela,diva," & 
             !is.na(play.scans$time),]$time <- "18:51"
play.scans[play.scans$date == "2016-12-11" & play.scans$play.group == ",diva,deja,bela," & 
             is.na(play.scans$time),]$time <- "19:21"
play.scans[play.scans$date == "2015-06-12" & play.scans$play.group == ",nerd,dork," & 
             !is.na(play.scans$time),]$time <- "18:53"
play.scans[play.scans$date == "2016-09-15" & play.scans$play.group == ",ring,quak," & 
             !is.na(play.scans$time),]$time <- "07:27"
play.scans$time <- as.POSIXct(play.scans$time, format = "%H:%M")

#Assign session number
play.scans$session <- NA
for(i in 1:nrow(play.scans)){
  #Select clan, date, time
  clan.i <- play.scans$clan[i]
  date.i <- play.scans$date[i]
  time.i <- play.scans$time[i]
  #Select hyenas
  hyenas.i <- play.scans$play.group[i]
  #Filter out unIDs
  hyenas.i <- gsub("/", ",", hyenas.i)
  hyenas.i <- gsub("\\?", "", hyenas.i)
  #Split character string and get rid of blanks
  hyenas.i <- str_split(hyenas.i, pattern = ",")[[1]]
  hyenas.i <- hyenas.i[hyenas.i != ""]
  #If it hasn't already been assigned to a session
  if(is.na(play.scans$session[i])){
    #Filter to sessions on same date and where time is within the session times
    sessions.i <- unique(filter(tblSessions, date == date.i, (time.i >= start), (time.i <= stop))$session)
    #Filter to hyenas in those sessions who are also in the play group
    hps.i <- unique(filter(tblHyenasPerSession, session %in% sessions.i & id %in% hyenas.i)$session)
    #If only one session identified, assign it
    if(length(hps.i) == 1){
      play.scans$session[i] <- hps.i
    }
  }
}
rm(clan.i)
rm(date.i)
rm(time.i)
rm(hyenas.i)
rm(sessions.i)
rm(hps.i)
rm(i)

#Assign session number to ones that didn't work with for-loop
play.scans[play.scans$date == "2016-08-11" & play.scans$clan == "serena.s" & 
             play.scans$play.group == ",fyre,dang," ,]$session <- "s20208"	
play.scans[play.scans$date == "2017-06-27" & play.scans$clan == "serena.s" & 
             (grepl("tngy", play.scans$play.group) | grepl("swt", play.scans$play.group)),]$session <- "s21703"	
play.scans[play.scans$date == "2017-06-27" & play.scans$clan == "serena.s" & 
             (grepl("prdg", play.scans$play.group) | grepl("sami", play.scans$play.group)),]$session <- "s21706"	
play.scans[play.scans$date == "2017-09-17" & play.scans$clan == "serena.n" & 
             play.scans$play.group == ",naga,huk," ,]$session <- "s22506"	

#Check remaining ones
summary(as.factor(play.scans$session))   #3 NAs
check <- filter(play.scans, is.na(session))
check$play.group   #all unIDs - okay
rm(check)


########## 3.3 Create play dataset by individual id instead of by play group ########## 

play.scans.all <- play.scans
rm(play.scans)

play.scans.all <- unique(play.scans.all)
play.list = list()
counter = 1
for(i in 1:nrow(play.scans.all)){
  #Select clan, date, time
  session.current <- play.scans.all$session[i]
  clan.current <- play.scans.all$clan[i]
  date.current <- play.scans.all$date[i]
  time.current <- play.scans.all$time[i]
  if(!is.na(clan.current) & !is.na(date.current) & !is.na(time.current)){
    #Split play.group character string
    hyenas.current <- filter(play.scans.all, clan == clan.current & date == date.current & 
                               time == time.current)$play.group %>% strsplit(split = ',') %>% .[[1]]
    if(length(hyenas.current) > 0){
      #Create dataframe of hyenas
      hyenas.new <- data.frame(session = session.current, clan = clan.current, 
                               date = date.current, time = time.current, id = hyenas.current)
    }
    if(length(hyenas.current) == 0){
      #Create dataframe with NA for hyenas
      hyenas.new <- data.frame(session = session.current, clan = clan.current, 
                               date = date.current, time = time.current, id = NA)
    }
    play.list[[counter]] <- hyenas.new
  }
  counter <- counter+1
}

play.scans.sep <- do.call(rbind, play.list)
play.scans <- play.scans.sep
play.scans <- filter(play.scans, id != "")
play.scans$behav <- "play"

rm(hyenas.new)
rm(play.list)
rm(play.scans.sep)
rm(clan.current)
rm(date.current)
rm(hyenas.current)
rm(time.current)
rm(counter)
rm(i)

#Final dataset
play.scans <- unique(play.scans)
summary(play.scans)


######################################################################
##### 4.0 Categorize saliva behavior data #####
######################################################################

########## 4.1 Assign behaviors to categories ##########

#Define categories for variables
#These variables are matched partially, so 'car' would match 'carry, 'cares', 'carrot', caring is sharing, etc'
play <- c('play')
feed <- c('fd', "chew", "stick")
nurse <- c('nu')
active <- c('arr', 'app', 'wan', 'paste', 'inv', 'st', 'wlk', "roll", 'grass', 'grt',
            'whoop', "alarm", "pres", "poop", 'romp', "dig", "lope", "run", "ll")
rest <- c('so', 'den', 'oos', 'poke', 'sit')

saliva.behav$behav_cat <- ''

#Order is important
##active - standing or moving
saliva.behav[saliva.behav$behavior %>% 
               grep(pattern = paste(active, collapse = '|'), x = .),]$behav_cat <- 'active'

##play - involvement in any play
saliva.behav[saliva.behav$behavior %>% 
               grep(pattern = paste(play, collapse = '|'), x = .),]$behav_cat <- 'play'

##rest - sacked out or sitting (including in den)
saliva.behav[saliva.behav$behavior %>% 
               grep(pattern = paste(rest, collapse = '|'), x = .),]$behav_cat <- 'rest'

##nurse
saliva.behav[saliva.behav$behavior %>% 
               grep(pattern = paste(nurse, collapse = '|'), x = .),]$behav_cat <- 'nurse'

##feed - feeding or chewing (incl mouth contact with saliva stick)
saliva.behav[saliva.behav$behavior %>% 
               grep(pattern = paste(feed, collapse = '|'), x = .),]$behav_cat <- 'feed'


########## 4.2 Fix errors and clean up data ##########

#Fix ones that didn't match appropriately
saliva.behav$behav_cat[saliva.behav$behavior == 'so npp from'] <- 'nurse'
saliva.behav$behav_cat[saliva.behav$behavior == 'carry scrap'] <- 'feed'         #carry = likely fed
saliva.behav$behav_cat[saliva.behav$behavior == 'emerge D carry scrap'] <- 'feed'
saliva.behav$behav_cat[saliva.behav$behavior == 'soc snf spot'] <- 'active'
saliva.behav$behav_cat[saliva.behav$behavior == 'brt soc snf spot'] <- 'active'
saliva.behav$behav_cat[saliva.behav$behavior == 'brt soc snf'] <- 'active'
saliva.behav$behav_cat[saliva.behav$behavior == 'lv'] <- 'active'
saliva.behav$behav_cat[saliva.behav$behavior == 'groan ov'] <- 'active'     
saliva.behav$behav_cat[saliva.behav$behavior == 'inv saliva stick'] <- 'active'     
saliva.behav$behav_cat[saliva.behav$behavior == 'oos'] <- '' 
saliva.behav$behav_cat[saliva.behav$behavior == 'grm phallus'] <- ''     
saliva.behav$behav_cat[saliva.behav$behavior == 'snf phallus by'] <- ''     
saliva.behav$behav_cat[saliva.behav$behavior == 'snf phallus '] <- ''     
saliva.behav$behav_cat[saliva.behav$behavior == 'snf phallus'] <- ''     
saliva.behav$behav_cat[saliva.behav$behavior == 'go oos tall grass'] <- ''
saliva.behav$behav_cat[saliva.behav$behavior == 'emerge tall grass'] <- ''

#Check categories for errors
unique(filter(saliva.behav, behav_cat == "nurse")$behavior)
unique(filter(saliva.behav, behav_cat == "feed")$behavior)
unique(filter(saliva.behav, behav_cat == "rest")$behavior)
unique(filter(saliva.behav, behav_cat == "play")$behavior)
unique(filter(saliva.behav, behav_cat == "active")$behavior)
unique(filter(saliva.behav, behav_cat == "")$behavior)

#Change all others to here
saliva.behav$behav_cat[saliva.behav$behav_cat == ''] <- 'here'

#Format dataset for combining later
saliva.behav <- saliva.behav[,c(2,3,5,4,14)]
colnames(saliva.behav)[4] <- "id"
colnames(saliva.behav)[5] <- "behav"

#Final dataset
saliva.behav <- unique(saliva.behav)
summary(saliva.behav)


######################################################################
##### 5.0 Create final datasets #####
######################################################################

########## 5.1 Format behavior data ##########

#Create aggressor and recipient datasets
aggressors <- tblAggression[,c(1,4,5,6)]
aggressors <- unique(aggressors)
colnames(aggressors)[4] <- "id"
aggressors$behav <- "aggressor"

recipients <- tblAggression[,c(1,4,5,7)]
recipients <- unique(recipients)
colnames(recipients)[4] <- "id"
recipients$behav <- "recipient"


########## 5.2 Combine saliva behavior data ##########

saliva.behavior <- rbind(aggressors, recipients, play.scans[,2:6], saliva.behav)
saliva.behavior$clan <- gsub(" ", "", saliva.behavior$clan)
saliva.behavior$id <- gsub(" ", "", saliva.behavior$id)
saliva.behavior$behav <- gsub(" ", "", saliva.behavior$behav)
saliva.behavior <- filter(saliva.behavior, !is.na(date) & !is.na(time) & id != "")     #remove 1

saliva.behavior <- unique(saliva.behavior)
summary(saliva.behavior)
summary(as.factor(saliva.behavior$behav))


########## 5.3 Check for feeding ##########

#Fix times
saliva.behavior$time <- as.POSIXct(format(saliva.behavior$time, format = "%H:%M"), format = "%H:%M")
saliva.cortisol$start_time <- as.POSIXct(format(saliva.cortisol$start_time, format = "%H:%M"), format = "%H:%M")

#Filter to hyenas in saliva dataset
saliva.behavior <- filter(saliva.behavior, date %in% saliva.cortisol$date & id %in% saliva.cortisol$hyena_id)
summary(as.factor(saliva.behavior$behav))

#Assign feeding behavior
saliva.cortisol$feeding <- ""
for(i in 1:nrow(saliva.cortisol)){
  #Select clan, date, time, id
  clan.i <- saliva.cortisol$clan[i]
  date.i <- saliva.cortisol$date[i]
  time.i <- saliva.cortisol$start_time[i]
  id.i <- saliva.cortisol$hyena_id[i]
  #Check for feeding that occurred 0-10 min before saliva sample
  feed.i <- filter(saliva.behavior, clan == clan.i & date == date.i & id == id.i & 
                     time <= time.i & time >= (time.i - 10*60))
  feed.i <- filter(feed.i, behav == "feed" | behav == "nurse")$behav
  if(length(feed.i) >= 1){
    saliva.cortisol$feeding[i] <- "feed"
  }
}
summary(as.factor(saliva.cortisol$feeding))

#Filter out anyone who fed or chewed on something in the last 10 minutes - only 6 samples
saliva.cortisol <- filter(saliva.cortisol, feeding != "feed")


########## 5.4 Save data ##########

save(file = "07.saliva_behavior.Rdata", list = c("saliva.behavior", "saliva.cortisol"))




