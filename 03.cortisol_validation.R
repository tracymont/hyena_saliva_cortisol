################################################################################
#####                    03: Cortisol ACTH validation                      #####
################################################################################

########## 1.1 Set working directory & download packages ##########

rm(list = ls())
setwd("~/Documents/R/Saliva")
options(stringsAsFactors = F)
library(viridis)
library(tidyverse)

#Set plot colors using viridis
viridis_2 <- viridis(7)[-c(1,3,5,6,7)]
viridis_3 <- viridis(7)[-c(2,4,6,7)]


########## 1.2 ACTH (biological) validation ##########

acth <- read.csv("00.raw_data/ACTH.blood.saliva.csv")
acth <- acth[,c(1:6)]
acth$min_inj <- as.numeric(acth$min_inj)
acth$conc <- as.numeric(acth$conc)
acth$cv <- as.numeric(acth$cv)
acth <- filter(acth, !is.na(conc))
acth$log_conc <- log(acth$conc)

#Align samples for comparison
plasma <- filter(acth, type == "plasma")
plasma <- arrange(plasma, hyena_id, min_inj)

saliva <- filter(acth, type == "saliva")
saliva <- arrange(saliva, hyena_id, min_inj)

#Find nearest saliva sample (require it to be within 5 minutes of plasma sample)
plasma$nearest.saliva <- NA
for(i in 1:nrow(plasma)){
  id.i <- plasma$hyena_id[i]
  time.i <- plasma$min_inj[i]
  saliva.i <- filter(saliva, hyena_id == id.i)
  poss.matches <- filter(saliva.i, (min_inj <= (time.i + 5)) & (min_inj >= (time.i - 5)))
  if(nrow(poss.matches) == 0){
    plasma$nearest.saliva[i] <- "none"
  }
  if(nrow(poss.matches) == 1){
    plasma$nearest.saliva[i] <- poss.matches$sample_id
  }
  if(nrow(poss.matches) > 1){
    plasma$nearest.saliva[i] <- "multiple"
  }
}
rm(id.i)
rm(time.i)
rm(saliva.i)
rm(poss.matches)

#Rename columns for clarity
colnames(plasma) <- paste(colnames(plasma), "plasma", sep = "_")
colnames(saliva) <- paste(colnames(saliva), "saliva", sep = "_")

#Correlation test
comparison <- left_join(plasma, saliva, by = c("nearest.saliva_plasma" = "sample_id_saliva"))
comparison <- filter(comparison, !is.na(conc_saliva))
cor.test(comparison$log_conc_saliva, comparison$log_conc_plasma)
# Pearson's product-moment correlation
# data:  comparison$log_conc_saliva and comparison$log_conc_plasma
# t = 5.7848, df = 23, p-value = 6.811e-06
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.5385274 0.8932561
# sample estimates:
#       cor 
# 0.7698441 


########## 1.3 Look at time differences between plasma and saliva ##########

#Find max plasma point
filter(acth, type == "plasma") %>% arrange(hyena_id, conc)
# hyena sample_id min_inj    conc    cv sample    log_conc      conc_z
#  brgr     11736      49 47.8700  1.92 plasma  3.86848900  1.73952714
#  gaza     11744      54 38.2700  0.81 plasma  3.64466630  1.12312785
#  mchl     11747      56 23.4600  4.34 plasma  3.15529684  0.17220354
plyr::round_any(mean(c(49, 54, 56)), 5)     #55 min is average plasma peak

#Find max saliva point
filter(acth, type == "saliva") %>% arrange(hyena_id, conc)
# hyena sample_id min_inj    conc    cv sample    log_conc      conc_z
# brgr        R9      69  5.4850 11.73 saliva  1.70201709   2.3980249
# gaza       G12      80  3.8670  8.10 saliva  1.35247901   1.4157676
# mchl        M9      96  2.6090  4.74 saliva  0.95896701   0.6520594
plyr::round_any(mean(c(69, 80, 96)), 5)     #80 min is average saliva peak 

#Find last saliva point - same as max point
filter(acth, type == "saliva") %>% arrange(hyena_id, min_inj)
# brgr        R9      69  5.4850 11.73 saliva  1.70201709   2.3980249
# gaza       G12      80  3.8670  8.10 saliva  1.35247901   1.4157676
# mchl        M9      96  2.6090  4.74 saliva  0.95896701   0.6520594
plyr::round_any(mean(c(69, 80, 96)), 5)     #80 min is also the last time point - could be peak but don't know


########## 1.4 ACTH (biological) validation plot ##########

#Average samples prior to ACTH injection to create a single "baseline" point
means <- filter(acth, min_inj <= 0)
means <- means %>% group_by(hyena_id, type) %>% summarise(mean = mean(conc), min_inj = 0)
means$sample_id <- NA
means$cv <- NA
means$log_conc <- log(means$mean)
means <- means[,c(1,5,4,3,6,2,7)]
colnames(means)[4] <- "conc"
means <- as.data.frame(means)

#Combine with rest of dataset
acth <- filter(acth, min_inj != 0)     #remove 1
acth <- rbind(acth, means)

#Standardize by sample type
acth$conc_z <- NA
acth[acth$type == "saliva",]$conc_z <- scale(acth[acth$type == "saliva",]$conc, center = T, scale = T)
acth[acth$type == "plasma",]$conc_z <- scale(acth[acth$type == "plasma",]$conc, center = T, scale = T)

#Make ACTH plot
colnames(acth)[1] <- "hyena"
colnames(acth)[6] <- "sample"
#Legend
plot.acth <- ggplot(data = acth, aes(y = conc_z, x = min_inj, col = hyena, linetype = sample)) + 
  geom_vline(xintercept = 0, color = "dark grey") + 
  geom_point(size = 2) + 
  geom_line(size = 1) + 
  scale_color_manual(values = viridis_3) +
  xlab('Time since ACTH injection (min)') + 
  ylab('Cortisol (standardized)')+
  theme_classic() +
  theme(legend.title = element_text(size = 14), legend.text = element_text(size = 12),
        axis.title = element_text(size = 20), axis.text = element_text(size = 16)) + 
  scale_x_continuous(breaks = seq(-40, 120, 20))
pdf('03.plot.acth.pdf', width = 7, height = 5)
plot.acth
dev.off()

#No legend
plot.acth <- ggplot(data = acth, aes(y = conc_z, x = min_inj, col = hyena, linetype = sample)) + 
  geom_vline(xintercept = 0, color = "dark grey") + 
  geom_point(size = 2) + 
  geom_line(size = 1) + 
  scale_color_manual(values = viridis_3) +
  xlab('Time since ACTH injection (min)') + 
  ylab('Cortisol (standardized)')+
  theme_classic() +
  theme(legend.title = element_text(size = 14), legend.text = element_text(size = 12),
        axis.title = element_text(size = 20), axis.text = element_text(size = 16), 
        legend.position = "none") + 
  scale_x_continuous(breaks = seq(-40, 120, 20))
pdf('03.plot.acth.pdf', width = 7, height = 5)
plot.acth
dev.off()


########## 1.5 ACTH plot depicting lag-time in salivary cortisol ##########

#Plot for lag-time
plot.acth <- ggplot(data = acth, aes(y = conc_z, x = min_inj, col = hyena, linetype = sample)) + 
  geom_vline(xintercept = 0, color = "dark grey") + 
  geom_vline(xintercept = 20, color = "dark grey", lty = 2) + 
  geom_point(size = 2) + 
  geom_line(size = 1) + 
  scale_color_manual(values = viridis_3) +
  xlab('Time since ACTH injection (min)') + 
  ylab('Cortisol (standardized)')+
  theme_classic() +
  theme(legend.title = element_text(size = 14), legend.text = element_text(size = 12),
        axis.title = element_text(size = 20), axis.text = element_text(size = 16), 
        legend.position = "none") + 
  scale_x_continuous(breaks = seq(-40, 120, 20))
pdf('03.plot.acth.increase.pdf', width = 7, height = 5)
plot.acth
dev.off()


########## 1.6 ACTH plot depicting behavioral sampling window ##########

#Plot for behavioral sampling window
plot.acth <- ggplot(data = filter(acth, sample == "saliva"), aes(y = conc_z, x = min_inj, col = hyena)) + 
  geom_rect(xmin = 15, xmax = 45, ymin = -2, ymax = 3, fill = "gray80", color = NA) + 
  geom_point(size = 2) + 
  geom_line(size = 1, linetype = "dashed") + 
  scale_color_manual(values = viridis_3) +
  xlab('Time since ACTH injection (min)') + 
  ylab('Salivary cortisol (standardized)')+
  theme_classic() +
  theme(legend.title = element_text(size = 14), legend.text = element_text(size = 12),
        axis.title = element_text(size = 20), axis.text = element_text(size = 16), 
        legend.position = "none") + 
  scale_x_continuous(breaks = seq(0, 120, 10)) + 
  coord_cartesian(xlim = c(0,60))
pdf('03.plot.acth.zoom.pdf', width = 5, height = 5)
plot.acth
dev.off()




