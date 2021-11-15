##################################################################################################
#Swaziland EMOD plotting
##################################################################################################
library(reshape2)
library(plyr)
library(rlang)
library(dplyr)
library(ggplot2)
library(GGally)
library(dplyr)
library(mgcv)
library(data.table)
library(tidyr)
library(stringr)

options(scipen=999)

##################################################################################################
#Bring in reportHIVbyageandgender Baseline and 100% ART scenario for Eswatini model with ART stages
##################################################################################################

#Bring in newest downloaded baseline transmission data + reportHIVbyageandgender from same sim folders
input_dir <- "C:\\Users\\aakullian\\Dropbox (IDM)\\GitHub\\EMOD_eswatini\\Output\\Baseline_transmission_ExtraData"
primary_dirs = list.files(input_dir)
transmission.data.frame <- data.frame("Year2"=seq(1980,2050,1))

#reporthivbyageandgender data
i=0
for(dir in primary_dirs) {
  i=i+1
  temp_table <- as.data.table(read.csv(file = paste(input_dir,dir,"ReportHIVByAgeAndGender.csv",sep="/"))) 
  temp_table$Year2 <- floor((temp_table$Year-0.5))
  temp_table$Uninfected.Population = temp_table$Population-temp_table$Infected
  trajectories_IR.1a <- aggregate(cbind(Newly.Infected,Transmitters) ~ Year2+Gender+Age+IsCircumcised+IP_Key.ARTstate, subset(temp_table, Age>10 & Year < 2051), FUN=sum) #sums number of new infections/transmissions in each year
  trajectories_IR.2 <- aggregate(cbind(Uninfected.Population,Infected,Population) ~ Year+Gender+Age+IsCircumcised+IP_Key.ARTstate, subset(temp_table, Age>10 & Year < 2051), FUN=sum) #sums at risk population in each year
  trajectories_IR.2$Year2 <- floor(trajectories_IR.2$Year-0.5)
  trajectories_IR.2 <- trajectories_IR.2[!duplicated(trajectories_IR.2[c("Year2","Gender", "Age", "IsCircumcised", "IP_Key.ARTstate")]),] #remove second instance of duplicate rows
  trajectories_IR.2 <- trajectories_IR.2[-match("Year",names(trajectories_IR.2))]
  trajectories_IRoverall <- merge(trajectories_IR.1a, trajectories_IR.2, by=c("Year2","Gender","Age","IsCircumcised","IP_Key.ARTstate"))
  trajectories_IRoverall$sim.id <- paste(dir)
  if (i == 1){
    transmission.data.frame <- trajectories_IRoverall
  }
  else{
    transmission.data.frame <- rbind(transmission.data.frame, trajectories_IRoverall)
  }
  print(paste("working on folder", i, dir,sep=" "))
}

reporthivbyageandgender.baseline <- transmission.data.frame

#transmission data
i=0
vars <-  c("SRC_ID", "SRC_GENDER", "DEST_ID", "DEST_GENDER", "SRC_STAGE","src_age_int","dest_age_int","Year2")
for(dir in primary_dirs) {
  i=i+1
  temp_table <- as.data.table(read.csv(file = paste(input_dir,dir,"TransmissionReport.csv",sep="/"))) 
  temp_table <- subset(temp_table, YEAR>=1980 & YEAR <2050)
  temp_table$Year2 <- floor(temp_table$YEAR)
  temp_table$src_age_int <- floor(temp_table$SRC_AGE / 365.25)
  temp_table$dest_age_int <- floor(temp_table$DEST_AGE / 365.25)
  temp_table <- subset(temp_table, dest_age_int > 14 & src_age_int > 14)
  temp_table <- temp_table[, vars, with = FALSE]
  temp_table$sim.id <- paste(dir)
  if (i == 1){
    transmission.data.frame <- temp_table
  }
  else{
    transmission.data.frame <- rbind(transmission.data.frame, temp_table)
  }
  print(paste("working on folder", i, dir,sep=" "))
}

transmissionreport.master.baseline <- transmission.data.frame

#100% ART

#Bring in newest downloaded baseline transmission data + reportHIVbyageandgender from same sim folders
input_dir <- "C:\\Users\\aakullian\\Dropbox (IDM)\\GitHub\\EMOD_eswatini\\Output\\ART100pct_transmission_ExtraData"
primary_dirs = list.files(input_dir)
transmission.data.frame <- data.frame("Year2"=seq(1980,2050,1))

#reporthivbyageandgender data
i=0
for(dir in primary_dirs) {
  i=i+1
  temp_table <- as.data.table(read.csv(file = paste(input_dir,dir,"ReportHIVByAgeAndGender.csv",sep="/"))) 
  temp_table$Year2 <- floor((temp_table$Year-0.5))
  temp_table$Uninfected.Population = temp_table$Population-temp_table$Infected
  trajectories_IR.1a <- aggregate(cbind(Newly.Infected,Transmitters) ~ Year2+Gender+Age+IsCircumcised+IP_Key.ARTstate, subset(temp_table, Age>10 & Year < 2051), FUN=sum) #sums number of new infections/transmissions in each year
  trajectories_IR.2 <- aggregate(cbind(Uninfected.Population,Infected,Population) ~ Year+Gender+Age+IsCircumcised+IP_Key.ARTstate, subset(temp_table, Age>10 & Year < 2051), FUN=sum) #sums at risk population in each year
  trajectories_IR.2$Year2 <- floor(trajectories_IR.2$Year-0.5)
  trajectories_IR.2 <- trajectories_IR.2[!duplicated(trajectories_IR.2[c("Year2","Gender", "Age", "IsCircumcised", "IP_Key.ARTstate")]),] #remove second instance of duplicate rows
  trajectories_IR.2 <- trajectories_IR.2[-match("Year",names(trajectories_IR.2))]
  trajectories_IRoverall <- merge(trajectories_IR.1a, trajectories_IR.2, by=c("Year2","Gender","Age","IsCircumcised","IP_Key.ARTstate"))
  trajectories_IRoverall$sim.id <- paste(dir)
  if (i == 1){
    transmission.data.frame <- trajectories_IRoverall
  }
  else{
    transmission.data.frame <- rbind(transmission.data.frame, trajectories_IRoverall)
  }
  print(paste("working on folder", i, dir,sep=" "))
}

reporthivbyageandgender.ART100pct <- transmission.data.frame

#transmission data
i=0
vars <-  c("SRC_ID", "SRC_GENDER", "DEST_ID", "DEST_GENDER", "SRC_STAGE","src_age_int","dest_age_int","Year2")
for(dir in primary_dirs) {
  i=i+1
  temp_table <- as.data.table(read.csv(file = paste(input_dir,dir,"TransmissionReport.csv",sep="/"))) 
  temp_table <- subset(temp_table, YEAR>=1980 & YEAR <2050)
  temp_table$Year2 <- floor(temp_table$YEAR)
  temp_table$src_age_int <- floor(temp_table$SRC_AGE / 365.25)
  temp_table$dest_age_int <- floor(temp_table$DEST_AGE / 365.25)
  temp_table <- subset(temp_table, dest_age_int > 14 & src_age_int > 14)
  temp_table <- temp_table[, vars, with = FALSE]
  temp_table$sim.id <- paste(dir)
  if (i == 1){
    transmission.data.frame <- temp_table
  }
  else{
    transmission.data.frame <- rbind(transmission.data.frame, temp_table)
  }
  print(paste("working on folder", i, dir,sep=" "))
}

transmissionreport.master.ART100pct <- transmission.data.frame

#combine into one reportHIV and one transmission datafile
reporthivbyageandgender.baseline$scenario = "baseline"
reporthivbyageandgender.ART100pct$scenario = "ART100pct"
transmissionreport.master.baseline$scenario = "baseline"
transmissionreport.master.ART100pct$scenario = "ART100pct"

reporthivbyageandgender <- rbind(reporthivbyageandgender.baseline, reporthivbyageandgender.ART100pct)
transmissionreport <- rbind(transmissionreport.master.baseline, transmissionreport.master.ART100pct)

#save an Rda data file for inputs: 
setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\GitHub\\EMOD_eswatini\\Output")
save(reporthivbyageandgender, transmissionreport, file = "inputfiles.rda")

#################################################################################################
#Bring in Model data
#################################################################################################

#Documentation on all variables: https://idmod.org/docs/hiv/software-report-transmission.html?searchText=TransmissionReport

#bring in full saved dataset (bring in just the first sim with all vars for now)
setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\GitHub\\EMOD_eswatini\\Output")
load(file = "inputfiles.rda")

#create variables
names(transmissionreport)
transmissionreport$age.at.acquisition.cat <- cut(transmissionreport$dest_age_int, c(15,25,35,200), right=F)
transmissionreport$age.at.transmission.cat <- cut(transmissionreport$src_age_int, c(15,25,35,200), right=F)
transmissionreport$age.at.acquisition.cat5 <- cut(transmissionreport$dest_age_int, c(15,20,25,30,35,40,45,200), right=F)
transmissionreport$age.at.transmission.cat5 <- cut(transmissionreport$src_age_int, c(15,20,25,30,35,40,45,200), right=F)
transmissionreport$AgeCat <- cut(transmissionreport$dest_age_int, c(15,20,25,30,35,40,45,50,55,60,100), right=F)
transmissionreport$yearcat3 <- cut(transmissionreport$Year2, c(1980,1990,2004,2016,2030,2050), right=F, dig.lab=4)

#################################################################################################
# Set up transmission data
#################################################################################################
#library(igraph)

#Create transmission network for plotting (choose subset)
head(transmissionreport)
table(transmissionreport$scenario)
network.data <- subset(transmissionreport) #all sims
date.of.acquisition <- data.frame("SRC_ID"=network.data$DEST_ID, "acq.year"=network.data$Year2, "sim.id"=network.data$sim.id, "scenario"=network.data$scenario) #source id used to merge in with source of transmitters
network.data <- left_join(network.data, date.of.acquisition, by=c("SRC_ID","sim.id","scenario"))
names(network.data)

network.data.m <- data.frame("sim.id"=network.data$sim.id, 
                             "src"=network.data$SRC_ID, 
                             "dest"=network.data$DEST_ID, 
                             "src.gender"=network.data$SRC_GENDER, 
                             "dest.gender"=network.data$DEST_GENDER, 
                             "dest.age.cat"=network.data$age.at.acquisition.cat5, 
                             "src.age.cat"=network.data$age.at.transmission.cat5,
                             "dest.age.int"=network.data$dest_age_int, 
                             "src.age.int"=network.data$src_age_int,
                             "transm.year.int"=floor(network.data$Year2),
                             "acq.year"=network.data$acq.year,
                             "transm.year.cat"=network.data$yearcat3,
                             "src.stage"=network.data$SRC_STAGE,
                             "scenario"=network.data$scenario,
                             "AgeCat"=network.data$AgeCat)
names(network.data.m)

nrow(network.data.m[is.na(network.data.m$acq.year)==T,]) #number of individuals who's infection was seeded in (see campaign.json)
#network.data.m$acq.year[is.na(network.data.m$acq.year)==T] <- 1982.095633 #fills in missing infection year with seed year (left out for now)
network.data.m$acq.age=floor(network.data.m$src.age.int-(network.data.m$transm.year.int-network.data.m$acq.year))
network.data.m$acq.year.cat <- cut(network.data.m$acq.year, c(1990,2004,2016,2030,2050), right=F, dig.lab=4)
network.data.m$transm.year.cat <- cut(network.data.m$transm.year.int, c(1990,2004,2016,2030,2050), right=F, dig.lab=4)
network.data.m$acq.age.cat <- cut(network.data.m$acq.age, c(15,20,25,30,35,40,45,200), right=F)
network.data.m$src.age.cat.all <- cut(network.data.m$src.age.int, seq(15,100,5), right=F)
network.data.m$src.age.25 <- cut(network.data.m$src.age.int, c(15,25,100), right=F)
table(network.data.m$src.age.25)
network.data.m$src.age.cat3 <- cut(network.data.m$src.age.int, c(15,25,35,200), right=F)
network.data.m$dest.age.cat3 <- cut(network.data.m$dest.age.int, c(15,25,35,200), right=F)
network.data.m$era <- "UTT era" #era defined by year of transmission
network.data.m$era[network.data.m$transm.year < 2016] <- "Pre-UTT"
network.data.m$era.acq <- "UTT era"
network.data.m$era.acq[network.data.m$acq.year < 2016] <- "Pre-UTT"

#network.data.m <- network.data.m[order(network.data.m$simtime),]
head(network.data.m,20)
summary(network.data.m$transm.year.int)

network.data.m$src.stage.f <- factor(network.data.m$src.stage,
                                     levels = c(1,2,3,4),
                                     labels = c("Acute","Latent","AIDS","on ART"))

network.data.m$TimeToTransmission <- network.data.m$transm.year.int - network.data.m$acq.year

#################################################################################################
# clean reportHIVbyageandgender file
#################################################################################################

reporthivbyageandgender$src.gender <- reporthivbyageandgender$Gender
reporthivbyageandgender$year <- reporthivbyageandgender$Year2
reporthivbyageandgender$NeverOnART <- 0
reporthivbyageandgender$NeverOnART[reporthivbyageandgender$IP_Key.ARTstate=="NeverOnART"] <- 1
table(reporthivbyageandgender$NeverOnART, reporthivbyageandgender$IP_Key.ARTstate)
reporthivbyageandgender$src.age.cat <- cut(reporthivbyageandgender$Age, c(15,20,25,30,35,40,45,200), right=F)
reporthivbyageandgender$AgeCat <- cut(reporthivbyageandgender$Age, c(15,20,25,30,35,40,45,50,55,60,200), right=F)

#################################################################################################
#Risk of HIV transmission by subgroups - number of transmissions per 100 py infected
#################################################################################################
labs <- c("1" = "Women", "0" = "Men")
labs.age <- c("[15,25)"="15-24","[25,35)"="25-34","[35,200)"="35+")

# # # Overall # # #
head(network.data.m)
transmission.table.year <- subset(network.data.m)  %>%
  group_by(scenario, sim.id, transm.year.int, src.stage.f) %>%
  summarize(transmissions=n()) 
transmission.table.year$year <- transmission.table.year$transm.year.int
head(transmission.table.year)

head(reporthivbyageandgender)
reporthivbyageandgender.table.incidence <- subset(reporthivbyageandgender)  %>%
  group_by(scenario, sim.id, year) %>%
  summarize(Uninfected.Population=sum(Uninfected.Population)) 

transm.risk.table.all <- merge(transmission.table.year, reporthivbyageandgender.table.incidence, by=c("scenario","year","sim.id"), all.x=T)
transm.risk.table.all$inc <- transm.risk.table.all$transmissions/transm.risk.table.all$Uninfected.Population

transm.risk.table.all.freq <- subset(transm.risk.table.all) %>%
  group_by(scenario, year, sim.id) %>%
  mutate(freq = transmissions / sum(transmissions)) %>%
  group_by(scenario, year, src.stage.f) %>%
  summarize(Median=median(freq), lb = quantile(freq, probs = 0.025), ub = quantile(freq, probs = 0.975)) 

transm.risk.table.all.ci <- subset(transm.risk.table.all)  %>%
  group_by(scenario, year, src.stage.f) %>%
  summarize(MedianInc=median(inc), lb = quantile(inc, probs = 0.025), ub = quantile(inc, probs = 0.975)) 

#Baseline - incidence
ggplot(subset(transm.risk.table.all.ci, scenario=="baseline"), aes(fill=src.stage.f, y=MedianInc*100, x=year)) + 
  geom_bar(position="stack", stat="identity") +
  scale_x_continuous(expand = c(0,0), breaks=seq(1980,2045,5)) +
  scale_y_continuous(breaks=seq(0,0.035,0.005)*100,expand = c(0,0),"new infections / 100 py") +
  coord_cartesian(xlim=c(1980, 2050)) +
  theme(text=element_text(size=24))+
  theme(strip.background = element_rect(colour="white", fill="white")) +
  theme(axis.text.x = element_text(angle = 90)) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)+
  theme(legend.title=element_blank(),axis.line = element_line(colour = "black"),panel.background = element_blank(),legend.position = "bottom")

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\Transmission in Era of UTT\\Figures")
ggsave("HIVincidence_by_src.stage_baseline_bar.jpg", height=10, width=10)

ggplot(subset(transm.risk.table.all.ci, scenario=="baseline"),aes(x=year)) + 
  geom_line(aes(y=MedianInc*100, color=src.stage.f), size=2) +
  #geom_line(aes(y=MedianPopInc*100), size=2, color="black") +
  geom_ribbon(aes(ymin=lb*100, ymax=ub*100, fill=src.stage.f), size=2, alpha=0.2) +
  #geom_ribbon(aes(ymin=pop.lb*100, ymax=pop.ub*100), fill="black", size=2, alpha=0.2) +
  #facet_grid(src.stage.f~scenario) +
  ggtitle("Incidence of HIV by source's stage of infection") +
  #scale_color_manual(values = c("blue", "red"),labels = c("Male", "Female"))+
  scale_x_continuous(expand = c(0,0), breaks=seq(1980,2045,5)) +
  scale_y_continuous(breaks=seq(0,0.05,0.005)*100,expand = c(0,0),"new infections / 100 py") +
  coord_cartesian(xlim=c(1980, 2050)) +
  theme(text=element_text(size=24))+
  theme(strip.background = element_rect(colour="white", fill="white")) +
  theme(axis.text.x = element_text(angle = 90)) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)+
  theme(legend.title=element_blank(),axis.line = element_line(colour = "black"),panel.background = element_blank(),legend.position = "bottom")

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\Transmission in Era of UTT\\Figures")
ggsave("HIVincidence_by_src.stage_baseline.jpg", height=10, width=10)

#Baseline - proportion of incidence bar stacked
ggplot(subset(transm.risk.table.all.freq, scenario=="baseline"), aes(fill=src.stage.f, y=Median*100, x=year)) + 
  geom_bar(position="stack", stat="identity") +
  scale_x_continuous(expand = c(0,0), breaks=seq(1980,2045,5)) +
  scale_y_continuous(breaks=seq(0,100,20),expand = c(0,0),"% of incidence") +
  coord_cartesian(xlim=c(1980, 2050)) +
  theme(text=element_text(size=24))+
  theme(strip.background = element_rect(colour="white", fill="white")) +
  theme(axis.text.x = element_text(angle = 90)) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)+
  theme(legend.title=element_blank(),axis.line = element_line(colour = "black"),panel.background = element_blank(),legend.position = "bottom")

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\Transmission in Era of UTT\\Figures")
ggsave("HIVincidence_percent_baseline_bar.jpg", height=10, width=10)

#Baseline - proportion of incidence
ggplot(subset(transm.risk.table.all.freq, scenario=="baseline"),aes(x=year)) + 
  geom_line(aes(y=Median*100, color=src.stage.f), size=2) +
  geom_ribbon(aes(ymin=lb*100, ymax=ub*100, fill=src.stage.f), size=2, alpha=0.2) +
  facet_grid(labeller=labeller(Gender = labs)) +
  ggtitle("% of HIV incidence by source's stage of infection") +
  scale_x_continuous(expand = c(0,0), breaks=seq(1980,2045,5)) +
  scale_y_continuous(breaks=seq(0,1,0.1)*100,expand = c(0,0),"% of total incidence") +
  #coord_cartesian(ylim=c(0, 0.027)*100, xlim=c(1980, 2050)) +
  theme(text=element_text(size=24))+
  theme(strip.background = element_rect(colour="white", fill="white")) +
  theme(axis.text.x = element_text(angle = 90)) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)+
  theme(legend.title=element_blank(),axis.line = element_line(colour = "black"),panel.background = element_blank(),legend.position = "bottom")

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\Transmission in Era of UTT\\Figures")
ggsave("HIVincidence_percent_baseline.jpg", height=10, width=10)

#ART100 - incidence bar stacked
ggplot(subset(transm.risk.table.all.ci, scenario=="ART100pct"), aes(fill=src.stage.f, y=MedianInc*100, x=year)) + 
  geom_bar(position="stack", stat="identity") +
  scale_x_continuous(expand = c(0,0), breaks=seq(1980,2045,5)) +
  scale_y_continuous(breaks=seq(0,0.035,0.005)*100,expand = c(0,0),"new infections / 100 py") +
  coord_cartesian(xlim=c(1980, 2050)) +
  theme(text=element_text(size=24))+
  theme(strip.background = element_rect(colour="white", fill="white")) +
  theme(axis.text.x = element_text(angle = 90)) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)+
  theme(legend.title=element_blank(),axis.line = element_line(colour = "black"),panel.background = element_blank(),legend.position = "bottom")

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\Transmission in Era of UTT\\Figures")
ggsave("HIVincidence_percent_ART100pct_bar.jpg", height=10, width=10)

#ART100 - incidence
ggplot(subset(transm.risk.table.all.ci, scenario=="ART100pct"),aes(x=year)) + 
  geom_line(aes(y=MedianInc*100, color=src.stage.f), size=2) +
  geom_ribbon(aes(ymin=lb*100, ymax=ub*100, fill=src.stage.f), size=2, alpha=0.2) +
  #facet_grid(src.stage.f~scenario) +
  ggtitle("Incidence of HIV by source's stage of infection") +
  #scale_color_manual(values = c("blue", "red"),labels = c("Male", "Female"))+
  scale_x_continuous(expand = c(0,0), breaks=seq(1980,2045,5)) +
  scale_y_continuous(breaks=seq(0,0.05,0.005)*100,expand = c(0,0),"new infections / 100 py") +
  coord_cartesian(xlim=c(1980, 2050)) +
  theme(text=element_text(size=24))+
  theme(strip.background = element_rect(colour="white", fill="white")) +
  theme(axis.text.x = element_text(angle = 90)) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)+
  theme(legend.title=element_blank(),axis.line = element_line(colour = "black"),panel.background = element_blank(),legend.position = "bottom")

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\Transmission in Era of UTT\\Figures")
ggsave("HIVincidence_by_src.stage_ART100pct.jpg", height=10, width=10)

#ART100 - proportion of incidence bar stacked
ggplot(subset(transm.risk.table.all.freq, scenario=="ART100pct"), aes(fill=src.stage.f, y=Median*100, x=year)) + 
  geom_bar(position="stack", stat="identity") +
  scale_x_continuous(expand = c(0,0), breaks=seq(1980,2045,5)) +
  scale_y_continuous(breaks=seq(0,100,20),expand = c(0,0),"% of incidence") +
  coord_cartesian(xlim=c(1980, 2050)) +
  theme(text=element_text(size=24))+
  theme(strip.background = element_rect(colour="white", fill="white")) +
  theme(axis.text.x = element_text(angle = 90)) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)+
  theme(legend.title=element_blank(),axis.line = element_line(colour = "black"),panel.background = element_blank(),legend.position = "bottom")

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\Transmission in Era of UTT\\Figures")
ggsave("HIVincidence_percent_ART100_bar.jpg", height=10, width=10)

#ART100 - proportion of incidence
ggplot(subset(transm.risk.table.all.freq, scenario=="ART100pct"),aes(x=year)) + 
  geom_line(aes(y=Median*100, color=src.stage.f), size=2) +
  geom_ribbon(aes(ymin=lb*100, ymax=ub*100, fill=src.stage.f), size=2, alpha=0.2) +
  facet_grid(labeller=labeller(Gender = labs)) +
  ggtitle("% of HIV incidence by source's stage of infection") +
  scale_x_continuous(expand = c(0,0), breaks=seq(1980,2045,5)) +
  scale_y_continuous(breaks=seq(0,1,0.1)*100,expand = c(0,0),"% of total incidence") +
  #coord_cartesian(ylim=c(0, 0.027)*100, xlim=c(1980, 2050)) +
  theme(text=element_text(size=24))+
  theme(strip.background = element_rect(colour="white", fill="white")) +
  theme(axis.text.x = element_text(angle = 90)) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)+
  theme(legend.title=element_blank(),axis.line = element_line(colour = "black"),panel.background = element_blank(),legend.position = "bottom")

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\Transmission in Era of UTT\\Figures")
ggsave("HIVincidence_percent_ART100.jpg", height=10, width=10)

# # # By age and gender of at risk of infection # # #
transmission.table.age.gender.year <- subset(network.data.m)  %>%
  group_by(scenario, sim.id, dest.gender, dest.age.cat3, transm.year.int, src.stage.f) %>%
  summarize(transmissions=n()) 
transmission.table.age.gender.year$year <- transmission.table.age.gender.year$transm.year.int
transmission.table.age.gender.year$Gender <- transmission.table.age.gender.year$dest.gender

reporthivbyageandgender.table.incidence <- subset(reporthivbyageandgender)  %>%
  group_by(scenario, sim.id, Gender, dest.age.cat3, year) %>%
  summarize(Uninfected.Population=sum(Uninfected.Population)) 

transm.risk.table <- merge(transmission.table.age.gender.year, reporthivbyageandgender.table.incidence, 
                           by=c("scenario","sim.id","Gender","dest.age.cat3","year"), all.x=T)
transm.risk.table$inc <- transm.risk.table$transmissions/transm.risk.table$Uninfected.Population
transm.risk.table.dest.age <- subset(transm.risk.table)  %>%
  group_by(scenario, dest.age.cat3, Gender, year, src.stage.f) %>%
  summarize(MedianInc=median(inc), lb = quantile(inc, probs = 0.025), ub = quantile(inc, probs = 0.975)) 

transm.risk.table.dest.age.freq <- subset(transm.risk.table) %>%
  group_by(scenario, dest.age.cat3, Gender, year, sim.id) %>%
  mutate(freq = transmissions / sum(transmissions)) %>%
  group_by(scenario, dest.age.cat3, Gender, year, src.stage.f) %>%
  summarize(Median=median(freq), lb = quantile(freq, probs = 0.025), ub = quantile(freq, probs = 0.975)) 

table(transm.risk.table.dest.age$scenario)
#Baseline
ggplot(subset(transm.risk.table.dest.age, scenario=="baseline" & year > 2014),aes(x=year, fill=src.stage.f, y=MedianInc*100)) + 
  geom_line(aes(y=MedianInc*100, color=src.stage.f), size=2) +
  #geom_bar(position="stack", stat="identity") +
  geom_ribbon(aes(ymin=lb*100, ymax=ub*100, fill=src.stage.f), size=2, alpha=0.2) +
  facet_grid(Gender~dest.age.cat3, labeller=labeller(Gender = labs, dest.age.cat3=labs.age)) +
  scale_x_continuous(expand = c(0,0), breaks=seq(1970,2045,5)) +
  scale_y_continuous(breaks=seq(0,0.02,0.002)*100,expand = c(0,0),"new infections / 100 py") +
  #coord_cartesian(ylim=c(0, 0.013)*100, xlim=c(2015, 2050)) +
  theme(text=element_text(size=24))+
  theme(strip.background = element_rect(colour="white", fill="white")) +
  theme(axis.text.x = element_text(angle = 90)) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)+
  theme(legend.title=element_blank(),axis.line = element_line(colour = "black"),panel.background = element_blank(),legend.position = "bottom")

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\Transmission in Era of UTT\\Figures")
ggsave("HIVincidence_by_src.stage_age_gender_baseline.jpg", height=10, width=15)

ggplot(subset(transm.risk.table.dest.age.freq, scenario=="baseline" & year > 2014),aes(x=year, y=Median*100)) + 
  geom_line(aes(y=Median*100, color=src.stage.f), size=2) +
  #geom_bar(position="stack", stat="identity") +
  geom_ribbon(aes(ymin=lb*100, ymax=ub*100, fill=src.stage.f), size=2, alpha=0.2) +
  facet_grid(Gender~dest.age.cat3, labeller=labeller(Gender = labs, dest.age.cat3=labs.age)) +
  scale_x_continuous(expand = c(0,0), breaks=seq(1970,2045,5)) +
  scale_y_continuous(expand = c(0,0),"percent of incidence") +
  #coord_cartesian(ylim=c(0, 0.013)*100, xlim=c(2015, 2050)) +
  theme(text=element_text(size=24))+
  theme(strip.background = element_rect(colour="white", fill="white")) +
  theme(axis.text.x = element_text(angle = 90)) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)+
  theme(legend.title=element_blank(),axis.line = element_line(colour = "black"),panel.background = element_blank(),legend.position = "bottom")

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\Transmission in Era of UTT\\Figures")
ggsave("HIVincidence_by_src.stage_age_gender_baseline.jpg", height=10, width=15)

# # # By gender (all ages) # # #
transmission.table.year <- subset(network.data.m)  %>%
  group_by(scenario, sim.id, dest.gender, transm.year.int, src.stage.f) %>%
  summarize(transmissions=n()) 
transmission.table.year$year <- transmission.table.year$transm.year.int
transmission.table.year$Gender <- transmission.table.year$dest.gender

reporthivbyageandgender.table.incidence <- subset(reporthivbyageandgender)  %>%
  group_by(scenario, sim.id, Gender, year) %>%
  summarize(Uninfected.Population=sum(Uninfected.Population)) 

transm.risk.table.both.genders <- merge(transmission.table.year, reporthivbyageandgender.table.incidence, by=c("scenario","Gender","year","sim.id"), all.x=T)
transm.risk.table.both.genders$inc <- transm.risk.table.both.genders$transmissions/transm.risk.table.both.genders$Uninfected.Population
transm.risk.table.both.genders <- subset(transm.risk.table.both.genders)  %>%
  group_by(scenario, Gender, year, src.stage.f) %>%
  summarize(MedianInc=median(inc), lb = quantile(inc, probs = 0.025), ub = quantile(inc, probs = 0.975)) 

#Baseline
ggplot(subset(transm.risk.table.both.genders, scenario=="baseline"),aes(x=year)) + 
  geom_line(aes(y=MedianInc*100, color=src.stage.f), size=2) +
  geom_ribbon(aes(ymin=lb*100, ymax=ub*100, fill=src.stage.f), size=2, alpha=0.2) +
  facet_grid(~Gender, labeller=labeller(Gender = labs)) +
  ggtitle("Incidence of HIV by source's stage of infection") +
  #scale_color_manual(values = c("blue", "red"),labels = c("Male", "Female"))+
  scale_x_continuous(expand = c(0,0), breaks=seq(1970,2045,5)) +
  scale_y_continuous(breaks=seq(0,0.05,0.01)*100,expand = c(0,0),"new infections / 100 py") +
  coord_cartesian(ylim=c(0, 0.035)*100, xlim=c(1980, 2050)) +
  theme(text=element_text(size=24))+
  theme(strip.background = element_rect(colour="white", fill="white")) +
  theme(axis.text.x = element_text(angle = 90)) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)+
  theme(legend.title=element_blank(),axis.line = element_line(colour = "black"),panel.background = element_blank(),legend.position = "bottom")

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\Transmission in Era of UTT\\Figures")
ggsave("HIVincidence_by_src.stage_gender_baseline.jpg", height=10, width=15)

#ART 100%
ggplot(subset(transm.risk.table.both.genders, scenario=="ART100pct"),aes(x=year)) + 
  geom_line(aes(y=MedianInc*100, color=src.stage.f), size=2) +
  geom_ribbon(aes(ymin=lb*100, ymax=ub*100, fill=src.stage.f), size=2, alpha=0.2) +
  facet_grid(~Gender, labeller=labeller(Gender = labs)) +
  ggtitle("Incidence of HIV by source's stage of infection") +
  #scale_color_manual(values = c("blue", "red"),labels = c("Male", "Female"))+
  scale_x_continuous(expand = c(0,0), breaks=seq(1970,2045,5)) +
  scale_y_continuous(breaks=seq(0,0.05,0.01)*100,expand = c(0,0),"new infections / 100 py") +
  coord_cartesian(ylim=c(0, 0.035)*100, xlim=c(1980, 2050)) +
  theme(text=element_text(size=24))+
  theme(strip.background = element_rect(colour="white", fill="white")) +
  theme(axis.text.x = element_text(angle = 90)) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)+
  theme(legend.title=element_blank(),axis.line = element_line(colour = "black"),panel.background = element_blank(),legend.position = "bottom")

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\Transmission in Era of UTT\\Figures")
ggsave("HIVincidence_by_src.stage_gender_ART100pct.jpg", height=10, width=15)

# # # By age-cohort effect # # #

table(network.data.m$AgeCat) #age at acquisition - to match repotHIVbyageandgender
transmission.cohort <- subset(network.data.m, dest.age.int<60)  %>%
  group_by(scenario, sim.id, dest.gender, transm.year.int, AgeCat, src.stage.f) %>%
  summarize(transmissions=n()) 
transmission.cohort$year <- transmission.cohort$transm.year.int
transmission.cohort$Gender <- transmission.cohort$dest.gender
summary(transmission.cohort$AgeCat)

table(reporthivbyageandgender$AgeCat)
reporthivbyageandgender.cohort <- subset(reporthivbyageandgender, Age<60)  %>%
  group_by(scenario, sim.id, AgeCat, Gender, year) %>%
  summarize(Uninfected.Population=sum(Uninfected.Population)) 

transm.risk.table.cohort <- merge(transmission.cohort, reporthivbyageandgender.cohort, by=c("scenario","Gender","year","AgeCat","sim.id"), all.x=T)
transm.risk.table.cohort$inc <- transm.risk.table.cohort$transmissions/transm.risk.table.cohort$Uninfected.Population

transm.risk.table.cohort.inc <- subset(transm.risk.table.cohort)  %>%
  group_by(scenario, Gender, AgeCat, year, src.stage.f) %>%
  summarize(MedianInc=median(inc), lb = quantile(inc, probs = 0.025), ub = quantile(inc, probs = 0.975)) 
transm.risk.table.cohort.inc$Age <- as.numeric(str_sub(transm.risk.table.cohort$AgeCat, 2, 3))
names(transm.risk.table.cohort.inc)
transm.risk.table.cohort.inc$Age <- as.numeric(str_sub(transm.risk.table.cohort.inc$AgeCat, 2, 3))

transm.risk.table.cohort.freq <- subset(transm.risk.table.cohort) %>%
  group_by(scenario, sim.id, Gender, AgeCat, year) %>%
  mutate(freq = transmissions / sum(transmissions)) %>%
  group_by(scenario, Gender, AgeCat, year, src.stage.f) %>%
  summarize(Median=median(freq), lb = quantile(freq, probs = 0.025), ub = quantile(freq, probs = 0.975)) 
transm.risk.table.cohort.freq$Age <- as.numeric(str_sub(transm.risk.table.cohort.freq$AgeCat, 2, 3))

max.stage <- transm.risk.table.cohort.freq %>% group_by(scenario,Gender,AgeCat,year) %>% slice(which.max(Median))

ggplot(max.stage, aes(year, Age)) + 
  geom_tile(aes(fill = src.stage.f),colour = "white") + 
  facet_grid(~Gender, labeller=labeller(Gender = labs)) +
  scale_fill_manual(values=c("red", "blue", "black","orange")) 

ggplot(subset(max.stage, scenario=="baseline"), aes(year,Age)) + 
  geom_tile(aes(fill = src.stage.f)) + 
  facet_grid(scenario~Gender, labeller=labeller(Gender = labs)) +
  scale_x_continuous(expand = c(0,0), breaks=seq(1970,2050,5)) +
  theme(text=element_text(size=24))+
  scale_y_continuous(breaks=seq(15,55,5),expand = c(0,0),"Age") +
  theme(strip.background = element_rect(colour="white", fill="white")) +
  theme(axis.text.x = element_text(angle = 90)) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)+
  theme(legend.title=element_blank(),axis.line = element_line(colour = "black"),
        panel.background = element_blank())

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\Transmission in Era of UTT\\Figures")
ggsave("Stage_most_common_incidence_baseline.jpg", height=8, width=15)

ggplot(subset(max.stage, scenario=="ART100pct"), aes(year,Age)) + 
  geom_tile(aes(fill = src.stage.f)) + 
  facet_grid(scenario~Gender, labeller=labeller(Gender = labs)) +
  scale_x_continuous(expand = c(0,0), breaks=seq(1970,2050,5)) +
  theme(text=element_text(size=24))+
  scale_y_continuous(breaks=seq(15,55,5),expand = c(0,0),"Age") +
  theme(strip.background = element_rect(colour="white", fill="white")) +
  theme(axis.text.x = element_text(angle = 90)) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)+
  theme(legend.title=element_blank(),axis.line = element_line(colour = "black"),
        panel.background = element_blank())

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\Transmission in Era of UTT\\Figures")
ggsave("Stage_most_common_incidence_ART100pct.jpg", height=8, width=15)

ggplot(subset(transm.risk.table.cohort.freq, scenario=="baseline" & src.stage.f=="Acute" & year>1990)) + 
  geom_raster(aes(y=Age, x=year, fill=Median*100))+ 
  scale_fill_gradientn(name="% of transmissions",
                       colours=c("white","red","darkred"),
                       limits=c(0,60)) +
  #geom_line(aes(y=MedianInc*100, color=src.stage.f), size=2) +
  #geom_ribbon(aes(ymin=lb*100, ymax=ub*100, fill=src.stage.f), size=2, alpha=0.2) +
  facet_grid(~Gender, labeller=labeller(Gender = labs)) +
  theme(legend.key.height = unit(2, "cm")) +
  #scale_color_manual(values = c("blue", "red"),labels = c("Male", "Female"))+
  #scale_x_continuous(expand = c(0,0), breaks=seq(1970,2045,5)) +
  #scale_y_continuous(breaks=seq(0,0.05,0.01)*100,expand = c(0,0),"new infections / 100 py") +
  #coord_cartesian(ylim=c(0, 0.035)*100, xlim=c(1980, 2050)) +
  theme(text=element_text(size=24))+
  theme(strip.background = element_rect(colour="white", fill="white")) +
  theme(axis.text.x = element_text(angle = 90)) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)+
  theme(legend.title=element_blank(),axis.line = element_line(colour = "black"),panel.background = element_blank())

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\Transmission in Era of UTT\\Figures")
ggsave("HIVincidence_by_acute_by_ageyear_freq_baseline.jpg", height=8, width=15)

ggplot(subset(transm.risk.table.cohort.inc, scenario=="baseline" & src.stage.f=="Acute" & year>2000)) + 
  geom_raster(aes(y=Age, x=year, fill=MedianInc*100))+ 
  scale_fill_gradientn(name="incidence attributed to stage",
                       colours=c("white","red","darkred"),
                       limits=c(0,1)) +
  #geom_line(aes(y=MedianInc*100, color=src.stage.f), size=2) +
  #geom_ribbon(aes(ymin=lb*100, ymax=ub*100, fill=src.stage.f), size=2, alpha=0.2) +
  facet_grid(~Gender, labeller=labeller(Gender = labs)) +
  theme(legend.key.height = unit(2, "cm")) +
  #scale_color_manual(values = c("blue", "red"),labels = c("Male", "Female"))+
  #scale_x_continuous(expand = c(0,0), breaks=seq(1970,2045,5)) +
  scale_y_continuous(breaks=seq(15,55,5),expand = c(0,0),"Age") +
  #coord_cartesian(ylim=c(0, 0.035)*100, xlim=c(1980, 2050)) +
  theme(text=element_text(size=24))+
  theme(strip.background = element_rect(colour="white", fill="white")) +
  theme(axis.text.x = element_text(angle = 90)) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)+
  theme(legend.title=element_blank(),axis.line = element_line(colour = "black"),panel.background = element_blank())

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\Transmission in Era of UTT\\Figures")
ggsave("HIVincidence_by_acute_by_ageyear_freq_baseline.jpg", height=8, width=15)


table(transm.risk.table.cohort.freq$src.stage.f)
ggplot(subset(transm.risk.table.cohort.freq, scenario=="baseline" & src.stage.f=="on ART" & year>2014)) + 
  geom_raster(aes(y=Age, x=year, fill=Median*100))+ 
  scale_fill_gradientn(name="% of transmissions",
                       colours=c("white","darkorchid1","darkorchid4"),
                       limits=c(0,70)) +
  #geom_line(aes(y=MedianInc*100, color=src.stage.f), size=2) +
  #geom_ribbon(aes(ymin=lb*100, ymax=ub*100, fill=src.stage.f), size=2, alpha=0.2) +
  facet_grid(~Gender, labeller=labeller(Gender = labs)) +
  theme(legend.key.height = unit(2, "cm")) +
  #scale_color_manual(values = c("blue", "red"),labels = c("Male", "Female"))+
  #scale_x_continuous(expand = c(0,0), breaks=seq(1970,2045,5)) +
  #scale_y_continuous(breaks=seq(0,0.05,0.01)*100,expand = c(0,0),"new infections / 100 py") +
  #coord_cartesian(ylim=c(0, 0.035)*100, xlim=c(1980, 2050)) +
  theme(text=element_text(size=24))+
  theme(strip.background = element_rect(colour="white", fill="white")) +
  theme(axis.text.x = element_text(angle = 90)) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)+
  theme(legend.title=element_blank(),axis.line = element_line(colour = "black"),panel.background = element_blank())

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\Transmission in Era of UTT\\Figures")
ggsave("HIVincidence_by_acute_by_ageyear_freq_baseline.jpg", height=8, width=15)

ggplot(subset(transm.risk.table.cohort.freq, scenario=="ART100pct" & src.stage.f=="Acute" & year>2014)) + 
  geom_raster(aes(y=Age, x=year, fill=Median*100))+ 
  scale_fill_gradientn(name="% of transmissions",
                       colours=c("white","red","darkred"),
                       limits=c(0,100)) +
  #geom_line(aes(y=MedianInc*100, color=src.stage.f), size=2) +
  #geom_ribbon(aes(ymin=lb*100, ymax=ub*100, fill=src.stage.f), size=2, alpha=0.2) +
  facet_grid(~Gender, labeller=labeller(Gender = labs)) +
  theme(legend.key.height = unit(2, "cm")) +
  #scale_color_manual(values = c("blue", "red"),labels = c("Male", "Female"))+
  #scale_x_continuous(expand = c(0,0), breaks=seq(1970,2045,5)) +
  #scale_y_continuous(breaks=seq(0,0.05,0.01)*100,expand = c(0,0),"new infections / 100 py") +
  #coord_cartesian(ylim=c(0, 0.035)*100, xlim=c(1980, 2050)) +
  theme(text=element_text(size=24))+
  theme(strip.background = element_rect(colour="white", fill="white")) +
  theme(axis.text.x = element_text(angle = 90)) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)+
  theme(legend.title=element_blank(),axis.line = element_line(colour = "black"),panel.background = element_blank())

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\Transmission in Era of UTT\\Figures")
ggsave("HIVincidence_by_acute_by_ageyear_freq_ART100pct.jpg", height=8, width=15)

### By Time from acquisition ###

network.data.m$TimeToTransmissionCat <- cut(network.data.m$TimeToTransmission, c(0,1,100), right=F, dig.lab=4)
table(network.data.m$TimeToTransmissionCat)

table(network.data.m$AgeCat) #age at acquisition - to match repotHIVbyageandgender
transmission.cohort <- subset(network.data.m, dest.age.int<60 & transm.year.int>2000)  %>%
  group_by(scenario, sim.id, dest.gender, transm.year.int, AgeCat, TimeToTransmissionCat) %>%
  summarize(transmissions=n()) 
transmission.cohort$year <- transmission.cohort$transm.year.int
transmission.cohort$Gender <- transmission.cohort$dest.gender
summary(transmission.cohort$AgeCat)

table(reporthivbyageandgender$AgeCat)
reporthivbyageandgender.cohort <- subset(reporthivbyageandgender, Age<60 & year > 2000)  %>%
  group_by(scenario, sim.id, AgeCat, Gender, year) %>%
  summarize(Uninfected.Population=sum(Uninfected.Population)) 

transm.risk.table.1yr <- merge(transmission.cohort, reporthivbyageandgender.cohort, by=c("scenario","Gender","year","AgeCat","sim.id"), all.x=T)
transm.risk.table.1yr$inc <- transm.risk.table.1yr$transmissions/transm.risk.table.1yr$Uninfected.Population

transm.risk.table.1yr.inc <- subset(transm.risk.table.1yr, TimeToTransmissionCat=="[0,1)" | TimeToTransmissionCat=="[1,100)")  %>%
  group_by(scenario, Gender, AgeCat, year, TimeToTransmissionCat) %>%
  summarize(MedianInc=median(inc), lb = quantile(inc, probs = 0.025), ub = quantile(inc, probs = 0.975)) 
transm.risk.table.1yr.inc$Age <- as.numeric(str_sub(transm.risk.table.1yr.inc$AgeCat, 2, 3))

transm.risk.table.1yr.freq <- subset(transm.risk.table.1yr, TimeToTransmissionCat=="[0,1)" | TimeToTransmissionCat=="[1,100)") %>%
  group_by(scenario, sim.id, Gender, AgeCat, year) %>%
  mutate(freq = transmissions / sum(transmissions)) %>%
  group_by(scenario, Gender, AgeCat, year, TimeToTransmissionCat) %>%
  summarize(Median=median(freq), lb = quantile(freq, probs = 0.025), ub = quantile(freq, probs = 0.975)) 
transm.risk.table.1yr.freq$Age <- as.numeric(str_sub(transm.risk.table.1yr.freq$AgeCat, 2, 3))

table(transm.risk.table.1yr.inc$scenario)
ggplot(subset(transm.risk.table.1yr.inc, scenario=="ART100pct" & year>2010 & TimeToTransmissionCat=="[0,1)")) + 
  geom_raster(aes(y=Age, x=year, fill=MedianInc*1000))+ 
  scale_fill_gradientn(name="incidence from new infections",
                       colours=c("darkblue", "blue", "white","red","darkred"),
                       limits=c(0,6)) +
  #geom_line(aes(y=MedianInc*100, color=src.stage.f), size=2) +
  #geom_ribbon(aes(ymin=lb*100, ymax=ub*100, fill=src.stage.f), size=2, alpha=0.2) +
  facet_grid(scenario~Gender, labeller=labeller(Gender = labs)) +
  theme(legend.key.height = unit(3, "cm")) +
  #scale_color_manual(values = c("blue", "red"),labels = c("Male", "Female"))+
  scale_x_continuous(expand = c(0,0), breaks=seq(2010,2050,2)) +
  scale_y_continuous(breaks=seq(15,55,5),"Age") +
  #coord_cartesian(ylim=c(0, 0.035)*100, xlim=c(1980, 2050)) +
  theme(text=element_text(size=24))+
  theme(strip.background = element_rect(colour="white", fill="white")) +
  theme(axis.text.x = element_text(angle = 90)) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)+
  theme(legend.title=element_blank(),axis.line = element_line(colour = "black"),panel.background = element_blank())

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\Transmission in Era of UTT\\Figures")
ggsave("HIVincidence_by_transmission1yr_since_infection_ART100pct.jpg", height=10, width=15)

head(transm.risk.table.1yr.freq)
ggplot(subset(transm.risk.table.1yr.freq, year>2015 & TimeToTransmissionCat=="[0,1)")) + 
  geom_raster(aes(y=Age, x=year, fill=Median*100))+ 
  scale_fill_gradientn(name="incidence from new infections",
                       colours=c("darkblue", "blue", "white","red","darkred"),
                       limits=c(0,100)) +
  #geom_line(aes(y=MedianInc*100, color=src.stage.f), size=2) +
  #geom_ribbon(aes(ymin=lb*100, ymax=ub*100, fill=src.stage.f), size=2, alpha=0.2) +
  facet_grid(scenario~Gender, labeller=labeller(Gender = labs)) +
  theme(legend.key.height = unit(6, "cm")) +
  #scale_color_manual(values = c("blue", "red"),labels = c("Male", "Female"))+
  scale_x_continuous(expand = c(0,0), breaks=seq(2016,2050,2)) +
  scale_y_continuous(breaks=seq(15,55,5),"Age") +
  #coord_cartesian(ylim=c(0, 0.035)*100, xlim=c(1980, 2050)) +
  theme(text=element_text(size=24))+
  theme(strip.background = element_rect(colour="white", fill="white")) +
  theme(axis.text.x = element_text(angle = 90)) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)+
  theme(legend.title=element_blank(),axis.line = element_line(colour = "black"),panel.background = element_blank())

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\Transmission in Era of UTT\\Figures")
ggsave("HIVincidence_freq_by_transmission1yr_since_infection.jpg", height=15, width=15)

ggplot(subset(transm.risk.table.1yr.inc, scenario=="baseline" & year==2020 & TimeToTransmissionCat=="[0,1)")) + 
  geom_line(aes(x=Age, y=MedianInc*100, color=Gender), size=2) +
  geom_ribbon(aes(x=Age, ymin=lb*100, ymax=ub*100, fill=Gender), size=2, alpha=0.2) +
  facet_grid(scenario~Gender, labeller=labeller(Gender = labs)) +
  theme(legend.key.height = unit(3, "cm")) +
  #scale_color_manual(values = c("blue", "red"),labels = c("Male", "Female"))+
  scale_x_continuous(expand = c(0,0), breaks=seq(2010,2050,2)) +
  scale_y_continuous(breaks=seq(15,55,5),"Age") +
  #coord_cartesian(ylim=c(0, 0.035)*100, xlim=c(1980, 2050)) +
  theme(text=element_text(size=24))+
  theme(strip.background = element_rect(colour="white", fill="white")) +
  theme(axis.text.x = element_text(angle = 90)) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)+
  theme(legend.title=element_blank(),axis.line = element_line(colour = "black"),panel.background = element_blank())


#################################################################################################
#Time to transmission by sex/age/year 
#################################################################################################
head(network.data.m)
network.data.m$TimeToTransmission <- network.data.m$transm.year.int - network.data.m$acq.year

summary(network.data.m$TimeToTransmission) #median 5.2 years to transmission (1.2-11.6 IQR)
summary(network.data.m$TimeToTransmission[network.data.m$src.gender==0]) #Men: median 4.7 years to transmission (1.0-10.3 IQR)
summary(network.data.m$TimeToTransmission[network.data.m$src.gender==1]) #Women: median 6 years to transmission (1.8-13.6 IQR)

#Proportion of HIV transmissions occuring soon after acqusision 1,2,3...50 years (transmission perspective)?
network.data.m$TimeToTransmissionCat <- cut(network.data.m$TimeToTransmission, c(0,1,5,10,50,100), right=F, dig.lab=4)
table(network.data.m$TimeToTransmissionCat)
table(network.data.m$era)

TimeToTransmissionProportions <- subset(network.data.m,  is.na(TimeToTransmissionCat)==F) %>%
  group_by(sim.id, era, src.age.int, src.gender, TimeToTransmissionCat) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  group_by(era, src.age.int, src.gender, TimeToTransmissionCat) %>%
  summarize(median = quantile(freq, probs = 0.5),
            lb = quantile(freq, probs = 0.025),
            ub = quantile(freq, probs = 0.975))
head(TimeToTransmissionProportions)

labs <- c("1" = "Women", "0" = "Men")
ggplot(subset(TimeToTransmissionProportions, TimeToTransmissionCat!="[50,100)"), aes(fill=forcats::fct_rev(TimeToTransmissionCat), y=median, x=src.age.int)) + 
  geom_bar(position="stack", stat="identity", width=1) +
  facet_grid(era~src.gender, labeller=labeller(src.gender = labs))+
  scale_x_continuous(expand = c(0,0), breaks=seq(15,49,5)) +  scale_y_continuous(expand = c(0,0), breaks=seq(0,1,0.1)) +
  coord_cartesian(ylim=c(0, 1),xlim=c(15, 50)) +
  theme_bw(base_size=16) +
  xlab("Age at HIV acquisition")+ ylab("Proportion") +
  scale_fill_manual(values = c("white", "blue3", "darkblue","black"),
                    labels = c("10+","5-9","1-4","<1"),name= "Years", guide = guide_legend(reverse = TRUE))+
  ggtitle("Proportion of secondary transmissions by time from index acquisition") +
  guides(fill = guide_legend(reverse=T, override.aes = list(colour = "black"))) +
  theme(legend.background = element_rect(fill="transparent"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        legend.position = "bottom")

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\Transmission in Era of UTT\\Figures")
ggsave("TimeToTransmission_Cat_AcquisitionPerspective.jpg", height=10, width=15)

labs <- c("1" = "Women", "0" = "Men")
head(TimeToTransmissionProportions)
ggplot(subset(TimeToTransmissionProportions, TimeToTransmissionCat=="[0,1)")) + 
  geom_line(aes(x=src.age.int, y=median, color=era)) +
  geom_ribbon(aes(x=src.age.int,ymin = lb, ymax = ub, fill = era, alpha=0.4)) +
  facet_grid(~src.gender, labeller=labeller(src.gender = labs))+
  scale_x_continuous(expand = c(0,0), breaks=seq(15,49,5)) + scale_y_continuous(expand = c(0,0), breaks=seq(0,1,0.1)) +
  coord_cartesian(ylim=c(0, 1),xlim=c(15, 50)) +
  theme_bw(base_size=16) +
  xlab("Age at HIV acquisition")+ ylab("Proportion") +
  ggtitle("Proportion of secondary transmissions by time from index acquisition") +
  guides(fill = guide_legend(reverse=T, override.aes = list(colour = "black"))) +
  theme(legend.background = element_rect(fill="transparent"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        legend.position = "bottom")

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\Transmission in Era of UTT\\Figures")
ggsave("TimeToTransmission_Cat_AcquisitionPerspective.jpg", height=10, width=15)

#Proportion of HIV transmissions occuring soon after acqusision 1,2,3...50 years (acquision perspective)?
network.data.m$TimeToTransmissionCat <- cut(network.data.m$TimeToTransmission, c(0,1,5,10,50,100), right=F, dig.lab=4)
table(network.data.m$TimeToTransmissionCat)
table(network.data.m$era)

TimeToTransmissionProportions.preUTT <- subset(network.data.m, era=="Pre-UTT" & is.na(TimeToTransmissionCat)==F) %>%
  group_by(dest.age.int, dest.gender, TimeToTransmissionCat) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))
TimeToTransmissionProportions.preUTT$era <- "pre-UTT"

TimeToTransmissionProportions.UTT <- subset(network.data.m, era=="UTT era" & is.na(TimeToTransmissionCat)==F) %>%
  group_by(dest.age.int, dest.gender, TimeToTransmissionCat) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))
TimeToTransmissionProportions.UTT$era <- "UTT era"

TimeToTransmissionProportions <- rbind(TimeToTransmissionProportions.preUTT,TimeToTransmissionProportions.UTT)

labs <- c("1" = "Women", "0" = "Men")
ggplot(subset(TimeToTransmissionProportions, era=="UTT era" & TimeToTransmissionCat!="[50,100)"), aes(fill=forcats::fct_rev(TimeToTransmissionCat), y=freq, x=dest.age.int)) + 
  geom_bar(position="stack", stat="identity", width=1) +
  facet_grid(~dest.gender, labeller=labeller(dest.gender = labs))+
  scale_x_continuous(expand = c(0,0), breaks=seq(15,49,5)) +
  scale_y_continuous(expand = c(0,0), breaks=seq(0,1,0.1)) +
  coord_cartesian(ylim=c(0, 1),xlim=c(15, 50)) +
  #geom_hline(yintercept=0.5, color='red', linetype=2) +
  theme_bw(base_size=16) +
  xlab("Age at HIV acquisition")+ ylab("Proportion") +
  scale_fill_manual(values = c("white", "blue3", "darkblue","black"),
                    labels = c("10+","5-9","1-4","<1"),name= "Years", guide = guide_legend(reverse = TRUE))+
  ggtitle("Proportion of secondary transmissions by time from index acquisition") +
  guides(fill = guide_legend(reverse=T, override.aes = list(colour = "black"))) +
  theme(legend.background = element_rect(fill="transparent"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        legend.position = "bottom")

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\Transmission in Era of UTT\\Figures")
ggsave("TimeToTransmission_Cat_AcquisitionPerspective.jpg", height=10, width=15)

#What is the perecent of transmissions from those on ART, initiated but off, never initiated (with a band showing 1 year from acquisition and acute transmission)
network.data.m$src.stage.f <- factor(network.data.m$src.stage,
                    levels = c(1,2,3,4),
                    labels = c("Acute","Latent","AIDS","on ART"))
table(network.data.m$src.stage.f)
TimeToTransmissionProportions <- subset(network.data.m, is.na(TimeToTransmissionCat)==F) %>%
  group_by(sim.id, era, src.age.cat, src.stage.f) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  group_by(era, src.age.cat, src.stage.f) %>%
  summarize(median = quantile(freq, probs = 0.5),
            lb = quantile(freq, probs = 0.025),
            ub = quantile(freq, probs = 0.975))
head(TimeToTransmissionProportions)

#Key: The stage of infection for the transmitting indivdiual. 
# 1 for untreated acute HIV, 2 for untreated latent HIV, 3 for untreated late/AIDS stage, 4 for on ART
head(TimeToTransmissionProportions)

ggplot(data=TimeToTransmissionProportions,aes(x=src.age.cat, y=median, color=src.stage.f, group=src.stage.f)) + 
  geom_point() +
  geom_ribbon(aes(ymin=lb, ymax=ub, fill=src.stage.f, group=src.stage.f), alpha=0.3) +
  facet_grid(~era, labeller=labeller(src.gender = labs))+
  #scale_x_discrete(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), breaks=seq(0,1,0.1)) +
  scale_x_discrete(labels=c("15-19","20-24","25-29","30-34","35-39","40-44","45+")) +
  #coord_cartesian(ylim=c(0, 1),xlim=c(15, 50)) +
  theme_bw(base_size=24) +
  xlab("Age at HIV transmission")+ ylab("Proportion") +
  #scale_color_manual(values = c("1"="Red","2"="Blue","3"="Purple","4"="Dark Green"),labels = c("Acute","Latent","AIDS","on ART"),name= "Stage")+
  #scale_fill_manual(values = c("1"="Red","2"="Blue","3"="Purple","4"="Dark Green"),labels = c("Acute","Latent","AIDS","on ART"),name= "Stage")+
  ggtitle("Proportion of transmissions by HIV stage of source") +
  guides(fill = guide_legend(reverse=T, override.aes = list(colour = "black"))) +
  theme(legend.background = element_rect(fill="transparent"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        legend.position = "bottom")

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\Transmission in Era of UTT\\Figures")
ggsave("Stage_of_HIV_transmission.jpg", height=10, width=15)

############################################################################################
#Time from HIV acquisition to transmission
############################################################################################

#What is the median time from acquisition to transmission by age,gender, at year of transmission).
TimeToTransmission <- subset(network.data.m)  %>%
  group_by(transm.year.cat, src.gender, src.age.cat3) %>%
  summarize(time.to.transmission = median(TimeToTransmission, na.rm = TRUE))

ggplot(data=subset(network.data.m), aes(x=src.age.cat, y=TimeToTransmission, fill=factor(src.gender))) +
  geom_boxplot(notch=FALSE, outlier.shape=NA,outlier.size = 0, coef=0, width=0.3, alpha=0.3)+
  xlab("Age at transmission")+ ylab("Years")+
  facet_wrap(~era) +
  theme_bw(base_size=16) +
  scale_fill_manual(values = c("blue3","red2"), labels=c("Men","Women")) +
  scale_y_continuous(expand = c(0,0), breaks=seq(0,50,1)) +
  coord_cartesian(ylim=c(0, 35)) +
  ggtitle("Time from acquisition to transmission") +
  guides(color=guide_legend(title="", nrow=1))+
  theme(legend.position="bottom") +
  #theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(strip.background = element_rect(colour="black", fill="white")) +
  theme(legend.background = element_rect(fill="transparent"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"))

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\Transmission in Era of UTT\\Figures")
ggsave("TimeToTransmission_TransmissionPerspective.jpg", height=8, width=12)

#What is the median time from acquisition to transmission by age,gender, by year/age of acquisition)?
head(network.data.m)
ggplot(data=subset(network.data.m, acq.age>14), aes(x=acq.age.cat, y=TimeToTransmission, fill=factor(src.gender))) +
  geom_boxplot(notch=FALSE, outlier.shape=NA,outlier.size = 0, coef=0, width=0.3, alpha=0.3)+
  xlab("Age at acquisition")+ ylab("Years")+
  facet_wrap(~era) +
  theme_bw(base_size=16) +
  scale_fill_manual(values = c("blue3","red2"), labels=c("Men","Women")) +
  scale_y_continuous(expand = c(0,0), breaks=seq(0,50,1)) +
  coord_cartesian(ylim=c(0, 25)) +
  ggtitle("Time from HIV acquisition until transmission (among transmitters)") +
  guides(color=guide_legend(title="", nrow=1))+
  theme(legend.position="bottom") +
  #theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(strip.background = element_rect(colour="black", fill="white")) +
  theme(legend.background = element_rect(fill="transparent"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"))

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\Transmission in Era of UTT\\Figures")
ggsave("TimeToTransmission_TransmissionPerspective_AgeYearAtAcquisition.jpg", height=8, width=12)

#Time to transmission density
d2 <- subset(network.data.m, acq.age>14) %>%
  group_by(era, src.gender) %>%
  summarize(median = quantile(TimeToTransmission, probs = 0.5))

labs <- c("1" = "Women", "0" = "Men")
head(network.data.m)
ggplot(subset(network.data.m, acq.age>14), aes(x=TimeToTransmission, fill=factor(src.gender), color=factor(src.gender))) +
  geom_density(alpha=0.3, adjust=3) +
  facet_grid(src.gender~era, labeller=labeller(src.gender = labs)) +
  scale_color_manual(values = c("blue3","red2"), labels=c("Men","Women")) +
  scale_fill_manual(values = c("blue3","red2"), labels=c("Men","Women")) +
  ggtitle("Distribution of time to transmission") +
  theme_bw(base_size=24) +
  xlab("Years")+ 
  scale_x_continuous(expand = c(0,0),breaks=seq(0,35,5)) +
  scale_y_continuous(expand = c(0,0)) +
  coord_cartesian(xlim=c(0, 40)) +
  geom_vline(data = d2, aes(xintercept = median, color=factor(src.gender))) +
  theme(panel.border = element_blank(), 
        strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\Transmission in Era of UTT\\Figures")
ggsave("TimeToTransmission_Density_Allages_with_median.jpg", height=10, width=15)

head(network.data.m)
network.data.m$acq.age.cat5 <- cut(network.data.m$acq.age, c(15,20,25,30,35,200), right=F)
table(network.data.m$acq.age.cat5)
labs.age <- c("[15,20)"="15-19","[20,25)"="20-24","[25,30)"="25-29","[30,35)"="30-34","[35,200)"="35+")

d3 <- subset(network.data.m, acq.age>14) %>%
  group_by(era, src.gender, acq.age.cat5) %>%
  summarize(median = quantile(TimeToTransmission, probs = 0.5))

#distributino of time from acquisition to transmission by age when the transmitter acquired HIV
head(network.data.m)
ggplot(subset(network.data.m, acq.age>14), aes(x=TimeToTransmission, fill=factor(src.gender))) +
  geom_density(alpha=0.3, adjust=3) +
  facet_grid(era~acq.age.cat5,labeller=labeller(acq.age.cat5 = labs.age)) +
  scale_color_manual(values = c("blue3","red2"), labels=c("Men","Women")) +
  scale_fill_manual(values = c("blue3","red2"), labels=c("Men","Women")) +
  ggtitle("Time from HIV acquisition to transmission by age of transmitter at HIV acquisition") +
  theme_bw(base_size=24) +
  geom_vline(data = d3, aes(xintercept = median, color=factor(src.gender))) +
  xlab("Years")+ 
  scale_x_continuous(expand = c(0,0),breaks=seq(0,40,5)) +
  coord_cartesian(xlim=c(0, 45)) +
  scale_y_continuous(expand = c(0,0)) +
  theme(panel.border = element_blank(), 
        strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank())

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\Transmission in Era of UTT\\Figures")
ggsave("TimeToTransmission_Density_age_gender.jpg", height=10, width=20)

#How has median time to transmission changed over time?
network.data.m$acq.year.int <- floor(network.data.m$acq.year)
network.data.m$trans.year.int <- floor(network.data.m$transm.year)
head(network.data.m)
d4 <- subset(network.data.m, acq.age>14) %>%
  group_by(trans.year.int, src.gender) %>%
  summarize(median = quantile(TimeToTransmission, probs = 0.5),
            p25 = quantile(TimeToTransmission, probs = 0.25),
            p75 = quantile(TimeToTransmission, probs = 0.75))

labs <- c("1" = "Women", "0" = "Men")
ggplot(d4) +
  #geom_line(aes(x=trans.year.int, y=median, color=factor(src.gender))) +
  #geom_line(aes(x=trans.year.int, y=p25, color=factor(src.gender))) +
  #geom_line(aes(x=trans.year.int, y=p75, color=factor(src.gender))) +
  #geom_point(position = position_dodge(width = 0.05)) +
  geom_crossbar(aes(x=trans.year.int, y=median, color=factor(src.gender),ymin = p25, ymax = p75), width = 1) +
  facet_grid(~src.gender, labeller=labeller(src.gender = labs)) +
  scale_color_manual(values = c("blue3","red2"), labels=c("Men","Women")) +
  #scale_fill_manual(values = c("blue3","red2"), labels=c("Men","Women")) +
  ggtitle("Time from HIV acquisition to transmission") +
  theme_bw(base_size=24) +
  xlab("Years")+ ylab("median (IQR)") +
  coord_cartesian(xlim=c(1990, 2030), ylim=c(0,20)) +
  scale_y_continuous(expand = c(0,0), breaks=seq(0,30,5)) +
  theme(axis.text.x=element_text(angle=90,hjust=1)) +
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        legend.position = "none",
        legend.title = element_blank())

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\Transmission in Era of UTT\\Figures")
ggsave("Med_TimeToTransmission_byyearoftransmission.jpg", height=10, width=18)

head(trajectory_ART_calib.master.plot)
ggplot(d4) +
  #geom_line(aes(x=trans.year.int, y=median, color=factor(src.gender))) +
  #geom_line(aes(x=trans.year.int, y=p25, color=factor(src.gender))) +
  #geom_line(aes(x=trans.year.int, y=p75, color=factor(src.gender))) +
  #geom_point(position = position_dodge(width = 0.05)) +
  geom_smooth(data=subset(trajectory_ART_calib.master.plot,scenario_f=='baseline'),
              aes(x=Year, y=ART_coverage*100), color="black",linetype=2, method="loess",span=0.1, se=F, size=1.2) +
  geom_crossbar(aes(x=trans.year.int, y=median, color=factor(src.gender),ymin = p25, ymax = p75), width = 1) +
  facet_grid(~src.gender, labeller=labeller(src.gender = labs)) +
  scale_color_manual(values = c("blue3","red2"), labels=c("Men","Women")) +
  #scale_fill_manual(values = c("blue3","red2"), labels=c("Men","Women")) +
  ggtitle("Time from HIV acquisition to transmission") +
  geom_hline(yintercept=1, linetype=2)+
  theme_bw(base_size=24) +
  xlab("Years")+ ylab("median (IQR)") +
  coord_cartesian(xlim=c(1990, 2035), ylim=c(0,25)) +
  scale_y_continuous(expand = c(0,0), breaks=seq(0,30,5)) +
  theme(axis.text.x=element_text(angle=90,hjust=1)) +
  theme(strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        legend.position = "none",
        legend.title = element_blank())

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\Transmission in Era of UTT\\Figures")
ggsave("Med_TimeToTransmission_byyearoftransmission_dashed.jpg", height=10, width=18)

#Heatmap of proportion of transmissions by time to transmission and age 

network.data.m$yr1transmission <- cut(network.data.m$TimeToTransmission, c(0,1,100), right=F, dig.lab=4)
network.data.m$transm.year.cat5 <- cut(network.data.m$transm.year.int, seq(1980,2050,5), right=F, dig.lab=4)
table(network.data.m$transm.year.cat5)
summary(network.data.m$transm.year.int)
table(as.numeric(network.data.m$src.age.cat))
head(network.data.m)

TimeToTransmissionProportions.oneyear <- subset(network.data.m, is.na(yr1transmission)==F) %>%
  group_by(sim.id, transm.year.int, src.age.cat, src.gender, yr1transmission) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  group_by(transm.year.int, src.age.cat, src.gender, yr1transmission) %>%
  summarize(median = quantile(freq, probs = 0.5),
            lb = quantile(freq, probs = 0.025),
            ub = quantile(freq, probs = 0.975))

head(TimeToTransmissionProportions.oneyear)
summary(TimeToTransmissionProportions.oneyear$median)

#2d proprtion of <1 year transmissions
labs <- c("1" = "Women", "0" = "Men")
limits = c(0,max(TimeToTransmissionProportions.oneyear$median)*100)
breaks = seq(0,max(TimeToTransmissionProportions.oneyear$median)*100,10)
TimeToTransmissionProportions.oneyear$src.age.cat <- as.numeric(TimeToTransmissionProportions.oneyear$src.age.cat)
table(TimeToTransmissionProportions.oneyear$src.age.cat)

head(TimeToTransmissionProportions.oneyear)
ggplot() +
  geom_tile(data=subset(TimeToTransmissionProportions.oneyear, yr1transmission=="[0,1)"), aes(y=src.age.cat, x=as.numeric(transm.year.int), fill=median*100))+ 
  scale_fill_gradientn(name="%",
                       colours=c("blue","yellow","red"),
                       limits=limits,
                       breaks=breaks) +
  facet_grid(~src.gender, labeller=labeller(src.gender = labs)) +
  geom_vline(xintercept=2016, color='black', linetype=2, size=2) +
  #stat_contour(data=subset(TimeToTransmissionProportions.oneyear, yr1transmission=="[0,1)"), aes(y=as.numeric(src.age.cat), x=as.numeric(transm.year.int), z=median*100),  
  #            color="black", size=0.3, linetype=1, binwidth=20) +
  theme(legend.position="bottom") +  xlab("Year")+ ylab("Age")+
  scale_x_continuous(expand=c(0,0), breaks=seq(1980,2050,5)) +
  scale_y_continuous(expand=c(0,0), breaks=seq(1,7,1),labels=c("15-19","20-24","25-29","30-34","35-39","40-44","45+")) +
  theme_bw(base_size=36) +
  guides(fill = guide_legend(keywidth = 2, keyheight = 1, nrow=1)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(legend.position="bottom") +
  theme(strip.background = element_rect(colour="black", fill="white"))

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\Transmission in Era of UTT\\Figures")
ggsave("oneyeartransmissionheatmap.jpg", height=15, width=25)

ggplot() +
  geom_line(data=subset(TimeToTransmissionProportions.oneyear, yr1transmission=="[0,1)"), aes(y=median*100, x=as.numeric(transm.year.int), color=src.age.cat))+ 
  facet_grid(src.age.cat~src.gender, labeller=labeller(src.gender = labs)) +
  geom_vline(xintercept=2016, color='black', linetype=2, size=2) +
  theme(legend.position="bottom") +  xlab("Year")+ ylab("Age")+
  scale_x_continuous(expand=c(0,0), breaks=seq(1980,2050,5)) +
  scale_y_continuous(expand=c(0,0), breaks=seq(1,7,1),labels=c("15-19","20-24","25-29","30-34","35-39","40-44","45+")) +
  theme_bw(base_size=36) +
  guides(fill = guide_legend(keywidth = 2, keyheight = 1, nrow=1)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(legend.position="bottom") +
  theme(strip.background = element_rect(colour="black", fill="white"))

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\Transmission in Era of UTT\\Figures")

ggsave("oneyeartransmissionheatmap.jpg", height=15, width=25)

#from acquiring partner perspective - % of transmissions wthin one year of acquisition

TimeToTransmissionProportions.oneyear <- subset(network.data.m, is.na(yr1transmission)==F) %>%
  group_by(sim.id, transm.year.int, dest.age.cat, dest.gender, yr1transmission) %>% #change to dest
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  group_by(transm.year.int, dest.age.cat, dest.gender, yr1transmission) %>%
  summarize(median = quantile(freq, probs = 0.5),
            lb = quantile(freq, probs = 0.025),
            ub = quantile(freq, probs = 0.975))

head(TimeToTransmissionProportions.oneyear)
summary(TimeToTransmissionProportions.oneyear$median)

#2d proprtion of <1 year transmissions
labs <- c("1" = "Women", "0" = "Men")
limits = c(0,max(TimeToTransmissionProportions.oneyear$median)*100)
breaks = seq(0,max(TimeToTransmissionProportions.oneyear$median)*100,10)
TimeToTransmissionProportions.oneyear$dest.age.cat <- as.numeric(TimeToTransmissionProportions.oneyear$dest.age.cat)

head(TimeToTransmissionProportions.oneyear)
ggplot() +
  geom_tile(data=subset(TimeToTransmissionProportions.oneyear, yr1transmission=="[0,1)"), aes(y=dest.age.cat, x=as.numeric(transm.year.int), fill=median*100))+ 
  scale_fill_gradientn(name="%",
                       colours=c("blue","yellow","red"),
                       limits=limits,
                       breaks=breaks) +
  facet_grid(~dest.gender, labeller=labeller(src.gender = labs)) +
  geom_vline(xintercept=2016, color='black', linetype=2, size=2) +
  #stat_contour(data=subset(TimeToTransmissionProportions.oneyear, yr1transmission=="[0,1)"), aes(y=as.numeric(dest.age.cat), x=as.numeric(transm.year.int), z=median*100),  
  #            color="black", size=0.3, linetype=1, binwidth=20) +
  theme(legend.position="bottom") +  xlab("Year")+ ylab("Age")+
  scale_x_continuous(expand=c(0,0), breaks=seq(1980,2050,5)) +
  scale_y_continuous(expand=c(0,0), breaks=seq(1,7,1),labels=c("15-19","20-24","25-29","30-34","35-39","40-44","45+")) +
  theme_bw(base_size=36) +
  guides(fill = guide_legend(keywidth = 2, keyheight = 1, nrow=1)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(legend.position="bottom") +
  theme(strip.background = element_rect(colour="black", fill="white"))

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\Transmission in Era of UTT\\Figures")
ggsave("oneyeartransmissionheatmap_acqperspective.jpg", height=15, width=25)

#acute infection (within 3 months)
TimeToTransmissionProportions.acute <- subset(network.data.m, is.na(src.stage)==F) %>%
  group_by(sim.id, transm.year.int, src.age.cat, src.gender, src.stage) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  group_by(transm.year.int, src.age.cat, src.gender, src.stage) %>%
  summarize(median = quantile(freq, probs = 0.5),
            lb = quantile(freq, probs = 0.025),
            ub = quantile(freq, probs = 0.975))

head(TimeToTransmissionProportions.acute)
summary(TimeToTransmissionProportions.acute$median)

#2d proprtion of <1 year transmissions
labs <- c("1" = "Women", "0" = "Men")
limits = c(0,max(TimeToTransmissionProportions.oneyear$median)*100) #use same bounds as 1 year ests
breaks = seq(0,max(TimeToTransmissionProportions.oneyear$median)*100,10)
TimeToTransmissionProportions.acute$src.age.cat <- as.numeric(TimeToTransmissionProportions.acute$src.age.cat)
table(TimeToTransmissionProportions.acute$src.age.cat)

head(TimeToTransmissionProportions.acute)
ggplot() +
  geom_tile(data=subset(TimeToTransmissionProportions.acute, src.stage==1), aes(y=src.age.cat, x=as.numeric(transm.year.int), fill=median*100))+ 
  scale_fill_gradientn(name="%",
                       colours=c("blue","yellow","red"),
                       limits=limits,
                       breaks=breaks) +
  facet_grid(~src.gender, labeller=labeller(src.gender = labs)) +
  geom_vline(xintercept=2016, color='black', linetype=2, size=2) +
  #stat_contour(data=subset(TimeToTransmissionProportions.oneyear, yr1transmission=="[0,1)"), aes(y=as.numeric(src.age.cat), x=as.numeric(transm.year.int), z=median*100),  
  #            color="black", size=0.3, linetype=1, binwidth=20) +
  theme(legend.position="bottom") +  xlab("Year")+ ylab("Age")+
  scale_x_continuous(expand=c(0,0), breaks=seq(1980,2050,5)) +
  scale_y_continuous(expand=c(0,0), breaks=seq(1,7,1),labels=c("15-19","20-24","25-29","30-34","35-39","40-44","45+")) +
  theme_bw(base_size=36) +
  guides(fill = guide_legend(keywidth = 2, keyheight = 1, nrow=1)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(legend.position="bottom") +
  theme(strip.background = element_rect(colour="black", fill="white"))

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\Transmission in Era of UTT\\Figures")
ggsave("acutetransmissionheatmap.jpg", height=15, width=25)

#Proportion of transmissions by time-cat

#Among those acquiring HIV (by age and gender) - how long had the transmitter been infected before infecting?
network.data.m %>%
  group_by(dest.gender, dest.age) %>%
  summarize(time.to.transmission = median(TimeToTransmission, na.rm = TRUE))

head(network.data.m)
ggplot(data=subset(network.data.m), aes(x=dest.age.cat3, y=TimeToTransmission, fill=factor(dest.gender))) +
  #geom_point(color=dest.age),size=1, alpha=0.5,position = position_dodge(width = 0.05)) +
  geom_boxplot(notch=FALSE, outlier.shape=NA,outlier.size = 0, coef=0, width=0.5)+
  xlab("Age of acquiring person")+ ylab("Years")+
  facet_wrap(~era) +
  theme_bw(base_size=16) +
  scale_fill_manual(values = c("white","grey"), labels=c("Men","Women")) +
  scale_y_continuous(expand = c(0,0), breaks=seq(0,100,1)) +
  coord_cartesian(ylim=c(0, 35)) +
  ggtitle("Time from acquisition to transmission (from acquiring partner perspective)") +
  #theme(legend.position=c(0.86,0.97)) +
  guides(color=guide_legend(title="", nrow=1))+
  theme(legend.position="bottom") +
  #theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(strip.background = element_rect(colour="black", fill="white")) +
  theme(legend.background = element_rect(fill="transparent"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"))

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\Transmission in Era of UTT\\Figures")
ggsave("TimeToTransmission_AcquisitionPerspective.jpg", height=8, width=12)

#################################################################################################
#Time from infection until transmission by age and sex year and scenario + Proportion acute / early
#################################################################################################

#CREATE DATA FRAME OF DESTINATION TRANSMISSIONS AS THE SOURCE OF FUTURE TRANSMISSIONS
source.df <- data.frame("sim.id"=transmissionreport$sim.id, "SRC_ID" = transmissionreport$DEST_ID, "src.gender"=transmissionreport$DEST_GENDER, "year.acq" = transmissionreport$year, "age.at.acq" = floor(transmissionreport$dest_age_int))
source.df <- source.df[order(source.df$sim.id, source.df$SRC_ID, source.df$year.acq),]

#merge in source infection date and age and keep those who did not transmit
transmissionreport.merge <- merge(transmissionreport, source.df, by=c("sim.id","SRC_ID")) #only keep matching records (missing ones will be those seeded in and transmissions from those infected at birth) 
transmissionreport.merge <- transmissionreport.merge[order(transmissionreport.merge$sim.id, transmissionreport.merge$SRC_ID, transmissionreport.merge$year),]
head(transmissionreport.merge,10)
transmissionreport.merge$yrs.to.transmit <- transmissionreport.merge$year - transmissionreport.merge$year.acq
transmissionreport.merge$yrs.to.transmit.cat <- cut(transmissionreport.merge$yrs.to.transmit, c(0,1,5,200), right=F)
transmissionreport.merge$transmission <- 1
transmissionreport.merge$acute <- 0
transmissionreport.merge$acute[transmissionreport.merge$SRC_STAGE==1] <- 1
transmissionreport.merge$early <- 0
transmissionreport.merge$early[transmissionreport.merge$yrs.to.transmit.cat=="[0,1)"] <- 1
head(transmissionreport.merge,100)
transmissionreport.merge$year.fl <- floor(transmissionreport.merge$year)

#aggregate the number of transmissions by sim.id and year
transmissionreport.acute <- aggregate(cbind(transmission,acute) ~ sim.id+year.fl+SRC_GENDER, transmissionreport.merge, sum) #number of transmissions by gender, year
transmissionreport.acute$stage <- "acute"
transmissionreport.acute$p.transmission <- transmissionreport.acute$acute / transmissionreport.acute$transmission
head(transmissionreport.acute)
transmissionreport.acute <- transmissionreport.acute[,c(1:3,6,7)]
transmissionreport.early <- aggregate(cbind(transmission,early) ~ sim.id+year.fl+SRC_GENDER, transmissionreport.merge, sum) #number of transmissions by gender, year
transmissionreport.early$stage <- "early"
transmissionreport.early$p.transmission <- transmissionreport.early$early / transmissionreport.early$transmission
transmissionreport.early <- transmissionreport.early[,c(1:3,6,7)]
transmissions.acute.early <- rbind(transmissionreport.early,transmissionreport.acute)
head(transmissions.acute.early,5)
table(transmissions.acute.early$sim.id)

#Proportion of transmissions by time-cat
labs <- c("1" = "Women", "0" = "Men")
head(transmissions.acute.early)
p.transmission.time.cat <- ggplot(transmissions.acute.early) +
  geom_smooth(aes(x=year.fl, y=p.transmission, color=stage),method="loess", span=0.3, se = F, size=1.2, linetype=1) +
  geom_line(aes(x=year.fl, y=p.transmission, color=stage, group=interaction(sim.id, stage)), alpha=0.02)+ 
  facet_grid(~SRC_GENDER, labeller=labeller(SRC_GENDER = labs)) +
  xlab("Year") +
  ylab("Proportion of transmissions") +
  scale_x_continuous(limits = c(1984, 2050), breaks=seq(1985,2045,5),expand=c(0,0)) +
  scale_y_continuous(limits = c(0, 0.6), breaks=seq(0,1,0.1),expand=c(0,0)) +
  # guides(fill = guide_legend(keywidth = 2, keyheight = 1)) +
  scale_color_manual(labels = c("acute","early"), values = c("red","blue")) +
  theme_bw(base_size=20) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(legend.position="bottom") +
  theme(legend.title=element_blank())+
  theme(strip.background = element_rect(colour="black", fill="white")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p.transmission.time.cat

loess <- loess(p.transmission ~ year.fl, subset(transmissions.acute.early, SRC_GENDER==0 & stage=="acute"), span=0.3)
predict(loess, data.frame(year.fl = 2030, se = F))
loess <- loess(p.transmission ~ year.fl, subset(transmissions.acute.early, SRC_GENDER==0 & stage=="early"), span=0.3)
predict(loess, data.frame(year.fl = 2030, se = F))

loess <- loess(p.transmission ~ year.fl, subset(transmissions.acute.early, SRC_GENDER==1 & stage=="acute"), span=0.3)
predict(loess, data.frame(year.fl = 2030, se = F))
loess <- loess(p.transmission ~ year.fl, subset(transmissions.acute.early, SRC_GENDER==1 & stage=="early"), span=0.3)
predict(loess, data.frame(year.fl = 2030, se = F))

###########################################################################################################
#Transmission pair formation
###########################################################################################################
head(transmissionreport,5)
table(transmissionreport$scenario)

#Set the gender of the destination:
transmissionreport$Gender[transmissionreport$SRC_GENDER==0] <- "Female"
transmissionreport$Gender[transmissionreport$SRC_GENDER==1] <- "Male"
transmissionreport$year.cat <- cut(transmissionreport$year, c(1980,1990,2000,2010,2020,2030,2040,2050), dig.lab=4, right=F)
transmissionreport$year.fl <- floor(transmissionreport$year)

#Proportion of all *acquisisions* in each year from each age/gender group
p.new.infection <- aggregate(acquisition ~ sim.id+year.fl+Gender, transmissionreport, FUN=sum)
p.new.infection.age.gender <- aggregate(acquisition ~ sim.id+year.fl+Gender+dest_age_int, transmissionreport, FUN=sum)
p.new.infection.age.gender.merge <- merge(p.new.infection.age.gender,p.new.infection, by=c("sim.id","year.fl","Gender"))
p.new.infection.age.gender.merge$prop <- p.new.infection.age.gender.merge$acquisition.x / p.new.infection.age.gender.merge$acquisition.y
head(p.new.infection)
head(p.new.infection.age.gender.merge)

#Calculate the median probability of transmissions by pair across all 250 sims
transmissionreport.n.median <- aggregate(prop ~ year.fl+Gender+dest_age_int, p.new.infection.age.gender.merge, median) #take the median of 250 sims
head(transmissionreport.n.median)
hiv.acquisitions <- transmissionreport.n.median
hiv.acquisitions$direction <- "acquisitions"

year1=2005
year2=2050
labs.gender <- c("1" = "Women", "0" = "Men")
p.age.dist.transmissions <- ggplot(subset(transmissionreport.n.median, year.fl >= year1 & year.fl <= year2)) +
  geom_smooth(aes(x=dest_age_int, y=prop*100, color=year.fl, group=year.fl),method="loess", span=0.8, se = F, size=1, linetype=1) +
  facet_grid(~Gender) +
  scale_colour_gradientn(colours=rainbow(2),guide = guide_colourbar(barwidth=30, barheight=1), breaks=seq(year1,year2,10)) +
  #scale_colour_gradientn(colours=rainbow(2),limits=c(year1, year2), breaks = c(year1,round((year1+year2)/2,2),year2), guide_colorbar(barheight = 30)) +
  xlab("Age") +
  ylab("Percent of new HIV acquisitions") +
  scale_x_continuous(limits = c(15, 75), breaks=seq(15,70,5),expand=c(0,0)) +
  scale_y_continuous(limits = c(0, 6),expand=c(0,0)) +
  theme_bw(base_size=20) +
  guides(fill = guide_legend(keywidth = 2, keyheight = 1)) +
  theme(legend.title=element_blank())+
  theme(legend.position="bottom") +
  theme(strip.background = element_rect(colour="black", fill="white")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p.age.dist.transmissions

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\BMGF Lancet Paper\\Final_Plots")
ggsave("swaziland_acquisitions_by_age_year_20052050.jpg", height=8, width=10)

#Proportion of all *transmissions* in each year from each age/gender group
p.new.infection <- aggregate(acquisition ~ sim.id+year.fl+SRC_GENDER, transmissionreport, FUN=sum)
p.new.infection.age.gender <- aggregate(acquisition ~ sim.id+year.fl+SRC_GENDER+src_age_int, transmissionreport, FUN=sum)
p.new.infection.age.gender.merge <- merge(p.new.infection.age.gender,p.new.infection, by=c("sim.id","year.fl","SRC_GENDER"))
p.new.infection.age.gender.merge$prop <- p.new.infection.age.gender.merge$acquisition.x / p.new.infection.age.gender.merge$acquisition.y
head(p.new.infection)
head(p.new.infection.age.gender.merge)

#Calculate the median probability of transmissions by pair across all 250 sims
transmissionreport.n.median <- aggregate(prop ~ year.fl+SRC_GENDER+src_age_int, p.new.infection.age.gender.merge, median) #take the median of 250 sims
head(transmissionreport.n.median)

hiv.transmissions <- transmissionreport.n.median
hiv.transmissions$direction <- "transmission"

year1=2005
year2=2050
labs.gender <- c("0" = "Male","1" = "Female")
p.age.dist.transmissions <- ggplot(subset(transmissionreport.n.median, year.fl >= year1 & year.fl <= year2)) +
  geom_smooth(aes(x=src_age_int, y=prop*100, color=year.fl, group=year.fl),method="loess", span=0.8, se = F, size=1, linetype=1) +
  facet_grid(~SRC_GENDER, labeller=labeller(SRC_GENDER = labs.gender)) +
  scale_colour_gradientn(colours=rainbow(2),guide = guide_colourbar(barwidth=30, barheight=1), breaks=seq(year1,year2,10)) +
  #scale_colour_gradientn(colours=rainbow(2),limits=c(year1, year2), breaks = c(year1,round((year1+year2)/2,2),year2), guide_colorbar(barheight = 30)) +
  xlab("Age") +
  ylab("Percent of new HIV transmissions") +
  scale_x_continuous(limits = c(15, 75), breaks=seq(15,70,5),expand=c(0,0)) +
  scale_y_continuous(limits = c(0, 6),expand=c(0,0)) +
  theme_bw(base_size=20) +
  guides(fill = guide_legend(keywidth = 2, keyheight = 1)) +
  theme(legend.title=element_blank())+
  theme(legend.position="bottom") +
  theme(strip.background = element_rect(colour="black", fill="white")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p.age.dist.transmissions

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\BMGF Lancet Paper\\Revised_Plots")
ggsave("swaziland_transmissions_by_age_year_20042050.jpg", height=8, width=10)


hiv.transmissions$age <- hiv.transmissions$src_age_int
hiv.transmissions$sex <- hiv.transmissions$SRC_GENDER
hiv.acquisitions$age <- hiv.acquisitions$dest_age_int
hiv.transmissions$sex <- "Female"
hiv.transmissions$sex[hiv.transmissions$SRC_GENDER==0] <- "Male"

head(hiv.transmissions)
head(hiv.acquisitions)
combined <- rbind(hiv.acquisitions[,c(1,4,5,6,7)], hiv.transmissions[,c(1,4,5,6,7)])
table(combined$sex)

year1=2005
year2=2050
labs.gender <- c("0" = "Male","1" = "Female")
p.age.combined <- ggplot(subset(combined, year.fl >= year1 & year.fl <= year2)) +
  geom_smooth(aes(x=age, y=prop*100, color=year.fl, group=year.fl),method="loess", span=0.8, se = F, size=1, linetype=1) +
  facet_grid(direction~sex) +
  scale_colour_gradientn(colours=rainbow(2),guide = guide_colourbar(barwidth=30, barheight=1), breaks=seq(year1,year2,10)) +
  xlab("Age") +
  ylab("Percent of new HIV acquisitions / transmissions") +
  scale_x_continuous(limits = c(15, 75), breaks=seq(15,70,5),expand=c(0,0)) +
  scale_y_continuous(limits = c(0, 6),expand=c(0,0)) +
  theme_bw(base_size=20) +
  guides(fill = guide_legend(keywidth = 2, keyheight = 1)) +
  theme(legend.title=element_blank())+
  theme(legend.position="bottom") +
  theme(strip.background = element_rect(colour="black", fill="white")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p.age.combined

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\BMGF Lancet Paper\\Revised_Plots")
ggsave("swaziland_transmissions_by_age_year_20042050.jpg", height=8, width=10)

agex.min = 15
agex.max = 49
agey.min = 15
agey.max = 54
limits = c(0,max(transmissionreport.n.median$prop))
breaks = seq(0, max(transmissionreport.n.median$prop), max(transmissionreport.n.median$prop)/6)
breaks=round(breaks, 3)

TransmissionMatrixPlot <-  ggplot() +
  geom_raster(data=subset(transmissionreport.n.median, Gender=="Male"), aes(y=dest_age_int, x=src_age_int, fill=prop*100))+ 
  scale_fill_gradientn(name="% of transmissions",
                       colours=c("darkblue","blue","yellow","darkred"),
                       limits=limits*100,
                       breaks=breaks*100) +
  # stat_contour(data=subset(transmissionreport.n.median, Gender=="Female" & year.cat=="[1980,2005)"), aes(y=dest_age_int, x=src_age_int, z=prop*100),
  #              color="black", size=0.1, linetype=1, binwidth=0.05) +
  facet_grid(~year) +
  geom_abline(intercept = 0, slope = 1, linetype=2) +
  theme(legend.position="bottom") +
  xlab("Age of source")+
  ylab("Age of destination")+
  scale_x_continuous(expand=c(0,0), breaks=seq(agex.min,agex.max,5), limits=c(agex.min,agex.max)) +
  scale_y_continuous(expand=c(0,0), breaks=seq(agey.min,agey.max,5), limits=c(agey.min,agey.max)) +
  theme_bw(base_size=14) +
  guides(fill = guide_legend(keywidth = 2, keyheight = 1, nrow=1)) +
  theme(legend.position="top") +
  theme(strip.background = element_rect(colour="black", fill="white")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
TransmissionMatrixPlot

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\BMGF Lancet Paper\\Final_Plots")
ggsave("swaziland_TransmissionByAgePairing_year_baseline_Male.jpg", height=5, width=17)
setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\BMGF Lancet Paper\\Final_Plots")
ggsave("transmissions.by.age.gender.year.jpg", height=8, width=8)

head(transmissionreport,100)
p.new.infection <- aggregate(acquisition ~ year+SRC_GENDER, transmissionreport, FUN=sum)
p.new.infection.age.gender <- aggregate(acquisition ~ year + src_age_int + SRC_GENDER, transmissionreport, FUN=sum)
p.new.infection.age.gender <- merge(p.new.infection.age.gender,p.new.infection, by=c("year", "SRC_GENDER"))
p.new.infection.age.gender$prop <- p.new.infection.age.gender$acquisition.x / p.new.infection.age.gender$acquisition.y
head(p.new.infection.age.gender)

year1=1985
year2=2050
labs <- c("1" = "Women", "0" = "Men")
p.age.dist.transmissions <- ggplot(subset(p.new.infection.age.gender, year >= year1 & year <= year2)) +
  geom_smooth(aes(x=src_age_int, y=prop, color=year, group=year),method="loess", span=0.6, se = F, size=1, linetype=1) +
  facet_grid(~SRC_GENDER, labeller=labeller(SRC_GENDER = labs)) +
  scale_colour_gradientn(colours=rainbow(2),limits=c(year1, year2), breaks = c(year1,round((year1+year2)/2,2),year2)) +
  xlab("Year") +
  ylab("Proportion of Transmissions") +
  scale_x_continuous(limits = c(15, 100), breaks=seq(15,80,10),expand=c(0,0)) +
  scale_y_continuous(limits = c(0, 0.06), breaks=seq(0,0.06,0.01),expand=c(0,0)) +
  theme_bw(base_size=20) +
  guides(fill = guide_legend(keywidth = 2, keyheight = 1)) +
  theme(legend.title=element_blank())+
  theme(legend.position="bottom") +
  theme(strip.background = element_rect(colour="black", fill="white")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p.age.dist.transmissions

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\BMGF Lancet Paper\\Final_Plots")
ggsave("transmissions.by.agecontinuous.gender.year.jpg", height=8, width=12)

#noART
scenario="noART"
limits = c(0,max(TransmissionMatrix.master.median$prop[TransmissionMatrix.master.median$scenario==scenario]*100))
breaks = seq(0, max(TransmissionMatrix.master.median$prop[TransmissionMatrix.master.median$scenario==scenario]*100), max(TransmissionMatrix.master.median$prop[TransmissionMatrix.master.median$scenario==scenario]*100)/6)
breaks=round(breaks, 2)
TransmissionMatrixPlot <-  ggplot() +
  geom_raster(data=subset(TransmissionMatrix.master.median, scenario=="noART"), aes(y=dest_age_int, x=src_age_int, fill=prop*100))+ 
  scale_fill_gradientn(name="Perecent of all transmissions",
                       colours=c("darkblue","blue","yellow","darkred"),
                       limits=limits,
                       breaks=breaks) +
  # stat_contour(data=subset(TransmissionMatrix.master.median, scenario==scenario), aes(y=dest_age_int, x=src_age_int, z=prop*100),  
  #              color="black", size=0.1, linetype=1, binwidth=0.1) +
  facet_grid(Gender~year.cat) +
  geom_abline(intercept = 0, slope = 1, linetype=2) +
  theme(legend.position="bottom") +
  xlab("Age of source")+
  ylab("Age of destination")+
  scale_x_continuous(expand=c(0,0), breaks=seq(agex.min,agex.max,5), limits=c(agex.min,agex.max)) +
  scale_y_continuous(expand=c(0,0), breaks=seq(agey.min,agey.max,5), limits=c(agey.min,agey.max)) +
  theme_bw(base_size=14) +
  guides(fill = guide_legend(keywidth = 2, keyheight = 1, nrow=1)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(legend.position="bottom") +
  theme(strip.background = element_rect(colour="black", fill="white"))
TransmissionMatrixPlot

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\BMGF Lancet Paper\\Final_Plots")
ggsave("swaziland_TransmissionByAgePairing_year_noART.jpg", height=10, width=18)

#noART
scenario="noARTnoVMMC"
limits = c(0,max(TransmissionMatrix.master.median$prop[TransmissionMatrix.master.median$scenario==scenario]*100))
breaks = seq(0, max(TransmissionMatrix.master.median$prop[TransmissionMatrix.master.median$scenario==scenario]*100), max(TransmissionMatrix.master.median$prop[TransmissionMatrix.master.median$scenario==scenario]*100)/6)
breaks=round(breaks, 2)
TransmissionMatrixPlot <-  ggplot() +
  geom_raster(data=subset(TransmissionMatrix.master.median, scenario=="noARTnoVMMC"), aes(y=dest_age_int, x=src_age_int, fill=prop*100))+ 
  scale_fill_gradientn(name="Perecent of all transmissions",
                       colours=c("darkblue","blue","yellow","darkred"),
                       limits=limits,
                       breaks=breaks) +
  # stat_contour(data=subset(TransmissionMatrix.master.median, scenario==scenario), aes(y=dest_age_int, x=src_age_int, z=prop*100),  
  #              color="black", size=0.1, linetype=1, binwidth=0.1) +
  facet_grid(Gender~year.cat) +
  geom_abline(intercept = 0, slope = 1, linetype=2) +
  theme(legend.position="bottom") +
  xlab("Age of source")+
  ylab("Age of destination")+
  scale_x_continuous(expand=c(0,0), breaks=seq(agex.min,agex.max,5), limits=c(agex.min,agex.max)) +
  scale_y_continuous(expand=c(0,0), breaks=seq(agey.min,agey.max,5), limits=c(agey.min,agey.max)) +
  theme_bw(base_size=14) +
  guides(fill = guide_legend(keywidth = 2, keyheight = 1, nrow=1)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(legend.position="bottom") +
  theme(strip.background = element_rect(colour="black", fill="white"))
TransmissionMatrixPlot

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\BMGF Lancet Paper\\Final_Plots")
ggsave("swaziland_TransmissionByAgePairing_year_noARTnoVMMC.jpg", height=10, width=18)

###############################################################################
#What is the average time until first transmission by age at acquisition and gender?
###############################################################################
names(transmissionreport.merge)
transmissionreport.merge <- transmissionreport.merge[order(transmissionreport.merge$SRC_ID, transmissionreport.merge$yrs.to.transmit),]
transmissionreport.merge$num <- ave(transmissionreport.merge$acquisition, transmissionreport.merge$SRC_ID, FUN=seq_along)
head(transmissionreport.merge,15)
time.to.first.transmission <- aggregate(yrs.to.transmit ~ src.gender + year.acq.int, subset(transmissionreport.merge, num==1), FUN=mean)
head(time.to.first.transmission)

time.to.first.transmission.plot <- ggplot() +
  geom_smooth(data = time.to.first.transmission, aes(x=year.acq.int, y=yrs.to.transmit, color=factor(src.gender), group=factor(src.gender)),method="loess", span=0.5, se = F, size=1.2, linetype=1) +
  geom_point(data = time.to.first.transmission, size=1.2, aes(x = year.acq.int, y=yrs.to.transmit, color=factor(src.gender), group=factor(src.gender)))+ 
  xlab("Year of HIV acquisition") +
  ylab("Mean time from acquisition to transmission (years)") +
  scale_x_continuous(limits = c(1985, 2030), breaks=seq(1985,2029,5),expand=c(0,0)) +
  scale_y_continuous(limits = c(0, 15), breaks=seq(0,15,1),expand=c(0,0)) +
  # guides(fill = guide_legend(keywidth = 2, keyheight = 1)) +
  scale_color_manual(labels = c("Men","Women"), values = c("blue","red")) +
  theme_bw(base_size=20) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(legend.position="bottom") +
  theme(legend.title=element_blank())+
  theme(strip.background = element_rect(colour="black", fill="white")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
time.to.first.transmission.plot

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\BMGF Lancet Paper\\Final_Plots")
ggsave("meantimefromacquisitiontotransmission.jpg", height=10, width=10)

names(transmissionreport.merge)
transmissionreport.merge <- transmissionreport.merge[order(transmissionreport.merge$SRC_ID, transmissionreport.merge$yrs.to.transmit),]
transmissionreport.merge$num <- ave(transmissionreport.merge$acquisition, transmissionreport.merge$SRC_ID, FUN=seq_along)
head(transmissionreport.merge,15)

time.to.first.transmission.plot <-
  ggplot(transmissionreport.merge, aes(x=factor(year.acq.int), y=yrs.to.transmit, fill=factor(src.gender)))+
  geom_boxplot(notch=FALSE, outlier.shape=NA, alpha=0.2,outlier.size = 0, coef=0)+
  xlab("Year of HIV acquisition") +
  ylab("Mean time from acquisition to transmission (years)") +
  facet_grid(~src.gender,labeller=labeller(src.gender = labs)) +
  #scale_x_continuous(limits = c(1985, 2030), breaks=seq(1985,2029,5),expand=c(0,0)) +
  scale_y_continuous(limits = c(0, 25), breaks=seq(0,25,5),expand=c(0,0)) +
  #guides(fill = guide_legend(keywidth = 2, keyheight = 1)) +
  #scale_color_manual(labels = c("Men","Women"), values = c("blue","red")) +
  scale_fill_manual(values = c("blue", "red")) +
  theme_bw(base_size=20) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(legend.position="bottom") +
  theme(legend.title=element_blank())+
  theme(strip.background = element_rect(colour="black", fill="white")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
time.to.first.transmission.plot

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\BMGF Lancet Paper\\Final_Plots")
ggsave("meantimefromacquisitiontotransmission.jpg", height=10, width=10)























