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

# baseline sweeps

##################################################################################################
#Bring in ART data
##################################################################################################

#condensed file
setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\HIV\\EMOD\\Swaziland_05June2018\\Scenarios_calibration_run3_v05Sep2018_newcampaign_FINAL\\ScenariosOct_FINAL")
reporthivbyageandgender.master <- read.csv("reporthivbyageandgender.master.final.condensed.csv")
table(reporthivbyageandgender.master$scenario)
table(reporthivbyageandgender.master$Agecat)
head(reporthivbyageandgender.master)
reporthivbyageandgender.master.final <- reporthivbyageandgender.master

names(reporthivbyageandgender.master.final)
table(reporthivbyageandgender.master.final$scenario)
table(reporthivbyageandgender.master.final$Agecat)

trajectory_ART_calib <- aggregate(cbind(On_ART, Infected) ~ Year+Agecat+Gender+sim.id+scenario,
                                  subset(reporthivbyageandgender.master.final, Year > 2000 & Year < 2050), FUN=sum)
trajectory_ART_calib$ART_coverage <- trajectory_ART_calib$On_ART / trajectory_ART_calib$Infected
trajectory_ART_calib.master<-trajectory_ART_calib
head(trajectory_ART_calib.master)
table(trajectory_ART_calib.master$Agecat)

#ART scale-up to 2016
trajectory_ART_calib.master.plot <- subset(trajectory_ART_calib.master, scenario=="baseline")

##################################################################################################
#Bring in reportHIVbyageandgender Baseline and 100% ART scenario for Eswatini model with ART stages
##################################################################################################

#read in reportHIVbyageandgender files
input_dir <- "C:\\Users\\aakullian\\Dropbox (IDM)\\GitHub\\EMOD_eswatini\\Output\\Baseline_transmission_ARTstate"
primary_dirs = list.files(input_dir)
transmission.data.frame <- data.frame("Year2"=seq(1980,2050,1))

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

names(temp_table)
head(transmission.data.frame)
table(transmission.data.frame$Age)
table(transmission.data.frame$Year2)

reporthivbyageandgender <- subset(transmission.data.frame) #create a dataframe to bring into transmission network specifics, Ages 15+
table(reporthivbyageandgender$Age)

reporthivbyageandgender$agecat3 <- cut(reporthivbyageandgender$Age, c(15,25,35,101), right=F)
table(reporthivbyageandgender$agecat3, reporthivbyageandgender$Age)
reporthivbyageandgender$yearcat3 <- cut(reporthivbyageandgender$Year2, c(1990,2004,2016,2030,2050), right=F)
reporthivbyageandgender$incidence <- reporthivbyageandgender$Newly.Infected / reporthivbyageandgender$Uninfected.Population
reporthivbyageandgender$riskoftransm <- reporthivbyageandgender$Transmitters / reporthivbyageandgender$Infected
head(reporthivbyageandgender)

reporthivbyageandgender.age.gender.year <- aggregate(cbind(Transmitters, Newly.Infected, Uninfected.Population,Infected,Population) ~ yearcat3+Gender+agecat3+sim.id, reporthivbyageandgender, FUN=sum)
reporthivbyageandgender.age.gender.year$incidence <- reporthivbyageandgender.age.gender.year$Newly.Infected / reporthivbyageandgender.age.gender.year$Uninfected.Population
reporthivbyageandgender.age.gender.year$riskoftransm <- reporthivbyageandgender.age.gender.year$Transmitters / reporthivbyageandgender.age.gender.year$Infected
head(reporthivbyageandgender.age.gender.year)

table(reporthivbyageandgender.age.gender.year$Gender)
table(reporthivbyageandgender.age.gender.year$agecat3)
head(reporthivbyageandgender.age.gender.year)
reporthivbyageandgender.age.gender.year.mean <- aggregate(riskoftransm ~ agecat3+Gender+yearcat3, reporthivbyageandgender.age.gender.year, FUN=mean)

#save an Rda with three data files: original report and two derived summary files
setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\GitHub\\EMOD_eswatini\\Output\\Baseline_transmission_ARTstate")
save(reporthivbyageandgender, reporthivbyageandgender.age.gender.year, reporthivbyageandgender.age.gender.year.mean, file = "reporthivbyageandgender.rda")

# # #

#read in 100% ART reportHIVbyageandgender files
input_dir <- "C:\\Users\\aakullian\\Dropbox (IDM)\\GitHub\\EMOD_eswatini\\Output\\ART100pct_transmission_ARTstate"
primary_dirs = list.files(input_dir)
transmission.data.frame <- data.frame("Year2"=seq(1980,2050,1))

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

reporthivbyageandgender.100pctART <- subset(transmission.data.frame) #create a dataframe to bring into transmission network specifics, Ages 15+
reporthivbyageandgender.100pctART$agecat3 <- cut(reporthivbyageandgender.100pctART$Age, c(15,25,35,101), right=F)
reporthivbyageandgender.100pctART$yearcat3 <- cut(reporthivbyageandgender.100pctART$Year2, c(1990,2004,2016,2030,2050), right=F)
reporthivbyageandgender.100pctART$incidence <- reporthivbyageandgender.100pctART$Newly.Infected / reporthivbyageandgender.100pctART$Uninfected.Population
reporthivbyageandgender.100pctART$riskoftransm <- reporthivbyageandgender.100pctART$Transmitters / reporthivbyageandgender.100pctART$Infected

reporthivbyageandgender.100pctART.age.gender.year <- aggregate(cbind(Transmitters, Newly.Infected, Uninfected.Population,Infected,Population) ~ yearcat3+Gender+agecat3+sim.id, reporthivbyageandgender.100pctART, FUN=sum)
reporthivbyageandgender.100pctART.age.gender.year$incidence <- reporthivbyageandgender.100pctART.age.gender.year$Newly.Infected / reporthivbyageandgender.100pctART.age.gender.year$Uninfected.Population
reporthivbyageandgender.100pctART.age.gender.year$riskoftransm <- reporthivbyageandgender.100pctART.age.gender.year$Transmitters / reporthivbyageandgender.100pctART.age.gender.year$Infected
reporthivbyageandgender.100pctART.age.gender.year.mean <- aggregate(riskoftransm ~ agecat3+Gender+yearcat3, reporthivbyageandgender.100pctART.age.gender.year, FUN=mean)

#save an Rda with three data files: original report and two derived summary files
setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\GitHub\\EMOD_eswatini\\Output\\ART100pct_transmission_ARTstate")
save(reporthivbyageandgender.100pctART, reporthivbyageandgender.100pctART.age.gender.year, reporthivbyageandgender.100pctART.age.gender.year.mean, file = "ART100pct_transmission_ARTstate.rda")

# # #

#Bring in newest downloaded baseline transmission data + reportHIVbyageandgender from same sim folders
input_dir <- "C:\\Users\\aakullian\\Dropbox (IDM)\\GitHub\\EMOD_eswatini\\Output\\00c54188-9a43-ea11-a2c3-c4346bcb1551"
test.data <- as.data.table(read.csv(file = paste(input_dir,"ReportHIVByAgeAndGender.csv",sep="/")))
test.tr.data <- as.data.table(read.csv(file = paste(input_dir,"TransmissionReport.csv",sep="/")))
transmission.table <- subset(test.tr.data)  %>%
  group_by(floor(YEAR),DEST_GENDER, SRC_STAGE) %>%
  summarize(transmissions=n()) 
demominator.table <- subset(test.data)  %>%
  group_by(Year,Gender) %>%
  summarize(Population=sum(Population), Infected=sum(Infected)) 

temp_table$Year2 <- floor((temp_table$Year-0.5))
temp_table$Uninfected.Population = temp_table$Population-temp_table$Infected


#################################################################################################
#Bring in Transmission data
#################################################################################################

#Documentation on all variables: https://idmod.org/docs/hiv/software-report-transmission.html?searchText=TransmissionReport

#bring in full saved dataset (bring in just the first sim with all vars for now)
setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\HIV\\EMOD\\Swaziland_05June2018\\Scenarios_calibration_run3_v05Sep2018_newcampaign_FINAL\\ScenariosOct_FINAL\\Baseline_TransmissionReport\\TransmissionReport\\")
transmissionreport.master <- readRDS("transmissionreport.master.rds")
head(transmissionreport.master)
#transmissionreport.master$sim.id <- 1 #placeholder for simid

#create variables
head(transmissionreport.master)
table(transmissionreport.master$sim.id)
transmissionreport <- transmissionreport.master[order(transmissionreport.master$sim.id, transmissionreport.master$SRC_ID, transmissionreport.master$YEAR),]
transmissionreport <- subset(transmissionreport, YEAR>=1980 & YEAR <2050)
summary(transmissionreport$YEAR)
transmissionreport$src_age_int <- floor(transmissionreport$SRC_AGE / 365.25)
transmissionreport$dest_age_int <- floor(transmissionreport$DEST_AGE / 365.25)
transmissionreport <- subset(transmissionreport, dest_age_int > 14 & src_age_int > 14)
transmissionreport$age.at.acquisition.cat <- cut(transmissionreport$dest_age_int, c(15,25,35,200), right=F)
transmissionreport$age.at.transmission.cat <- cut(transmissionreport$src_age_int, c(15,25,35,200), right=F)
transmissionreport$age.at.acquisition.cat5 <- cut(transmissionreport$dest_age_int, c(15,20,25,30,35,40,45,200), right=F)
transmissionreport$age.at.transmission.cat5 <- cut(transmissionreport$src_age_int, c(15,20,25,30,35,40,45,200), right=F)
transmissionreport$yearcat3 <- cut(transmissionreport$YEAR, c(1980,1990,2004,2016,2030,2050), right=F, dig.lab=4)

#################################################################################################
# Set up transmission data
#################################################################################################
library(igraph)

#Create transmission network for plotting (choose subset)
head(transmissionreport)
#network.data <- subset(transmissionreport, sim.id=="TransmissionReport_TPI0000_REP0001.csv") #first sim
network.data <- subset(transmissionreport) #all sims
table(network.data$sim.id)
date.of.acquisition <- data.frame("SRC_ID"=network.data$DEST_ID, "acq.year"=network.data$YEAR, "sim.id"=network.data$sim.id) #source id used to merge in with source of transmitters
summary(date.of.acquisition$acq.year)

network.data <- left_join(network.data, date.of.acquisition, by=c("SRC_ID","sim.id"))
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
                             "transm.year"=network.data$YEAR,
                             "transm.year.int"=floor(network.data$YEAR),
                             "acq.year"=network.data$acq.year,
                             "transm.year.cat"=network.data$yearcat3,
                             "simtime"=network.data$SIM_TIME,
                             "src.stage"=network.data$SRC_STAGE)
head(network.data.m)

nrow(network.data.m[is.na(network.data.m$acq.year)==T,]) #number of individuals who's infection was seeded in (see campaign.json)
#network.data.m$acq.year[is.na(network.data.m$acq.year)==T] <- 1982.095633 #fills in missing infection year with seed year (left out for now)
network.data.m$acq.age=floor(network.data.m$src.age.int-(network.data.m$transm.year-network.data.m$acq.year))
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

network.data.m <- network.data.m[order(network.data.m$simtime),]
head(network.data.m,20)
summary(network.data.m$transm.year)

#################################################################################################
# Bring in reportHIVbyageandgender rda files
#################################################################################################

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\GitHub\\EMOD_eswatini\\Output\\Baseline_transmission_ARTstate")
load(file="reporthivbyageandgender.rda") #load report HIVbyageandgenderfiles scenario 1 (baseline)
setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\GitHub\\EMOD_eswatini\\Output\\ART100pct_transmission_ARTstate")
load(file="ART100pct_transmission_ARTstate.rda") #load report HIVbyageandgenderfiles scenario 2

reporthivbyageandgender$scenario <- "1.baseline"
reporthivbyageandgender.100pctART$scenario <- "2.ART100pct"

#merge the two scenarios into one file:

reporthivbyageandgender <- rbind(reporthivbyageandgender, reporthivbyageandgender.100pctART)

reporthivbyageandgender$src.gender <- reporthivbyageandgender$Gender
reporthivbyageandgender$year <- reporthivbyageandgender$Year2
reporthivbyageandgender$NeverOnART <- 0
reporthivbyageandgender$NeverOnART[reporthivbyageandgender$IP_Key.ARTstate=="NeverOnART"] <- 1
table(reporthivbyageandgender$NeverOnART, reporthivbyageandgender$IP_Key.ARTstate)
head(reporthivbyageandgender)

#################################################################################################
# Final analyses / plots
#################################################################################################

#1) How has the percent of the infected pool who have never initiated ART changed over time?

#get percent of infected population by ART state

pct.neveronART <- subset(reporthivbyageandgender)  %>%
  group_by(scenario, sim.id, src.gender, year, IP_Key.ARTstate) %>% 
  summarize(count = sum(Infected)) %>% group_by(scenario, sim.id, src.gender, year) %>% 
  mutate(percent = count/sum(count)) %>% group_by(scenario, src.gender, year, IP_Key.ARTstate) %>%
  summarize(median = quantile(percent, probs = 0.5, na.rm=T),
            lb = quantile(percent, probs = 0.025, na.rm=T),
            ub = quantile(percent, probs = 0.975, na.rm=T))
head(pct.neveronART)  

pct.neveronART.age <- subset(reporthivbyageandgender)  %>%
  group_by(scenario,sim.id, src.gender, Age , year, IP_Key.ARTstate) %>% 
  summarize(count = sum(Infected)) %>% group_by(scenario,sim.id, src.gender, Age , year) %>% 
  mutate(percent = count/sum(count)) %>% group_by(scenario,src.gender, Age , year, IP_Key.ARTstate) %>%
  summarize(median = quantile(percent, probs = 0.5, na.rm=T),
            lb = quantile(percent, probs = 0.025, na.rm=T),
            ub = quantile(percent, probs = 0.975, na.rm=T))
head(pct.neveronART.age)  

ggplot(subset(pct.neveronART, IP_Key.ARTstate=="NeverOnART")) + 
  geom_line(aes(y=median, x=year, color=factor(src.gender)), size=2) +
  #geom_ribbon(aes(ymin=lb, ymax=ub, x=year, fill=factor(src.gender), group=factor(src.gender))) +
  facet_grid(~scenario) +
  scale_x_continuous(expand = c(0,0), breaks=seq(1980,2040,10)) +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw(base_size=24) +
  xlab("Year")+ ylab("proportion") + ggtitle("Proportion of PLHIV who have never initiated ART") +
  scale_color_manual(values = c("blue", "red"),labels = c("Male", "Female"))+
  theme(strip.background = element_rect(colour="black", fill="white")) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        legend.title=element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        legend.position = "bottom")

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\Transmission in Era of UTT\\Figures")
ggsave("NeverInitiatedART.jpg", height=10, width=16)

ggplot(subset(pct.neveronART.age, IP_Key.ARTstate=="NeverOnART" & scenario=="2.ART100pct" & Age<50)) + 
  geom_line(aes(y=1-median, x=year, color=factor(src.gender)), size=2) +
  geom_ribbon(aes(ymin=1-lb, ymax=1-ub, x=year, fill=factor(src.gender)), alpha=0.2) +
  scale_x_continuous(expand = c(0,0), breaks=seq(1980,2040,10)) + coord_cartesian(xlim=c(2010, 2050), ylim=c(0.5,1)) +
  scale_y_continuous(expand = c(0,0), breaks=seq(0,1,0.05)) +
  facet_grid(~Age) +
  theme_bw(base_size=24) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Year")+ ylab("proportion") +
  scale_color_manual(values = c("blue", "red"),labels = c("Male", "Female"))+
  scale_fill_manual(values = c("blue", "red"),labels = c("Male", "Female"))+
  theme(strip.background = element_rect(colour="black", fill="white")) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        legend.title=element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        legend.position = "bottom")

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\Transmission in Era of UTT\\Figures")
ggsave("InitiatedART_age_100pctART_v2.jpg", height=8, width=15)

ggplot(subset(pct.neveronART.age, IP_Key.ARTstate=="NeverOnART" & scenario=="1.baseline" & Age<50)) + 
  geom_line(aes(y=1-median, x=year, color=factor(src.gender)), size=2) +
  geom_ribbon(aes(ymin=1-lb, ymax=1-ub, x=year, fill=factor(src.gender)), alpha=0.2) +
  scale_x_continuous(expand = c(0,0), breaks=seq(1980,2040,10)) + coord_cartesian(xlim=c(2010, 2050), ylim=c(0,1)) +
  scale_y_continuous(expand = c(0,0), breaks=seq(0,1,0.05)) +
  facet_grid(~Age) +
  theme_bw(base_size=24) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Year")+ ylab("proportion") +
  scale_color_manual(values = c("blue", "red"),labels = c("Male", "Female"))+
  scale_fill_manual(values = c("blue", "red"),labels = c("Male", "Female"))+
  theme(strip.background = element_rect(colour="black", fill="white")) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        legend.title=element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        legend.position = "bottom")

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\Transmission in Era of UTT\\Figures")
ggsave("InitiatedART_age_baseline.jpg", height=8, width=15)

ggplot(subset(pct.neveronART.age, IP_Key.ARTstate=="NeverOnART" & Age<50)) + 
  geom_line(aes(y=1-median, x=year, color=factor(src.gender), linetype=scenario), size=2) +
  #geom_ribbon(aes(ymin=lb, ymax=ub, x=year, fill=factor(src.gender)), alpha=0.2) +
  scale_x_continuous(expand = c(0,0), breaks=seq(1980,2040,10)) + coord_cartesian(xlim=c(2009, 2050), ylim=c(0,1)) +
  scale_y_continuous(expand = c(0,0), breaks=seq(0,1,0.05)) +
  facet_grid(~Age) +
  theme_bw(base_size=24) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Year")+ ylab("proportion") +
  scale_color_manual(values = c("blue", "red"),labels = c("Male", "Female"))+
  scale_fill_manual(values = c("blue", "red"),labels = c("Male", "Female"))+
  theme(strip.background = element_rect(colour="black", fill="white")) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        legend.title=element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        legend.position = "bottom")

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\Transmission in Era of UTT\\Figures")
ggsave("NeverInitiatedART_age_baseline.jpg", height=8, width=15)


# # #

#2) How has the percent of the infected pool who are new infections changed over time?
reporthivbyageandgender.table.gender.year <- subset(reporthivbyageandgender)  %>%
  group_by(sim.id, scenario, src.gender, year) %>%
  summarize(Infected=sum(Infected), Population=sum(Population), NewlyInfected=sum(Newly.Infected), Prevalence=sum(Infected)/sum(Population), NewInfRatio=sum(Newly.Infected)/sum(Infected))  %>%
  group_by(scenario, src.gender, year) %>%
  summarize(medianInfected=median(Infected), lb.Infected=quantile(Infected, probs = 0.025), ub.Infected=quantile(Infected, probs = 0.975),
            medianPopulation=median(Population),
            medianPrevalence=median(Prevalence), 
            medianNewInfRatio=median(NewInfRatio),lb.NewInfRatio=quantile(NewInfRatio, probs = 0.025, na.rm=T),ub.NewInfRatio=quantile(NewInfRatio, probs = 0.975, na.rm=T))

reporthivbyageandgender.table.gender.year.age <- subset(reporthivbyageandgender)  %>%
  group_by(sim.id, scenario, src.gender, Age, year) %>%
  summarize(Infected=sum(Infected), Population=sum(Population), NewlyInfected=sum(Newly.Infected), Prevalence=sum(Infected)/sum(Population), NewInfRatio=sum(Newly.Infected)/sum(Infected))  %>%
  group_by(scenario, src.gender, Age, year) %>%
  summarize(medianInfected=median(Infected), lb.Infected=quantile(Infected, probs = 0.025), ub.Infected=quantile(Infected, probs = 0.975),
            medianPopulation=median(Population),
            medianPrevalence=median(Prevalence), 
            medianNewInfRatio=median(NewInfRatio),lb.NewInfRatio=quantile(NewInfRatio, probs = 0.025, na.rm=T),ub.NewInfRatio=quantile(NewInfRatio, probs = 0.975, na.rm=T))

labs <- c("1" = "Women", "0" = "Men")
head(transm.risk.table.gender)
ggplot(subset(reporthivbyageandgender.table.gender.year)) + 
  geom_line(aes(y=medianNewInfRatio, x=year, color=factor(src.gender)), size=2, linetype=1) +
  geom_ribbon(aes(ymin=lb.NewInfRatio, ymax=ub.NewInfRatio, x=year, fill=factor(src.gender)), alpha=0.3) +
  facet_grid(~scenario) +
  scale_x_continuous(expand = c(0,0), breaks=seq(1980,2050,5)) +
  scale_y_continuous(expand = c(0,0), breaks=seq(0,1,0.05)) +
  coord_cartesian(xlim=c(1990, 2048), ylim=c(0,0.3)) +
  theme_bw(base_size=24) +
  xlab("Year")+ ylab("Proportion") +
  scale_color_manual(values = c("blue", "red"),labels = c("Male", "Female")) + 
  scale_fill_manual(values = c("blue", "red"),labels = c("Male", "Female"))+
  ggtitle("Proportion of prevalent infections that are new infections") +
  geom_vline(xintercept = c(2004,2016), linetype=2) +
  theme(strip.background = element_rect(colour="black", fill="white")) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        plot.title = element_text(hjust = 0.5), legend.title=element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        legend.position = "bottom")

ggplot(subset(reporthivbyageandgender.table.gender.year.age, Age<50 & scenario=="1.baseline")) + 
  geom_line(aes(y=medianNewInfRatio, x=year, color=factor(src.gender)), size=2, linetype=1) +
  geom_ribbon(aes(ymin=lb.NewInfRatio, ymax=ub.NewInfRatio, x=year, fill=factor(src.gender)), alpha=0.3) +
  scale_x_continuous(expand = c(0,0), breaks=seq(1980,2040,10)) + scale_y_continuous(expand = c(0,0), breaks=seq(0,1,0.05)) +
  coord_cartesian(xlim=c(1990, 2048), ylim=c(0,0.75)) +
  facet_grid(~Age) +
  theme_bw(base_size=24) +
  xlab("Year")+ ylab("Proportion") +
  scale_color_manual(values = c("blue", "red"),labels = c("Male", "Female")) + 
  scale_fill_manual(values = c("blue", "red"),labels = c("Male", "Female"))+
  geom_vline(xintercept = c(2016), linetype=2) +
  theme(axis.text.x=element_text(angle=90,hjust=1)) +
  theme(strip.background = element_rect(colour="black", fill="white")) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        plot.title = element_text(hjust = 0.5), legend.title=element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        legend.position = "bottom")

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\Transmission in Era of UTT\\Figures")
ggsave("InfectionsThatAreNew_baseline_age.jpg", height=8, width=15)

ggplot(subset(reporthivbyageandgender.table.gender.year.age, Age<50)) + 
  geom_line(aes(y=medianNewInfRatio, x=year, color=factor(src.gender), linetype=scenario), size=2) +
  #geom_ribbon(aes(ymin=lb.NewInfRatio, ymax=ub.NewInfRatio, x=year, fill=factor(src.gender)), alpha=0.3) +
  scale_x_continuous(expand = c(0,0), breaks=seq(1980,2040,10)) + scale_y_continuous(expand = c(0,0), breaks=seq(0,1,0.05)) +
  coord_cartesian(xlim=c(1990, 2048), ylim=c(0,0.75)) +
  facet_grid(~Age) +
  theme_bw(base_size=24) +
  xlab("Year")+ ylab("Proportion") +
  scale_color_manual(values = c("blue", "red"),labels = c("Male", "Female")) + 
  scale_fill_manual(values = c("blue", "red"),labels = c("Male", "Female"))+
  geom_vline(xintercept = c(2016), linetype=2) +
  theme(axis.text.x=element_text(angle=90,hjust=1)) +
  theme(strip.background = element_rect(colour="black", fill="white")) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        plot.title = element_text(hjust = 0.5), legend.title=element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        legend.position = "bottom")

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\Transmission in Era of UTT\\Figures")
ggsave("InfectionsThatAreNew_baseline_age.jpg", height=8, width=15)

table(reporthivbyageandgender.table.gender.year.age$scenario)
ggplot(subset(reporthivbyageandgender.table.gender.year.age, Age<50 & scenario=="2.ART100pct")) + 
  geom_line(aes(y=medianNewInfRatio, x=year, color=factor(src.gender)), size=2, linetype=1) +
  geom_ribbon(aes(ymin=lb.NewInfRatio, ymax=ub.NewInfRatio, x=year, fill=factor(src.gender)), alpha=0.3) +
  scale_x_continuous(expand = c(0,0), breaks=seq(1980,2040,10)) + scale_y_continuous(expand = c(0,0), breaks=seq(0,1,0.05)) +
  coord_cartesian(xlim=c(2009, 2048), ylim=c(0,0.75)) +
  facet_grid(~Age) +
  theme_bw(base_size=24) +
  xlab("Year")+ ylab("Proportion") +
  scale_color_manual(values = c("blue", "red"),labels = c("Male", "Female")) + 
  scale_fill_manual(values = c("blue", "red"),labels = c("Male", "Female"))+
  geom_vline(xintercept = c(2016), linetype=2) +
  theme(axis.text.x=element_text(angle=90,hjust=1)) +
  theme(strip.background = element_rect(colour="black", fill="white")) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        plot.title = element_text(hjust = 0.5), legend.title=element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        legend.position = "bottom")

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\Transmission in Era of UTT\\Figures")
ggsave("InfectionsThatAreNew_100pctART_age.jpg", height=8, width=15)

#By gender over time
transmission.table.gender.year <- subset(network.data.m)  %>%
  group_by(sim.id, src.gender, transm.year.int) %>%
  summarize(transmissions=n()) %>%
  group_by(src.gender, transm.year.int) %>%
  summarize(medianTransmissions=median(transmissions),lb.Transmissions = quantile(transmissions, probs = 0.025),
            ub.Transmissions = quantile(transmissions, probs = 0.975))  
transmission.table.gender.year$year <- transmission.table.gender.year$transm.year.int
summary(transmission.table.gender.year$year)

#By Age, gender, year
head(network.data.m)
transmission.table.age.gender.year <- subset(network.data.m)  %>%
  group_by(sim.id, src.gender, src.age.cat, transm.year.int) %>%
  summarize(transmissions=n()) %>%
  group_by(src.gender, src.age.cat, transm.year.int) %>%
  summarize(medianTransmissions=median(transmissions))  
transmission.table.age.gender.year$year <- transmission.table.age.gender.year$transm.year.int
head(transmission.table.age.gender.year)

# # #

# 3) How has the # transmissions per infected changed over time?
#bring in reportHIVbyageandgender for persontime at risk of transmission (infection years)

transm.risk.table.gender <- merge(transmission.table.gender.year, reporthivbyageandgender.table.gender.year, by=c("src.gender","year"), all.y=T)
transm.risk.table.gender$tr.risk <- transm.risk.table.gender$medianTransmissions / transm.risk.table.gender$medianInfected
transm.risk.table.gender$tr.risk.lb <- transm.risk.table.gender$lb.Transmissions / transm.risk.table.gender$lb.Infected
transm.risk.table.gender$tr.risk.ub <- transm.risk.table.gender$ub.Transmissions / transm.risk.table.gender$ub.Infected
table(transm.risk.table.gender$medianNeverOnARTRatio)

labs <- c("1" = "Women", "0" = "Men")
ggplot(subset(transm.risk.table.gender)) + 
  geom_line(aes(y=tr.risk, x=year, color=factor(src.gender), group=src.gender), size=2) +
  #geom_ribbon(aes(ymin=tr.risk.lb, ymax=tr.risk.ub, x=year, fill=factor(src.gender), group=src.gender)) +
  #geom_line(aes(y=medianPrevalence, x=year, color=factor(src.gender), group=src.gender), size=2, linetype=2) +
  scale_x_continuous(expand = c(0,0), breaks=seq(1980,2050,5)) +
  scale_y_continuous(expand = c(0,0), breaks=seq(0,0.6,0.05)) + 
  coord_cartesian(xlim=c(1990, 2050), ylim=c(0, 0.6)) +
  theme_bw(base_size=16) +
  xlab("Year")+ ylab("Number of transmissions per infected individual") +
  scale_color_manual(values = c("blue", "red"),labels = c("Male", "Female"))+
  ggtitle("Risk of secondary transmission by sex") +
  theme(strip.background = element_rect(colour="black", fill="white")) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        legend.title=element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        legend.position = "bottom")

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\Transmission in Era of UTT\\Figures")
ggsave("Riskoftransmissionbyyear_withprev.jpg", height=12, width=12)


# # #


#bring in reportHIVbyageandgender for persontime at risk of transmission (infection years)
reporthivbyageandgender$src.gender <- reporthivbyageandgender$Gender
reporthivbyageandgender$year <- reporthivbyageandgender$Year2
reporthivbyageandgender$src.age.cat <- cut(reporthivbyageandgender$Age, c(15,20,25,30,35,40,45,200), right=F)
table(reporthivbyageandgender$Age, reporthivbyageandgender$src.age.cat)
head(reporthivbyageandgender)

reporthivbyageandgender.table.age.gender.year <- subset(reporthivbyageandgender)  %>%
  group_by(sim.id, src.age.cat, src.gender, year) %>%
  summarize(Infected=sum(Infected), Population=sum(Population), Prevalence=sum(Infected)/sum(Population)) %>%
  group_by(src.age.cat, src.gender, year) %>%
  summarize(medianInfected=median(Infected), medianPopulation=median(Population),medianPrevalence=median(Prevalence))  

summary(reporthivbyageandgender.table.age.gender.year$year)
tapply(reporthivbyageandgender.table.age.gender.year$medianPrevalence, reporthivbyageandgender.table.age.gender.year$src.age.cat, summary)
head(reporthivbyageandgender.table.age.gender.year)
head(transmission.table.age.gender.year)

transm.risk.table <- merge(transmission.table.age.gender.year[,c(1,2,4,5)], reporthivbyageandgender.table.age.gender.year, by=c("src.gender","src.age.cat","year"), all.y=T)
transm.risk.table$tr.risk <- transm.risk.table$medianTransmissions / transm.risk.table$medianInfected
head(transm.risk.table)
hist(transm.risk.table$tr.risk)
summary(transm.risk.table$year)

#plot
summary(transm.risk.table$year)
labs <- c("1" = "", "0" = "")
scl=1.6
ggplot(subset(transm.risk.table),mapping=aes(x=year,color=factor(src.gender))) + 
  geom_line(mapping=aes(y=tr.risk)) +
  geom_line(mapping=aes(y=medianPrevalence*scl), color="black", linetype=2) +
  facet_grid(src.gender~src.age.cat, labeller=labeller(src.gender = labs)) +
  ggtitle("Number of transmissions per infected individual per year") +
  scale_color_manual(values = c("blue", "red"),labels = c("Male", "Female"))+
  scale_x_continuous(expand = c(0,0), breaks=seq(1970,2050,10)) +
  scale_y_continuous(breaks=seq(0,1,0.1),expand = c(0,0),"transmission rate (per person-year)", sec.axis = sec_axis(~./scl , name="HIV prevalence", breaks=seq(0,0.6,0.05))) +
  coord_cartesian(ylim=c(0, 1)) +
  theme_minimal() +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
  theme(strip.background = element_rect(colour="white", fill="white")) +
  theme(axis.text.x = element_text(angle = 90)) +
  theme(strip.placement = "outside",panel.border = element_blank(), panel.grid.major = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.title=element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        legend.position = "bottom")

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\Transmission in Era of UTT\\Figures")
ggsave("Riskoftransmissionbyyear_byagegender.jpg", height=7, width=12)

#Overlay with proportion of infected individuals recenlty infected (incidence / prevalence ratio)

head(reporthivbyageandgender)
reporthivbyageandgender.table.age.gender.year <- subset(reporthivbyageandgender)  %>%
  group_by(sim.id, src.age.cat, src.gender, year) %>%
  summarize(NewlyInfected=sum(Newly.Infected), Infected=sum(Infected), Population=sum(Population), Prevalence=sum(Infected)/sum(Population), NewInfRatio=sum(Newly.Infected)/sum(Infected)) %>%
  group_by(src.age.cat, src.gender, year) %>%
  summarize(medianInfected=median(Infected), medianPopulation=median(Population),medianPrevalence=median(Prevalence),medianNewInfRatio=median(NewInfRatio))  

transm.risk.table <- merge(transmission.table.age.gender.year[,c(1,2,4,5)], reporthivbyageandgender.table.age.gender.year, by=c("src.gender","src.age.cat","year"), all.y=T)
transm.risk.table$tr.risk <- transm.risk.table$medianTransmissions / transm.risk.table$medianInfected
head(transm.risk.table)
hist(transm.risk.table$tr.risk)
summary(transm.risk.table$year)

ggplot(subset(transm.risk.table),mapping=aes(x=year,color=factor(src.gender))) + 
  geom_line(mapping=aes(y=tr.risk)) +
  geom_line(mapping=aes(y=medianNewInfRatio*scl), color="black", linetype=2) +
  facet_grid(src.gender~src.age.cat, labeller=labeller(src.gender = labs)) +
  ggtitle("Number of transmissions per infected individual per year") +
  scale_color_manual(values = c("blue", "red"),labels = c("Male", "Female"))+
  scale_x_continuous(expand = c(0,0), breaks=seq(1970,2040,10)) +
  scale_y_continuous(breaks=seq(0,1,0.1),expand = c(0,0),"transmission rate (per person-year)") +
  coord_cartesian(ylim=c(0, 1), xlim=c(1985, 2048)) +
  theme_minimal() +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
  theme(strip.background = element_rect(colour="white", fill="white")) +
  theme(axis.text.x = element_text(angle = 90)) +
  theme(strip.placement = "outside",panel.border = element_blank(), panel.grid.major = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.title=element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        legend.position = "bottom")

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\Transmission in Era of UTT\\Figures")
ggsave("Riskoftransmissionbyyear_byagegender_newinfprevratio.jpg", height=7, width=12)

# # #

# 2) who contributes most to transmissions? Age and gender contributions over time
head(network.data.m)
transmission.table <- subset(network.data.m, src.age.int>14 & src.age.int < 80 & transm.year.int>1989 & transm.year.int <2030)  %>%
  group_by(sim.id, era, src.gender, src.age.cat.all) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  group_by(era, src.gender, src.age.cat.all) %>%
  summarize(median = quantile(freq, probs = 0.5),
            lb = quantile(freq, probs = 0.025),
            ub = quantile(freq, probs = 0.975))
head(transmission.table)

labs <- c("1" = "Women", "0" = "Men")
ggplot(subset(transmission.table)) + 
  geom_point(aes(y=median, x=src.age.cat.all, color=factor(src.gender), group=src.gender), size=4) +
  geom_ribbon(aes(ymin=lb, ymax=ub, x=src.age.cat.all, fill=factor(src.gender), group=src.gender), alpha=0.3, size=2) +
  scale_x_discrete(labels=c("15-19","20-24","25-29","30-34","35-39","40-44","45-49","50-54","55-59","60-64","65-69","70-74","75-79")) +
  scale_y_continuous(expand = c(0,0), breaks=seq(0,0.3,0.05), limits=c(0,0.3)) +
  facet_grid(~era)+
  theme_bw(base_size=24) +
  xlab("Year")+ ylab("Proportion of transmissions") +
  scale_color_manual(values = c("blue", "red"),labels = c("Male", "Female"))+
  scale_fill_manual(values = c("blue", "red"),labels = c("Male", "Female"))+
  theme(strip.background = element_rect(colour="black", fill="white")) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        legend.title=element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        legend.position = "bottom")

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\Transmission in Era of UTT\\Figures")
ggsave("TransmAgeDistributionByEra.jpg", height=14, width=24)

#% of transmissions > / < 25 years
head(network.data.m)
transmission.table <- subset(network.data.m, src.age.int>14 & src.age.int < 80 & transm.year.int>1989 & transm.year.int <2030)  %>%
  group_by(sim.id, era, src.gender, src.age.25) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n)) %>%
  group_by(era, src.gender, src.age.25) %>%
  summarize(median = quantile(freq, probs = 0.5),
            lb = quantile(freq, probs = 0.025),
            ub = quantile(freq, probs = 0.975))
transmission.table

#####################################################################################################
#Analysis
#####################################################################################################
library(networkD3)

#pre-UTT era
#Plots are stratified on UTT era - chose one and then plot
table(network.data.m$dest.age.cat3)
transmission.table.simple <- subset(network.data.m, era=="Pre-UTT")  %>%
  group_by(sim.id, src.gender, dest.gender, src.age.cat3, dest.age.cat3) %>%
  summarize(transmissions=n()) %>%
  group_by(src.gender, dest.gender, src.age.cat3, dest.age.cat3) %>%
  summarize(meanTransmissions=mean(transmissions))  
head(transmission.table.simple)

##
transmission.table.simple$source.gender.age <- paste(transmission.table.simple$src.gender, transmission.table.simple$src.age.cat3)
transmission.table.simple$destination.gender.age <- paste(transmission.table.simple$dest.gender, transmission.table.simple$dest.age.cat3)

head(transmission.table.simple)
##
flow.data <- data.frame("orig_reg" = transmission.table.simple$source.gender.age, "dest_reg" = transmission.table.simple$destination.gender.age, "Transmissions" = transmission.table.simple$meanTransmissions)
table(flow.data$orig_reg)

flow.data$source <- as.numeric(flow.data$orig_reg)-1 #must be zero indexed so subtract 1
flow.data$target <- as.numeric(flow.data$dest_reg)-1 #must be zero indexed so subtract 1
head(flow.data)

plot.data = data.frame("source"=flow.data$source, "target"=flow.data$target, "value"=flow.data$Transmissions)
nodes <- data.frame("name"=c("Male 15-24","Male 25-34","Male 35+",
                             "Female 15-24","Female 25-34","Female 35+"))
plot.data.m <- plot.data[plot.data$source<3,] #plot for m -> f
plot.data.f <- plot.data[plot.data$source>2,] #plot for f -> m

sankeyNetwork(Links = plot.data.m, Nodes = nodes,
              Source = "source", Target = "target",
              Value = "value", NodeID = "name", nodeWidth = 30, fontSize = 24)

sankeyNetwork(Links = plot.data.f, Nodes = nodes,
              Source = "source", Target = "target",
              Value = "value", NodeID = "name", nodeWidth = 30, fontSize = 24)

#UTT era
table(network.data.m$dest.age.cat3)
transmission.table.simple <- subset(network.data.m, era=="UTT era")  %>%
  group_by(sim.id, src.gender, dest.gender, src.age.cat3, dest.age.cat3) %>%
  summarize(transmissions=n()) %>%
  group_by(src.gender, dest.gender, src.age.cat3, dest.age.cat3) %>%
  summarize(meanTransmissions=mean(transmissions))  
head(transmission.table.simple)

##
transmission.table.simple$source.gender.age <- paste(transmission.table.simple$src.gender, transmission.table.simple$src.age.cat3)
transmission.table.simple$destination.gender.age <- paste(transmission.table.simple$dest.gender, transmission.table.simple$dest.age.cat3)

head(transmission.table.simple)
##
flow.data <- data.frame("orig_reg" = transmission.table.simple$source.gender.age, "dest_reg" = transmission.table.simple$destination.gender.age, "Transmissions" = transmission.table.simple$meanTransmissions)
table(flow.data$orig_reg)

flow.data$source <- as.numeric(flow.data$orig_reg)-1 #must be zero indexed so subtract 1
flow.data$target <- as.numeric(flow.data$dest_reg)-1 #must be zero indexed so subtract 1
head(flow.data)

plot.data = data.frame("source"=flow.data$source, "target"=flow.data$target, "value"=flow.data$Transmissions)
nodes <- data.frame("name"=c("Male 15-24","Male 25-34","Male 35+",
                             "Female 15-24","Female 25-34","Female 35+"))
plot.data.m <- plot.data[plot.data$source<3,] #plot for m -> f
plot.data.f <- plot.data[plot.data$source>2,] #plot for f -> m

sankeyNetwork(Links = plot.data.m, Nodes = nodes,
              Source = "source", Target = "target",
              Value = "value", NodeID = "name", nodeWidth = 30, fontSize = 24)

sankeyNetwork(Links = plot.data.f, Nodes = nodes,
              Source = "source", Target = "target",
              Value = "value", NodeID = "name", nodeWidth = 30, fontSize = 24)

#################################################################################################
#Risk of HIV transmission - number of transmissions per 100 py infected
#################################################################################################

#risk of transmission by stage of disease and time since infection
table(network.data.m$dest.age.cat3)
transmission.table.age.gender.year <- subset(network.data.m)  %>%
  group_by(sim.id, dest.gender, dest.age.cat3, transm.year.int, src.stage.f) %>%
  summarize(transmissions=n()) %>%
  group_by(dest.gender, dest.age.cat3, transm.year.int, src.stage.f) %>%
  summarize(medianTransmissions=median(transmissions))  
transmission.table.age.gender.year$year <- transmission.table.age.gender.year$transm.year.int
transmission.table.age.gender.year$Gender <- transmission.table.age.gender.year$dest.gender
head(transmission.table.age.gender.year)

head(reporthivbyageandgender)
reporthivbyageandgender$dest.age.cat3 <- cut(reporthivbyageandgender$Age, c(15,25,35,200), right=F)
reporthivbyageandgender.table.incidence <- subset(reporthivbyageandgender, scenario=="1.baseline")  %>%
  group_by(sim.id, dest.age.cat3, Gender, year) %>%
  summarize(Newly.Infected=sum(Newly.Infected), Uninfected.Population=sum(Uninfected.Population), Incidence=sum(Newly.Infected)/sum(Uninfected.Population)) %>%
  group_by(dest.age.cat3, Gender, year) %>%
  summarize(medianNewly.Infected=median(Newly.Infected), medianUninfected.Population=median(Uninfected.Population),medianIncidence=median(Incidence))  

names(transmission.table.age.gender.year)
names(reporthivbyageandgender.table.incidence)
transm.risk.table <- merge(transmission.table.age.gender.year[,c(2,4,5,6,7)], reporthivbyageandgender.table.incidence, by=c("Gender","dest.age.cat3","year"), all.x=T)
head(transm.risk.table)

#include across all ages
names(network.data.m)
transmission.table.gender.year <- subset(network.data.m)  %>%
  group_by(sim.id, dest.gender, transm.year.int, src.stage.f) %>%
  summarize(transmissions=n()) %>%
  group_by(dest.gender, transm.year.int, src.stage.f) %>%
  summarize(medianTransmissions=median(transmissions))  
transmission.table.gender.year$year <- transmission.table.gender.year$transm.year.int
transmission.table.gender.year$Gender <- transmission.table.gender.year$dest.gender
head(transmission.table.gender.year)

transmission.table.sum <- subset(network.data.m)  %>% group_by(sim.id, dest.gender, transm.year.int) %>%
  summarize(transmissions=n()) %>%  group_by(dest.gender, transm.year.int) %>%
  summarize(medianTransmissions=median(transmissions))  
transmission.table.sum$year <- transmission.table.sum$transm.year.int
transmission.table.sum$Gender <- transmission.table.sum$dest.gender

head(reporthivbyageandgender)
reporthivbyageandgender.table.incidence <- subset(reporthivbyageandgender, scenario=="1.baseline")  %>%
  group_by(sim.id, Gender, year) %>%
  summarize(Uninfected.Population=sum(Uninfected.Population),NewlyInfected.Population=sum(Newly.Infected)) %>%
  group_by(Gender, year) %>%
  summarize(medianUninfected.Population=median(Uninfected.Population),medianNewlyInfected.Population=median(NewlyInfected.Population))  
head(reporthivbyageandgender.table.incidence)
head(transmission.table.gender.year)

transm.risk.table.all <- merge(transmission.table.gender.year[,c(3:6)], reporthivbyageandgender.table.incidence, by=c("Gender","year"), all.x=T)
transm.risk.table.sum <- merge(transmission.table.sum[,c(3:5)], reporthivbyageandgender.table.incidence, by=c("Gender","year"), all.x=T)

head(transm.risk.table.all)

ggplot(subset(transm.risk.table.all),mapping=aes(x=year,color=src.stage.f)) + 
  geom_line(mapping=aes(y=medianTransmissions/medianUninfected.Population), size=2) +
  geom_line(data=transm.risk.table.sum, mapping=aes(y=medianTransmissions/medianUninfected.Population), size=2, color="black") +
  facet_grid(~Gender, labeller=labeller(Gender = labs)) +
  ggtitle("Incidence of HIV by source's stage of infection") +
  #scale_color_manual(values = c("blue", "red"),labels = c("Male", "Female"))+
  scale_x_continuous(expand = c(0,0), breaks=seq(1970,2040,10)) +
  scale_y_continuous(breaks=seq(0,1,0.01),expand = c(0,0),"transmission rate (per person-years)") +
  #coord_cartesian(ylim=c(0, 0.06), xlim=c(1980, 2050)) +
  theme_minimal() +
  #annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
  #annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
  theme(strip.background = element_rect(colour="white", fill="white")) +
  theme(axis.text.x = element_text(angle = 90)) +
  theme(strip.placement = "outside",panel.border = element_blank(), panel.grid.major = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.title=element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        legend.position = "bottom")

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\Transmission in Era of UTT\\Figures")
ggsave("HIVincidence_by_src.stage.jpg", height=12, width=12)

#################################################################################################
#Time to transmission by sex/age/year 
#################################################################################################
head(network.data.m)
network.data.m$TimeToTransmission <- network.data.m$transm.year - network.data.m$acq.year

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























