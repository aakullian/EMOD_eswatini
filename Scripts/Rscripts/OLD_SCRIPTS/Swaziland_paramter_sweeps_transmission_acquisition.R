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
#Bring in reportHIVbyageandgender Baseline Eswatini model with ART stages
##################################################################################################

#read in reportHIVbyageandgender files
input_dir <- "C:\\Users\\aakullian\\Dropbox (IDM)\\GitHub\\EMOD_eswatini\\Baseline_transmission_ARTstate"
primary_dirs = list.files(input_dir)
transmission.data.frame <- data.frame("Year2"=seq(1990,2050,1))

i=0
for(dir in primary_dirs) {
  i=i+1
  temp_table <- as.data.table(read.csv(file = paste(input_dir,dir,"ReportHIVByAgeAndGender.csv",sep="/")))
  temp_table$Year2 <- floor((temp_table$Year-0.5))
  temp_table$Uninfected.Population = temp_table$Population-temp_table$Infected
  trajectories_IR.1a <- aggregate(cbind(Newly.Infected,Transmitters) ~ Year2+Gender+Age+IsCircumcised+IP_Key.ARTstate, subset(temp_table, Age>10 & Year > 1990 & Year < 2051), FUN=sum) #sums number of new infections/transmissions in each year
  trajectories_IR.2 <- aggregate(cbind(Uninfected.Population,Infected,Population) ~ Year+Gender+Age+IsCircumcised+IP_Key.ARTstate, subset(temp_table, Age>10 & Year > 1990 & Year < 2051), FUN=sum) #sums at risk population in each year
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
reporthivbyageandgender.age.gender.year.mean <- aggregate(riskoftransm ~ agecat3+Gender+yearcat3, transmission.data.frame.age.gender.year, FUN=mean)

#inc_values6.agg.m$IP_Key.ARTstate_f <- ordered(inc_values6.agg.m$IP_Key.ARTstate, levels=c("NeverOnART","OnART","OffART"))
labs <- c("1" = "Women", "0" = "Men")
labs.age <- c("[15,25)"="15-24","[25,35)"="25-34","[35,101)"="35+")
ggplot(data=subset(reporthivbyageandgender)) +
  #geom_line(aes(x=Year2, y=riskoftransm,group=interaction(Gender,agecat3,sim.id), color=agecat3, alpha=0.04),size=1,alpha=0.02) +
  geom_smooth(aes(x=Year2, y=riskoftransm, color=agecat3, group=interaction(agecat3, Gender), linetype = factor(Gender)),method="loess", span=0.2, se=F,size=2) +
  xlab("Year")+
  ylab("Proportion")+
  facet_grid(~agecat3,labeller=labeller(agecat3 = labs.age)) +
  theme_bw(base_size=16) +
  scale_x_continuous(breaks = seq(1990,2035,10),limits=c(1990,2035), expand = c(0,0)) +
  scale_y_continuous(breaks = seq(0,0.5,0.1),limits=c(0,0.5), expand = c(0,0)) +
  scale_color_manual(labels = c("15-24","25-34", "35+"), values = c("red", "purple", "blue")) +
  ggtitle("Proportion of PLHIV who transmit per year by age") +
  theme(legend.position="bottom") +
  #guides(color=guide_legend(title="ART status", nrow=1))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(strip.background = element_rect(colour="black", fill="white")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"))

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\Transmission in Era of UTT\\Figures")
ggsave("ProportionWhoTransmit_overtime.jpg", height=8, width=12)

table(transmission.data.frame.age.gender.year$yearcat3)
ggplot(data=subset(transmission.data.frame.age.gender.year)) +
  geom_point(aes(x=yearcat3, y=riskoftransm, group=interaction(sim.id), color=agecat3),size=1, alpha=0.5,position = position_dodge(width = 0.05)) +
  geom_point(data=transmission.data.frame.age.gender.year.mean,aes(x=yearcat3, y=riskoftransm, color=agecat3),shape=95, size=10) +
  #geom_smooth(aes(x=Year2, y=riskoftransm, color=agecat3, group=interaction(agecat3, Gender), linetype = factor(Gender)),method="loess", span=0.2, se=F,size=2) +
  xlab("Year")+
  ylab("Proportion")+
  facet_grid(~Gender,labeller=labeller(Gender = labs)) +
  theme_bw(base_size=16) +
  scale_x_discrete(labels=c("[1.99e+03,2e+03)"="1990-2003","[2e+03,2.02e+03)"="2004-2015","[2.02e+03,2.03e+03)"="2016-2029","[2.03e+03,2.05e+03)"="2030-2050")) +
  scale_y_continuous(breaks = seq(0,0.35,0.05),limits=c(0,0.35), expand = c(0,0)) +
  scale_color_manual(labels = c("15-24","25-34", "35+"), values = c("red", "purple", "blue")) +
  ggtitle("Proportion of PLHIV who transmit per year") +
  theme(legend.position=c(0.86,0.97)) +
  guides(color=guide_legend(title="", nrow=1))+
  #theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(strip.background = element_rect(colour="black", fill="white")) +
  theme(legend.background = element_rect(fill="transparent"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"))

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\Transmission in Era of UTT\\Figures")
ggsave("ProportionWhoTransmit_bar.jpg", height=8, width=12)

#################################################################################################
#Bring in Transmission data
#################################################################################################

#bring in full saved dataset
setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\HIV\\EMOD\\Swaziland_05June2018\\Scenarios_calibration_run3_v05Sep2018_newcampaign_FINAL\\ScenariosOct_FINAL\\Baseline_TransmissionReport\\TransmissionReport\\")
transmissionreport.master <- readRDS("transmissionreport.master.rds")

#create variables
head(transmissionreport.master)
table(transmissionreport.master$sim.id)
transmissionreport <- transmissionreport.master[order(transmissionreport.master$sim.id, transmissionreport.master$SRC_ID, transmissionreport.master$YEAR),]
transmissionreport <- subset(transmissionreport, YEAR>=1990 & YEAR <2050)
summary(transmissionreport$YEAR)
transmissionreport$src_age_int <- floor(transmissionreport$SRC_AGE / 365.25)
transmissionreport$dest_age_int <- floor(transmissionreport$DEST_AGE / 365.25)
transmissionreport <- subset(transmissionreport, dest_age_int > 14 & src_age_int > 14)
transmissionreport$age.at.acquisition.cat <- cut(transmissionreport$dest_age_int, c(15,25,35,200), right=F)
transmissionreport$age.at.transmission.cat <- cut(transmissionreport$src_age_int, c(15,25,35,200), right=F)
transmissionreport$age.at.acquisition.cat5 <- cut(transmissionreport$dest_age_int, c(15,20,25,30,35,40,45,200), right=F)
transmissionreport$age.at.transmission.cat5 <- cut(transmissionreport$src_age_int, c(15,20,25,30,35,40,45,200), right=F)
transmissionreport$yearcat3 <- cut(transmissionreport$YEAR, c(1990,2004,2016,2030,2050), right=F, dig.lab=4)

head(transmissionreport)
table(transmissionreport$SRC_STAGE)

#################################################################################################
# Transmisison network - to plot
#################################################################################################
library(igraph)

#Create transmission network for plotting (choose subset)
head(transmissionreport)
network.data <- subset(transmissionreport, sim.id=="TransmissionReport_TPI0000_REP0001.csv") #first sim
network.data <- subset(transmissionreport) #all sims
date.of.acquisition <- data.frame("SRC_ID"=network.data$DEST_ID, "acq.year"=network.data$YEAR, "sim.id"=network.data$sim.id) #source id used to merge in with source of transmitters
summary(date.of.acquisition$acq.year)
sum(is.na(date.of.acquisition$acq.year))

network.data <- merge(network.data, date.of.acquisition, by="SRC_ID")
network.data.m <- data.frame("src"=network.data$SRC_ID, "dest"=network.data$DEST_ID, 
                           "src.gender"=network.data$SRC_GENDER, 
                           "dest.gender"=network.data$DEST_GENDER, 
                           "dest.age"=network.data$age.at.acquisition.cat5, 
                           "src.age"=network.data$age.at.transmission.cat5,
                           "dest.age.int"=network.data$dest_age_int, 
                           "src.age.int"=network.data$src_age_int,
                           "transm.year"=network.data$YEAR,
                           "transm.year.int"=floor(network.data$YEAR),
                           "acq.year"=network.data$acq.year,
                           "acq.age"=floor(network.data$src_age_int-(network.data$YEAR-network.data$acq.year)),
                           "transm.year.cat"=network.data$yearcat3)

network.data.m$acq.year.cat <- cut(network.data.m$acq.year, c(1990,2004,2016,2030,2050), right=F, dig.lab=4)
network.data.m$acq.age.cat <- cut(network.data.m$acq.age, c(15,20,25,30,35,40,45,200), right=F)
network.data.m$src.age.cat3 <- cut(network.data.m$src.age.int, c(15,25,35,200), right=F)
network.data.m$dest.age.cat3 <- cut(network.data.m$dest.age.int, c(15,25,35,200), right=F)
network.data.m$era <- "UTT era"
network.data.m$era[network.data.m$transm.year < 2016] <- "Pre-UTT"
network.data.m$era.acq <- "UTT era"
network.data.m$era.acq[network.data.m$acq.year < 2016] <- "Pre-UTT"

head(network.data.m,20)

#################################################################################################
# Transmisison flows
#################################################################################################
# Libraries
library(networkD3)

#Plots are stratified on UTT era - chose one and then plot
table(network.data.m$dest.age.cat3)
transmission.table.simple <- subset(network.data.m, era=="Pre-UTT")  %>%
  group_by(src.gender, dest.gender, src.age.cat3, dest.age.cat3) %>%
  summarize(transmissions=n())

table(network.data.m$dest.age.cat3)
transmission.table.simple <- subset(network.data.m, era=="UTT era")  %>%
  group_by(src.gender, dest.gender, src.age.cat3, dest.age.cat3) %>%
  summarize(transmissions=n())

##
transmission.table.simple$source.gender.age <- paste(transmission.table.simple$src.gender, transmission.table.simple$src.age.cat3)
transmission.table.simple$destination.gender.age <- paste(transmission.table.simple$dest.gender, transmission.table.simple$dest.age.cat3)

head(transmission.table.simple)
##
flow.data <- data.frame("orig_reg" = transmission.table.simple$source.gender.age, "dest_reg" = transmission.table.simple$destination.gender.age, "Transmissions" = transmission.table.simple$transmissions)
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
head(network.data.m)
transmission.table.age.gender.year <- subset(network.data.m)  %>%
  group_by(src.gender, src.age, transm.year.int) %>%
  summarize(transmissions=n())

#bring in reportHIVbyageandgender for persontime at risk of transmission (infection years)
reporthivbyageandgender$src.age <- cut(reporthivbyageandgender$Age, c(15,20,25,30,35,40,45,200), right=F)
table(reporthivbyageandgender$Age, reporthivbyageandgender$src.age)
table(transmission.table.age.gender.year$src.age)
table(reporthivbyageandgender$Age)
head(reporthivbyageandgender)
head(transmission.table.age.gender.year)

#################################################################################################
#Time to transmission by sex/age/year 
#################################################################################################
network.data.m$TimeToTransmission <- network.data.m$transm.year - network.data.m$acq.year
summary(network.data.m$TimeToTransmission) #median 5.2 years to transmission (1.2-11.6 IQR)
summary(network.data.m$TimeToTransmission[network.data.m$src.gender==0]) #Men: median 4.7 years to transmission (1.0-10.3 IQR)
summary(network.data.m$TimeToTransmission[network.data.m$src.gender==1]) #Women: median 6 years to transmission (1.8-13.6 IQR)

#Proportion of HIV transmissions occuring soon after acqusision 1,2,3...50 years?
network.data.m$TimeToTransmissionCat <- cut(network.data.m$TimeToTransmission, c(0,1,10,50,100), right=F, dig.lab=4)
table(network.data.m$TimeToTransmissionCat)
table(network.data.m$era)

TimeToTransmissionProportions.preUTT <- subset(network.data.m, era=="Pre-UTT") %>%
  group_by(src.age.int, src.gender, TimeToTransmissionCat) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))
TimeToTransmissionProportions.preUTT$era <- "pre-UTT"

TimeToTransmissionProportions.UTT <- subset(network.data.m, era=="UTT era") %>%
  group_by(src.age.int, src.gender, TimeToTransmissionCat) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))
TimeToTransmissionProportions.UTT$era <- "UTT era"

TimeToTransmissionProportions <- rbind(TimeToTransmissionProportions.preUTT,TimeToTransmissionProportions.UTT)

labs <- c("1" = "Women", "0" = "Men")
ggplot(subset(TimeToTransmissionProportions, TimeToTransmissionCat!="[50,100)"), aes(fill=forcats::fct_rev(TimeToTransmissionCat), y=freq, x=src.age.int)) + 
  geom_bar(position="stack", stat="identity", width=1) +
  facet_grid(era~src.gender, labeller=labeller(src.gender = labs))+
  scale_x_continuous(expand = c(0,0), breaks=seq(15,49,5)) +
  scale_y_continuous(expand = c(0,0), breaks=seq(0,1,0.1)) +
  coord_cartesian(ylim=c(0, 1),xlim=c(15, 50)) +
  geom_hline(yintercept=0.5, color='red', linetype=2) +
  theme_bw(base_size=16) +
  xlab("Age at HIV transmission")+ ylab("Probability") +
  scale_fill_manual(values = c("white", "blue", "darkblue"),labels = c("10+", "1-9", "<1"),name= "Years", guide = guide_legend(reverse = TRUE))+
  ggtitle("Time from acquisition to transmission") +
  guides(fill = guide_legend(reverse=T)) +
  theme(strip.background = element_rect(colour="black", fill="white")) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        legend.position = "bottom")

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\Transmission in Era of UTT\\Figures")
ggsave("TimeToTransmission_Cat_TransmissionPerspective.jpg", height=12, width=12)

#Proportion of HIV transmissions occuring soon after acqusision 1,2,3...50 years (acquision perspective)?
network.data.m$TimeToTransmissionCat <- cut(network.data.m$TimeToTransmission, c(0,1,5,10,50,100), right=F, dig.lab=4)
table(network.data.m$TimeToTransmissionCat)
table(network.data.m$era)

TimeToTransmissionProportions.preUTT <- subset(network.data.m, era=="Pre-UTT") %>%
  group_by(dest.age.int, dest.gender, TimeToTransmissionCat) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))
TimeToTransmissionProportions.preUTT$era <- "pre-UTT"

TimeToTransmissionProportions.UTT <- subset(network.data.m, era=="UTT era") %>%
  group_by(dest.age.int, dest.gender, TimeToTransmissionCat) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))
TimeToTransmissionProportions.UTT$era <- "UTT era"

TimeToTransmissionProportions <- rbind(TimeToTransmissionProportions.preUTT,TimeToTransmissionProportions.UTT)

labs <- c("1" = "Women", "0" = "Men")
ggplot(subset(TimeToTransmissionProportions, era=="UTT era" & TimeToTransmissionCat!="[50,100)"), aes(fill=forcats::fct_rev(TimeToTransmissionCat), y=freq, x=dest.age.int)) + 
  geom_area(position="stack") +
  facet_grid(~dest.gender, labeller=labeller(dest.gender = labs))+
  scale_x_continuous(expand = c(0,0), breaks=seq(15,49,5)) +
  scale_y_continuous(expand = c(0,0), breaks=seq(0,1,0.1)) +
  coord_cartesian(ylim=c(0, 1),xlim=c(15, 50)) +
  #geom_hline(yintercept=0.5, color='red', linetype=2) +
  theme_bw(base_size=16) +
  xlab("Age at HIV acquisition")+ ylab("Proportion") +
  scale_fill_manual(values = c("lightblue", "blue3", "darkblue","black"),
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


#What is the median time from acquisition to transmission by age,gender, at year of transmission).
TimeToTransmission <- subset(network.data.m)  %>%
  group_by(transm.year.cat, src.gender, src.age) %>%
  summarize(time.to.transmission = median(TimeToTransmission, na.rm = TRUE))

ggplot(data=subset(network.data.m), aes(x=src.age, y=TimeToTransmission, fill=factor(src.gender))) +
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

#Look at the distribution of time from HIV acquisition to tranmsission across all transmissions
labs <- c("1" = "Women", "0" = "Men")
head(network.data.m)
ggplot(subset(network.data.m, acq.year<2016), aes(x=acq.age, y=TimeToTransmission, color=factor(src.gender), alpha=0.01)) +
  geom_point(size=3) +
  facet_grid(~src.gender, labeller=labeller(src.gender = labs)) +
  scale_color_manual(values = c("red2","blue3")) +
  ggtitle("Time from HIV acquisition to transmission") +
  xlab("Age at acquisition")+ ylab("Years")+
  scale_x_continuous(expand = c(0,0), breaks=seq(15,100,5)) +
  scale_y_continuous(expand = c(0,0), breaks=seq(0,60,5)) +
  coord_cartesian(xlim=c(14.5, 60)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        legend.position = "none")

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\Transmission in Era of UTT\\Figures")
ggsave("TimeToTransmission_AcquisitionPerspective_Points.jpg", height=8, width=12)

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

head(network.data.m)
ggplot(subset(network.data.m, acq.age>14), aes(x=TimeToTransmission, fill=factor(src.gender))) +
  geom_density(alpha=0.3, adjust=3) +
  facet_grid(era~acq.age.cat5,labeller=labeller(acq.age.cat5 = labs.age)) +
  scale_color_manual(values = c("blue3","red2"), labels=c("Men","Women")) +
  scale_fill_manual(values = c("blue3","red2"), labels=c("Men","Women")) +
  ggtitle("Time from HIV acquisition to transmission") +
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

head(network.data.m)
ggplot(network.data.m, aes(x=dest.age.int, y=TimeToTransmission)) +
  scale_fill_gradientn(name="Incidence rate (per 100 py)",
                       colours=c("darkblue","blue","yellow","darkred"),
                       limits=limits,
                       breaks=breaks) +
  facet_grid(era~dest.gender) +
  stat_contour(data=subset(pred.inc.age.year), aes(y=Age, x=Year, z=inc*100),
               color="black", size=0.3, linetype=1, binwidth=0.5) +
  theme(legend.position="bottom") +
  xlab("age")+
  ylab("time to transmissions")+
  scale_x_continuous(expand=c(0,0), breaks=seq(15,49,5)) +
  scale_y_continuous(expand=c(0,0), breaks=c(1,10,50)) +
  coord_cartesian(xlim=c(15, 50), ylim=c(0,50)) +
  theme_bw(base_size=16) +
  guides(fill = guide_legend(keywidth = 2, keyheight = 1, nrow=1)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(legend.position="bottom") +
  theme(strip.background = element_rect(colour="black", fill="white"))

#Proportion of transmissions by time-cat

#Among those acquiring HIV (by age and gender) - how long had the transmitter been infected before infecting?
network.data.m %>%
  group_by(dest.gender, dest.age) %>%
  summarize(time.to.transmission = median(TimeToTransmission, na.rm = TRUE))

head(network.data.m)
ggplot(data=subset(network.data.m), aes(x=dest.age, y=TimeToTransmission, fill=factor(dest.gender))) +
  #geom_point(color=dest.age),size=1, alpha=0.5,position = position_dodge(width = 0.05)) +
  geom_boxplot(notch=FALSE, outlier.shape=NA,outlier.size = 0, coef=0, width=0.5)+
  xlab("Age of acquiring person")+ ylab("Years")+
  facet_wrap(~era) +
  theme_bw(base_size=16) +
  scale_fill_manual(values = c("white","grey"), labels=c("Men","Women")) +
  scale_y_continuous(expand = c(0,0), breaks=seq(0,100,1)) +
  coord_cartesian(ylim=c(0, 35)) +
  ggtitle("Time from acquisition to transmission (from acquiring patner perspective)") +
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























