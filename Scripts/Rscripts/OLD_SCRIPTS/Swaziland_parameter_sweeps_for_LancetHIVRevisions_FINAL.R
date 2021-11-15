##################################################################################################
#Swaziland EMOD plotting
##################################################################################################
library(reshape2)
library(plyr)
library(dplyr)
library(ggplot2)
library(GGally)
library(dplyr)
library(mgcv)
library(data.table)
library(tidyr)
library(matrixStats)
library(stringr)

options(scipen=999)

# baseline sweeps

##################################################################################################
#Acute Paramter Sweep (includes a 2-dim sweep over acute duration and acute stage multiplier)
##################################################################################################
input_dir <- "C:/Users/aakullian/Documents/GitHub/EMOD_eswatini/ParamterSweepOutput_Baseline/acuteness"
primary_dirs = list.files(input_dir)
dir = "Acute_Stage_Infectivity_Multiplier-10--Acute_Duration_In_Months-1"
sub = "045538e9-acf8-e911-a2c3-c4346bcb1551"

inc_values <- data.frame("year"=seq(1980,2056,1))
median_inc <- data.frame("year"=seq(1980,2056,1))
i = 0

for(dir in primary_dirs) {
  sub_folders = list.files(paste(input_dir,dir,sep = "/"))
  inc_values <- data.frame("year"=seq(1980,2056,1))
  for(sub in sub_folders){
    i=i+1
    temp_table <- as.data.table(read.csv(file = paste(input_dir,dir,sub,"ReportHIVByAgeAndGender.csv",sep="/")))
    temp_table$Year2 <- floor((temp_table$Year-0.5))
    temp_table$Uninfected.Population = temp_table$Population-temp_table$Infected
    trajectories_IR.1a <- aggregate(Newly.Infected ~ Year2, subset(temp_table, Age>10 & Age<50), FUN=sum) #sums number of new infections in each year
    trajectories_IR.2 <- aggregate(Uninfected.Population ~ Year, subset(temp_table, Age>10 & Age<50), FUN=sum)
    trajectories_IR.2$Year2 <- floor(trajectories_IR.2$Year-0.5)
    trajectories_IR.2 <- trajectories_IR.2[!duplicated(trajectories_IR.2[c("Year2")]),] #remove second instance of duplicate rows
    trajectories_IR.2 <- trajectories_IR.2[-match("Year",names(trajectories_IR.2))]
    trajectories_IRoverall <- merge(trajectories_IR.1a, trajectories_IR.2, by=c("Year2"))
    trajectories_IRoverall$incidence <- trajectories_IRoverall$Newly.Infected / trajectories_IRoverall$Uninfected.Population
    inc_values$newCol1 <- trajectories_IRoverall$incidence
    colnames(inc_values)[ncol(inc_values)] <- paste(sub)
    print(paste("working on sub-folder",sub,"sim",i,sep=" "))
  }
  median_inc$val1=rowQuantiles(as.matrix(inc_values[,2:ncol(inc_values)]), probs = 0.5)
  colnames(median_inc)[ncol(median_inc)] <- paste(str_extract_all(dir,"\\(?[0-9,.]+\\)?")[[1]][1],str_extract_all(dir,"\\(?[0-9,.]+\\)?")[[1]][2], sep=",") 
  print(paste("working on sub-folder",dir,sep=" "))
}

head(median_inc)
rtest <- gather(median_inc, parameter, incidence, `10,1`:`30,5`, factor_key=TRUE)
acute_stage <- data.frame(within(rtest, parameter<-data.frame(do.call('rbind', strsplit(as.character(parameter), ',', fixed=TRUE)))))
acute_stage <- data.frame("year"=acute_stage$year,"incidence"=acute_stage$incidence,"acute_multi"=acute_stage$parameter[1], "acute_dur"=acute_stage$parameter[2])
acute_stage$acute_multi <- acute_stage$X1
acute_stage$acute_dur <- acute_stage$X2
acute_stage <- acute_stage[,c(1:2,5:6)]
head(acute_stage)

acute_stage$acute_multi <- as.numeric(as.character(acute_stage$acute_multi))
summary(acute_stage$acute_multi)
acute_stage$acute_dur <- as.numeric(as.character(acute_stage$acute_dur))
summary(acute_stage$acute_dur)

head(acute_stage)

ggplot() +
  geom_line(data=subset(acute_stage, acute_dur==3), aes(x=year, y=incidence*100, color=as.factor(acute_multi), group=acute_multi),size=2) +
  xlab("Year")+
  ylab("Incidence (per 100 py)")+
  theme_bw(base_size=16) +
  scale_x_continuous(breaks = seq(1980,2050,10),limits=c(1980,2051), expand = c(0,0)) +
  ggtitle("Sensitivity analysis of acute stage infectivity multiplier") +
  theme(legend.position="bottom") +
  guides(color=guide_legend(title="Multiplier on acute stage infectivity"))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(strip.background = element_rect(colour="black", fill="white")) +
  #guides(col = guide_legend(nrow=3)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"))

setwd("C:\\Users\\aakullian\\Documents\\GitHub\\EMOD_eswatini\\LancetHIVParameterSweeps\\Figures")
ggsave("inc_trends_acutemult_acutedur3mo.jpg", height=8, width=8)

ggplot() +
  geom_line(data=subset(acute_stage, acute_multi==25), aes(x=year, y=incidence*100, color=as.factor(acute_dur), group=acute_dur),size=2) +
  xlab("Year")+
  ylab("Incidence (per 100 py)")+
  theme_bw(base_size=16) +
  scale_x_continuous(breaks = seq(1980,2050,10),limits=c(1980,2051), expand = c(0,0)) +
  ggtitle("Sensitivity analysis of acute stage duration") +
  theme(legend.position="bottom") +
  guides(color=guide_legend(title="Acute stage duration in months"))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(strip.background = element_rect(colour="black", fill="white")) +
  #guides(col = guide_legend(nrow=3)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"))

setwd("C:\\Users\\aakullian\\Documents\\GitHub\\EMOD_eswatini\\LancetHIVParameterSweeps\\Figures")
ggsave("inc_trends_acutedur_acutemulti26.jpg", height=8, width=8)

ggplot(data=subset(acute_stage, year==2016 | year==2020 | year == 2030 | year==2050)) +
  geom_raster(aes(y=as.numeric(as.character(acute_multi)), x=as.numeric(as.character(acute_dur)), fill=incidence*100))+ 
  geom_point(aes(x=3, y=25), shape=3) +
  scale_fill_gradient(name="incidence per 100 py", low="blue", high="yellow",breaks=c(0.6,0.8,1,1.2,1.4,1.6)-0.1)+
  stat_contour(aes(y=acute_multi, x=acute_dur, z=incidence*100),
               color="black", size=0.1, linetype=1, binwidth=0.1) +
  theme(legend.position="bottom") +
  facet_wrap(~year)+
  xlab("acute phase duration (months)") +
  ylab("acute phase infectivity")+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(keywidth = 2, keyheight = 1, nrow=1)) +
  theme(legend.position="bottom") +
  theme(strip.background = element_rect(colour="black", fill="white")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

setwd("C:\\Users\\aakullian\\Documents\\GitHub\\EMOD_eswatini\\LancetHIVParameterSweeps\\Figures")
ggsave("inc_raster_acute_dur_multi.jpg", height=8, width=8)

##################################################################################################
#Paramter sweep of ART reduce acquire & delay from infection to ART initiation 
##################################################################################################

# Set working directories
input_dir <- "C:/Users/aakullian/Documents/GitHub/EMOD_eswatini/ParamterSweepOutput_Baseline/delay_and_supression"
primary_dirs = list.files(input_dir)
inc_values2 <- data.frame("year"=seq(1980,2056,1))
median_inc2 <- data.frame("year"=seq(1980,2056,1))

dir="Delay_Period_Mean-60--ART_Viral_Suppression_Multiplier-0.0"

for(dir in primary_dirs) {
  i=0
  sub_folders = list.files(paste(input_dir,dir,sep = "/"))
  inc_values2 <- data.frame("year"=seq(1980,2056,1))
  for(sub in sub_folders){
    i=i+1
    temp_table <- as.data.table(read.csv(file = paste(input_dir,dir,sub,"ReportHIVByAgeAndGender.csv",sep="/")))
    temp_table$Year2 <- floor((temp_table$Year-0.5))
    temp_table$Uninfected.Population = temp_table$Population-temp_table$Infected
    trajectories_IR.1a <- aggregate(Newly.Infected ~ Year2, subset(temp_table, Age>10 & Age<50), FUN=sum) #sums number of new infections in each year
    trajectories_IR.2 <- aggregate(Uninfected.Population ~ Year, subset(temp_table, Age>10 & Age<50), FUN=sum)
    trajectories_IR.2$Year2 <- floor(trajectories_IR.2$Year-0.5)
    trajectories_IR.2 <- trajectories_IR.2[!duplicated(trajectories_IR.2[c("Year2")]),] #remove second instance of duplicate rows
    trajectories_IR.2 <- trajectories_IR.2[-match("Year",names(trajectories_IR.2))]
    trajectories_IRoverall <- merge(trajectories_IR.1a, trajectories_IR.2, by=c("Year2"))
    trajectories_IRoverall$incidence <- trajectories_IRoverall$Newly.Infected / trajectories_IRoverall$Uninfected.Population
    inc_values2$newCol1 <- trajectories_IRoverall$incidence
    colnames(inc_values2)[ncol(inc_values2)] <- paste(sub)
    print(paste("working on folder", dir, "sub-folder",sub,"sim",i,sep=" "))
  }
  median_inc2$val1=rowQuantiles(as.matrix(inc_values2[,2:ncol(inc_values2)]), probs = 0.5)
  colnames(median_inc2)[ncol(median_inc2)] <- paste(str_extract_all(dir,"\\(?[0-9,.]+\\)?")[[1]][1],str_extract_all(dir,"\\(?[0-9,.]+\\)?")[[1]][2], sep=",") 
  print(paste("working on sub-folder",dir,sep=" "))
}

rtest <- gather(median_inc2, parameter, incidence, `0,0.0`:`60,0.2`, factor_key=TRUE)
head(rtest,100)
ARTeffect <- data.frame(within(rtest, parameter<-data.frame(do.call('rbind', strsplit(as.character(parameter), ',', fixed=TRUE)))))
ARTeffect <- data.frame("year"=ARTeffect$year,"incidence"=ARTeffect$incidence,"timetoart"=ARTeffect$parameter[1], "artefficacy"=ARTeffect$parameter[2])
ARTeffect$timetoart <- ARTeffect$X1
ARTeffect$artefficacy <- ARTeffect$X2
ARTeffect <- ARTeffect[,c(1:2,5:6)]
class(ARTeffect$timetoart)
class(ARTeffect$artefficacy)
head(ARTeffect)

ARTeffect$timetoart <- as.numeric(as.character(ARTeffect$timetoart))
summary(ARTeffect$timetoart)
ARTeffect$artefficacy <- as.numeric(as.character(ARTeffect$artefficacy))
summary(ARTeffect$artefficacy)

ggplot(data=subset(ARTeffect, timetoart==180)) +
  geom_line(aes(x=year, y=incidence*100, color=as.factor(artefficacy),group=artefficacy),size=2) +
  xlab("Year")+
  ylab("Incidence (per 100 py)")+
  theme_bw(base_size=16) +
  scale_x_continuous(breaks = seq(1980,2051,10),limits=c(1980,2051), expand = c(0,0)) +
  ggtitle("Sensitivity analysis of ART efficacy") +
  theme(legend.position="bottom") +
  guides(fill=guide_legend(title="ART efficacy"))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(strip.background = element_rect(colour="black", fill="white")) +
  #guides(col = guide_legend(nrow=2)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"))

setwd("C:\\Users\\aakullian\\Documents\\GitHub\\EMOD_eswatini\\LancetHIVParameterSweeps\\Figures")
ggsave("inc_trends_artefficacy_timetoart180.jpg", height=8, width=8)

ggplot() +
  geom_line(data=subset(ARTeffect, artefficacy==0.08), aes(x=year, y=incidence*100, color=factor(timetoart), group=timetoart),size=2) +
  xlab("Year")+
  ylab("Incidence (per 100 py)")+
  theme_bw(base_size=16) +
  scale_x_continuous(breaks = seq(1980,2050,10),limits=c(1980,2051), expand = c(0,0)) +
  ggtitle("Sensitivity analysis of time from infection to ART intitation") +
  theme(legend.position="bottom") +
  guides(color=guide_legend(title="Days from infection to treatment", nrow=1))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(strip.background = element_rect(colour="black", fill="white")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"))


setwd("C:\\Users\\aakullian\\Documents\\GitHub\\EMOD_eswatini\\LancetHIVParameterSweeps\\Figures")
ggsave("inc_trends_timetoart_artefficacy0.08.jpg", height=8, width=8)

ggplot(data=subset(ARTeffect, year==2016 | year==2020 | year == 2030 | year==2050)) +
  geom_raster(aes(y=timetoart, x=1-artefficacy, fill=incidence*100))+ 
  geom_point(aes(x=0.92, y=180), shape=3) +
  scale_fill_gradient(name="incidence per 100 py", low="blue", high="yellow",breaks=c(0.6,0.8,1,1.2,1.4,1.6))+
  stat_contour(aes(y=timetoart, x=1-artefficacy, z=incidence*100),
               color="black", size=0.1, linetype=1, binwidth=0.1) +
  theme(legend.position="bottom") +
  facet_wrap(~year)+
  xlab("ART efficacy") +
  ylab("Time from infection to ART initiation (months)")+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(keywidth = 2, keyheight = 1, nrow=1)) +
  theme(legend.position="bottom") +
  theme(strip.background = element_rect(colour="black", fill="white")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

setwd("C:\\Users\\aakullian\\Documents\\GitHub\\EMOD_eswatini\\LancetHIVParameterSweeps\\Figures")
ggsave("inc_raster_arteffic_timetoart.jpg", height=8, width=8)

#save incidence output so it doesn't have to be run again
setwd("C:\\Users\\aakullian\\Documents\\GitHub\\EMOD_eswatini\\LancetHIVParameterSweeps")
save.image(file="parameter_sweep_baseline.RData")
load("parameter_sweep_baseline.RData")



# 100 % ART sweeps

##################################################################################################
#Acute Paramter Sweep (includes a 2-dim sweep over acute duration and acute stage multiplier)
##################################################################################################
input_dir <- "C:/Users/aakullian/Documents/GitHub/EMOD_eswatini/ParamterSweepOutput_100pctART/acuteness"
primary_dirs = list.files(input_dir)

inc_values3 <- data.frame("year"=seq(1980,2056,1))
median_inc3 <- data.frame("year"=seq(1980,2056,1))
i = 0

dir="Acute_Stage_Infectivity_Multiplier-10--Acute_Duration_In_Months-1"
sub="0b0c6cf3-45fb-e911-a2c3-c4346bcb1551"

for(dir in primary_dirs) {
  sub_folders = list.files(paste(input_dir,dir,sep = "/"))
  inc_values3 <- data.frame("year"=seq(1980,2056,1))
  for(sub in sub_folders){
    i=i+1
    temp_table <- as.data.table(read.csv(file = paste(input_dir,dir,sub,"ReportHIVByAgeAndGender.csv",sep="/")))
    temp_table$Year2 <- floor((temp_table$Year-0.5))
    temp_table$Uninfected.Population = temp_table$Population-temp_table$Infected
    trajectories_IR.1a <- aggregate(Newly.Infected ~ Year2, subset(temp_table, Age>10 & Age<50), FUN=sum) #sums number of new infections in each year
    trajectories_IR.2 <- aggregate(Uninfected.Population ~ Year, subset(temp_table, Age>10 & Age<50), FUN=sum)
    trajectories_IR.2$Year2 <- floor(trajectories_IR.2$Year-0.5)
    trajectories_IR.2 <- trajectories_IR.2[!duplicated(trajectories_IR.2[c("Year2")]),] #remove second instance of duplicate rows
    trajectories_IR.2 <- trajectories_IR.2[-match("Year",names(trajectories_IR.2))]
    trajectories_IRoverall <- merge(trajectories_IR.1a, trajectories_IR.2, by=c("Year2"))
    trajectories_IRoverall$incidence <- trajectories_IRoverall$Newly.Infected / trajectories_IRoverall$Uninfected.Population
    inc_values3$newCol1 <- trajectories_IRoverall$incidence
    colnames(inc_values3)[ncol(inc_values3)] <- paste(sub)
    print(paste("working on sub-folder",sub,"sim",i,sep=" "))
  }
  median_inc3$val1=rowQuantiles(as.matrix(inc_values3[,2:ncol(inc_values3)]), probs = 0.5)
  colnames(median_inc3)[ncol(median_inc3)] <- paste(str_extract_all(dir,"\\(?[0-9,.]+\\)?")[[1]][1],str_extract_all(dir,"\\(?[0-9,.]+\\)?")[[1]][2], sep=",") 
  print(paste("working on", dir, "sub-folder",dir,sep=" "))
}

head(median_inc3)
rtest <- gather(median_inc3, parameter, incidence, `10,1`:`30,5`, factor_key=TRUE)
acute_stage <- data.frame(within(rtest, parameter<-data.frame(do.call('rbind', strsplit(as.character(parameter), ',', fixed=TRUE)))))
acute_stage <- data.frame("year"=acute_stage$year,"incidence"=acute_stage$incidence,"acute_multi"=acute_stage$parameter[1], "acute_dur"=acute_stage$parameter[2])
acute_stage$acute_multi <- acute_stage$X1
acute_stage$acute_dur <- acute_stage$X2
acute_stage <- acute_stage[,c(1:2,5:6)]
head(acute_stage)

acute_stage$acute_multi <- as.numeric(as.character(acute_stage$acute_multi))
summary(acute_stage$acute_multi)
acute_stage$acute_dur <- as.numeric(as.character(acute_stage$acute_dur))
summary(acute_stage$acute_dur)

head(acute_stage)

ggplot() +
  geom_line(data=subset(acute_stage, acute_dur==3), aes(x=year, y=incidence*100, color=as.factor(acute_multi), group=acute_multi),size=2) +
  xlab("Year")+
  ylab("Incidence (per 100 py)")+
  theme_bw(base_size=16) +
  scale_x_continuous(breaks = seq(1980,2050,10),limits=c(1980,2051), expand = c(0,0)) +
  scale_y_continuous(breaks = seq(0,4,1),limits=c(0,4), expand = c(0,0)) +
  ggtitle("Sensitivity analysis of acute stage infectivity multiplier") +
  theme(legend.position="bottom") +
  guides(color=guide_legend(title="Multiplier on acute stage infectivity"))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(strip.background = element_rect(colour="black", fill="white")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"))

setwd("C:\\Users\\aakullian\\Documents\\GitHub\\EMOD_eswatini\\LancetHIVParameterSweeps\\Figures")
ggsave("inc_trends_100ART_acutemult_acutedur3mo.jpg", height=8, width=8)

ggplot() +
  geom_line(data=subset(acute_stage, acute_multi==25), aes(x=year, y=incidence*100, color=as.factor(acute_dur), group=acute_dur),size=2) +
  xlab("Year")+
  ylab("Incidence (per 100 py)")+
  theme_bw(base_size=16) +
  scale_x_continuous(breaks = seq(1980,2050,10),limits=c(1980,2051), expand = c(0,0)) +
  scale_y_continuous(breaks = seq(0,4,1),limits=c(0,4), expand = c(0,0)) +
  ggtitle("Sensitivity analysis of acute stage duration") +
  theme(legend.position="bottom") +
  guides(color=guide_legend(title="Acute stage duration in months"))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(strip.background = element_rect(colour="black", fill="white")) +
  #guides(col = guide_legend(nrow=3)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"))

setwd("C:\\Users\\aakullian\\Documents\\GitHub\\EMOD_eswatini\\LancetHIVParameterSweeps\\Figures")
ggsave("inc_trends_100ART_acutedur_acutemulti26.jpg", height=8, width=8)

summary(acute_stage$incidence[acute_stage$year==2020])
ggplot(data=subset(acute_stage, year==2016 | year==2020 | year == 2030 | year==2050)) +
  geom_raster(aes(y=as.numeric(as.character(acute_multi)), x=as.numeric(as.character(acute_dur)), fill=incidence*100))+ 
  geom_point(aes(x=3, y=25), shape=3) +
  scale_fill_gradient(name="incidence per 100 py", low="blue", high="yellow",breaks=seq(0.2,1.1,0.15))+
  stat_contour(aes(y=acute_multi, x=acute_dur, z=incidence*100),
               color="black", size=0.1, linetype=1, binwidth=0.1) +
  theme(legend.position="bottom") +
  facet_wrap(~year)+
  xlab("acute phase duration (months)") +
  ylab("acute phase infectivity")+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(keywidth = 2, keyheight = 1, nrow=1)) +
  theme(legend.position="bottom") +
  theme(strip.background = element_rect(colour="black", fill="white")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

setwd("C:\\Users\\aakullian\\Documents\\GitHub\\EMOD_eswatini\\LancetHIVParameterSweeps\\Figures")
ggsave("inc_raster100ART_acute_dur_multi.jpg", height=8, width=8)

##################################################################################################
#Parameter sweep of ART reduce acquire & delay from infection to ART initiation 
##################################################################################################

# Set working directories
input_dir <- "C:/Users/aakullian/Documents/GitHub/EMOD_eswatini/ParamterSweepOutput_100pctART/delay_and_supression"
primary_dirs = list.files(input_dir)
inc_values4 <- data.frame("year"=seq(1980,2056,1))
median_inc4 <- data.frame("year"=seq(1980,2056,1))

for(dir in primary_dirs) {
  i=0
  sub_folders = list.files(paste(input_dir,dir,sep = "/"))
  inc_values4 <- data.frame("year"=seq(1980,2056,1))
  for(sub in sub_folders){
    i=i+1
    temp_table <- as.data.table(read.csv(file = paste(input_dir,dir,sub,"ReportHIVByAgeAndGender.csv",sep="/")))
    temp_table$Year2 <- floor((temp_table$Year-0.5))
    temp_table$Uninfected.Population = temp_table$Population-temp_table$Infected
    trajectories_IR.1a <- aggregate(Newly.Infected ~ Year2, subset(temp_table, Age>10 & Age<50), FUN=sum) #sums number of new infections in each year
    trajectories_IR.2 <- aggregate(Uninfected.Population ~ Year, subset(temp_table, Age>10 & Age<50), FUN=sum)
    trajectories_IR.2$Year2 <- floor(trajectories_IR.2$Year-0.5)
    trajectories_IR.2 <- trajectories_IR.2[!duplicated(trajectories_IR.2[c("Year2")]),] #remove second instance of duplicate rows
    trajectories_IR.2 <- trajectories_IR.2[-match("Year",names(trajectories_IR.2))]
    trajectories_IRoverall <- merge(trajectories_IR.1a, trajectories_IR.2, by=c("Year2"))
    trajectories_IRoverall$incidence <- trajectories_IRoverall$Newly.Infected / trajectories_IRoverall$Uninfected.Population
    inc_values4$newCol1 <- trajectories_IRoverall$incidence
    colnames(inc_values4)[ncol(inc_values4)] <- paste(sub)
    print(paste("working on folder", dir, "sub-folder",sub,"sim",i,sep=" "))
  }
  median_inc4$val1=rowQuantiles(as.matrix(inc_values4[,2:ncol(inc_values4)]), probs = 0.5)
  colnames(median_inc4)[ncol(median_inc4)] <- paste(str_extract_all(dir,"\\(?[0-9,.]+\\)?")[[1]][1],str_extract_all(dir,"\\(?[0-9,.]+\\)?")[[1]][2], sep=",") 
  print(paste("working on sub-folder",dir,sep=" "))
}

rtest <- gather(median_inc4, parameter, incidence, `0,0.0`:`60,0.2`, factor_key=TRUE)
head(rtest,100)
ARTeffect <- data.frame(within(rtest, parameter<-data.frame(do.call('rbind', strsplit(as.character(parameter), ',', fixed=TRUE)))))
ARTeffect <- data.frame("year"=ARTeffect$year,"incidence"=ARTeffect$incidence,"timetoart"=ARTeffect$parameter[1], "artefficacy"=ARTeffect$parameter[2])
ARTeffect$timetoart <- ARTeffect$X1
ARTeffect$artefficacy <- ARTeffect$X2
ARTeffect <- ARTeffect[,c(1:2,5:6)]
class(ARTeffect$timetoart)
class(ARTeffect$artefficacy)
head(ARTeffect)

ARTeffect$timetoart <- as.numeric(as.character(ARTeffect$timetoart))
summary(ARTeffect$timetoart)
ARTeffect$artefficacy <- as.numeric(as.character(ARTeffect$artefficacy))
summary(ARTeffect$artefficacy)

ggplot(data=subset(ARTeffect, timetoart==180)) +
  geom_line(aes(x=year, y=incidence*100, color=as.factor(artefficacy),group=artefficacy),size=2) +
  xlab("Year")+
  ylab("Incidence (per 100 py)")+
  theme_bw(base_size=16) +
  scale_x_continuous(breaks = seq(1980,2051,10),limits=c(1980,2051), expand = c(0,0)) +
  ggtitle("Sensitivity analysis of ART efficacy") +
  theme(legend.position="bottom") +
  guides(fill=guide_legend(title="ART efficacy"))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(strip.background = element_rect(colour="black", fill="white")) +
  #guides(col = guide_legend(nrow=2)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"))

setwd("C:\\Users\\aakullian\\Documents\\GitHub\\EMOD_eswatini\\LancetHIVParameterSweeps\\Figures")
ggsave("inc_trends_artefficacy_timetoart180.jpg", height=8, width=8)

ggplot() +
  geom_line(data=subset(ARTeffect, artefficacy==0.08), aes(x=year, y=incidence*100, color=factor(timetoart), group=timetoart),size=2) +
  xlab("Year")+
  ylab("Incidence (per 100 py)")+
  theme_bw(base_size=16) +
  scale_x_continuous(breaks = seq(1980,2050,10),limits=c(1980,2051), expand = c(0,0)) +
  ggtitle("Sensitivity analysis of time from infection to ART intitation") +
  theme(legend.position="bottom") +
  guides(color=guide_legend(title="Days from infection to treatment", nrow=1))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(strip.background = element_rect(colour="black", fill="white")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"))


setwd("C:\\Users\\aakullian\\Documents\\GitHub\\EMOD_eswatini\\LancetHIVParameterSweeps\\Figures")
ggsave("inc_trends_timetoart_artefficacy0.08.jpg", height=8, width=8)

ggplot(data=subset(ARTeffect, year==2016 | year==2020 | year == 2030 | year==2050)) +
  geom_raster(aes(y=timetoart, x=1-artefficacy, fill=incidence*100))+ 
  geom_point(aes(x=0.92, y=180), shape=3) +
  scale_fill_gradient(name="incidence per 100 py", low="blue", high="yellow",breaks=c(0.6,0.8,1,1.2,1.4,1.6))+
  stat_contour(aes(y=timetoart, x=1-artefficacy, z=incidence*100),
               color="black", size=0.1, linetype=1, binwidth=0.1) +
  theme(legend.position="bottom") +
  facet_wrap(~year)+
  xlab("ART efficacy") +
  ylab("Time from infection to ART initiation (months)")+
  theme_bw(base_size=14) +
  guides(fill = guide_legend(keywidth = 2, keyheight = 1, nrow=1)) +
  theme(legend.position="bottom") +
  theme(strip.background = element_rect(colour="black", fill="white")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

setwd("C:\\Users\\aakullian\\Documents\\GitHub\\EMOD_eswatini\\LancetHIVParameterSweeps\\Figures")
ggsave("inc_raster_arteffic_timetoart.jpg", height=8, width=8)

#save incidence output so it doesn't have to be run again
setwd("C:\\Users\\aakullian\\Documents\\GitHub\\EMOD_eswatini\\LancetHIVParameterSweeps")
save.image(file="parameter_sweep_baseline.RData")
load("parameter_sweep_baseline.RData")









##################################################################################################
#Paramters = Granich et al 2009
##################################################################################################

# Set working directories
input_dir <- "C:\\Users\\aakullian\\Dropbox (IDM)\\GitHub\\EMOD_eswatini\\GranichParmSweeps"
primary_dirs = list.files(input_dir)
inc_values5 <- data.frame("Year2"=seq(1980,2056,1))

i=0
for(dir in primary_dirs) {
  i=i+1
  temp_table <- as.data.table(read.csv(file = paste(input_dir,dir,"ReportHIVByAgeAndGender.csv",sep="/")))
  temp_table$Year2 <- floor((temp_table$Year-0.5))
  temp_table$Uninfected.Population = temp_table$Population-temp_table$Infected
  trajectories_IR.1a <- aggregate(Newly.Infected ~ Year2, subset(temp_table, Age>10 & Age<50), FUN=sum) #sums number of new infections in each year
  trajectories_IR.2 <- aggregate(Uninfected.Population ~ Year, subset(temp_table, Age>10 & Age<50), FUN=sum)
  trajectories_IR.2$Year2 <- floor(trajectories_IR.2$Year-0.5)
  trajectories_IR.2 <- trajectories_IR.2[!duplicated(trajectories_IR.2[c("Year2")]),] #remove second instance of duplicate rows
  trajectories_IR.2 <- trajectories_IR.2[-match("Year",names(trajectories_IR.2))]
  trajectories_IRoverall <- merge(trajectories_IR.1a, trajectories_IR.2, by=c("Year2"))
  trajectories_IRoverall$incidence <- trajectories_IRoverall$Newly.Infected / trajectories_IRoverall$Uninfected.Population
  trajectories_IRoverall$sim.id <- paste(dir)
  if (i == 1){
    inc_values5$incidence=trajectories_IRoverall$incidence
    inc_values5$sim.id=trajectories_IRoverall$sim.id
  }
  else{
    inc_values5 <- rbind(inc_values5, trajectories_IRoverall[,c(1,4,5)])
  }
  print(paste("working on folder", dir,sep=" "))
}

head(trajectories_IRoverall)
head(inc_values5)

ggplot(data=subset(inc_values5)) +
  geom_line(aes(x=Year2, y=incidence*100,group=sim.id),color="lightblue",size=2) +
  geom_smooth(aes(x=Year2, y=incidence*100),method="loess", span=0.1, se = T, size=1, color="blue", linetype=1) +
  #geom_line(data=subset(ARTeffect, timetoart==180 & artefficacy==0.08),aes(x=year, y=incidence*100, color=as.factor(1-artefficacy),group=1-artefficacy),size=2) +
  #geom_line(data=subset(ARTeffect, timetoart==180 & artefficacy==0.08),aes(x=year, y=incidence*100, color=as.factor(1-artefficacy),group=1-artefficacy), color="black", linetype=2,size=2) +
  geom_hline(yintercept=0.1, linetype="dashed", color = "red", size=1) +
  xlab("Year")+
  ylab("Incidence (per 100 py)")+
  theme_bw(base_size=10) +
  scale_x_continuous(breaks = seq(1980,2051,10),limits=c(1980,2051), expand = c(0,0)) +
  scale_y_continuous(breaks = seq(0,5,1),limits=c(0,5), expand = c(0,0)) +
  ggtitle("Incidence under Granich et al., 2009 paramterization") +
  theme(legend.position="bottom") +
  guides(color=guide_legend(title="ART efficacy", nrow=1))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(strip.background = element_rect(colour="black", fill="white")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"))

setwd("C:\\Users\\aakullian\\Documents\\GitHub\\EMOD_eswatini\\LancetHIVParameterSweeps\\Figures")
ggsave("inc_trends_granich2009_estimates.jpg", height=8, width=8)


#Sexual risk behavior report

##################################################################################################
#ART stages in 100% ART scenario
##################################################################################################

library(rjson)

#Read in the Accessibility_and_Risk_IP_Overlay

Accessibility_and_Risk_IP_Overlay <- fromJSON(file = "C:\\Users\\aakullian\\Dropbox (IDM)\\GitHub\\EMOD_eswatini\\Calibration\\SeparateInput\\Accessibility_and_Risk_IP_Overlay.json")
Accessibility_and_Risk_IP_Overlay_json <- toJSON(Accessibility_and_Risk_IP_Overlay)

# Set working directories
input_dir <- "C:\\Users\\aakullian\\Dropbox (IDM)\\GitHub\\EMOD_eswatini\\Calibration\\SeparateInput\\ART100pctTransmissionCounter"
primary_dirs = list.files(input_dir)
dir="01d11cc3-ffd6-e811-a2bd-c4346bcb1555"

for(dir in primary_dirs) {
  write(Accessibility_and_Risk_IP_Overlay_json, file = paste(input_dir,dir,"Accessibility_and_Risk_IP_Overlay.json",sep="\\"))
}

#read in reportHIVbyageandgender files
input_dir <- "C:\\Users\\aakullian\\Dropbox (IDM)\\GitHub\\EMOD_eswatini\\Baseline_transmission_ARTstate"
primary_dirs = list.files(input_dir)
inc_values6 <- data.frame("Year2"=seq(1980,2056,1))

i=0
for(dir in primary_dirs) {
  i=i+1
  temp_table <- as.data.table(read.csv(file = paste(input_dir,dir,"ReportHIVByAgeAndGender.csv",sep="/")))
  temp_table$Year2 <- floor((temp_table$Year-0.5))
  temp_table$Uninfected.Population = temp_table$Population-temp_table$Infected
  trajectories_IR.1a <- aggregate(cbind(Newly.Infected,Transmitters) ~ Year2+IP_Key.ARTstate, subset(temp_table, Age>10 & Age<50), FUN=sum) #sums number of new infections in each year
  trajectories_IR.2 <- aggregate(cbind(Uninfected.Population,Infected,Population) ~ Year, subset(temp_table, Age>10 & Age<50), FUN=sum)
  trajectories_IR.2$Year2 <- floor(trajectories_IR.2$Year-0.5)
  trajectories_IR.2 <- trajectories_IR.2[!duplicated(trajectories_IR.2[c("Year2")]),] #remove second instance of duplicate rows
  trajectories_IR.2 <- trajectories_IR.2[-match("Year",names(trajectories_IR.2))]
  trajectories_IRoverall <- merge(trajectories_IR.1a, trajectories_IR.2, by=c("Year2"))
  trajectories_IRoverall$incidence <- trajectories_IRoverall$Newly.Infected / trajectories_IRoverall$Uninfected.Population
  trajectories_IRoverall$sim.id <- paste(dir)
  if (i == 1){
    inc_values6 <- trajectories_IRoverall
  }
  else{
    inc_values6 <- rbind(inc_values6, trajectories_IRoverall)
  }
  print(paste("working on folder", i, dir,sep=" "))
}

head(temp_table)
head(inc_values6,100)

inc_values6.agg <- aggregate(Transmitters ~ Year2+sim.id, inc_values6, FUN=sum)
inc_values6.agg.m <- merge(inc_values6, inc_values6.agg, by=c("Year2","sim.id"))
inc_values6.agg.m$IP_Key.ARTstate_P <- inc_values6.agg.m$Transmitters.x / inc_values6.agg.m$Transmitters.y
inc_values6.agg.m$IP_Key.ARTstate_INC <- inc_values6.agg.m$Transmitters.x / inc_values6.agg.m$Infected

summary(inc_values6.agg.m$IP_Key.ARTstate_P)
summary(inc_values6.agg.m$IP_Key.ARTstate_INC)
names(inc_values6.agg.m)
head(subset(inc_values6.agg.m[c(1,3:8,10:12)],IP_Key.ARTstate_INC>1),100)

#mean value of all sims by year
table(inc_values6.agg.m$IP_Key.ARTstate_P)
inc_values6.agg.mean <- aggregate(IP_Key.ARTstate_P ~ IP_Key.ARTstate+Year2, inc_values6.agg.m, FUN=mean)
head(subset(inc_values6.agg.mean, Year2==2030 | Year2==2050))

table(inc_values6.agg.m$IP_Key.ARTstate)
inc_values6.agg.m$IP_Key.ARTstate_f <- ordered(inc_values6.agg.m$IP_Key.ARTstate, levels=c("NeverOnART","OnART","OffART"))
ggplot(data=subset(inc_values6.agg.m)) +
  geom_line(aes(x=Year2, y=IP_Key.ARTstate_P,group=interaction(IP_Key.ARTstate_f,sim.id), color=IP_Key.ARTstate_f, alpha=0.02),size=1,alpha=0.02) +
  geom_smooth(aes(x=Year2, y=IP_Key.ARTstate_P, color=IP_Key.ARTstate_f, group=IP_Key.ARTstate_f),method="loess", span=0.2, se=F,size=2) +
  xlab("Year")+
  ylab("Proportion of transmissions")+
  theme_bw(base_size=16) +
  scale_x_continuous(breaks = seq(1980,2051,10),limits=c(1980,2051), expand = c(0,0)) +
  scale_y_continuous(breaks = seq(0,1,0.1),limits=c(0,1), expand = c(0,0)) +
  scale_color_manual(labels = c("Never initiated ART","On ART", "Dropped out from ART"), values = c("red", "blue", "purple")) +
  ggtitle("Transmissions by ART status of transmitter") +
  theme(legend.position="bottom") +
  guides(color=guide_legend(title="ART status", nrow=1))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(strip.background = element_rect(colour="black", fill="white")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"))

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\GitHub\\EMOD_eswatini\\LancetHIVParameterSweeps\\Figures")
ggsave("TransmittersByARTstate_proportions_baseline.jpg", height=8, width=8)

#art scale-up
i=0
for(dir in primary_dirs) {
  i=i+1
  temp_table <- as.data.table(read.csv(file = paste(input_dir,dir,"ReportHIVByAgeAndGender.csv",sep="/")))
  art_temp <- aggregate(cbind(On_ART,Infected) ~ Year, subset(temp_table, Age>10 & Age<50), FUN=sum) #sums number of new infections in each year
  art_temp$ARTcoverage <- art_temp$On_ART / art_temp$Infected
  art_temp$sim.id <- paste(dir)
  if (i == 1){
    art_values1 <- art_temp
  }
  else{
    art_values1 <- rbind(art_values1, art_temp)
  }
  print(paste("working on folder", i, dir,sep=" "))
}

head(art_values1)

ggplot(data=subset(inc_values6.agg.m)) +
  geom_line(aes(x=Year2, y=IP_Key.ARTstate_P,group=interaction(IP_Key.ARTstate_f,sim.id), color=IP_Key.ARTstate_f, alpha=0.02),size=1,alpha=0.02) +
  geom_smooth(aes(x=Year2, y=IP_Key.ARTstate_P, color=IP_Key.ARTstate_f, group=IP_Key.ARTstate_f),method="loess", span=0.2, se=F,size=2) +
  geom_smooth(data=art_values1, aes(x=Year, y=ARTcoverage, color="black"),method="loess", color="black", span=0.2, se=F,size=1,linetype=2) +
  xlab("Year")+
  ylab("Proportion of transmissions")+
  theme_bw(base_size=16) +
  scale_x_continuous(breaks = seq(1980,2051,10),limits=c(1980,2051), expand = c(0,0)) +
  scale_y_continuous(breaks = seq(0,1,0.1),limits=c(0,1), expand = c(0,0)) +
  scale_color_manual(labels = c("Never initiated ART","On ART", "Dropped out from ART", "ART coverage"), values = c("red", "blue", "purple","black")) +
  ggtitle("Transmissions by ART status of transmitter") +
  theme(legend.position="bottom") +
  guides(color=guide_legend(title="ART status", nrow=1))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(strip.background = element_rect(colour="black", fill="white")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"))

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\GitHub\\EMOD_eswatini\\LancetHIVParameterSweeps\\Figures")
ggsave("TransmittersByARTstate_proportions_ARTdata_baseline.jpg", height=8, width=8)
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
head(transmissionreport)
transmissionreport$src_age_int <- floor(transmissionreport$SRC_AGE / 365.25)
transmissionreport$dest_age_int <- floor(transmissionreport$DEST_AGE / 365.25)
head(transmissionreport)
transmissionreport$acquisition <- 1
transmissionreport <- subset(transmissionreport, dest_age_int > 14 & src_age_int > 14)
#transmissionreport$age.at.acquisition.cat <- cut(transmissionreport$dest_age_int, c(15,25,35,200), right=F)
#transmissionreport$age.at.transmission.cat <- cut(transmissionreport$src_age_int, c(15,25,35,200), right=F)
transmissionreport$year <- transmissionreport$YEAR
names(transmissionreport)

#################################################################################################
#Descriptive stats
#################################################################################################
# #baseline transmission report from one run:
# wd13 <- "C:\\Users\\aakullian\\Dropbox (IDM)\\HIV\\EMOD\\Swaziland_05June2018\\Scenarios_calibration_run3_v05Sep2018_newcampaign_FINAL\\ScenariosOct_FINAL\\Baseline_TransmissionReport\\TransmissionReport\\"
# setwd(paste0(wd13))
# files <- list.files(full.names = F)
# head(files)
# f <- paste0(wd13, files)
# length(f)
# i = 2
# 
# transmissionreport.master <- read.csv(f[i])
# transmissionreport.master$scenario <- "baseline"
# transmissionreport.master$sim.id <- paste0(files[i])
# names(transmissionreport)
# table(transmissionreport$DEST_INFECTED) #all are infected
# table(transmissionreport$SRC_VIRAL_LOAD)
# library(tableone)
# 
# #CREATE DATA FRAME OF acquisitions (the destination of transmissions)
# acquisitions <- data.frame("sim.id"=transmissionreport$sim.id, "ID" = transmissionreport$DEST_ID, "dest.gender"=transmissionreport$DEST_GENDER, "year.acq" = floor(transmissionreport$year), "age.at.acq" = floor(transmissionreport$dest_age_int),
#                            "current_relationship_count" = transmissionreport$DEST_current_relationship_count,
#                            "lifetime_relationship_count" = transmissionreport$DEST_lifetime_relationship_count,
#                            "relationships_in_last_6_months" = transmissionreport$DEST_relationships_in_last_6_months,
#                            "circumcised" = transmissionreport$DEST_CIRCUMSIZED,
#                            "sti"=transmissionreport$DEST_STI,
#                            "sourcestage"=transmissionreport$SRC_STAGE,
#                            "sourceCD4"=transmissionreport$SRC_CD4)
# acquisitions <- acquisitions[order(acquisitions$sim.id, acquisitions$ID, acquisitions$year.acq),]
# acquisitions$Gender[acquisitions$dest.gender==1] <- "Female"
# acquisitions$Gender[acquisitions$dest.gender==0] <- "Male"
# head(acquisitions)
# 
# ## Create a variable list
# vars <- c("year.acq","age.at.acq","Gender","current_relationship_count","lifetime_relationship_count",
#           "relationships_in_last_6_months","sti","circumcised","sourcestage","sourceCD4")
# 
# catVars <- c("sti","circumcised","sourcestage")
# 
# ## Create Table 1 stratified by trt
# tableOne <- CreateTableOne(vars = vars, factorVars = catVars, strata = c("Gender"), data = acquisitions)
# tableOne

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























