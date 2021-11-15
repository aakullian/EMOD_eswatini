##################################################################################################
#set working directories for scenarios
##################################################################################################
wd1 <- "C:\\Users\\aakullian\\Dropbox (IDM)\\HIV\\EMOD\\Swaziland_05June2018\\Scenarios_calibration_run3_v05Sep2018_newcampaign_FINAL\\ScenariosOct_FINAL\\baseline\\ReportHIVByAgeAndGender\\"
wd2 <- "C:\\Users\\aakullian\\Dropbox (IDM)\\HIV\\EMOD\\Swaziland_05June2018\\Scenarios_calibration_run3_v05Sep2018_newcampaign_FINAL\\ScenariosOct_FINAL\\noARTnoVMMC\\ReportHIVByAgeAndGender\\"
wd3 <- "C:\\Users\\aakullian\\Dropbox (IDM)\\HIV\\EMOD\\Swaziland_05June2018\\Scenarios_calibration_run3_v05Sep2018_newcampaign_FINAL\\ScenariosOct_FINAL\\noART\\ReportHIVByAgeAndGender\\"
wd4 <- "C:\\Users\\aakullian\\Dropbox (IDM)\\HIV\\EMOD\\Swaziland_05June2018\\Scenarios_calibration_run3_v05Sep2018_newcampaign_FINAL\\ScenariosOct_FINAL\\agetarg909090\\ReportHIVByAgeAndGender\\"
wd5 <- "C:\\Users\\aakullian\\Dropbox (IDM)\\HIV\\EMOD\\Swaziland_05June2018\\Scenarios_calibration_run3_v05Sep2018_newcampaign_FINAL\\ScenariosOct_FINAL\\hist909090\\ReportHIVByAgeAndGender\\"
wd6 <- "C:\\Users\\aakullian\\Dropbox (IDM)\\HIV\\EMOD\\Swaziland_05June2018\\Scenarios_calibration_run3_v05Sep2018_newcampaign_FINAL\\ScenariosOct_FINAL\\ART100pct\\ReportHIVByAgeAndGender\\"
wd9 <- "C:\\Users\\aakullian\\Dropbox (IDM)\\HIV\\EMOD\\Swaziland_05June2018\\Scenarios_calibration_run3_v05Sep2018_newcampaign_FINAL\\ScenariosOct_FINAL\\HIGHRISK75\\ReportHIVByAgeAndGender\\"
wd10 <- "C:\\Users\\aakullian\\Dropbox (IDM)\\HIV\\EMOD\\Swaziland_05June2018\\Scenarios_calibration_run3_v05Sep2018_newcampaign_FINAL\\ScenariosOct_FINAL\\HIGHRISK100_909090_riskgroup_3\\ReportHIVByAgeAndGender\\"
wd11 <- "C:\\Users\\aakullian\\Dropbox (IDM)\\HIV\\EMOD\\Swaziland_05June2018\\Scenarios_calibration_run3_v05Sep2018_newcampaign_FINAL\\ScenariosOct_FINAL\\LOWRISK100_909090_riskgroup_3\\ReportHIVByAgeAndGender\\"
wd12 <- "C:\\Users\\aakullian\\Dropbox (IDM)\\HIV\\EMOD\\Swaziland_05June2018\\Scenarios_calibration_run3_v05Sep2018_newcampaign_FINAL\\ScenariosOct_FINAL\\baseline_riskgroup\\ReportHIVByAgeAndGender\\"

#Bring in 250 parameter sets
setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\HIV\\EMOD\\Swaziland_05June2018\\Scenarios_calibration_run3_v05Sep2018_newcampaign_FINAL")
infile.parms <- "ReEvaluatedLikelihood_full_MASTER_new_r"
infile.parms <- read.csv(paste0('./', infile.parms, '.csv'))
names(infile.parms)

stack(lapply(infile.parms, quantile, prob = c(0.25,0.5,0.75), names = T))

infile.parms %>%
  summarise_all(funs(list(quantile(., probs = c(0.5,0.25,0.75) )))) %>%
  gather("Column", "Median p25 p75")

#reshape to create ggplots of histograms
ggplot(gather(infile.parms), aes(value)) + 
  geom_histogram(bins = 10, fill="white", col="black") + 
  facet_wrap(~key, scales = 'free_x',ncol=5) +
  xlab("Parameter value") +
  ylab("Count") +
  theme(strip.background = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\BMGF Lancet Paper\\Final_Plots")
ggsave("parameter_histogram.jpg", height=12, width=12)

#Bring in incidence reference data
setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\HIV\\EMOD\\Swaziland_05June2018\\Data\\Swaziland_point_estimates")
infile4 <- "Swaziland_Incidence_Data"
incidence.data <- read.csv(paste0('./', infile4, '.csv'))

#Bring in prevalence reference data
setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\HIV\\EMOD\\Swaziland_05June2018\\Data\\Swaziland_point_estimates")
infileprev <- "SWAZILAND_nationalprevalence"
prevalence <- read.csv(paste0('./', infileprev, '.csv'))
prevalence$Age <- prevalence$start.age
infileprev <- "SWAZILAND_nationalprevalence_bothgenders"
prevalence.bothgenders <- read.csv(paste0('./', infileprev, '.csv'))
head(prevalence.bothgenders)

#Bring in ART coverage data
setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\HIV\\EMOD\\Swaziland_05June2018\\Data\\Swaziland_point_estimates")
infile1 <- "SWAZILAND_calibration_nationalARTprevalence"
artdata <- read.csv(paste0('./', infile1, '.csv'))
artdata$Agecat <- artdata$AgeBin
artdata$ART_coverage <- artdata$NationalARTPrevalence
labs <- c("1" = "Women", "0" = "Men")

###########################################
#Load scenarios
###########################################

#Long file
setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\HIV\\EMOD\\Swaziland_05June2018\\Scenarios_calibration_run3_v05Sep2018_newcampaign_FINAL\\ScenariosOct_FINAL")
reporthivbyageandgender.master.final <- read.csv("reporthivbyageandgender.master.final.csv")

#condensed file
setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\HIV\\EMOD\\Swaziland_05June2018\\Scenarios_calibration_run3_v05Sep2018_newcampaign_FINAL\\ScenariosOct_FINAL")
reporthivbyageandgender.master <- read.csv("reporthivbyageandgender.master.final.condensed.csv")
table(reporthivbyageandgender.master$scenario)
table(reporthivbyageandgender.master$Agecat)
head(reporthivbyageandgender.master)

# #bring in risk scenarios
# reporthivbyageandgender.master.final.risk <- read.csv("reporthivbyageandgender.master.final.risk.csv")
# head(reporthivbyageandgender.master.final.risk)

#Calculate Incidence overall from condensed
table(reporthivbyageandgender.master$Agecat)
reporthivbyageandgender.master.final <- reporthivbyageandgender.master
reporthivbyageandgender.master.final$Year2 <- floor((reporthivbyageandgender.master.final$Year-0.5))
reporthivbyageandgender.master.final$Uninfected.Population = reporthivbyageandgender.master.final$Population-reporthivbyageandgender.master.final$Infected
trajectories_IR.1a <- aggregate(Newly.Infected ~ Year2+sim.id+scenario, subset(reporthivbyageandgender.master.final), FUN=sum) #sums number of new infections in each year
trajectories_IR.2 <- aggregate(Uninfected.Population ~ Year+sim.id+scenario, subset(reporthivbyageandgender.master.final), FUN=sum)
trajectories_IR.2$Year2 <- floor(trajectories_IR.2$Year-0.5)
trajectories_IR.2 <- trajectories_IR.2[!duplicated(trajectories_IR.2[c("Year2","sim.id","scenario")]),] #remove second instance of duplicate rows
trajectories_IR.2 <- trajectories_IR.2[-match("Year",names(trajectories_IR.2))]
trajectories_IRoverall <- merge(trajectories_IR.1a, trajectories_IR.2, by=c("Year2","sim.id","scenario"))
trajectories_IRoverall$incidence <- trajectories_IRoverall$Newly.Infected / trajectories_IRoverall$Uninfected.Population
summary(trajectories_IRoverall$incidence[trajectories_IRoverall$scenario=="baseline"])

#add in noVMMC with ART plot

#noVMMC (with ART) #Added this in upon reviewer request
library(data.table)
#head(trajectories_IRoverall) # <- make to look like this
dir = "C:\\Users\\aakullian\\Dropbox (IDM)\\GitHub\\EMOD_eswatini\\noVMMCwithART"
sub_folders = list.files(paste(dir,sep = "/"))
sub <- "00056499-03fd-e911-a2c3-c4346bcb1551"
temp_table <- as.data.table(read.csv(file = "C:\\Users\\aakullian\\Dropbox (IDM)\\GitHub\\EMOD_eswatini\\ReportHIVByAgeAndGender.csv"))

i = 0
for(sub in sub_folders){
  i = i+1
  temp_table <- as.data.table(read.csv(file = paste(dir,sub,"ReportHIVByAgeAndGender.csv",sep="/")))
  temp_table <- subset(temp_table, Age>10 & Age <50)
  temp_table$Year2 <- floor((temp_table$Year-0.5))
  temp_table$Uninfected.Population = temp_table$Population-temp_table$Infected
  trajectories_IR.1a <- aggregate(Newly.Infected ~ Year2, subset(temp_table), FUN=sum) #sums number of new infections in each year
  trajectories_IR.2 <- aggregate(Uninfected.Population ~ Year, subset(temp_table), FUN=sum)
  trajectories_IR.2$Year2 <- floor(trajectories_IR.2$Year-0.5)
  trajectories_IR.2 <- trajectories_IR.2[!duplicated(trajectories_IR.2[c("Year2")]),] #remove second instance of duplicate rows
  trajectories_IR.2 <- trajectories_IR.2[-match("Year",names(trajectories_IR.2))]
  trajectories_IRoverall.novmmc.art <- merge(trajectories_IR.1a, trajectories_IR.2, by=c("Year2"))
  trajectories_IRoverall.novmmc.art$incidence <- trajectories_IRoverall.novmmc.art$Newly.Infected / trajectories_IRoverall.novmmc.art$Uninfected.Population
  trajectories_IRoverall.novmmc.art <- trajectories_IRoverall.novmmc.art[,c(1,4)]
  trajectories_IRoverall.novmmc.art$sim.id <- paste(sub)
  trajectories_IRoverall.novmmc.art$scenario <- "noVMMCwithART"
  if (i==1){
    trajectories_IRoverall.novmmc.art.temp <- trajectories_IRoverall.novmmc.art
  } else {
    trajectories_IRoverall.novmmc.art.temp <- rbind(trajectories_IRoverall.novmmc.art.temp, trajectories_IRoverall.novmmc.art)
    print(paste0("Why hello there, I'm working on noVMMCwithART folder",i))
  }
}

names(trajectories_IRoverall.novmmc.art.temp)

summary(trajectories_IRoverall.novmmc.art.temp$incidence)

ggplot() +
  geom_line(data=subset(trajectories_IRoverall.novmmc.art.temp, Year2 < year), aes(x=Year2, y=incidence*100, group=sim.id)) +
  #geom_smooth(data=subset(trajectories_IRoverall.novmmc.art.temp, Year2 < year), aes(x=Year2, y=incidence*100),method="loess", span=span, se = F, size=2) +
  xlab("Year")+
  ylab("Incidence (per 100 py)")+
  theme_bw(base_size=24) +
  guides(fill = guide_legend(keywidth = 2, keyheight = 1)) +
  geom_hline(yintercept=0.1, linetype=2) +
  scale_y_continuous(breaks = seq(0,6,0.5),limits=c(0,5), expand = c(0,0)) +
  scale_x_continuous(breaks = seq(1975,year,5), limits=c(1975, year),expand = c(0,0)) +
  theme(legend.position="bottom") +
  theme(legend.title=element_blank())+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(strip.background = element_rect(colour="black", fill="white")) +
  guides(col = guide_legend(nrow=3)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"))


head(trajectories_IRoverall[,c(1,2,3,6)])
trajectories_IRoverall.merge <- rbind(trajectories_IRoverall[,c(1,2,3,6)], trajectories_IRoverall.novmmc.art.temp)
table(trajectories_IRoverall.merge$scenario)
head(trajectories_IRoverall.merge[trajectories_IRoverall.merge$scenario=="noVMMCwithART",],100)
summary(trajectories_IRoverall.merge$incidence[trajectories_IRoverall.merge$scenario=="noVMMCwithART"])

year = 2051
span = 0.2
trajectories_IRoverall.merge$scenario_f <- ordered(trajectories_IRoverall.merge$scenario, levels=c("noVMMCwithART","ART100pct","noARTnoVMMC","noART","909090","baseline"))
table(trajectories_IRoverall.merge$scenario_f)
ggplot() +
  geom_smooth(data=subset(trajectories_IRoverall.merge, Year2 < year & (scenario_f=="baseline"|scenario_f=="noARTnoVMMC"|scenario_f=="noART"|scenario_f=="noVMMCwithART")), aes(x=Year2, y=incidence*100, color=scenario_f),method="loess", span=span, se = F, size=2) +
  xlab("Year")+
  ylab("Incidence (per 100 py)")+
  theme_bw(base_size=24) +
  scale_color_manual(labels = c("No VMMC with ART","No VMMC without ART", "VMMC without ART", "VMMC with ART (baseline)"), values = c("yellow","red", "orange", "blue")) +
  guides(fill = guide_legend(keywidth = 2, keyheight = 1)) +
  geom_hline(yintercept=0.1, linetype=2) +
  scale_y_continuous(breaks = seq(0,6,0.5),limits=c(0,5), expand = c(0,0)) +
  scale_x_continuous(breaks = seq(1975,year,5), limits=c(1975, year),expand = c(0,0)) +
  theme(legend.position="bottom") +
  theme(legend.title=element_blank())+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(strip.background = element_rect(colour="black", fill="white")) +
  guides(col = guide_legend(nrow=2),byrow=T) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"))

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\BMGF Lancet Paper\\Revised_Plots")
ggsave("swaziland.incidence.6scenarios_FINAL.jpg", height=8, width=8)

###########################################
#Scenarios
###########################################

#baseline
setwd(paste0(wd1))
files <- list.files(full.names = F)
head(files)
f <- paste0(wd1, files)
length(f)

for (i in seq(1,250,1)){
  reporthivbyageandgender <- read.csv(f[i])
  reporthivbyageandgender$scenario <- "baseline"
  reporthivbyageandgender$sim.id <- paste0(files[i])
  if (i==1){
    reporthivbyageandgender.master <- reporthivbyageandgender
  } else {
    reporthivbyageandgender.master <- rbind(reporthivbyageandgender.master, reporthivbyageandgender)
    print(paste0("Why hello there, I'm working on baseline folder",i))
  }
}

reporthivbyageandgender.master.baseline <- reporthivbyageandgender.master
names(reporthivbyageandgender.master.baseline)
reporthivbyageandgender.master.baseline <- reporthivbyageandgender.master.baseline[ , -which(names(reporthivbyageandgender.master.baseline) %in% c("Newly.Tested.Positive","Newly.Tested.Negative","Tested.Past.Year.or.On_ART", "Tested.Ever","Diagnosed","NodeId","HasIntervention.Traditional_MC.","Non_Program_MMC", "Program_VMMC"))]

#noARTnoVMMC
setwd(paste0(wd2))
files <- list.files(full.names = F)
head(files)
f <- paste0(wd2, files)
length(f)

for (i in seq(1,250,1)){
  reporthivbyageandgender <- read.csv(f[i])
  reporthivbyageandgender$scenario <- "noARTnoVMMC"
  reporthivbyageandgender$sim.id <- paste0(files[i])
  #reporthivbyageandgender <- reporthivbyageandgender[,c(1,3,7,8,9,10,11,12,22,23)]
  if (i==1){
    reporthivbyageandgender.master <- reporthivbyageandgender
  } else {
    reporthivbyageandgender.master <- rbind(reporthivbyageandgender.master, reporthivbyageandgender)
    print(paste0("Why hello there, I'm working on noARTnoVMMC folder",i))
  }
}

reporthivbyageandgender.master.noartnovmmc <- reporthivbyageandgender.master
names(reporthivbyageandgender.master.noartnovmmc)
reporthivbyageandgender.master.noartnovmmc <- reporthivbyageandgender.master.noartnovmmc[ , -which(names(reporthivbyageandgender.master.noartnovmmc) %in% c("Newly.Tested.Positive","Newly.Tested.Negative","Tested.Past.Year.or.On_ART", "Tested.Ever","Diagnosed","NodeId","HasIntervention.Traditional_MC.","Non_Program_MMC", "Program_VMMC"))]

#noART
setwd(paste0(wd3))
files <- list.files(full.names = F)
head(files)
f <- paste0(wd3, files)
length(f)

for (i in seq(1,250,1)){
  reporthivbyageandgender <- read.csv(f[i])
  reporthivbyageandgender$scenario <- "noART"
  reporthivbyageandgender$sim.id <- paste0(files[i])
  #reporthivbyageandgender <- reporthivbyageandgender[,c(1,3,7,8,9,10,11,12,22,23)]
  if (i==1){
    reporthivbyageandgender.master <- reporthivbyageandgender
  } else {
    reporthivbyageandgender.master <- rbind(reporthivbyageandgender.master, reporthivbyageandgender)
    print(paste0("Why hello there, I'm working on noART folder",i))
  }
}

reporthivbyageandgender.master.noart <- reporthivbyageandgender.master
names(reporthivbyageandgender.master.noart)
reporthivbyageandgender.master.noart <- reporthivbyageandgender.master.noart[ , -which(names(reporthivbyageandgender.master.noart) %in% c("Newly.Tested.Positive","Newly.Tested.Negative","Tested.Past.Year.or.On_ART", "Tested.Ever","Diagnosed","NodeId","HasIntervention.Traditional_MC.","Non_Program_MMC", "Program_VMMC"))]


#909090
setwd(paste0(wd4))
files <- list.files(full.names = F)
head(files)
f <- paste0(wd4, files)
length(f)

for (i in seq(1,250,1)){
  reporthivbyageandgender <- read.csv(f[i])
  reporthivbyageandgender$scenario <- "909090"
  reporthivbyageandgender$sim.id <- paste0(files[i])
  if (i==1){
    reporthivbyageandgender.master <- reporthivbyageandgender
  } else {
    reporthivbyageandgender.master <- rbind(reporthivbyageandgender.master, reporthivbyageandgender)
    print(paste0("Why hello there, I'm working on 909090agetarg folder",i))
  }
}

reporthivbyageandgender.master.909090agetarg <- reporthivbyageandgender.master
names(reporthivbyageandgender.master.909090agetarg)
reporthivbyageandgender.master.909090agetarg <- reporthivbyageandgender.master.909090agetarg[ , -which(names(reporthivbyageandgender.master.909090agetarg) %in% c("Newly.Tested.Positive","Newly.Tested.Negative","Tested.Past.Year.or.On_ART", "Tested.Ever","Diagnosed","NodeId","HasIntervention.Traditional_MC.","Non_Program_MMC", "Program_VMMC"))]

#ART100pct
setwd(paste0(wd6))
files <- list.files(full.names = F)
head(files)
f <- paste0(wd6, files)
length(f)

for (i in seq(1,250,1)){
  reporthivbyageandgender <- read.csv(f[i])
  reporthivbyageandgender$scenario <- "ART100pct"
  reporthivbyageandgender$sim.id <- paste0(files[i])
  #reporthivbyageandgender <- reporthivbyageandgender[,c(1,3,7,8,9,10,11,12,22,23)]
  if (i==1){
    reporthivbyageandgender.master <- reporthivbyageandgender
  } else {
    reporthivbyageandgender.master <- rbind(reporthivbyageandgender.master, reporthivbyageandgender)
    print(paste0("Why hello there, I'm working on ART100pct folder",i))
  }
}

reporthivbyageandgender.master.ART100pct <- reporthivbyageandgender.master
names(reporthivbyageandgender.master.ART100pct)
reporthivbyageandgender.master.ART100pct <- reporthivbyageandgender.master.ART100pct[ , -which(names(reporthivbyageandgender.master.ART100pct) %in% c("Newly.Tested.Positive","Newly.Tested.Negative","Tested.Past.Year.or.On_ART", "Tested.Ever","Diagnosed","NodeId","HasIntervention.Traditional_MC.","Non_Program_MMC", "Program_VMMC"))]

#Stitch together reporthivbyageandgender files
reporthivbyageandgender.master.final <- rbind(reporthivbyageandgender.master.baseline,
                                              reporthivbyageandgender.master.noartnovmmc,
                                              reporthivbyageandgender.master.noart,
                                              reporthivbyageandgender.master.909090agetarg,
                                              reporthivbyageandgender.master.ART100pct)


table(reporthivbyageandgender.master.final$scenario)

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\HIV\\EMOD\\Swaziland_05June2018\\Scenarios_calibration_run3_v05Sep2018_newcampaign_FINAL\\ScenariosOct_FINAL")
write.csv(reporthivbyageandgender.master.final, "reporthivbyageandgender.master.final.csv")

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\HIV\\EMOD\\Swaziland_05June2018\\Scenarios_calibration_run3_v05Sep2018_newcampaign_FINAL\\ScenariosOct_FINAL")
reporthivbyageandgender.master.final <- read.csv("reporthivbyageandgender.master.final.csv")

###########################################
#Bring in risk scenarios
###########################################

# #909090
# setwd(paste0(wd4))
# files <- list.files(full.names = F)
# head(files)
# f <- paste0(wd4, files)
# length(f)
# 
# for (i in seq(1,250,1)){
#   reporthivbyageandgender <- read.csv(f[i])
#   reporthivbyageandgender$scenario <- "909090"
#   reporthivbyageandgender$sim.id <- paste0(files[i])
#   reporthivbyageandgender <- reporthivbyageandgender[,c(1,3,7,8,9,10,11,23,24)]
#   if (i==1){
#     reporthivbyageandgender.master <- reporthivbyageandgender
#   } else {
#     reporthivbyageandgender.master <- rbind(reporthivbyageandgender.master, reporthivbyageandgender)
#     print(paste0("Why hello there, I'm working on 909090 folder",i))
#   }
# }
# 
# reporthivbyageandgender.master.909090 <- reporthivbyageandgender.master
# names(reporthivbyageandgender.master.909090)
# 
# setwd(paste0(wd10))
# files <- list.files(full.names = F)
# head(files)
# f <- paste0(wd10, files)
# length(f)
# 
# for (i in seq(1,250,1)){
#   reporthivbyageandgender <- read.csv(f[i])
#   reporthivbyageandgender$scenario <- "highrisk100"
#   reporthivbyageandgender$sim.id <- paste0(files[i])
#   reporthivbyageandgender <- reporthivbyageandgender[,c(1,3,7,8,9,10,11,23,24)]
#   if (i==1){
#     reporthivbyageandgender.master <- reporthivbyageandgender
#   } else {
#     reporthivbyageandgender.master <- rbind(reporthivbyageandgender.master, reporthivbyageandgender)
#     print(paste0("Why hello there, I'm working on highrisk100 folder",i))
#   }
# }
# 
# reporthivbyageandgender.master.highrisk100 <- reporthivbyageandgender.master
# 
# setwd(paste0(wd11))
# files <- list.files(full.names = F)
# head(files)
# f <- paste0(wd11, files)
# length(f)
# 
# for (i in seq(1,250,1)){
#   reporthivbyageandgender <- read.csv(f[i])
#   reporthivbyageandgender$scenario <- "lowrisk100"
#   reporthivbyageandgender$sim.id <- paste0(files[i])
#   reporthivbyageandgender <- reporthivbyageandgender[,c(1,3,7,8,9,10,11,23,24)]
#   if (i==1){
#     reporthivbyageandgender.master <- reporthivbyageandgender
#   } else {
#     reporthivbyageandgender.master <- rbind(reporthivbyageandgender.master, reporthivbyageandgender)
#     print(paste0("Why hello there, I'm working on lowrisk100 folder",i))
#   }
# }
# 
# reporthivbyageandgender.master.lowrisk100 <- reporthivbyageandgender.master
# 
# reporthivbyageandgender.master.final.risk <- rbind(reporthivbyageandgender.master.909090,
#                                                    reporthivbyageandgender.master.highrisk100,
#                                                    reporthivbyageandgender.master.lowrisk100)
# 
# table(reporthivbyageandgender.master.final.risk$scenario)
# 
# setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\HIV\\EMOD\\Swaziland_05June2018\\Scenarios_calibration_run3_v05Sep2018_newcampaign_FINAL\\ScenariosOct_FINAL")
# write.csv(reporthivbyageandgender.master.final.risk, "reporthivbyageandgender.master.final.risk.csv")
# 
# setwd(paste0(wd12))
# files <- list.files(full.names = F)
# head(files)
# f <- paste0(wd12, files)
# length(f)
# for (i in seq(1,250,1)){
#   reporthivbyageandgender <- read.csv(f[i])
#   reporthivbyageandgender <- subset(reporthivbyageandgender, Age >10 & Age < 50)
#   reporthivbyageandgender$sim.id <- paste0(files[i])
#   reporthivbyageandgender <- reporthivbyageandgender[,c(1,3,4,8,9,20)]
#   if (i==1){
#     reporthivbyageandgender.master <- reporthivbyageandgender
#   } else {
#     reporthivbyageandgender.master <- rbind(reporthivbyageandgender.master, reporthivbyageandgender)
#     print(paste0("Why hello there, I'm working on baselineriskgroup folder",i))
#   }
# }
# reporthivbyageandgender.master.baselineriskgroup <- reporthivbyageandgender.master 
# head(reporthivbyageandgender.master.baselineriskgroup)

###########################################
#create condensed file
###########################################

# setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\HIV\\EMOD\\Swaziland_05June2018\\Scenarios_calibration_run3_v05Sep2018_newcampaign_FINAL\\ScenariosOct_FINAL")
# 
# #bring in scale-up scenarios
# reporthivbyageandgender.master.final <- read.csv("reporthivbyageandgender.master.final.csv")
# reporthivbyageandgender.master <- subset(reporthivbyageandgender.master.final, Age>14 & Age < 50)
# reporthivbyageandgender.master$Agecat<-cut(reporthivbyageandgender.master$Age, c(15,25,35,50), right=F)
# reporthivbyageandgender.master <- aggregate(cbind(Population,Infected,Newly.Infected, On_ART, Died_from_HIV) ~ Year+Agecat+Gender+sim.id+scenario, reporthivbyageandgender.master, FUN=sum)
# names(reporthivbyageandgender.master)
# table(reporthivbyageandgender.master$scenario)
# reporthivbyageandgender.master <- subset(reporthivbyageandgender.master, scenario!="hist909090")
# 
# write.csv(reporthivbyageandgender.master, "reporthivbyageandgender.master.final.condensed.csv")
# 

###########################################
#Calculate incidence 
###########################################

#Calculate Incidence overall
reporthivbyageandgender.master.final$Year2 <- floor((reporthivbyageandgender.master.final$Year-0.5))
reporthivbyageandgender.master.final$Uninfected.Population = reporthivbyageandgender.master.final$Population-reporthivbyageandgender.master.final$Infected
trajectories_IR.1a <- aggregate(Newly.Infected ~ Year2+sim.id+scenario, subset(reporthivbyageandgender.master.final, Age>10 & Age<50), FUN=sum) #sums number of new infections in each year
trajectories_IR.2 <- aggregate(Uninfected.Population ~ Year+sim.id+scenario, subset(reporthivbyageandgender.master.final, Age>10 & Age<50), FUN=sum)
trajectories_IR.2$Year2 <- floor(trajectories_IR.2$Year-0.5)
trajectories_IR.2 <- trajectories_IR.2[!duplicated(trajectories_IR.2[c("Year2","sim.id","scenario")]),] #remove second instance of duplicate rows
trajectories_IR.2 <- trajectories_IR.2[-match("Year",names(trajectories_IR.2))]
trajectories_IRoverall <- merge(trajectories_IR.1a, trajectories_IR.2, by=c("Year2","sim.id","scenario"))
trajectories_IRoverall$incidence <- trajectories_IRoverall$Newly.Infected / trajectories_IRoverall$Uninfected.Population
summary(trajectories_IRoverall$incidence[trajectories_IRoverall$scenario=="baseline"])

#Calculate incidence by gender and append to overall
trajectories_IR.1a <- aggregate(Newly.Infected ~ Year2+Gender+sim.id+scenario, subset(reporthivbyageandgender.master.final, Age>10 & Age<50), FUN=sum) #sums number of new infections in each year
trajectories_IR.2 <- aggregate(Uninfected.Population ~ Year+Gender+sim.id+scenario, subset(reporthivbyageandgender.master.final, Age>10 & Age<50), FUN=sum)
trajectories_IR.2$Year2 <- floor(trajectories_IR.2$Year-0.5)
trajectories_IR.2 <- trajectories_IR.2[!duplicated(trajectories_IR.2[c("Year2","Gender","sim.id","scenario")]),] #remove second instance of duplicate rows
trajectories_IR.2 <- trajectories_IR.2[-match("Year",names(trajectories_IR.2))]
trajectories_IR <- merge(trajectories_IR.1a, trajectories_IR.2, by=c("Year2","Gender","sim.id","scenario"))
trajectories_IR$incidence <- trajectories_IR$Newly.Infected / trajectories_IR$Uninfected.Population
summary(trajectories_IR$incidence[trajectories_IR$scenario=="baseline" & trajectories_IR$Gender==0])
summary(trajectories_IR$incidence[trajectories_IR$scenario=="baseline" & trajectories_IR$Gender==1])

#combine gender-specific and overall IR trajectories
trajectories_IRoverall$Gender = 2
names(trajectories_IR)
names(trajectories_IRoverall)
trajectories_IR_comb <- rbind(trajectories_IR, trajectories_IRoverall)
names(trajectories_IR_comb)
#trajectories_IR_comb <- subset(trajectories_IR_comb, Gender < 2)
summary(trajectories_IR_comb$incidence[trajectories_IRoverall$scenario=="baseline" & trajectories_IR_comb$Gender==0])
summary(trajectories_IR_comb$incidence[trajectories_IRoverall$scenario=="baseline" & trajectories_IR_comb$Gender==1])

#get smoothed values
table(trajectories_IRoverall$scenario)
year=2050
baseline.smooth <- subset(trajectories_IRoverall, Year2 < year & scenario=="noARTnoVMMC")
smooth_vals = as.data.frame(cbind("incidence"= predict(loess(incidence~Year2,baseline.smooth, span=0.2)), "year" = baseline.smooth$Year2))
smooth_vals <- smooth_vals[!duplicated(smooth_vals[c("year","incidence")]),] #remove second instance of duplicate rows

year=2050
baseline.smooth.m <- subset(trajectories_IR_comb, Year2 < year & scenario=="baseline" & Gender==0)
baseline.smooth.f <- subset(trajectories_IR_comb, Year2 < year & scenario=="baseline" & Gender==1)
smooth_vals.m = as.data.frame(cbind("incidence"= predict(loess(incidence~Year2,baseline.smooth.m, span=0.2)), "year" = baseline.smooth.m$Year2))
smooth_vals.m <- smooth_vals.m[!duplicated(smooth_vals.m[c("year","incidence")]),] #remove second instance of duplicate rows
smooth_vals.f = as.data.frame(cbind("incidence"= predict(loess(incidence~Year2,baseline.smooth.f, span=0.2)), "year" = baseline.smooth.f$Year2))
smooth_vals.f <- smooth_vals.f[!duplicated(smooth_vals.f[c("year","incidence")]),] #remove second instance of duplicate rows

#Calibration plot
year=2050
swaziland.inc <- ggplot(data=subset(trajectories_IRoverall, Year2 < year & scenario=="baseline")) +
  geom_point(data=subset(trajectories_IRoverall,  Year2 < year & scenario=="baseline"), size=1.2, color = "grey", aes(x=Year2, y=incidence*100))+
  geom_line(data=subset(trajectories_IRoverall,  Year2 < year & scenario=="baseline"), color="grey",aes(x=Year2, y=incidence*100, group=sim.id)) +
  geom_smooth(aes(x=Year2, y=incidence*100),method="loess", span=0.1, se = T, size=1, color="blue", linetype=1) +
  geom_point(data = subset(incidence.data, Gender==2 & Year==2011), size=2, color = "black", aes(x=Year, y=Incidence))+
  geom_point(data = subset(incidence.data, Gender==2 & Year==2016), size=2, color = "black", aes(x=Year, y=Incidence))+
  geom_errorbar(data = subset(incidence.data, Gender==2 & Year==2011), aes(x=Year, ymin=lb, ymax=ub), color="black", width=2, size=1) +
  geom_errorbar(data = subset(incidence.data, Gender==2 & Year==2016), aes(x=Year, ymin=lb, ymax=ub), color="black", width=2, size=1) +
  xlab("Year")+
  ylab("Incidence (per 100 py)")+
  theme_bw(base_size=16) +
  guides(fill = guide_legend(keywidth = 2, keyheight = 1)) +
  scale_y_continuous(breaks = seq(0,5,1),limits=c(0,5),expand = c(0,0)) +
  scale_x_continuous(breaks = seq(1975,year,5),limits=c(1975,year+1),expand = c(0,0)) +
  theme(legend.position="bottom") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(strip.background = element_rect(colour="black", fill="white")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"))
swaziland.inc

#get smoothed values
year=2050
baseline.smooth <- subset(trajectories_IRoverall, Year2 < year & scenario=="baseline")
smooth_vals = as.data.frame(cbind("incidence"= predict(loess(incidence~Year2,baseline.smooth, span=0.2)), "year" = baseline.smooth$Year2))
smooth_vals <- smooth_vals[!duplicated(smooth_vals[c("year","incidence")]),] #remove second instance of duplicate rows

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\BMGF Lancet Paper\\Final_Plots")
ggsave("swaziland.incidence.calib_FINAL_gatesreview.jpg", height=8, width=8)

year = 2051
span = 0.2
table(trajectories_IRoverall$scenario)
trajectories_IRoverall$scenario_f <- ordered(trajectories_IRoverall$scenario, levels=c("ART100pct","noARTnoVMMC","noART","909090","baseline"))
swaziland.inc <- ggplot() +
  geom_smooth(data=subset(trajectories_IRoverall, Year2 < year), aes(x=Year2, y=incidence*100, color=scenario_f),method="loess", span=span, se = F, size=2) +
  xlab("Year")+
  ylab("Incidence (per 100 py)")+
  theme_bw(base_size=16) +
  scale_color_manual(labels = c("100% ART initiation","No interventions", "VMMC only","Age-targeted 90-90-90", "Baseline"), values = c("green","red", "orange", "purple", "blue")) +
  guides(fill = guide_legend(keywidth = 2, keyheight = 1)) +
  geom_hline(yintercept=0.1, linetype=2) +
  scale_y_continuous(breaks = seq(0,6,0.5),limits=c(0,6), expand = c(0,0)) +
  scale_x_continuous(breaks = seq(1975,year,5), limits=c(1975, year),expand = c(0,0)) +
  coord_cartesian(ylim=c(0, 5)) +
  theme(legend.position="bottom") +
  theme(legend.title=element_blank())+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(strip.background = element_rect(colour="black", fill="white")) +
  guides(col = guide_legend(nrow=3)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"))
swaziland.inc

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\BMGF Lancet Paper\\Submission\\Lancet HIV\\Resubmission2\\Figures")
ggsave("swaziland.incidence.5scenarios_FINAL.pdf", height=8, width=8)

###########################################
#Calculate incidence by gender
###########################################
table(trajectories_IR_comb$Gender)
labs <- c("1" = "Women", "0" = "Men")
#Calibration plot
swaziland.inc <- ggplot(data=subset(trajectories_IRoverall,  Year2 < year & scenario=="baseline")) +
  geom_line(data=subset(trajectories_IRoverall,  Year2 < year & scenario=="baseline"), color="grey",aes(x=Year2, y=incidence*100, group=sim.id)) +
  facet_grid(~Gender, labeller=labeller(Gender = labs)) +
  stat_smooth(aes(x=Year2, y=incidence*100),method="loess", span=0.2, se = T, size=1, color="blue") +
  geom_point(data = subset(incidence.data, Gender < 2), size=2, color = "black", aes(x=Year, y=Incidence))+
  geom_errorbar(data = subset(incidence.data, Gender < 2), aes(x=Year, ymin=lb, ymax=ub), color="black", width=2, size=1) +
  xlab("Year")+
  ylab("Incidence (per 100 py")+
  theme_bw(base_size=16) +
  guides(fill = guide_legend(keywidth = 2, keyheight = 1)) +
  scale_y_continuous(breaks = seq(0,6,1),expand = c(0,0)) +
  coord_cartesian(ylim=c(0, 6)) +
  theme(legend.position="bottom") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(strip.background = element_rect(colour="black", fill="white")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"))
swaziland.inc

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\BMGF Lancet Paper\\Submission\\Lancet HIV\\Resubmission2\\Figures")
ggsave("swaziland.incidence.bygender.calib_FINAL.jpg", height=8, width=12)

labs.3 <- c("2"="Combined", "1" = "Women", "0" = "Men")
table(trajectories_IR_comb$Gender)

swaziland.inc <- ggplot(data=subset(trajectories_IR_comb,  Year2 < year & scenario=="baseline")) +
  geom_line(data=subset(trajectories_IR_comb,  Year2 < year & scenario=="baseline"), color="grey",aes(x=Year2, y=incidence*100, group=sim.id)) +
  facet_grid(~Gender, labeller=labeller(Gender = labs.3)) +
  stat_smooth(aes(x=Year2, y=incidence*100),method="loess", span=0.2, se = T, size=1, color="blue") +
  geom_point(data = subset(incidence.data), size=2, color = "black", aes(x=Year, y=Incidence))+
  geom_errorbar(data = subset(incidence.data), aes(x=Year, ymin=lb, ymax=ub), color="black", width=2, size=1) +
  xlab("Year")+
  ylab("Incidence (per 100 py)")+
  theme_bw(base_size=16) +
  guides(fill = guide_legend(keywidth = 2, keyheight = 1)) +
  scale_y_continuous(breaks = seq(0,6,1),expand = c(0,0)) +
  theme(legend.position="bottom") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(strip.background = element_rect(colour="black", fill="white")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"))
swaziland.inc

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\BMGF Lancet Paper\\Submission\\Lancet HIV\\Resubmission2\\Figures")
ggsave("swaziland.incidence.bygender_andoverall.calib_FINAL.pdf", height=5, width=12)

###########################################
#Incidence scenarios by gender
###########################################

#Incidenc by (4) scenarios
year = 2051
span = 0.2
table(trajectories_IR_comb$scenario)
trajectories_IR_comb$scenario_f <- ordered(trajectories_IR_comb$scenario, levels=c("ART100pct","noARTnoVMMC","noART","909090","baseline"))
swaziland.inc <- ggplot() +
  geom_smooth(data=subset(trajectories_IR_comb, Year2 < year), aes(x=Year2, y=incidence*100, color=scenario_f),method="loess", span=span, se = F, size=2) +
  xlab("Year")+
  ylab("Incidence (per 100 py)")+
  facet_grid(~Gender, labeller=labeller(Gender = labs)) +
  theme_bw(base_size=24) +
  scale_color_manual(labels = c("100% ART initiation","No interventions", "VMMC only","Age-targeted 90-90-90", "Baseline"), values = c("green","red", "orange", "purple", "blue")) +
  guides(fill = guide_legend(keywidth = 2, keyheight = 1)) +
  scale_y_continuous(breaks = seq(0,6,0.5),limits=c(0,5), expand = c(0,0)) +
  scale_x_continuous(breaks = seq(1975,year,5), limits=c(1975, year),expand = c(0,0)) +
  theme(legend.position="bottom") +
  theme(legend.title=element_blank())+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(strip.background = element_rect(colour="black", fill="white")) +
  guides(col = guide_legend(nrow=3)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"))
swaziland.inc

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\BMGF Lancet Paper\\Submission\\Lancet HIV\\Resubmission2\\Figures")
ggsave("swaziland.incidence.bygender.4scenarios.pdf", height=8, width=12)

###########################################
#Incidence scenarios by age and gender
###########################################
reporthivbyageandgender.master$Year2 <- floor((reporthivbyageandgender.master$Year-0.5))
reporthivbyageandgender.master$Uninfected.Population <- reporthivbyageandgender.master$Population - reporthivbyageandgender.master$Infected
trajectories_IR.1a <- aggregate(Newly.Infected ~ Year2+Gender+Agecat+sim.id+scenario, subset(reporthivbyageandgender.master), FUN=sum) #sums number of new infections in each year
trajectories_IR.2 <- aggregate(Uninfected.Population ~ Year+Gender+Agecat+sim.id+scenario, subset(reporthivbyageandgender.master), FUN=sum)
trajectories_IR.2$Year2 <- floor(trajectories_IR.2$Year-0.5)
trajectories_IR.2 <- trajectories_IR.2[!duplicated(trajectories_IR.2[c("Year2","Gender","Agecat", "sim.id","scenario")]),] #remove second instance of duplicate rows
trajectories_IR.2 <- trajectories_IR.2[-match("Year",names(trajectories_IR.2))]
trajectories_IR <- merge(trajectories_IR.1a, trajectories_IR.2, by=c("Year2","Gender","Agecat","sim.id","scenario"))
trajectories_IR$incidence <- trajectories_IR$Newly.Infected / trajectories_IR$Uninfected.Population

#Incidence by scenario
year1=2015 
year2=2050
table(trajectories_IR$Agecat)
labs.age <- c("[15,25)" = "15-24", "[25,35)" = "25-34", "[35,50)" = "35-49")
year = 2051
span = 0.2
trajectories_IR$scenario_f <- ordered(trajectories_IR$scenario, levels=c("ART100pct","noARTnoVMMC","noART","909090","baseline"))
swaziland.inc <- ggplot() +
  geom_smooth(data=subset(trajectories_IR, Year2 < year & scenario_f != "noARTnoVMMC" & scenario_f!="noART"), aes(x=Year2, y=incidence*100, color=scenario_f),method="loess", span=span, se = F, size=2) +
  xlab("Year")+
  ylab("Incidence (per 100 py)")+
  facet_grid(Gender~Agecat, labeller=labeller(Gender = labs, Agecat = labs.age)) +
  theme_bw(base_size=24) +
  scale_color_manual(labels = c("100% ART initiation","Age-targeted 90-90-90", "Status-quo"), values = c("green", "purple", "blue")) +
  guides(fill = guide_legend(keywidth = 2, keyheight = 1)) +
  scale_y_continuous(breaks = seq(0,2.5,0.5),limits=c(0,2.6), expand = c(0,0)) +
  scale_x_continuous(breaks = seq(year1,year2,5), limits=c(year1, year2+1),expand = c(0,0)) +
  theme(legend.position="bottom") +
  theme(legend.title=element_blank())+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(strip.background = element_rect(colour="black", fill="white")) +
  guides(col = guide_legend(nrow=1)) +
  #geom_hline(yintercept=0.1, linetype=2) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"))
swaziland.inc

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\BMGF Lancet Paper\\Submission\\Lancet HIV\\Resubmission2\\Figures")
ggsave("swaziland.incidence.byageandgender.2scenarios.pdf", height=12, width=12)

#############################################
#Cumulative Infections & Number of infections averted by scenario
#############################################
reporthivbyageandgender.master$Year2 <- floor((reporthivbyageandgender.master$Year-0.5))
reporthivbyageandgender.master$Uninfected.Population <- reporthivbyageandgender.master$Population - reporthivbyageandgender.master$Infected
trajectories_IR.1a <- aggregate(Newly.Infected ~ Year2+sim.id+scenario, subset(reporthivbyageandgender.master), FUN=sum) #sums number of new infections in each year
trajectories_IR.2 <- aggregate(Uninfected.Population ~ Year+sim.id+scenario, subset(reporthivbyageandgender.master), FUN=sum)
trajectories_IR.2$Year2 <- floor(trajectories_IR.2$Year-0.5)
trajectories_IR.2 <- trajectories_IR.2[!duplicated(trajectories_IR.2[c("Year2","sim.id","scenario")]),] #remove second instance of duplicate rows
trajectories_IR.2 <- trajectories_IR.2[-match("Year",names(trajectories_IR.2))]
trajectories_IR <- merge(trajectories_IR.1a, trajectories_IR.2, by=c("Year2","sim.id","scenario"))
trajectories_IR$incidence <- trajectories_IR$Newly.Infected / trajectories_IR$Uninfected.Population

head(trajectories_IR,100)

#Add up the cumulative number of infections by year for each scenario.
#Combine scenarios by row (year, sim.id)
#Calculate the diference in cumulative infections between scenarios as infections averted (or additional infections averted)
#Number of infections averted
head(reporthivbyageandgender.master)
population.scale.factor = 681472  #15-49 year olds from population calibration ingest form
ref.year= 2014
reporthivbyageandgender.master.baseline <- subset(reporthivbyageandgender.master, scenario=="baseline")
pop <- aggregate(Population ~ Year+sim.id, reporthivbyageandgender.master.baseline, FUN=sum) 
head(pop,100)
sim.pop.2014 <- mean(pop$Population[pop$Year==2014])
scale.factor <- population.scale.factor / sim.pop.2014

#calculate cumulative incidence by year (by two strata)
trajectories_IR <- trajectories_IR[order(trajectories_IR$scenario, 
                                         trajectories_IR$sim.id, 
                                         trajectories_IR$Year2),]
trajectories_IR$csum <- ave(trajectories_IR$Newly.Infected, trajectories_IR$scenario, trajectories_IR$sim.id, FUN=cumsum)
head(trajectories_IR,100)

#convert long to wide
w <- reshape(trajectories_IR, 
             timevar = "scenario",
             idvar = c("sim.id", "Year2"),
             direction = "wide")
head(w)

#subtract cumulative infections from noARTnoVMMC to calculate infections averted
w$inf.av.909090 <- w$csum.noARTnoVMMC - w$csum.909090
w$inf.av.baseline <- w$csum.noARTnoVMMC - w$csum.baseline
w$inf.av.100pctART <- w$csum.noARTnoVMMC - w$csum.ART100pct
w$inf.av.noART <- w$csum.noARTnoVMMC - w$csum.noART
w$inf.av.noARTnoVMMC <- 0

#get smoothed values
year=2050
head(w)
smooth_vals = as.data.frame(cbind("infav"= predict(loess(inf.av.909090~Year2,w, span=0.2)), "year" = w$Year2))
smooth_vals <- smooth_vals[!duplicated(smooth_vals[c("year","infav")]),] #remove second instance of duplicate rows

year=2050
scaled = 1000
swaziland.inc <- ggplot(data=subset(w,  Year2 < year)) +
  stat_smooth(aes(x=Year2, y=(inf.av.100pctART*scale.factor)/scaled),method="loess", span=0.2, se = F, size=1, color="green") +
  stat_smooth(aes(x=Year2, y=(inf.av.baseline*scale.factor)/scaled),method="loess", span=0.2, se = F, size=1, color="blue") +
  stat_smooth(aes(x=Year2, y=(inf.av.909090*scale.factor)/scaled),method="loess", span=0.2, se = F, size=1, color="purple") +
  stat_smooth(aes(x=Year2, y=(inf.av.noART*scale.factor)/scaled),method="loess", span=0.2, se = F, size=1, color="orange") +
  xlab("Year")+
  ylab("Cumulative infections averted (1,000's)")+
  theme_bw(base_size=16) +
  guides(fill = guide_legend(keywidth = 2, keyheight = 1)) +
  scale_y_continuous(breaks = seq(0,700,50), expand = c(0,0)) +
  scale_x_continuous(breaks = seq(2005,year,10), limits=c(2000,year), expand = c(0,0)) +
  theme(legend.position="bottom") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(strip.background = element_rect(colour="black", fill="white")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"))
swaziland.inc

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\BMGF Lancet Paper\\Submission\\Lancet HIV\\Resubmission2\\Figures")
ggsave("swaziland.infections.averted.scenarios.pdf", height=8, width=8)

table(reporthivbyageandgender.master$scenario)
reporthivbyageandgender.master$Year2 <- floor((reporthivbyageandgender.master$Year-0.5))
reporthivbyageandgender.master$Uninfected.Population <- reporthivbyageandgender.master$Population - reporthivbyageandgender.master$Infected
trajectories_IR.1a <- aggregate(Newly.Infected ~ Year2+Gender+Agecat+sim.id+scenario, subset(reporthivbyageandgender.master), FUN=sum) #sums number of new infections in each year
trajectories_IR.2 <- aggregate(Uninfected.Population ~ Year+Gender+Agecat+sim.id+scenario, subset(reporthivbyageandgender.master), FUN=sum)
trajectories_IR.2$Year2 <- floor(trajectories_IR.2$Year-0.5)
trajectories_IR.2 <- trajectories_IR.2[!duplicated(trajectories_IR.2[c("Year2","Gender","Agecat", "sim.id","scenario")]),] #remove second instance of duplicate rows
trajectories_IR.2 <- trajectories_IR.2[-match("Year",names(trajectories_IR.2))]
trajectories_IR <- merge(trajectories_IR.1a, trajectories_IR.2, by=c("Year2","Gender","Agecat","sim.id","scenario"))
trajectories_IR$incidence <- trajectories_IR$Newly.Infected / trajectories_IR$Uninfected.Population
head(trajectories_IR,10)

#Add up the cumulative number of infections by year for each scenario.
#Combine scenarios by row (year, sim.id)
#Calculate the diference in cumulative infections between scenarios as infections averted (or additional infections averted)
#Number of infections averted

population.scale.factor = 681472
ref.year= 2014
head(reporthivbyageandgender.master)
reporthivbyageandgender.master.baseline <- subset(reporthivbyageandgender.master, scenario=="baseline")
table(reporthivbyageandgender.master.baseline$scenario)
pop <- aggregate(Population ~ Year+sim.id, reporthivbyageandgender.master.baseline, FUN=sum) 
head(pop)
sim.pop.2014 <- mean(pop$Population[pop$Year==2014])
scale.factor <- population.scale.factor / sim.pop.2014

#calculate cumulative incidence by year (by two strata)
trajectories_IR <- trajectories_IR[order(trajectories_IR$scenario, trajectories_IR$sim.id, trajectories_IR$Gender, 
                                         trajectories_IR$Agecat, trajectories_IR$Year2),]

trajectories_IR$csum <- ave(trajectories_IR$Newly.Infected, trajectories_IR$scenario, 
                            trajectories_IR$sim.id, trajectories_IR$Gender, 
                            trajectories_IR$Agecat, FUN=cumsum)

#convert long to wide
w <- reshape(trajectories_IR, 
             timevar = "scenario",
             idvar = c("sim.id","Gender", "Agecat", "Year2"),
             direction = "wide")
head(w)

#subtract cumulative infections from noARTnoVMMC to calculate infections averted
w$inf.av.909090 <- w$csum.baseline - w$csum.909090
w$inf.av.baseline <- w$csum.baseline - w$csum.baseline
w$inf.av.100pctART <- w$csum.baseline - w$csum.ART100pct
summary(w$inf.av.909090)

table(w$Gender)
table(w$Agecat)

labs.age <- c("[15,25)" = "15-24", "[25,35)" = "25-34", "[35,50)" = "35-49")
year=2050
swaziland.inc <- ggplot(data=subset(w,  Year2 < year)) +
  # geom_line(data=subset(trajectories_IRoverall,  Year2 < year & scenario=="baseline"), 
  #           color="blue",aes(x=Year2, y=csum, group=sim.id)) +
  # geom_line(data=subset(trajectories_IRoverall,  Year2 < year & scenario=="noARTnoVMMC"), 
  #           color="red",aes(x=Year2, y=csum, group=sim.id)) +
  facet_grid(Gender~Agecat, labeller=labeller(Gender = labs, Agecat = labs.age)) +
  stat_smooth(aes(x=Year2, y=(inf.av.100pctART*scale.factor)/1000),method="loess", span=0.2, se = F, size=2, color="green") +
  stat_smooth(aes(x=Year2, y=(inf.av.baseline*scale.factor)/1000),method="loess", span=0.2, se = F, size=2, color="blue") +
  stat_smooth(aes(x=Year2, y=(inf.av.909090*scale.factor)/1000),method="loess", span=0.2, se = F, size=2, color="purple") +
  xlab("Year")+
  ylab("Cumulative infections averted (1000's)")+
  theme_bw(base_size=24) +
  guides(fill = guide_legend(keywidth = 2, keyheight = 1)) +
  scale_y_continuous(breaks = seq(-5,25,5),limits=c(-5,30), expand = c(0,0)) +
  scale_x_continuous(breaks = seq(2015,year,5), limits=c(2015,year),expand = c(0,0)) +
  theme(legend.position="bottom") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(strip.background = element_rect(colour="black", fill="white")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"))
swaziland.inc

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\BMGF Lancet Paper\\Final_Plots")
ggsave("swaziland.infections.averted.scenarios.byageandgender.jpg", height=12, width=12)

#get smoothed values
year=2050
head(w)
w.f.1524 <- subset(w, Gender==1 & Agecat=="[15,25)")
smooth_vals = as.data.frame(cbind("infav"= predict(loess(inf.av.909090~Year2,w.f.1524, span=0.2)), "year" = w.f.1524$Year2))
smooth_vals <- smooth_vals[!duplicated(smooth_vals[c("year","infav")]),] #remove second instance of duplicate rows

#############################################
#Prevalence overall 15-49
#############################################
#Prevalence 15-49 
reporthivbyageandgender.master2 <- subset(reporthivbyageandgender.master.final, scenario=="baseline")
table(reporthivbyageandgender.master2$scenario)
head(reporthivbyageandgender.master2)
trajectories_infections <- aggregate(Infected ~ Year+sim.id, subset(reporthivbyageandgender.master2, Age>10 & Age<50), FUN=sum) #sums prev infections in each year
trajectories_population <- aggregate(Population ~ Year+sim.id, subset(reporthivbyageandgender.master2, Age>10 & Age<50), FUN=sum)
trajectories_prev <- merge(trajectories_infections, trajectories_population, by=c("Year","sim.id"))
trajectories_prev$prevalence <- trajectories_prev$Infected / trajectories_prev$Population
head(prevalence.bothgenders)
head(trajectories_prev)

prevalence.bothgenders.c <- subset(prevalence.bothgenders, is.na(prevalence.bothgenders$Gender==T))
year=2050
swaziland.prev <- ggplot(data=subset(trajectories_prev, Year < year)) +
  geom_line(size=1.2, aes(x=Year, y=prevalence*100, group=sim.id), alpha=0.01,color="grey")+ 
  geom_smooth(aes(x=Year, y=prevalence*100),color="blue", method="loess", span=0.1, se = F, size=2) +
  geom_point(data = prevalence.bothgenders.c, size=2, color = "black", aes(x=Year, y=NationalPrevalence*100)) + 
  geom_errorbar(data = prevalence.bothgenders.c, aes(x=Year, ymin=lb*100, ymax=ub*100), color="black", width=1, size=1) +
  xlab("Year")+
  ylab("Prevalence")+
  theme_bw(base_size=16) +
  guides(fill = guide_legend(keywidth = 2, keyheight = 1)) +
  scale_x_continuous(breaks = seq(1980,2020,5)) +
  theme(legend.position="bottom") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(strip.background = element_rect(colour="black", fill="white"))
swaziland.prev

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\BMGF Lancet Paper\\Final_Plots")
ggsave("swaziland.prevalence.calib.gender_baseline_age20_49.jpg", height=8, width=15)

#############################################
#Prevalence by age and gender
#############################################
#Prevalence 15-49 by gender
reporthivbyageandgender.master2 <- subset(reporthivbyageandgender.master.final, scenario=="baseline")
table(reporthivbyageandgender.master2$scenario)
head(reporthivbyageandgender.master2)
trajectories_infections <- aggregate(Infected ~ Year+Gender+sim.id+scenario, subset(reporthivbyageandgender.master2, Age>15 & Age<50), FUN=sum) #sums prev infections in each year
trajectories_population <- aggregate(Population ~ Year+Gender+sim.id+scenario, subset(reporthivbyageandgender.master2, Age>15 & Age<50), FUN=sum)
trajectories_prev <- merge(trajectories_infections, trajectories_population, by=c("Year","Gender","sim.id","scenario"))
trajectories_prev$prevalence <- trajectories_prev$Infected / trajectories_prev$Population
head(prevalence.bothgenders)
head(trajectories_prev)

trajectories_infections.age <- aggregate(Infected ~ Year+Gender+Age+sim.id+scenario, subset(reporthivbyageandgender.master2, Age>10 & Age<100), FUN=sum) #sums prev infections in each year
trajectories_population.age <- aggregate(Population ~ Year+Gender+Age+sim.id+scenario, subset(reporthivbyageandgender.master2, Age>10 & Age<100), FUN=sum)
trajectories_prev.age <- merge(trajectories_infections.age, trajectories_population.age, by=c("Year","Age","Gender","sim.id","scenario"))
trajectories_prev.age$prevalence <- trajectories_prev.age$Infected / trajectories_prev.age$Population

prevalence.bothgenders.c <- prevalence.bothgenders[c(1:6),]
labs <- c("1" = "Women", "0" = "Men")
year=2050
swaziland.prev <- ggplot(data=subset(trajectories_prev, Year < year)) +
  geom_line(size=1.2, aes(x=Year, y=prevalence*100, group=sim.id), alpha=0.01,color="blue")+ 
  geom_smooth(aes(x=Year, y=prevalence*100),color="blue", method="loess", span=0.1, se = F, size=2) +
  geom_point(data = prevalence.bothgenders.c, size=2, color = "black", aes(x=Year, y=NationalPrevalence*100)) + 
  geom_errorbar(data = prevalence.bothgenders.c, aes(x=Year, ymin=lb*100, ymax=ub*100), color="black", width=1, size=1) +
  geom_vline(xintercept=2004, linetype=2) +
  geom_vline(xintercept=2020, linetype=2) +
  facet_grid(~Gender, labeller=labeller(Gender = labs)) +
  xlab("Year")+
  ylab("Prevalence")+
  theme_bw(base_size=16) +
  guides(fill = guide_legend(keywidth = 2, keyheight = 1)) +
  scale_x_continuous(breaks = seq(1980,year,5)) +
  theme(legend.position="bottom") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(strip.background = element_rect(colour="black", fill="white")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
swaziland.prev

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\BMGF Lancet Paper\\Revised_Plots")
ggsave("swaziland.prevalence.calib.gender_baseline.png", height=8, width=15)

#HIV prev by age and gender
labs <- c("1" = "Women", "0" = "Men")
labs.age <- c("15"="15-19","20"="20-24","25"="25-29","30"="30-34","35"="35-39",
              "40"="40-44","45"="45-49","50"="50-54","55"="55-59","60"="60-64")
swaziland.prev <- ggplot(data=trajectories_prev.age) +
  geom_point(size=1.2, aes(x=Year, y=prevalence*100), alpha=0.01,color="blue")+ 
  geom_smooth(aes(x=Year, y=prevalence*100),color="blue", method="loess", span=0.03, se = F, size=2) +
  geom_point(data = prevalence, size=2, color = "black", aes(x=Year, y=NationalPrevalence*100)) + 
  geom_errorbar(data = prevalence, aes(x=Year, ymin=lb*100, ymax=ub*100), color="black", width=2, size=1) +
  facet_grid(Gender~Age, labeller=labeller(Gender = labs, Age=labs.age)) +
  xlab("Year")+
  ylab("Prevalence")+
  theme_bw(base_size=16) +
  guides(fill = guide_legend(keywidth = 2, keyheight = 1)) +
  #scale_y_continuous(breaks = seq(0,3,0.5),limits=c(0,3)) +
  scale_x_continuous(breaks = seq(1980,2050,10)) +
  theme(legend.position="bottom") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(strip.background = element_rect(colour="black", fill="white"))
swaziland.prev

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\BMGF Lancet Paper\\Final_Plots")
ggsave("swaziland.prevalence.calib.ageandgender_baseline.jpg", height=8, width=15)

#Swaziland prevalence over time by age and gender
labs <- c("1" = "Women", "0" = "Men")
head(prevalence)
head(trajectories_prev.age_year)
trajectories_prev.age_year <- subset(trajectories_prev.age, Year==2016|Year==2011)
table(trajectories_prev.age_year$Year)
swaziland.prev <- ggplot() +
  geom_point(data=subset(prevalence), size=2, aes(x=Age, y=NationalPrevalence*100, color=factor(Year))) + 
  geom_line(data=subset(prevalence), size=2, aes(x=Age, y=NationalPrevalence*100, color=factor(Year))) + 
  geom_errorbar(data=subset(prevalence), aes(x=Age, ymin=lb*100, ymax=ub*100,  color=factor(Year)), width=2, size=1) +
  geom_line(data=trajectories_prev.age_year,aes(x=Age, y=prevalence*100, color=factor(Year), group=sim.id), size=1, alpha=0.05) +
  scale_color_manual(values = c("blue","purple", "red")) +
  facet_grid(~Gender, labeller=labeller(Gender = labs)) +
  xlab("Year")+
  ylab("Prevalence")+
  theme_bw(base_size=20) +
  guides(fill = guide_legend(keywidth = 2, keyheight = 1)) +
  theme(legend.title=element_blank())+
  scale_y_continuous(breaks = seq(0,60,10),limits=c(0,60),expand=c(0,0)) +
  #scale_x_continuous(breaks = seq(15,60,5),expand=c(0,0)) +
  theme(legend.position="bottom") +
  theme(strip.background = element_rect(colour="black", fill="white")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
swaziland.prev

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\BMGF Lancet Paper\\Final_Plots")
ggsave("swaziland.prevalence.calib.ageandgender_baseline_byage.jpg", height=6, width=10)

##############################################
#ART coverage
################################################
names(reporthivbyageandgender.master.final)
table(reporthivbyageandgender.master.final$scenario)
reporthivbyageandgender.master.final$Agecat<-cut(reporthivbyageandgender.master.final$Age, c(15,25,35,45), right=F)
table(reporthivbyageandgender.master.final$Agecat)

trajectory_ART_calib <- aggregate(cbind(On_ART, Infected) ~ Year+Agecat+Gender+sim.id+scenario,
                                  subset(reporthivbyageandgender.master.final, Year > 2000 & Year < 2050), FUN=sum)
trajectory_ART_calib$ART_coverage <- trajectory_ART_calib$On_ART / trajectory_ART_calib$Infected
trajectory_ART_calib.master<-trajectory_ART_calib
head(trajectory_ART_calib.master)
table(trajectory_ART_calib.master$Agecat)

#ART scale-up to 2016
trajectory_ART_calib.master.plot <- subset(trajectory_ART_calib.master,  (scenario=="baseline" | scenario=="909090" | scenario=="ART100pct" ))
head(trajectory_ART_calib.master.plot)
table(trajectory_ART_calib.master.plot$Agecat)
table(trajectory_ART_calib.master.plot$scenario)
year=2049
labs.scenario <- c("baseline"="scenario 1", "909090"="scenario 4", "ART100pct"="scenario 5")
trajectory_ART_calib.master.plot$scenario_f = factor(trajectory_ART_calib.master.plot$scenario, levels=c('baseline','909090','ART100pct','noART','noARTnoVMMC'))
table(trajectory_ART_calib.master.plot$scenario_f)

swaziland.ART.calib <- ggplot() +
  geom_smooth(data=trajectory_ART_calib.master.plot,
              aes(x=Year, y=ART_coverage*100, group = factor(Agecat), color=factor(Agecat)), method="loess",span=0.1, se=F, size=1.2) +
  geom_smooth(data=trajectory_ART_calib.master.plot,
              aes(x=Year, y=ART_coverage*100), color="black",linetype=2, method="loess",span=0.1, se=F, size=1.2) +
  geom_point(data=subset(artdata), aes(x=Year, y=ART_coverage*100, color=factor(Agecat)), size=2) +
  geom_errorbar(data = subset(artdata), aes(x=Year, ymin=lb*100, ymax=ub*100, color=factor(Agecat)), width=1, size=1.2) +
  facet_grid(scenario_f~Gender, labeller=labeller(Gender = labs, scenario_f = labs.scenario)) +
  theme(legend.position="bottom") +
  #geom_hline(aes(yintercept=81)) +
  ylab("ART Coverage (%)")+
  theme_bw(base_size=18) +
  scale_color_manual(labels = c("15-24", "25-34", "35-44"), values = c("blue", "purple", "red")) +
  guides(fill = guide_legend(keywidth = 2, keyheight = 1)) +
  scale_y_continuous(limits=c(0,100), breaks = seq(0,90,10), expand=c(0,0)) +
  scale_x_continuous(limits = c(2000,year), breaks = seq(2000,year+1,5), expand=c(0,0)) +
  theme(legend.position="bottom") +
  theme(legend.direction='vertical') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(strip.background = element_rect(colour="black", fill="white")) +
  guides(colour = guide_legend(nrow = 1, title="")) +
  theme(strip.background = element_rect(colour="black", fill="white")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"))
swaziland.ART.calib

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\BMGF Lancet Paper\\Submission\\Lancet HIV\\Resubmission2\\Figures")
ggsave("swaziland_ART_scaleup_by_scenario.pdf", height=13, width=9)

################################################
#ART by risk scenario
################################################
#reporthivbyageandgender.master <- reporthivbyageandgender.master.final.risk
# names(reporthivbyageandgender.master)
# table(reporthivbyageandgender.master$scenario)
# reporthivbyageandgender.master$Agecat<-cut(reporthivbyageandgender.master$Age, c(15,25,35,45), right=F)
trajectory_ART_calib <- aggregate(cbind(On_ART, Infected) ~ Year+Agecat+Gender+sim.id+scenario,
                                  subset(reporthivbyageandgender.master, Age>14 & Age<50 & Year > 2000 & Year < 2050), FUN=sum)
trajectory_ART_calib$ART_coverage <- trajectory_ART_calib$On_ART / trajectory_ART_calib$Infected
trajectory_ART_calib.master<-trajectory_ART_calib

#ART scale-up to 2016
year=2049
swaziland.ART.calib <- ggplot() +
  geom_smooth(data=trajectory_ART_calib.master,
              aes(x=Year, y=ART_coverage*100, group = factor(Agecat), color=factor(Agecat)), method="loess",span=0.1, se=F, size=1.2) +
  geom_smooth(data=trajectory_ART_calib.master,
              aes(x=Year, y=ART_coverage*100), color="black",linetype=2, method="loess",span=0.1, se=F, size=1.2) +
  geom_point(data=subset(artdata), aes(x=Year, y=ART_coverage*100,
                                       color=factor(Agecat)), size=2) +
  geom_errorbar(data = subset(artdata), aes(x=Year, ymin=lb*100, ymax=ub*100, color=factor(Agecat)), width=1, size=1.2) +
  facet_grid(scenario~Gender, labeller=labeller(Gender = labs)) +
  theme(legend.position="bottom") +
  #geom_hline(aes(yintercept=81)) +
  ylab("ART Coverage %")+
  theme_bw(base_size=20) +
  scale_color_manual(labels = c("15-24", "25-34", "35-49"), values = c("blue", "purple", "red")) +
  guides(fill = guide_legend(keywidth = 2, keyheight = 1)) +
  scale_y_continuous(limits=c(0,100), breaks = seq(0,100,10), expand=c(0,0)) +
  scale_x_continuous(limits = c(2000,year), breaks = seq(2000,year+1,5), expand=c(0,0)) +
  theme(legend.position="bottom") +
  theme(legend.direction='vertical') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(strip.background = element_rect(colour="black", fill="white")) +
  guides(colour = guide_legend(nrow = 1, title=""))
swaziland.ART.calib

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\BMGF Lancet Paper\\Final_Plots")
ggsave("swaziland_ART_scaleup_by_risk_scenario.jpg", height=12, width=9)

################################################
#Numbers on ART by scenario
#################################################
head(trajectory_ART_calib.master,200)
table(trajectory_ART_calib.master$Agecat)
trajectory_ART_calib.master$Year2 <- floor((trajectory_ART_calib.master$Year-0.5))
trajectory_numbersonART <- aggregate(On_ART ~ Year2+Agecat+Gender+sim.id+scenario,
                                     subset(trajectory_ART_calib.master), FUN=sum)
table(trajectory_numbersonART$scenario)
trajectory_numbersonART$scenario_f = factor(trajectory_numbersonART$scenario, levels=c('ART100pct','909090',"baseline"))

table(trajectory_numbersonART$scenario_f)
table(trajectory_numbersonART$Agecat)
head(trajectory_numbersonART)
year=2048
swaziland.ART.calib <- ggplot() +
  geom_smooth(data=subset(trajectory_numbersonART, Agecat=="[15,25)"),
              aes(x=Year2, y=On_ART*scale.factor, group = factor(scenario_f), color=factor(scenario_f)), method="loess",span=0.2, se=F, size=1.2) +
  geom_smooth(data=subset(trajectory_numbersonART, Agecat=="[25,35)"),
              aes(x=Year2, y=On_ART*scale.factor, group = factor(scenario_f), color=factor(scenario_f)), method="loess",span=0.2, se=F, size=1.2, linetype=2) +
  geom_smooth(data=subset(trajectory_numbersonART, Agecat=="[35,50)"),
              aes(x=Year2, y=On_ART*scale.factor, group = factor(scenario_f), color=factor(scenario_f)), method="loess",span=0.2, se=F, size=1.2, linetype=3) +
  facet_grid(~Gender, labeller=labeller(Gender = labs)) +
  theme(legend.position="bottom") +
  #geom_hline(aes(yintercept=81)) +
  xlab("Year")+
  ylab("Number on ART")+
  theme_bw(base_size=20) +
  scale_color_manual(labels = c("ART100pct", "909090", "baseline"), values = c("green", "purple", "blue")) +
  guides(fill = guide_legend(keywidth = 2, keyheight = 1)) +
  #scale_y_continuous(limits=c(0,100), breaks = seq(0,100,10), expand=c(0,0)) +
  scale_x_continuous(limits = c(2000,year), breaks = seq(2000,year,5), expand=c(0,0)) +
  theme(legend.position="bottom") +
  theme(legend.direction='vertical') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(strip.background = element_rect(colour="black", fill="white")) +
  guides(colour = guide_legend(nrow = 1, title=""))
swaziland.ART.calib

labs.age <- c("[15,25)"="15-24","[25,35)"="25-34","[35,50)"="35-49")
table(trajectory_numbersonART$scenario)
test <- subset(trajectory_numbersonART, scenario=="baseline" | scenario=="ART100pct" | scenario=="909090")
table(test$scenario)
trajectory_numbersonART$scenario_f = factor(trajectory_numbersonART$scenario, levels=c('ART100pct','909090',"baseline"))
swaziland.ART.calib <- ggplot() +
  geom_smooth(data=subset(test),
              aes(x=Year2, y=On_ART*scale.factor, group = factor(scenario), color=factor(scenario)), method="loess",span=0.2, se=F, size=1.2) +
  facet_grid(Agecat~Gender, labeller=labeller(Gender = labs, Agecat=labs.age), scales = "free") +
  theme(legend.position="bottom") +
  #geom_hline(aes(yintercept=81)) +
  xlab("Year")+
  ylab("Number on ART")+
  theme_bw(base_size=20) +
  scale_color_manual(labels = c("Age-targeted 90-90-90 (scenario 4)", "100% ART initation (scenario 5)",  "Status-quo (scenario 1)"), values = c("purple", "green", "blue")) +
  guides(fill = guide_legend(keywidth = 2, keyheight = 1)) +
  #scale_y_continuous(limits=c(0,100), breaks = seq(0,100,10), expand=c(0,0)) +
  scale_x_continuous(limits = c(2000,year), breaks = seq(2000,year,5), expand=c(0,0)) +
  theme(legend.position="bottom",legend.direction='vertical') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(strip.background = element_rect(colour="black", fill="white")) +
  guides(colour = guide_legend(nrow = 3, title="")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"))
swaziland.ART.calib

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\BMGF Lancet Paper\\Submission\\Lancet HIV\\Resubmission2\\Figures")
ggsave("swaziland_numbersonART.byage.gender2_scaleup_by_scenario_new.pdf", height=13, width=9)

#by gender
head(trajectory_ART_calib.master,200)
trajectory_ART_calib.master$Year2 <- floor((trajectory_ART_calib.master$Year-0.5))
trajectory_numbersonART <- aggregate(On_ART ~ Year2+Gender+sim.id+scenario,
                                     subset(trajectory_ART_calib.master), FUN=sum)
trajectory_numbersonART$scenario_f = factor(trajectory_numbersonART$scenario, levels=c('baseline','909090',"ART100pct"))
table(trajectory_numbersonART$scenario_f)
year=2048
swaziland.ART.calib <- ggplot() +
  geom_smooth(data=subset(trajectory_numbersonART, scenario_f=="909090"),
              aes(x=Year2, y=On_ART*scale.factor), method="loess",span=0.2, se=F, size=1.2, color="purple") +
  geom_smooth(data=subset(trajectory_numbersonART, scenario_f=="ART100pct"),
              aes(x=Year2, y=On_ART*scale.factor), method="loess",span=0.2, se=F, size=1.2, color="green") +
  geom_smooth(data=subset(trajectory_numbersonART, scenario_f=="baseline"),
              aes(x=Year2, y=On_ART*scale.factor), method="loess",span=0.2, se=F, size=1.2, color="blue") +
  facet_grid(~Gender, labeller=labeller(Gender = labs)) +
  theme(legend.position="bottom") +
  ylab("Number on ART")+
  theme_bw(base_size=20) +
  guides(fill = guide_legend(keywidth = 2, keyheight = 1)) +
  scale_x_continuous(limits = c(2000,year), breaks = seq(2000,year,5), expand=c(0,0)) +
  theme(legend.position="bottom") +
  theme(legend.direction='vertical') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(strip.background = element_rect(colour="black", fill="white")) +
  guides(colour = guide_legend(nrow = 1, title=""))
swaziland.ART.calib

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\BMGF Lancet Paper\\Final_Plots")
ggsave("swaziland_numbersonART.bygender2_scaleup_by_scenario.jpg", height=8, width=10)

#################################################
#Risk group sizes
##################################################
head(reporthivbyageandgender.master.baselineriskgroup)
reporthivbyageandgender.master <- reporthivbyageandgender.master.baselineriskgroup
trajectories.pop <- aggregate(Population ~ Year+Gender+sim.id, subset(reporthivbyageandgender.master), FUN=sum) #sums prev infections in each year
trajectories_population.age.risk <- aggregate(Population ~ Year+IP_Key.Risk+Gender+sim.id, subset(reporthivbyageandgender.master), FUN=sum)
trajectories_prev.age <- merge(trajectories.pop, trajectories_population.age.risk, by=c("sim.id","Gender","Year"))
trajectories_prev.age$riskproportion <- trajectories_prev.age$Population.y/trajectories_prev.age$Population.x
head(trajectories_prev.age,100)

trajectories_prev.age$riskgroup_f = factor(trajectories_prev.age$IP_Key.Risk, levels=c('LOW','MEDIUM',"HIGH"))

year=2050
#labs.age <- c("[15,25)" = "15-24", "[25,35)" = "25-34", "[35,45)" = "35-44")
swaziland.risk.size <- ggplot() +
  geom_smooth(data=subset(trajectories_prev.age),
              aes(x=Year, y=riskproportion, group = factor(riskgroup_f), color=factor(riskgroup_f)), method="loess",span=0.2, se=F, size=1.2) +
  facet_grid(~Gender, labeller=labeller(Gender = labs)) +
  theme(legend.position="bottom") +
  ylab("Population size")+
  theme_bw(base_size=20) +
  #scale_color_manual(labels = c("15-24", "25-34", "35-44"), values = c("blue", "purple", "red")) +
  guides(fill = guide_legend(keywidth = 2, keyheight = 1)) +
  scale_y_continuous(limits=c(0,1), breaks = seq(0,1,.10), expand=c(0,0)) +
  scale_x_continuous(limits = c(1980,year), breaks = seq(1980,year,5), expand=c(0,0)) +
  theme(legend.position="bottom") +
  theme(legend.direction='vertical') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(strip.background = element_rect(colour="black", fill="white")) +
  guides(colour = guide_legend(nrow = 1, title=""))
swaziland.risk.size

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\BMGF Lancet Paper\\Final_Plots")
ggsave("Riskgroup.sizes.jpg", height=8, width=10)

##################################################
#Cumulative person years on ART by scenario
#################################################

#by age and gender
trajectory_ART_calib.master$Year2 <- floor((trajectory_ART_calib.master$Year-0.5))
trajectory_numbersonART <- aggregate(On_ART ~ Year2+Agecat+Gender+sim.id+scenario,
                                     subset(trajectory_ART_calib.master), FUN=sum)
trajectory_numbersonART$scenario_f = factor(trajectory_numbersonART$scenario, levels=c('baseline','909090',"ART100pct"))
table(trajectory_numbersonART$scenario_f)
head(trajectory_numbersonART)

test <- trajectory_numbersonART[order(trajectory_numbersonART$scenario, trajectory_numbersonART$sim.id, trajectory_numbersonART$Gender, 
                                      trajectory_numbersonART$Agecat, trajectory_numbersonART$Year2),]

#calculate cumulative incidence by year 
test$csum <- ave(test$On_ART, test$scenario, 
                 test$sim.id, test$Gender, 
                 test$Agecat, FUN=cumsum)

head(test)

year=2049
swaziland.ART.calib <- ggplot() +
  geom_smooth(data=subset(test, scenario_f=="baseline"),
              aes(x=Year2, y=csum*scale.factor, group = factor(Agecat), color=factor(Agecat)), method="loess",span=0.2, se=F, size=1.2) +
  geom_smooth(data=subset(test, scenario_f=="909090"),
              aes(x=Year2, y=csum*scale.factor, group = factor(Agecat), color=factor(Agecat)), method="loess",span=0.2, se=F, size=1.2, linetype=2) +
  geom_smooth(data=subset(test, scenario_f=="ART100pct"),
              aes(x=Year2, y=csum*scale.factor, group = factor(Agecat), color=factor(Agecat)), method="loess",span=0.2, se=F, size=1.2, linetype=3) +
  facet_grid(~Gender, labeller=labeller(Gender = labs)) +
  theme(legend.position="bottom") +
  #geom_hline(aes(yintercept=81)) +
  ylab("Number on ART")+
  theme_bw(base_size=20) +
  scale_color_manual(labels = c("15-24", "25-34", "35-44"), values = c("blue", "purple", "red")) +
  guides(fill = guide_legend(keywidth = 2, keyheight = 1)) +
  #scale_y_continuous(limits=c(0,100), breaks = seq(0,100,10), expand=c(0,0)) +
  scale_x_continuous(limits = c(2000,year+1), breaks = seq(2000,year,5), expand=c(0,0)) +
  theme(legend.position="bottom") +
  theme(legend.direction='vertical') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(strip.background = element_rect(colour="black", fill="white")) +
  guides(colour = guide_legend(nrow = 1, title=""))
swaziland.ART.calib

trajectory_ART_calib.master$Year2 <- floor((trajectory_ART_calib.master$Year-0.5))
trajectory_numbersonART <- aggregate(On_ART ~ Year2+Gender+sim.id+scenario,
                                     subset(trajectory_ART_calib.master), FUN=sum)
trajectory_numbersonART$scenario_f = factor(trajectory_numbersonART$scenario, levels=c('baseline','909090',"ART100pct"))
table(trajectory_numbersonART$scenario_f)
head(trajectory_numbersonART)

test <- trajectory_numbersonART[order(trajectory_numbersonART$scenario, trajectory_numbersonART$sim.id, trajectory_numbersonART$Gender, 
                                      trajectory_numbersonART$Year2),]

#calculate cumulative incidence by year 
test$csum <- ave(test$On_ART, test$scenario, 
                 test$sim.id, test$Gender, 
                 FUN=cumsum)

head(test)

year=2049
swaziland.ART.calib <- ggplot() +
  geom_smooth(data=subset(test, scenario_f=="baseline"),
              aes(x=Year2, y=csum*scale.factor), method="loess",span=0.2, se=F, size=1.2) +
  geom_smooth(data=subset(test, scenario_f=="909090"),
              aes(x=Year2, y=csum*scale.factor), method="loess",span=0.2, se=F, size=1.2, linetype=2) +
  geom_smooth(data=subset(test, scenario_f=="ART100pct"),
              aes(x=Year2, y=csum*scale.factor), method="loess",span=0.2, se=F, size=1.2, linetype=3) +
  facet_grid(~Gender, labeller=labeller(Gender = labs)) +
  theme(legend.position="bottom") +
  #geom_hline(aes(yintercept=81)) +
  ylab("Number on ART")+
  theme_bw(base_size=20) +
  scale_color_manual(labels = c("15-24", "25-34", "35-44"), values = c("blue", "purple", "red")) +
  guides(fill = guide_legend(keywidth = 2, keyheight = 1)) +
  #scale_y_continuous(limits=c(0,100), breaks = seq(0,100,10), expand=c(0,0)) +
  scale_x_continuous(limits = c(2000,year), breaks = seq(2000,year,5), expand=c(0,0)) +
  theme(legend.position="bottom") +
  theme(legend.direction='vertical') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(strip.background = element_rect(colour="black", fill="white")) +
  guides(colour = guide_legend(nrow = 1, title=""))
swaziland.ART.calib

trajectory_ART_calib.master$Year2 <- floor((trajectory_ART_calib.master$Year-0.5))
trajectory_numbersonART <- aggregate(On_ART ~ Year2+sim.id+scenario,
                                     subset(trajectory_ART_calib.master), FUN=sum)
trajectory_numbersonART$scenario_f = factor(trajectory_numbersonART$scenario, levels=c('baseline','909090',"ART100pct"))
table(trajectory_numbersonART$scenario_f)
head(trajectory_numbersonART)

test <- trajectory_numbersonART[order(trajectory_numbersonART$scenario, trajectory_numbersonART$sim.id,
                                      trajectory_numbersonART$Year2),]

#calculate cumulative incidence by year 
test$csum <- ave(test$On_ART, test$scenario, 
                 test$sim.id, 
                 FUN=cumsum)

head(test)

year=2049
swaziland.ART.calib <- ggplot() +
  geom_smooth(data=subset(test, scenario_f=="baseline"),
              aes(x=Year2, y=csum*scale.factor), method="loess",span=0.2, se=F, size=1.2) +
  geom_smooth(data=subset(test, scenario_f=="909090"),
              aes(x=Year2, y=csum*scale.factor), method="loess",span=0.2, se=F, size=1.2, linetype=2) +
  geom_smooth(data=subset(test, scenario_f=="ART100pct"),
              aes(x=Year2, y=csum*scale.factor), method="loess",span=0.2, se=F, size=1.2, linetype=3) +
  #facet_grid(~Gender, labeller=labeller(Gender = labs)) +
  theme(legend.position="bottom") +
  #geom_hline(aes(yintercept=81)) +
  ylab("Year")+
  ylab("Number on ART")+
  theme_bw(base_size=20) +
  scale_color_manual(labels = c("15-24", "25-34", "35-44"), values = c("blue", "purple", "red")) +
  guides(fill = guide_legend(keywidth = 2, keyheight = 1)) +
  #scale_y_continuous(trans='log10',limits=c(0,100), breaks = seq(0,100,10), expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  scale_x_continuous(limits = c(2000,year+1), breaks = seq(2000,year,5), expand=c(0,0)) +
  theme(legend.position="bottom") +
  theme(legend.direction='vertical') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(strip.background = element_rect(colour="black", fill="white")) +
  guides(colour = guide_legend(nrow = 1, title=""))
swaziland.ART.calib

###########################################
#Calculate mortality
###########################################

#Calculate Incidence and plot with calibration points
reporthivbyageandgender.master$Year2 <- floor((reporthivbyageandgender.master$Year-0.5))
trajectories_MR.1a <- aggregate(Died_from_HIV ~ Year2+sim.id+scenario, subset(reporthivbyageandgender.master), FUN=sum) #sums number of new deaths in each year
trajectories_MR.2 <- aggregate(Population ~ Year+sim.id+scenario, subset(reporthivbyageandgender.master), FUN=sum) #population at midyear
trajectories_MR.2$Year2 <- floor(trajectories_MR.2$Year-0.5)
trajectories_MR.2 <- trajectories_MR.2[!duplicated(trajectories_MR.2[c("Year2","sim.id","scenario")]),] #remove second instance of duplicate rows
trajectories_MR.2 <- trajectories_MR.2[-match("Year",names(trajectories_MR.2))]
trajectories_MRoverall <- merge(trajectories_MR.1a, trajectories_MR.2, by=c("Year2","sim.id","scenario"))
trajectories_MRoverall$mortality <- trajectories_MRoverall$Died_from_HIV / trajectories_MRoverall$Population
head(trajectories_MRoverall)

#Plot baseline mortality with calibration points
year=2050
#Calibration plot
swaziland.mortality <- ggplot(data=subset(trajectories_MRoverall, Year2 < year & scenario=="baseline")) +
  geom_point(data=subset(trajectories_MRoverall,  Year2 < year & scenario=="baseline"), size=1.2, color = "grey", aes(x=Year2, y=mortality*100))+
  geom_line(data=subset(trajectories_MRoverall,  Year2 < year & scenario=="baseline"), color="grey",aes(x=Year2, y=mortality*100, group=sim.id)) +
  geom_smooth(aes(x=Year2, y=mortality*100),method="loess", span=0.1, se = T, size=1, color="blue", linetype=1) +
  xlab("Year")+
  ylab("Mortality (per 100 py)")+
  theme_bw(base_size=16) +
  guides(fill = guide_legend(keywidth = 2, keyheight = 1)) +
  scale_y_continuous(breaks = seq(0,2,0.5),limits=c(0,2),expand = c(0,0)) +
  scale_x_continuous(breaks = seq(1975,year,5),limits=c(1975,year+1),expand = c(0,0)) +
  theme(legend.position="bottom") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(strip.background = element_rect(colour="black", fill="white")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"))
swaziland.mortality

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\BMGF Lancet Paper\\Final_Plots")
ggsave("swaziland.mortality.baseline.jpg", height=8, width=8)

year1 = 1975
year2 = 2050
head(trajectories_MRoverall)
trajectories_MRoverall$scenario_f <- ordered(trajectories_MRoverall$scenario, levels=c("ART100pct","noARTnoVMMC","noART","909090","baseline"))
swaziland.mortality.scenario <- ggplot() +
  geom_smooth(data=subset(trajectories_MRoverall, Year2 < year), aes(x=Year2, y=mortality*100, color=scenario_f),method="loess", span=span, se = F, size=2) +
  xlab("Year")+
  ylab("Mortality (per 100 py)")+
  #facet_grid(Gender~Agecat, labeller=labeller(Gender = labs, Agecat = labs.age)) +
  theme_bw(base_size=24) +
  scale_color_manual(labels = c("100% ART uptake","No interventions", "VMMC only","Age-targeted 90-90-90", "Baseline"), values = c("green","red", "orange", "purple", "blue")) +
  guides(fill = guide_legend(keywidth = 2, keyheight = 1)) +
  scale_y_continuous(breaks = seq(0,3,0.5),limits=c(0,3), expand = c(0,0)) +
  scale_x_continuous(breaks = seq(year1,year2,5), limits=c(year1, year2+1),expand = c(0,0)) +
  theme(legend.position="bottom") +
  theme(legend.title=element_blank())+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(strip.background = element_rect(colour="black", fill="white")) +
  guides(col = guide_legend(nrow=1)) +
  #geom_hline(yintercept=0.1, linetype=2) +
  guides(col = guide_legend(nrow=3)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"))
swaziland.mortality.scenario

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\BMGF Lancet Paper\\Submission\\Lancet HIV\\Resubmission2\\Figures")
ggsave("swaziland.mortality.by.scenario.pdf", height=8, width=8)

#combine incidence and mortality by scenarios
year = 2051
span = 0.25
names(trajectories_MRoverall)
trajectories_MRoverall$metric <- "Mortality"
trajectories_MRoverall$value <- trajectories_MRoverall$mortality
names(trajectories_IRoverall)
trajectories_IRoverall$metric <- "Incidence"
trajectories_IRoverall$value <- trajectories_IRoverall$incidence
trajectories_MR_IR <- rbind(trajectories_IRoverall[,c(1,2,3,8,9,10)], trajectories_MRoverall[,c(1,2,3,7,8,9)])
trajectories_MR_IR$epidemiccontrol <- NA
trajectories_MR_IR$epidemiccontrol[trajectories_MR_IR$metric=="Incidence"] <- 0.1
trajectories_MR_IR$epidemiccontrol[trajectories_MR_IR$metric=="Mortality"] <- NA
summary(trajectories_MR_IR$epidemiccontrol)

trajectories_MR_IR$scenario_f <- ordered(trajectories_MR_IR$scenario, levels=c("ART100pct","noARTnoVMMC","noART","909090","baseline"))

ggplot() +
  geom_smooth(data=subset(trajectories_MR_IR, Year2 < year), aes(x=Year2, y=value*100, color=scenario_f),method="loess", span=span, se = F, size=2) +
  geom_line(data=subset(trajectories_MR_IR,Year2<year), aes(x=Year2, y=epidemiccontrol),color="black", size=1, linetype=2) +
  facet_grid(~metric) +
  xlab("Year")+
  ylab("Rate (per 100 py)")+
  theme_bw(base_size=16) +
  scale_color_manual(labels = c("100% ART initiation","No interventions", "VMMC only","Age-targeted 90-90-90", "Status-quo"), values = c("green","red", "orange", "purple", "blue")) +
  guides(fill = guide_legend(keywidth = 2, keyheight = 1)) +
  #geom_hline(yintercept=trajectories_MR_IR$epidemiccontrol, linetype=2) +
  scale_y_continuous(breaks = seq(0,5,0.5), expand = c(0,0)) +
  coord_cartesian(ylim=c(0, 4.5)) +
  scale_x_continuous(breaks = seq(1980,year,5), limits=c(1979, year+0.5),expand = c(0,0)) +
  theme(legend.position="bottom") +
  theme(legend.title=element_blank())+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(strip.background = element_rect(colour="black", fill="white")) +
  guides(col = guide_legend(nrow=1)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"))

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\BMGF Lancet Paper\\Submission\\Lancet HIV\\Resubmission2\\Figures")
ggsave("swaziland.incidence.mortality.by.scenario.pdf", height=10, width=17)

#Crude mortality numbers

table(trajectories_MRoverall$scenario)
names(trajectories_MRoverall)
swaziland.mort <- ggplot() +
  geom_smooth(data=subset(trajectories_MRoverall, Year2 < year), aes(x=Year2, y=Died_from_HIV*scale.factor, color=scenario_f),method="loess", span=span, se = F, size=2) +
  xlab("Year")+
  ylab("Annual deaths from HIV/AIDS")+
  theme_bw(base_size=24) +
  guides(fill = guide_legend(keywidth = 2, keyheight = 1)) +
  scale_color_manual(labels = c("100% ART uptake","No interventions", "VMMC only","Age-targeted 90-90-90", "Baseline"), values = c("green","red", "orange", "purple", "blue")) +
  scale_y_continuous(breaks = seq(0,20000,5000),limits=c(0,22000),expand = c(0,0)) +
  scale_x_continuous(breaks = seq(1980,year,5), limits=c(1980, year+1),expand = c(0,0)) +
  theme(legend.position="bottom") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(strip.background = element_rect(colour="black", fill="white")) +
  guides(col = guide_legend(nrow=3)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.title=element_blank(),
        panel.border = element_rect(colour = "black"))
swaziland.mort

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\BMGF Lancet Paper\\Final_Plots")
ggsave("swaziland.deaths.5scenarios.jpg", height=8, width=8)

########################################################
#Sensitivity of allocating ART to different risk groups
########################################################

reporthivbyageandgender.master <- reporthivbyageandgender.master.final.risk

reporthivbyageandgender.master$Year2 <- floor((reporthivbyageandgender.master$Year-0.5))
reporthivbyageandgender.master$Uninfected.Population = reporthivbyageandgender.master$Population-reporthivbyageandgender.master$Infected
trajectories_IR.1a <- aggregate(Newly.Infected ~ Year2+sim.id+scenario, subset(reporthivbyageandgender.master), FUN=sum) #sums number of new infections in each year
trajectories_IR.2 <- aggregate(Uninfected.Population ~ Year+sim.id+scenario, subset(reporthivbyageandgender.master), FUN=sum)
trajectories_IR.2$Year2 <- floor(trajectories_IR.2$Year-0.5)
trajectories_IR.2 <- trajectories_IR.2[!duplicated(trajectories_IR.2[c("Year2","sim.id","scenario")]),] #remove second instance of duplicate rows
trajectories_IR.2 <- trajectories_IR.2[-match("Year",names(trajectories_IR.2))]
trajectories_IRoverall <- merge(trajectories_IR.1a, trajectories_IR.2, by=c("Year2","sim.id","scenario"))
trajectories_IRoverall$incidence <- trajectories_IRoverall$Newly.Infected / trajectories_IRoverall$Uninfected.Population

table(trajectories_IRoverall$scenario)

#Plot
year = 2050
span = 0.7
alpha = 0.02
swaziland.inc.risk <- ggplot(trajectories_IRoverall) +
  geom_line(data=subset(trajectories_IRoverall, Year2 < year & scenario=="909090"),
            aes(x=Year2, y=incidence*100, group=sim.id), alpha=alpha, color="slateblue3") +
  geom_line(data=subset(trajectories_IRoverall, Year2 < year & scenario=="lowrisk100"),
            aes(x=Year2, y=incidence*100, group=sim.id), alpha=alpha, color="mediumpurple1") +
  geom_line(data=subset(trajectories_IRoverall, Year2 < year & scenario=="highrisk100"),
            aes(x=Year2, y=incidence*100, group=sim.id), alpha=alpha, color="slateblue4") +
  geom_smooth(data=subset(trajectories_IRoverall, Year2 < year & scenario=="lowrisk100"),
              aes(x=Year2, y=incidence*100),method="loess", span=span, se = F, size=2, color="mediumpurple1") +
  geom_smooth(data=subset(trajectories_IRoverall, Year2 < year & scenario=="909090"),
              aes(x=Year2, y=incidence*100),method="loess", span=span, se = F, size=2, color="slateblue3") +
  geom_smooth(data=subset(trajectories_IRoverall, Year2 < year & scenario=="highrisk100"),
              aes(x=Year2, y=incidence*100),method="loess", span=span, se = F, size=2, color="slateblue4") +
  xlab("Year")+
  ylab("Incidence (per 100 py)")+
  theme_bw(base_size=24) +
  guides(fill = guide_legend(keywidth = 2, keyheight = 1)) +
  scale_y_continuous(breaks = seq(0,1.2,0.1),limits=c(0,1.25),expand = c(0,0)) +
  scale_x_continuous(breaks = seq(1980,year,5), limits=c(2015, year+1),expand = c(0,0)) +
  theme(legend.position="bottom") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(strip.background = element_rect(colour="black", fill="white")) 
swaziland.inc.risk

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\BMGF Lancet Paper\\Final_Plots")
ggsave("swaziland.incidence.risk.scenarios_4.jpg", height=8, width=8)

#By Gender

table(reporthivbyageandgender.master$scenario)
reporthivbyageandgender.master$Year2 <- floor((reporthivbyageandgender.master$Year-0.5))
reporthivbyageandgender.master$Uninfected.Population <- reporthivbyageandgender.master$Population - reporthivbyageandgender.master$Infected
trajectories_IR.1a <- aggregate(Newly.Infected ~ Year2+Gender+sim.id+scenario, subset(reporthivbyageandgender.master, Age>10 & Age<50), FUN=sum) #sums number of new infections in each year
trajectories_IR.2 <- aggregate(Uninfected.Population ~ Year+Gender+sim.id+scenario, subset(reporthivbyageandgender.master, Age>10 & Age<50), FUN=sum)
trajectories_IR.2$Year2 <- floor(trajectories_IR.2$Year-0.5)
trajectories_IR.2 <- trajectories_IR.2[!duplicated(trajectories_IR.2[c("Year2","Gender","sim.id","scenario")]),] #remove second instance of duplicate rows
trajectories_IR.2 <- trajectories_IR.2[-match("Year",names(trajectories_IR.2))]
trajectories_IR <- merge(trajectories_IR.1a, trajectories_IR.2, by=c("Year2","Gender","sim.id","scenario"))
trajectories_IR$incidence <- trajectories_IR$Newly.Infected / trajectories_IR$Uninfected.Population
head(trajectories_IR)

trajectories_IRoverall$Gender = 2
trajectories_IR_comb <- rbind(trajectories_IR, trajectories_IRoverall)
trajectories_IR_comb <- subset(trajectories_IR_comb, Gender < 2)
incidence.data <- subset(incidence.data, Gender < 2)

#Plot
year = 2050
span = 0.7
alpha = 0.02
swaziland.inc.risk <- ggplot(trajectories_IR) +
  geom_smooth(data=subset(trajectories_IR, Year2 < year & scenario=="lowrisk100"),
              aes(x=Year2, y=incidence*100),method="loess", span=span, se = F, size=2, color="mediumpurple1") +
  geom_smooth(data=subset(trajectories_IR, Year2 < year & scenario=="909090"),
              aes(x=Year2, y=incidence*100),method="loess", span=span, se = F, size=2, color="slateblue3") +
  geom_smooth(data=subset(trajectories_IR, Year2 < year & scenario=="highrisk100"),
              aes(x=Year2, y=incidence*100),method="loess", span=span, se = F, size=2, color="slateblue4") +
  facet_grid(~Gender, labeller=labeller(Gender = labs)) +
  xlab("Year")+
  ylab("Incidence (per 100 py)")+
  theme_bw(base_size=24) +
  guides(fill = guide_legend(keywidth = 2, keyheight = 1)) +
  scale_y_continuous(breaks = seq(0,2,0.2),limits=c(0,2),expand = c(0,0)) +
  scale_x_continuous(breaks = seq(2010,year,5), limits=c(2015, year+1),expand = c(0,0)) +
  theme(legend.position="bottom") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(strip.background = element_rect(colour="black", fill="white")) 
swaziland.inc.risk

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\BMGF Lancet Paper\\Final_Plots")
ggsave("swaziland.incidence.riskbygender.scenarios_4.jpg", height=8, width=12)

#By Age and Gender
reporthivbyageandgender.master$Agecat<-cut(reporthivbyageandgender.master$Age, c(15,25,35,45), right=F)
table(reporthivbyageandgender.master$scenario)
names(reporthivbyageandgender.master)
reporthivbyageandgender.master$Year2 <- floor((reporthivbyageandgender.master$Year-0.5))
reporthivbyageandgender.master$Uninfected.Population <- reporthivbyageandgender.master$Population - reporthivbyageandgender.master$Infected
trajectories_IR.1a <- aggregate(Newly.Infected ~ Year2+Gender+Agecat+sim.id+scenario, subset(reporthivbyageandgender.master, Age>10 & Age<50), FUN=sum) #sums number of new infections in each year
trajectories_IR.2 <- aggregate(Uninfected.Population ~ Year+Gender+Agecat+sim.id+scenario, subset(reporthivbyageandgender.master, Age>10 & Age<50), FUN=sum)
trajectories_IR.2$Year2 <- floor(trajectories_IR.2$Year-0.5)
trajectories_IR.2 <- trajectories_IR.2[!duplicated(trajectories_IR.2[c("Year2","Gender","Agecat","sim.id","scenario")]),] #remove second instance of duplicate rows
trajectories_IR.2 <- trajectories_IR.2[-match("Year",names(trajectories_IR.2))]
trajectories_IR <- merge(trajectories_IR.1a, trajectories_IR.2, by=c("Year2","Gender","Agecat","sim.id","scenario"))
head(trajectories_IR)
trajectories_IR$incidence <- trajectories_IR$Newly.Infected / trajectories_IR$Uninfected.Population
# trajectories_IRoverall$Gender = 2
# trajectories_IR_comb <- rbind(trajectories_IR, trajectories_IRoverall)
# trajectories_IR_comb <- subset(trajectories_IR_comb, Gender < 2)
# incidence.data <- subset(incidence.data, Gender < 2)

#Plot
year = 2050
span = 0.7
alpha = 0.02
labs.age <- c("[15,25)" = "15-24", "[25,35)" = "25-34", "[35,45)" = "35-44")
swaziland.inc.risk <- ggplot(trajectories_IR) +
  # geom_line(data=subset(trajectories_IR, Year2 < year & scenario=="highrisk100"),
  #           aes(x=Year2, y=incidence*100, group=sim.id), alpha=alpha, color="firebrick4") +
  # geom_line(data=subset(trajectories_IR, Year2 < year & scenario=="baseline"),
  #           aes(x=Year2, y=incidence*100, group=sim.id), alpha=alpha, color="salmon") +
  # geom_line(data=subset(trajectories_IR, Year2 < year & scenario=="lowrisk100"),
  #           aes(x=Year2, y=incidence*100, group=sim.id), alpha=alpha, color="rosybrown1") +
  geom_smooth(data=subset(trajectories_IR, Year2 < year & scenario=="lowrisk100"),
              aes(x=Year2, y=incidence*100),method="loess", span=span, se = F, size=2, color="mediumpurple1") +
  geom_smooth(data=subset(trajectories_IR, Year2 < year & scenario=="909090"),
              aes(x=Year2, y=incidence*100),method="loess", span=span, se = F, size=2, color="slateblue3") +
  geom_smooth(data=subset(trajectories_IR, Year2 < year & scenario=="highrisk100"),
              aes(x=Year2, y=incidence*100),method="loess", span=span, se = F, size=2, color="slateblue4") +
  facet_grid(Agecat~Gender, labeller=labeller(Gender = labs, Agecat= labs.age)) +
  
  xlab("Year")+
  ylab("Incidence (per 100 py)")+
  theme_bw(base_size=24) +
  guides(fill = guide_legend(keywidth = 2, keyheight = 1)) +
  scale_y_continuous(breaks = seq(0,2,0.5),limits=c(0,2.1),expand = c(0,0)) +
  scale_x_continuous(breaks = seq(2015,year,5), limits=c(2015, year+1),expand = c(0,0)) +
  theme(legend.position="bottom") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"))
swaziland.inc.risk

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\BMGF Lancet Paper\\Final_Plots")
ggsave("swaziland.incidence.risk.scenarios.by.ageandgender_4.jpg", height=12, width=8)

#########################################################
#Transmitters
#########################################################
table(reporthivbyageandgender.master$scenario)
summary(reporthivbyageandgender.master$Transmitters)
reporthivbyageandgender.master$Agecat<-cut(reporthivbyageandgender.master$Age, c(15,25,35,45), right=F)
trajectories.transmitters <- aggregate(Transmitters ~ Year+Gender+Agecat+sim.id+scenario, 
                                       subset(reporthivbyageandgender.master, Age>10 & Age <50 & scenario != "hist909090" & scenario!="noARTnoVMMC"), FUN=sum) #sums prev infections in each year
head(trajectories.transmitters,10)

year=2050
labs <- c("1" = "Women", "0" = "Men")
labs.age <- c("[15,25)" = "15-24", "[25,35)" = "25-34", "[35,45)" = "35-44")
swaziland.transmitters <- ggplot() +
  geom_smooth(data=subset(trajectories.transmitters),
              aes(x=Year, y=Transmitters, group = factor(Agecat), color=factor(Agecat)), method="loess",span=0.2, se=F, size=1.2) +
  facet_grid(scenario~Gender, labeller=labeller(Gender = labs)) +
  theme(legend.position="bottom") +
  ylab("Population size")+
  theme_bw(base_size=20) +
  #scale_color_manual(labels = c("15-24", "25-34", "35-44"), values = c("blue", "purple", "red")) +
  guides(fill = guide_legend(keywidth = 2, keyheight = 1)) +
  scale_y_continuous(expand=c(0,0)) +
  scale_x_continuous(limits = c(1980,year), breaks = seq(1980,year,5), expand=c(0,0)) +
  theme(legend.position="bottom") +
  theme(legend.direction='vertical') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(strip.background = element_rect(colour="black", fill="white")) +
  guides(colour = guide_legend(nrow = 1, title=""))
swaziland.transmitters

######################################################################################
#Get the values
######################################################################################
##########################################################################
#Incidence point estimates at future time-points
#uses the aggregated reporthivbyageandgender by age-group
##########################################################################
setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\BMGF Lancet Paper\\Final_Plots")
summary(reporthivbyageandgender.master$Agecat)
reporthivbyageandgender.master$Uninfected.Population = reporthivbyageandgender.master$Population-reporthivbyageandgender.master$Infected
reporthivbyageandgender.master$Year2 <- floor(reporthivbyageandgender.master$Year-0.5)
head(reporthivbyageandgender.master)

#sum new infections over a one year period for each end year (2010, 2016, 2030, and 2050)
new.infections <- aggregate(Newly.Infected ~ sim.id+scenario+Year2, subset(reporthivbyageandgender.master), FUN=sum)
pop.uninfected <- aggregate(Uninfected.Population ~ sim.id+scenario+Year, subset(reporthivbyageandgender.master), FUN=sum) 
pop.uninfected <- pop.uninfected[!duplicated(pop.uninfected[c("Year","sim.id","scenario")]),] #remove second instance of duplicate rows
pop.uninfected$Year2 <- pop.uninfected$Year
incidence <- merge(new.infections,pop.uninfected,by=c("Year2","sim.id","scenario"))
incidence$inc <- incidence$Newly.Infected / incidence$Uninfected.Population
incidence <- incidence[incidence$Year==2010 | incidence$Year==2016 | incidence$Year==2030 | incidence$Year==2050, c(2:7)]
incidence.w <- reshape(incidence, timevar = "Year", idvar = c("sim.id","scenario"), direction = "wide")
incidence.w$rr.2010 <- incidence.w$inc.2010 / incidence.w$inc.2010
incidence.w$rr.2016 <- incidence.w$inc.2016 / incidence.w$inc.2010
incidence.w$rr.2030 <- incidence.w$inc.2030 / incidence.w$inc.2010
incidence.w$rr.2050 <- incidence.w$inc.2050 / incidence.w$inc.2010
head(incidence.w2)
incidence.w <- incidence.w[,c(1,2,5,8,11,14:18)]
incidence.w2 <- reshape(incidence.w, timevar = "scenario", idvar = c("sim.id"), direction = "wide")

quants <- c(0.5,0.025,0.975)
incidence.w2.q <- apply(incidence.w2[2:ncol(incidence.w2)], 2 ,quantile ,probs = quants)

#flip rows / columns

##########################################################################
#set years for time in which to calculate incidence. ignore start years, only there for consistent headers
y1.start = 2010
y2.start = 2016
y3.start = 2016
y1.end = 2016
y2.end = 2030
y3.end = 2050

#sum new infections over a one year period for each end year (2010, 2016, 2030, and 2050)
new.infections.0 <- aggregate(Newly.Infected ~ sim.id+scenario, subset(reporthivbyageandgender.master, Year==y1.start+0.5 | Year==y1.start+1), FUN=sum)
new.infections.1 <- aggregate(Newly.Infected ~ sim.id+scenario, subset(reporthivbyageandgender.master, Year==y1.end+0.5 | Year==y1.end+1), FUN=sum)
new.infections.2 <- aggregate(Newly.Infected ~ sim.id+scenario, subset(reporthivbyageandgender.master, Year==y2.end+0.5 | Year==y2.end+1), FUN=sum)
new.infections.3 <- aggregate(Newly.Infected ~ sim.id+scenario, subset(reporthivbyageandgender.master, Year==y3.end+0.5 | Year==y3.end+1), FUN=sum)
new.infections.0.g <- aggregate(Newly.Infected ~ Gender+sim.id+scenario, subset(reporthivbyageandgender.master, Year==y1.start+0.5 | Year==y1.start+1), FUN=sum)
new.infections.1.g <- aggregate(Newly.Infected ~ Gender+sim.id+scenario, subset(reporthivbyageandgender.master, Year==y1.end+0.5 | Year==y1.end+1), FUN=sum)
new.infections.2.g <- aggregate(Newly.Infected ~ Gender+sim.id+scenario, subset(reporthivbyageandgender.master, Year==y2.end+0.5 | Year==y2.end+1), FUN=sum)
new.infections.3.g <- aggregate(Newly.Infected ~ Gender+sim.id+scenario, subset(reporthivbyageandgender.master, Year==y3.end+0.5 | Year==y3.end+1), FUN=sum)

#sum uninfected pop at midyear for the same years
pop.uninfected.0 <- aggregate(Uninfected.Population ~ sim.id+scenario, subset(reporthivbyageandgender.master,Year==y1.start+0.5), FUN=sum) 
pop.uninfected.1 <- aggregate(Uninfected.Population ~ sim.id+scenario, subset(reporthivbyageandgender.master,Year==y1.end+0.5), FUN=sum) 
pop.uninfected.2 <- aggregate(Uninfected.Population ~ sim.id+scenario, subset(reporthivbyageandgender.master,Year==y2.end+0.5), FUN=sum) 
pop.uninfected.3 <- aggregate(Uninfected.Population ~ sim.id+scenario, subset(reporthivbyageandgender.master,Year==y3.end+0.5), FUN=sum) 
pop.uninfected.0.g <- aggregate(Uninfected.Population ~ Gender+sim.id+scenario, subset(reporthivbyageandgender.master,Year==y1.start+0.5), FUN=sum) 
pop.uninfected.1.g <- aggregate(Uninfected.Population ~ Gender+sim.id+scenario, subset(reporthivbyageandgender.master,Year==y1.end+0.5), FUN=sum) 
pop.uninfected.2.g <- aggregate(Uninfected.Population ~ Gender+sim.id+scenario, subset(reporthivbyageandgender.master,Year==y2.end+0.5), FUN=sum) 
pop.uninfected.3.g <- aggregate(Uninfected.Population ~ Gender+sim.id+scenario, subset(reporthivbyageandgender.master,Year==y3.end+0.5), FUN=sum) 

#merge new infections with uninfected denominator
cumulative.incidence.0 <- merge(new.infections.0, pop.uninfected.0, by=c("sim.id","scenario"))
cumulative.incidence.1 <- merge(new.infections.1, pop.uninfected.1, by=c("sim.id","scenario"))
cumulative.incidence.2 <- merge(new.infections.2, pop.uninfected.2, by=c("sim.id","scenario"))
cumulative.incidence.3 <- merge(new.infections.3, pop.uninfected.3, by=c("sim.id","scenario"))
cumulative.incidence.0.g <- merge(new.infections.0.g, pop.uninfected.0.g, by=c("Gender","sim.id","scenario"))
cumulative.incidence.1.g <- merge(new.infections.1.g, pop.uninfected.1.g, by=c("Gender","sim.id","scenario"))
cumulative.incidence.2.g <- merge(new.infections.2.g, pop.uninfected.2.g, by=c("Gender","sim.id","scenario"))
cumulative.incidence.3.g <- merge(new.infections.3.g, pop.uninfected.3.g, by=c("Gender","sim.id","scenario"))

#calculate incidence
cumulative.incidence.0$c.inc <- cumulative.incidence.0$Newly.Infected/cumulative.incidence.0$Uninfected.Population
cumulative.incidence.1$c.inc <- cumulative.incidence.1$Newly.Infected/cumulative.incidence.1$Uninfected.Population
cumulative.incidence.2$c.inc <- cumulative.incidence.2$Newly.Infected/cumulative.incidence.2$Uninfected.Population
cumulative.incidence.3$c.inc <- cumulative.incidence.3$Newly.Infected/cumulative.incidence.3$Uninfected.Population
cumulative.incidence.0.g$c.inc <- cumulative.incidence.0.g$Newly.Infected/cumulative.incidence.0.g$Uninfected.Population
cumulative.incidence.1.g$c.inc <- cumulative.incidence.1.g$Newly.Infected/cumulative.incidence.1.g$Uninfected.Population
cumulative.incidence.2.g$c.inc <- cumulative.incidence.2.g$Newly.Infected/cumulative.incidence.2.g$Uninfected.Population
cumulative.incidence.3.g$c.inc <- cumulative.incidence.3.g$Newly.Infected/cumulative.incidence.3.g$Uninfected.Population
cumulative.incidence.0$start <- y1.start
cumulative.incidence.0$end <- y1.start
cumulative.incidence.1$start <- y1.start
cumulative.incidence.1$end <- y1.end
cumulative.incidence.2$start <- y2.start
cumulative.incidence.2$end <- y2.end
cumulative.incidence.3$start <- y3.start
cumulative.incidence.3$end <- y3.end
cumulative.incidence.0.g$start <- y1.start
cumulative.incidence.0.g$end <- y1.start
cumulative.incidence.1.g$start <- y1.start
cumulative.incidence.1.g$end <- y1.end
cumulative.incidence.2.g$start <- y2.start
cumulative.incidence.2.g$end <- y2.end
cumulative.incidence.3.g$start <- y3.start
cumulative.incidence.3.g$end <- y3.end

#reshame so that each column is scenario- and metric specific
w.0 <- reshape(cumulative.incidence.0, timevar = "scenario", idvar = "sim.id", direction = "wide")
w.1 <- reshape(cumulative.incidence.1, timevar = "scenario", idvar = "sim.id", direction = "wide")
w.2 <- reshape(cumulative.incidence.2, timevar = "scenario", idvar = "sim.id", direction = "wide")
w.3 <- reshape(cumulative.incidence.3, timevar = "scenario", idvar = "sim.id", direction = "wide")
w.0.g <- reshape(cumulative.incidence.0.g, timevar = "scenario", idvar = c("sim.id","Gender") , direction = "wide")
w.1.g <- reshape(cumulative.incidence.1.g, timevar = "scenario", idvar = c("sim.id","Gender") , direction = "wide")
w.2.g <- reshape(cumulative.incidence.2.g, timevar = "scenario", idvar = c("sim.id","Gender") , direction = "wide")
w.3.g <- reshape(cumulative.incidence.3.g, timevar = "scenario", idvar = c("sim.id","Gender") , direction = "wide")

master.c.inc <- rbind(w.0,w.1,w.2,w.3)
master.c.inc$Gender <- 2
master.c.inc.g <- rbind(w.0.g,w.1.g,w.2.g,w.3.g)
master.c.inc.c <- rbind(master.c.inc, master.c.inc.g)
table(master.c.inc.c$end.noART)

head(master.c.inc.c)
names(master.c.inc.c)
test <- master.c.inc.c[,c(1,4,9,14,19,21,24,27)]
test$year <- test$end.noART
head(test)
test <- test[,c(1:4,6:8)]
test.w <- reshape(test, timevar = "year", idvar = c("sim.id","Gender") , direction = "wide")


#calculate incidence per 100 py
master.c.inc.c$rinc.baseline <- master.c.inc.c$c.inc.baseline*100
master.c.inc.c$rinc.noartnovmmc <- master.c.inc.c$c.inc.noARTnoVMMC*100
master.c.inc.c$rinc.noart <- master.c.inc.c$c.inc.noART*100
master.c.inc.c$rinc.909090 <- master.c.inc.c$c.inc.909090*100
master.c.inc.c$rinc.ART100pct <- master.c.inc.c$c.inc.ART100pct*100
head(master.c.inc.c)

#calculate relative risk over time (ref to 2010)

t0 <- quantile(master.c.inc.c$rinc.baseline[master.c.inc.c$end.noART==y1.start & master.c.inc.c$Gender==2], probs = c(0.5,0.025,0.975)) # quartile
t1 <- quantile(master.c.inc.c$rinc.baseline[master.c.inc.c$end.noART==y1.end & master.c.inc.c$Gender==2], probs = c(0.5,0.025,0.975)) # quartile
t2 <- quantile(master.c.inc.c$rinc.baseline[master.c.inc.c$end.noART==y2.end & master.c.inc.c$Gender==2], probs = c(0.5,0.025,0.975)) # quartile
t3 <- quantile(master.c.inc.c$rinc.baseline[master.c.inc.c$end.noART==y3.end & master.c.inc.c$Gender==2], probs = c(0.5,0.025,0.975)) # quartile
t0.m <- quantile(master.c.inc.c$rinc.baseline[master.c.inc.c$end.noART==y1.start & master.c.inc.c$Gender==0], probs = c(0.5,0.025,0.975)) # quartile
t1.m <- quantile(master.c.inc.c$rinc.baseline[master.c.inc.c$end.noART==y1.end & master.c.inc.c$Gender==0], probs = c(0.5,0.025,0.975)) # quartile
t2.m <- quantile(master.c.inc.c$rinc.baseline[master.c.inc.c$end.noART==y2.end & master.c.inc.c$Gender==0], probs = c(0.5,0.025,0.975)) # quartile
t3.m <- quantile(master.c.inc.c$rinc.baseline[master.c.inc.c$end.noART==y3.end & master.c.inc.c$Gender==0], probs = c(0.5,0.025,0.975)) # quartile
t0.f <- quantile(master.c.inc.c$rinc.baseline[master.c.inc.c$end.noART==y1.start & master.c.inc.c$Gender==1], probs = c(0.5,0.025,0.975)) # quartile
t1.f <- quantile(master.c.inc.c$rinc.baseline[master.c.inc.c$end.noART==y1.end & master.c.inc.c$Gender==1], probs = c(0.5,0.025,0.975)) # quartile
t2.f <- quantile(master.c.inc.c$rinc.baseline[master.c.inc.c$end.noART==y2.end & master.c.inc.c$Gender==1], probs = c(0.5,0.025,0.975)) # quartile
t3.f <- quantile(master.c.inc.c$rinc.baseline[master.c.inc.c$end.noART==y3.end & master.c.inc.c$Gender==1], probs = c(0.5,0.025,0.975)) # quartile
baseline <- as.data.frame(rbind(t0,t1,t2,t3,t0.f, t1.f,t2.f,t3.f,t0.m,t1.m,t2.m,t3.m))
baseline$startyear <- c(y1.start,y1.start,y2.start,y3.start,y1.start,y1.start,y2.start,y3.start,y1.start,y1.start,y2.start,y3.start)
baseline$endyear <- c(y1.start,y1.end,y2.end,y3.end,y1.start,y1.end,y2.end,y3.end,y1.start,y1.end,y2.end,y3.end)
baseline$scenario <- "ART+VMMC"
t0 <- quantile(master.c.inc.c$rinc.noartnovmmc[master.c.inc.c$end.noART==y1.start & master.c.inc.c$Gender==2], probs = c(0.5,0.025,0.975)) # quartile
t1 <- quantile(master.c.inc.c$rinc.noartnovmmc[master.c.inc.c$end.noART==y1.end & master.c.inc.c$Gender==2], probs = c(0.5,0.025,0.975)) # quartile
t2 <- quantile(master.c.inc.c$rinc.noartnovmmc[master.c.inc.c$end.noART==y2.end & master.c.inc.c$Gender==2], probs = c(0.5,0.025,0.975)) # quartile
t3 <- quantile(master.c.inc.c$rinc.noartnovmmc[master.c.inc.c$end.noART==y3.end & master.c.inc.c$Gender==2], probs = c(0.5,0.025,0.975)) # quartile
t0.m <- quantile(master.c.inc.c$rinc.noartnovmmc[master.c.inc.c$end.noART==y1.start & master.c.inc.c$Gender==0], probs = c(0.5,0.025,0.975)) # quartile
t1.m <- quantile(master.c.inc.c$rinc.noartnovmmc[master.c.inc.c$end.noART==y1.end & master.c.inc.c$Gender==0], probs = c(0.5,0.025,0.975)) # quartile
t2.m <- quantile(master.c.inc.c$rinc.noartnovmmc[master.c.inc.c$end.noART==y2.end & master.c.inc.c$Gender==0], probs = c(0.5,0.025,0.975)) # quartile
t3.m <- quantile(master.c.inc.c$rinc.noartnovmmc[master.c.inc.c$end.noART==y3.end & master.c.inc.c$Gender==0], probs = c(0.5,0.025,0.975)) # quartile
t0.f <- quantile(master.c.inc.c$rinc.noartnovmmc[master.c.inc.c$end.noART==y1.start & master.c.inc.c$Gender==1], probs = c(0.5,0.025,0.975)) # quartile
t1.f <- quantile(master.c.inc.c$rinc.noartnovmmc[master.c.inc.c$end.noART==y1.end & master.c.inc.c$Gender==1], probs = c(0.5,0.025,0.975)) # quartile
t2.f <- quantile(master.c.inc.c$rinc.noartnovmmc[master.c.inc.c$end.noART==y2.end & master.c.inc.c$Gender==1], probs = c(0.5,0.025,0.975)) # quartile
t3.f <- quantile(master.c.inc.c$rinc.noartnovmmc[master.c.inc.c$end.noART==y3.end & master.c.inc.c$Gender==1], probs = c(0.5,0.025,0.975)) # quartile
novmmcnoart <- as.data.frame(rbind(t0,t1,t2,t3,t0.f, t1.f,t2.f,t3.f,t0.m,t1.m,t2.m,t3.m))
novmmcnoart$startyear <- c(y1.start,y1.start,y2.start,y3.start,y1.start,y1.start,y2.start,y3.start,y1.start,y1.start,y2.start,y3.start)
novmmcnoart$endyear <- c(y1.start,y1.end,y2.end,y3.end,y1.start,y1.end,y2.end,y3.end,y1.start,y1.end,y2.end,y3.end)
novmmcnoart$scenario <- "noVMMCnoART"
t0 <- quantile(master.c.inc.c$rinc.noart[master.c.inc.c$end.noART==y1.start & master.c.inc.c$Gender==2], probs = c(0.5,0.025,0.975)) # quartile
t1 <- quantile(master.c.inc.c$rinc.noart[master.c.inc.c$end.noART==y1.end & master.c.inc.c$Gender==2], probs = c(0.5,0.025,0.975)) # quartile
t2 <- quantile(master.c.inc.c$rinc.noart[master.c.inc.c$end.noART==y2.end & master.c.inc.c$Gender==2], probs = c(0.5,0.025,0.975)) # quartile
t3 <- quantile(master.c.inc.c$rinc.noart[master.c.inc.c$end.noART==y3.end & master.c.inc.c$Gender==2], probs = c(0.5,0.025,0.975)) # quartile
t0.m <- quantile(master.c.inc.c$rinc.noart[master.c.inc.c$end.noART==y1.start & master.c.inc.c$Gender==0], probs = c(0.5,0.025,0.975)) # quartile
t1.m <- quantile(master.c.inc.c$rinc.noart[master.c.inc.c$end.noART==y1.end & master.c.inc.c$Gender==0], probs = c(0.5,0.025,0.975)) # quartile
t2.m <- quantile(master.c.inc.c$rinc.noart[master.c.inc.c$end.noART==y2.end & master.c.inc.c$Gender==0], probs = c(0.5,0.025,0.975)) # quartile
t3.m <- quantile(master.c.inc.c$rinc.noart[master.c.inc.c$end.noART==y3.end & master.c.inc.c$Gender==0], probs = c(0.5,0.025,0.975)) # quartile
t0.f <- quantile(master.c.inc.c$rinc.noart[master.c.inc.c$end.noART==y1.start & master.c.inc.c$Gender==1], probs = c(0.5,0.025,0.975)) # quartile
t1.f <- quantile(master.c.inc.c$rinc.noart[master.c.inc.c$end.noART==y1.end & master.c.inc.c$Gender==1], probs = c(0.5,0.025,0.975)) # quartile
t2.f <- quantile(master.c.inc.c$rinc.noart[master.c.inc.c$end.noART==y2.end & master.c.inc.c$Gender==1], probs = c(0.5,0.025,0.975)) # quartile
t3.f <- quantile(master.c.inc.c$rinc.noart[master.c.inc.c$end.noART==y3.end & master.c.inc.c$Gender==1], probs = c(0.5,0.025,0.975)) # quartile
noart <- as.data.frame(rbind(t0,t1,t2,t3,t0.f, t1.f,t2.f,t3.f,t0.m,t1.m,t2.m,t3.m))
noart$startyear <- c(y1.start,y1.start,y2.start,y3.start,y1.start,y1.start,y2.start,y3.start,y1.start,y1.start,y2.start,y3.start)
noart$endyear <- c(y1.start,y1.end,y2.end,y3.end,y1.start,y1.end,y2.end,y3.end,y1.start,y1.end,y2.end,y3.end)
noart$scenario <- "VMMConly"
t0 <- quantile(master.c.inc.c$rinc.909090[master.c.inc.c$end.noART==y1.start & master.c.inc.c$Gender==2], probs = c(0.5,0.025,0.975)) # quartile
t1 <- quantile(master.c.inc.c$rinc.909090[master.c.inc.c$end.noART==y1.end & master.c.inc.c$Gender==2], probs = c(0.5,0.025,0.975)) # quartile
t2 <- quantile(master.c.inc.c$rinc.909090[master.c.inc.c$end.noART==y2.end & master.c.inc.c$Gender==2], probs = c(0.5,0.025,0.975)) # quartile
t3 <- quantile(master.c.inc.c$rinc.909090[master.c.inc.c$end.noART==y3.end & master.c.inc.c$Gender==2], probs = c(0.5,0.025,0.975)) # quartile
t0.m <- quantile(master.c.inc.c$rinc.909090[master.c.inc.c$end.noART==y1.start & master.c.inc.c$Gender==0], probs = c(0.5,0.025,0.975)) # quartile
t1.m <- quantile(master.c.inc.c$rinc.909090[master.c.inc.c$end.noART==y1.end & master.c.inc.c$Gender==0], probs = c(0.5,0.025,0.975)) # quartile
t2.m <- quantile(master.c.inc.c$rinc.909090[master.c.inc.c$end.noART==y2.end & master.c.inc.c$Gender==0], probs = c(0.5,0.025,0.975)) # quartile
t3.m <- quantile(master.c.inc.c$rinc.909090[master.c.inc.c$end.noART==y3.end & master.c.inc.c$Gender==0], probs = c(0.5,0.025,0.975)) # quartile
t0.f <- quantile(master.c.inc.c$rinc.909090[master.c.inc.c$end.noART==y1.start & master.c.inc.c$Gender==1], probs = c(0.5,0.025,0.975)) # quartile
t1.f <- quantile(master.c.inc.c$rinc.909090[master.c.inc.c$end.noART==y1.end & master.c.inc.c$Gender==1], probs = c(0.5,0.025,0.975)) # quartile
t2.f <- quantile(master.c.inc.c$rinc.909090[master.c.inc.c$end.noART==y2.end & master.c.inc.c$Gender==1], probs = c(0.5,0.025,0.975)) # quartile
t3.f <- quantile(master.c.inc.c$rinc.909090[master.c.inc.c$end.noART==y3.end & master.c.inc.c$Gender==1], probs = c(0.5,0.025,0.975)) # quartile
age909090 <- as.data.frame(rbind(t0,t1,t2,t3,t0.f, t1.f,t2.f,t3.f,t0.m,t1.m,t2.m,t3.m))
age909090$startyear <- c(y1.start,y1.start,y2.start,y3.start,y1.start,y1.start,y2.start,y3.start,y1.start,y1.start,y2.start,y3.start)
age909090$endyear <- c(y1.start,y1.end,y2.end,y3.end,y1.start,y1.end,y2.end,y3.end,y1.start,y1.end,y2.end,y3.end)
age909090$scenario <- "909090"
t0 <- quantile(master.c.inc.c$rinc.ART100pct[master.c.inc.c$end.noART==y1.start & master.c.inc.c$Gender==2], probs = c(0.5,0.025,0.975)) # quartile
t1 <- quantile(master.c.inc.c$rinc.ART100pct[master.c.inc.c$end.noART==y1.end & master.c.inc.c$Gender==2], probs = c(0.5,0.025,0.975)) # quartile
t2 <- quantile(master.c.inc.c$rinc.ART100pct[master.c.inc.c$end.noART==y2.end & master.c.inc.c$Gender==2], probs = c(0.5,0.025,0.975)) # quartile
t3 <- quantile(master.c.inc.c$rinc.ART100pct[master.c.inc.c$end.noART==y3.end & master.c.inc.c$Gender==2], probs = c(0.5,0.025,0.975)) # quartile
t0.m <- quantile(master.c.inc.c$rinc.ART100pct[master.c.inc.c$end.noART==y1.start & master.c.inc.c$Gender==0], probs = c(0.5,0.025,0.975)) # quartile
t1.m <- quantile(master.c.inc.c$rinc.ART100pct[master.c.inc.c$end.noART==y1.end & master.c.inc.c$Gender==0], probs = c(0.5,0.025,0.975)) # quartile
t2.m <- quantile(master.c.inc.c$rinc.ART100pct[master.c.inc.c$end.noART==y2.end & master.c.inc.c$Gender==0], probs = c(0.5,0.025,0.975)) # quartile
t3.m <- quantile(master.c.inc.c$rinc.ART100pct[master.c.inc.c$end.noART==y3.end & master.c.inc.c$Gender==0], probs = c(0.5,0.025,0.975)) # quartile
t0.f <- quantile(master.c.inc.c$rinc.ART100pct[master.c.inc.c$end.noART==y1.start & master.c.inc.c$Gender==1], probs = c(0.5,0.025,0.975)) # quartile
t1.f <- quantile(master.c.inc.c$rinc.ART100pct[master.c.inc.c$end.noART==y1.end & master.c.inc.c$Gender==1], probs = c(0.5,0.025,0.975)) # quartile
t2.f <- quantile(master.c.inc.c$rinc.ART100pct[master.c.inc.c$end.noART==y2.end & master.c.inc.c$Gender==1], probs = c(0.5,0.025,0.975)) # quartile
t3.f <- quantile(master.c.inc.c$rinc.ART100pct[master.c.inc.c$end.noART==y3.end & master.c.inc.c$Gender==1], probs = c(0.5,0.025,0.975)) # quartile
ART100pct <- as.data.frame(rbind(t0,t1,t2,t3,t0.f, t1.f,t2.f,t3.f,t0.m,t1.m,t2.m,t3.m))
ART100pct$startyear <- c(y1.start,y1.start,y2.start,y3.start,y1.start,y1.start,y2.start,y3.start,y1.start,y1.start,y2.start,y3.start)
ART100pct$endyear <- c(y1.start,y1.end,y2.end,y3.end,y1.start,y1.end,y2.end,y3.end,y1.start,y1.end,y2.end,y3.end)
ART100pct$scenario <- "ART100pct"

results.final.incidence <- rbind(baseline,novmmcnoart, noart, age909090 ,ART100pct)
results.final.incidence$Gender <- c("All","All","All","All","Female","Female","Female","Female","Male","Male","Male","Male","All","All","All","All","Female","Female","Female","Female","Male","Male","Male","Male","All","All","All","All","Female","Female","Female","Female","Male","Male","Male","Male","All","All","All","All","Female","Female","Female","Female","Male","Male","Male","Male","All","All","All","All","Female","Female","Female","Female","Male","Male","Male","Male")
##########################################################################
#Reduction in cumulative incidence by age and gender
##########################################################################
reporthivbyageandgender.master$Uninfected.Population = reporthivbyageandgender.master$Population-reporthivbyageandgender.master$Infected

#set years for time in which to calculate incidence. ignore start years, only there for consistent headers
y1.start = 2010
y2.start = 2016
y3.start = 2016
y1.end = 2016
y2.end = 2030
y3.end = 2050

new.infections.1 <- aggregate(Newly.Infected ~ sim.id+scenario+Agecat, subset(reporthivbyageandgender.master, Year>y1.start & Year<y1.end+1.5), FUN=sum)
new.infections.2 <- aggregate(Newly.Infected ~ sim.id+scenario+Agecat, subset(reporthivbyageandgender.master, Year>y2.start & Year<y2.end+1.5), FUN=sum)
new.infections.3 <- aggregate(Newly.Infected ~ sim.id+scenario+Agecat, subset(reporthivbyageandgender.master, Year>y3.start & Year<y3.end+1.5), FUN=sum)
new.infections.1.g <- aggregate(Newly.Infected ~ Gender+sim.id+scenario+Agecat, subset(reporthivbyageandgender.master, Year>y1.start & Year<y1.end+1.5), FUN=sum)
new.infections.2.g <- aggregate(Newly.Infected ~ Gender+sim.id+scenario+Agecat, subset(reporthivbyageandgender.master, Year>y2.start & Year<y2.end+1.5), FUN=sum)
new.infections.3.g <- aggregate(Newly.Infected ~ Gender+sim.id+scenario+Agecat, subset(reporthivbyageandgender.master, Year>y3.start & Year<y3.end+1.5), FUN=sum)

pop.uninfected.1 <- aggregate(Uninfected.Population ~ sim.id+scenario+Agecat, subset(reporthivbyageandgender.master,Year>y1.start-0.5 & Year<y1.end+1), FUN=sum) 
pop.uninfected.2 <- aggregate(Uninfected.Population ~ sim.id+scenario+Agecat, subset(reporthivbyageandgender.master,Year>y2.start-0.5 & Year<y2.end+1), FUN=sum) 
pop.uninfected.3 <- aggregate(Uninfected.Population ~ sim.id+scenario+Agecat, subset(reporthivbyageandgender.master,Year>y3.start-0.5 & Year<y3.end+1), FUN=sum) 
pop.uninfected.1.g <- aggregate(Uninfected.Population ~ Gender+sim.id+scenario+Agecat, subset(reporthivbyageandgender.master,Year>y1.start-0.5 & Year<y1.end+1), FUN=sum) 
pop.uninfected.2.g <- aggregate(Uninfected.Population ~ Gender+sim.id+scenario+Agecat, subset(reporthivbyageandgender.master,Year>y2.start-0.5 & Year<y2.end+1), FUN=sum) 
pop.uninfected.3.g <- aggregate(Uninfected.Population ~ Gender+sim.id+scenario+Agecat, subset(reporthivbyageandgender.master,Year>y3.start-0.5 & Year<y3.end+1), FUN=sum) 

cumulative.incidence.1 <- merge(new.infections.1, pop.uninfected.1, by=c("Agecat","sim.id","scenario"))
cumulative.incidence.2 <- merge(new.infections.2, pop.uninfected.2, by=c("Agecat","sim.id","scenario"))
cumulative.incidence.3 <- merge(new.infections.3, pop.uninfected.3, by=c("Agecat","sim.id","scenario"))
cumulative.incidence.1.g <- merge(new.infections.1.g, pop.uninfected.1.g, by=c("Agecat","Gender","sim.id","scenario"))
cumulative.incidence.2.g <- merge(new.infections.2.g, pop.uninfected.2.g, by=c("Agecat","Gender","sim.id","scenario"))
cumulative.incidence.3.g <- merge(new.infections.3.g, pop.uninfected.3.g, by=c("Agecat","Gender","sim.id","scenario"))

cumulative.incidence.1$c.inc <- cumulative.incidence.1$Newly.Infected/cumulative.incidence.1$Uninfected.Population
cumulative.incidence.2$c.inc <- cumulative.incidence.2$Newly.Infected/cumulative.incidence.2$Uninfected.Population
cumulative.incidence.3$c.inc <- cumulative.incidence.3$Newly.Infected/cumulative.incidence.3$Uninfected.Population
cumulative.incidence.1.g$c.inc <- cumulative.incidence.1.g$Newly.Infected/cumulative.incidence.1.g$Uninfected.Population
cumulative.incidence.2.g$c.inc <- cumulative.incidence.2.g$Newly.Infected/cumulative.incidence.2.g$Uninfected.Population
cumulative.incidence.3.g$c.inc <- cumulative.incidence.3.g$Newly.Infected/cumulative.incidence.3.g$Uninfected.Population
cumulative.incidence.1$start <- y1.start
cumulative.incidence.1$end <- y1.end
cumulative.incidence.2$start <- y2.start
cumulative.incidence.2$end <- y2.end
cumulative.incidence.3$start <- y3.start
cumulative.incidence.3$end <- y3.end
cumulative.incidence.1.g$start <- y1.start
cumulative.incidence.1.g$end <- y1.end
cumulative.incidence.2.g$start <- y2.start
cumulative.incidence.2.g$end <- y2.end
cumulative.incidence.3.g$start <- y3.start
cumulative.incidence.3.g$end <- y3.end

w.1 <- reshape(cumulative.incidence.1, timevar = "scenario", idvar = c("Agecat","sim.id"), direction = "wide")
w.2 <- reshape(cumulative.incidence.2, timevar = "scenario", idvar = c("Agecat","sim.id"), direction = "wide")
w.3 <- reshape(cumulative.incidence.3, timevar = "scenario", idvar = c("Agecat","sim.id"), direction = "wide")
w.1.g <- reshape(cumulative.incidence.1.g, timevar = "scenario", idvar = c("sim.id","Gender","Agecat") , direction = "wide")
w.2.g <- reshape(cumulative.incidence.2.g, timevar = "scenario", idvar = c("sim.id","Gender","Agecat") , direction = "wide")
w.3.g <- reshape(cumulative.incidence.3.g, timevar = "scenario", idvar = c("sim.id","Gender","Agecat") , direction = "wide")

master.c.inc <- rbind(w.1,w.2,w.3)
master.c.inc$Gender <- 2
master.c.inc.g <- rbind(w.1.g,w.2.g,w.3.g)
master.c.inc.c <- rbind(master.c.inc, master.c.inc.g)
names(master.c.inc.c)

#calculate ratio (AR%)
master.c.inc.c$rinc.909090 <- (1 - (master.c.inc.c$c.inc.909090 / master.c.inc.c$c.inc.baseline))*100
master.c.inc.c$rinc.ART100pct <- (1 - (master.c.inc.c$c.inc.ART100pct / master.c.inc.c$c.inc.baseline))*100
head(master.c.inc.c)

t1 <- quantile(master.c.inc.c$rinc.909090[master.c.inc.c$end.noART==y1.end & master.c.inc.c$Gender==2 & master.c.inc.c$Agecat=="[15,25)"], probs = c(0.5,0.025,0.975)) # quartile
t2 <- quantile(master.c.inc.c$rinc.909090[master.c.inc.c$end.noART==y2.end & master.c.inc.c$Gender==2 & master.c.inc.c$Agecat=="[15,25)"], probs = c(0.5,0.025,0.975)) # quartile
t3 <- quantile(master.c.inc.c$rinc.909090[master.c.inc.c$end.noART==y3.end & master.c.inc.c$Gender==2 & master.c.inc.c$Agecat=="[15,25)"], probs = c(0.5,0.025,0.975)) # quartile
t1.m <- quantile(master.c.inc.c$rinc.909090[master.c.inc.c$end.noART==y1.end & master.c.inc.c$Gender==0 & master.c.inc.c$Agecat=="[15,25)"], probs = c(0.5,0.025,0.975)) # quartile
t2.m <- quantile(master.c.inc.c$rinc.909090[master.c.inc.c$end.noART==y2.end & master.c.inc.c$Gender==0 & master.c.inc.c$Agecat=="[15,25)"], probs = c(0.5,0.025,0.975)) # quartile
t3.m <- quantile(master.c.inc.c$rinc.909090[master.c.inc.c$end.noART==y3.end & master.c.inc.c$Gender==0 & master.c.inc.c$Agecat=="[15,25)"], probs = c(0.5,0.025,0.975)) # quartile
t1.f <- quantile(master.c.inc.c$rinc.909090[master.c.inc.c$end.noART==y1.end & master.c.inc.c$Gender==1 & master.c.inc.c$Agecat=="[15,25)"], probs = c(0.5,0.025,0.975)) # quartile
t2.f <- quantile(master.c.inc.c$rinc.909090[master.c.inc.c$end.noART==y2.end & master.c.inc.c$Gender==1 & master.c.inc.c$Agecat=="[15,25)"], probs = c(0.5,0.025,0.975)) # quartile
t3.f <- quantile(master.c.inc.c$rinc.909090[master.c.inc.c$end.noART==y3.end & master.c.inc.c$Gender==1 & master.c.inc.c$Agecat=="[15,25)"], probs = c(0.5,0.025,0.975)) # quartile
age909090 <- as.data.frame(rbind(t1,t2,t3,t1.f,t2.f,t3.f,t1.m,t2.m,t3.m))
age909090$startyear <- c(y1.start,y2.start,y3.start,y1.start,y2.start,y3.start,y1.start,y2.start,y3.start)
age909090$endyear <- c(y1.end,y2.end,y3.end,y1.end,y2.end,y3.end,y1.end,y2.end,y3.end)
age909090$scenario <- "909090"
age909090$Agecat <- "15-24"
t1 <- quantile(master.c.inc.c$rinc.ART100pct[master.c.inc.c$end.noART==y1.end & master.c.inc.c$Gender==2 & master.c.inc.c$Agecat=="[15,25)"], probs = c(0.5,0.025,0.975)) # quartile
t2 <- quantile(master.c.inc.c$rinc.ART100pct[master.c.inc.c$end.noART==y2.end & master.c.inc.c$Gender==2 & master.c.inc.c$Agecat=="[15,25)"], probs = c(0.5,0.025,0.975)) # quartile
t3 <- quantile(master.c.inc.c$rinc.ART100pct[master.c.inc.c$end.noART==y3.end & master.c.inc.c$Gender==2 & master.c.inc.c$Agecat=="[15,25)"], probs = c(0.5,0.025,0.975)) # quartile
t1.m <- quantile(master.c.inc.c$rinc.ART100pct[master.c.inc.c$end.noART==y1.end & master.c.inc.c$Gender==0 & master.c.inc.c$Agecat=="[15,25)"], probs = c(0.5,0.025,0.975)) # quartile
t2.m <- quantile(master.c.inc.c$rinc.ART100pct[master.c.inc.c$end.noART==y2.end & master.c.inc.c$Gender==0 & master.c.inc.c$Agecat=="[15,25)"], probs = c(0.5,0.025,0.975)) # quartile
t3.m <- quantile(master.c.inc.c$rinc.ART100pct[master.c.inc.c$end.noART==y3.end & master.c.inc.c$Gender==0 & master.c.inc.c$Agecat=="[15,25)"], probs = c(0.5,0.025,0.975)) # quartile
t1.f <- quantile(master.c.inc.c$rinc.ART100pct[master.c.inc.c$end.noART==y1.end & master.c.inc.c$Gender==1 & master.c.inc.c$Agecat=="[15,25)"], probs = c(0.5,0.025,0.975)) # quartile
t2.f <- quantile(master.c.inc.c$rinc.ART100pct[master.c.inc.c$end.noART==y2.end & master.c.inc.c$Gender==1 & master.c.inc.c$Agecat=="[15,25)"], probs = c(0.5,0.025,0.975)) # quartile
t3.f <- quantile(master.c.inc.c$rinc.ART100pct[master.c.inc.c$end.noART==y3.end & master.c.inc.c$Gender==1 & master.c.inc.c$Agecat=="[15,25)"], probs = c(0.5,0.025,0.975)) # quartile
ART100pct <- as.data.frame(rbind(t1,t2,t3,t1.f,t2.f,t3.f,t1.m,t2.m,t3.m))
ART100pct$startyear <- c(y1.start,y2.start,y3.start,y1.start,y2.start,y3.start,y1.start,y2.start,y3.start)
ART100pct$endyear <- c(y1.end,y2.end,y3.end,y1.end,y2.end,y3.end,y1.end,y2.end,y3.end)
ART100pct$scenario <- "ART100pct"
ART100pct$Agecat <- "15-24"

results.final.incidence.1524 <- rbind(age909090 ,ART100pct)
results.final.incidence.1524$Gender <- c("All","All","All","Female","Female","Female","Male","Male","Male","All","All","All","Female","Female","Female","Male","Male","Male")

t1 <- quantile(master.c.inc.c$rinc.909090[master.c.inc.c$end.noART==y1.end & master.c.inc.c$Gender==2 & master.c.inc.c$Agecat=="[25,35)"], probs = c(0.5,0.025,0.975)) # quartile
t2 <- quantile(master.c.inc.c$rinc.909090[master.c.inc.c$end.noART==y2.end & master.c.inc.c$Gender==2 & master.c.inc.c$Agecat=="[25,35)"], probs = c(0.5,0.025,0.975)) # quartile
t3 <- quantile(master.c.inc.c$rinc.909090[master.c.inc.c$end.noART==y3.end & master.c.inc.c$Gender==2 & master.c.inc.c$Agecat=="[25,35)"], probs = c(0.5,0.025,0.975)) # quartile
t1.m <- quantile(master.c.inc.c$rinc.909090[master.c.inc.c$end.noART==y1.end & master.c.inc.c$Gender==0 & master.c.inc.c$Agecat=="[25,35)"], probs = c(0.5,0.025,0.975)) # quartile
t2.m <- quantile(master.c.inc.c$rinc.909090[master.c.inc.c$end.noART==y2.end & master.c.inc.c$Gender==0 & master.c.inc.c$Agecat=="[25,35)"], probs = c(0.5,0.025,0.975)) # quartile
t3.m <- quantile(master.c.inc.c$rinc.909090[master.c.inc.c$end.noART==y3.end & master.c.inc.c$Gender==0 & master.c.inc.c$Agecat=="[25,35)"], probs = c(0.5,0.025,0.975)) # quartile
t1.f <- quantile(master.c.inc.c$rinc.909090[master.c.inc.c$end.noART==y1.end & master.c.inc.c$Gender==1 & master.c.inc.c$Agecat=="[25,35)"], probs = c(0.5,0.025,0.975)) # quartile
t2.f <- quantile(master.c.inc.c$rinc.909090[master.c.inc.c$end.noART==y2.end & master.c.inc.c$Gender==1 & master.c.inc.c$Agecat=="[25,35)"], probs = c(0.5,0.025,0.975)) # quartile
t3.f <- quantile(master.c.inc.c$rinc.909090[master.c.inc.c$end.noART==y3.end & master.c.inc.c$Gender==1 & master.c.inc.c$Agecat=="[25,35)"], probs = c(0.5,0.025,0.975)) # quartile
age909090 <- as.data.frame(rbind(t1,t2,t3,t1.f,t2.f,t3.f,t1.m,t2.m,t3.m))
age909090$startyear <- c(y1.start,y2.start,y3.start,y1.start,y2.start,y3.start,y1.start,y2.start,y3.start)
age909090$endyear <- c(y1.end,y2.end,y3.end,y1.end,y2.end,y3.end,y1.end,y2.end,y3.end)
age909090$scenario <- "909090"
age909090$Agecat <- "25-34"
t1 <- quantile(master.c.inc.c$rinc.ART100pct[master.c.inc.c$end.noART==y1.end & master.c.inc.c$Gender==2 & master.c.inc.c$Agecat=="[25,35)"], probs = c(0.5,0.025,0.975)) # quartile
t2 <- quantile(master.c.inc.c$rinc.ART100pct[master.c.inc.c$end.noART==y2.end & master.c.inc.c$Gender==2 & master.c.inc.c$Agecat=="[25,35)"], probs = c(0.5,0.025,0.975)) # quartile
t3 <- quantile(master.c.inc.c$rinc.ART100pct[master.c.inc.c$end.noART==y3.end & master.c.inc.c$Gender==2 & master.c.inc.c$Agecat=="[25,35)"], probs = c(0.5,0.025,0.975)) # quartile
t1.m <- quantile(master.c.inc.c$rinc.ART100pct[master.c.inc.c$end.noART==y1.end & master.c.inc.c$Gender==0 & master.c.inc.c$Agecat=="[25,35)"], probs = c(0.5,0.025,0.975)) # quartile
t2.m <- quantile(master.c.inc.c$rinc.ART100pct[master.c.inc.c$end.noART==y2.end & master.c.inc.c$Gender==0 & master.c.inc.c$Agecat=="[25,35)"], probs = c(0.5,0.025,0.975)) # quartile
t3.m <- quantile(master.c.inc.c$rinc.ART100pct[master.c.inc.c$end.noART==y3.end & master.c.inc.c$Gender==0 & master.c.inc.c$Agecat=="[25,35)"], probs = c(0.5,0.025,0.975)) # quartile
t1.f <- quantile(master.c.inc.c$rinc.ART100pct[master.c.inc.c$end.noART==y1.end & master.c.inc.c$Gender==1 & master.c.inc.c$Agecat=="[25,35)"], probs = c(0.5,0.025,0.975)) # quartile
t2.f <- quantile(master.c.inc.c$rinc.ART100pct[master.c.inc.c$end.noART==y2.end & master.c.inc.c$Gender==1 & master.c.inc.c$Agecat=="[25,35)"], probs = c(0.5,0.025,0.975)) # quartile
t3.f <- quantile(master.c.inc.c$rinc.ART100pct[master.c.inc.c$end.noART==y3.end & master.c.inc.c$Gender==1 & master.c.inc.c$Agecat=="[25,35)"], probs = c(0.5,0.025,0.975)) # quartile
ART100pct <- as.data.frame(rbind(t1,t2,t3,t1.f,t2.f,t3.f,t1.m,t2.m,t3.m))
ART100pct$startyear <- c(y1.start,y2.start,y3.start,y1.start,y2.start,y3.start,y1.start,y2.start,y3.start)
ART100pct$endyear <- c(y1.end,y2.end,y3.end,y1.end,y2.end,y3.end,y1.end,y2.end,y3.end)
ART100pct$scenario <- "ART100pct"
ART100pct$Agecat <- "25-34"

results.final.incidence.2534 <- rbind(age909090 ,ART100pct)
results.final.incidence.2534$Gender <- c("All","All","All","Female","Female","Female","Male","Male","Male","All","All","All","Female","Female","Female","Male","Male","Male")

t1 <- quantile(master.c.inc.c$rinc.909090[master.c.inc.c$end.noART==y1.end & master.c.inc.c$Gender==2 & master.c.inc.c$Agecat=="[35,50)"], probs = c(0.5,0.025,0.975)) # quartile
t2 <- quantile(master.c.inc.c$rinc.909090[master.c.inc.c$end.noART==y2.end & master.c.inc.c$Gender==2 & master.c.inc.c$Agecat=="[35,50)"], probs = c(0.5,0.025,0.975)) # quartile
t3 <- quantile(master.c.inc.c$rinc.909090[master.c.inc.c$end.noART==y3.end & master.c.inc.c$Gender==2 & master.c.inc.c$Agecat=="[35,50)"], probs = c(0.5,0.025,0.975)) # quartile
t1.m <- quantile(master.c.inc.c$rinc.909090[master.c.inc.c$end.noART==y1.end & master.c.inc.c$Gender==0 & master.c.inc.c$Agecat=="[35,50)"], probs = c(0.5,0.025,0.975)) # quartile
t2.m <- quantile(master.c.inc.c$rinc.909090[master.c.inc.c$end.noART==y2.end & master.c.inc.c$Gender==0 & master.c.inc.c$Agecat=="[35,50)"], probs = c(0.5,0.025,0.975)) # quartile
t3.m <- quantile(master.c.inc.c$rinc.909090[master.c.inc.c$end.noART==y3.end & master.c.inc.c$Gender==0 & master.c.inc.c$Agecat=="[35,50)"], probs = c(0.5,0.025,0.975)) # quartile
t1.f <- quantile(master.c.inc.c$rinc.909090[master.c.inc.c$end.noART==y1.end & master.c.inc.c$Gender==1 & master.c.inc.c$Agecat=="[35,50)"], probs = c(0.5,0.025,0.975)) # quartile
t2.f <- quantile(master.c.inc.c$rinc.909090[master.c.inc.c$end.noART==y2.end & master.c.inc.c$Gender==1 & master.c.inc.c$Agecat=="[35,50)"], probs = c(0.5,0.025,0.975)) # quartile
t3.f <- quantile(master.c.inc.c$rinc.909090[master.c.inc.c$end.noART==y3.end & master.c.inc.c$Gender==1 & master.c.inc.c$Agecat=="[35,50)"], probs = c(0.5,0.025,0.975)) # quartile
age909090 <- as.data.frame(rbind(t1,t2,t3,t1.f,t2.f,t3.f,t1.m,t2.m,t3.m))
age909090$startyear <- c(y1.start,y2.start,y3.start,y1.start,y2.start,y3.start,y1.start,y2.start,y3.start)
age909090$endyear <- c(y1.end,y2.end,y3.end,y1.end,y2.end,y3.end,y1.end,y2.end,y3.end)
age909090$scenario <- "909090"
age909090$Agecat <- "35-49"
t1 <- quantile(master.c.inc.c$rinc.ART100pct[master.c.inc.c$end.noART==y1.end & master.c.inc.c$Gender==2 & master.c.inc.c$Agecat=="[35,50)"], probs = c(0.5,0.025,0.975)) # quartile
t2 <- quantile(master.c.inc.c$rinc.ART100pct[master.c.inc.c$end.noART==y2.end & master.c.inc.c$Gender==2 & master.c.inc.c$Agecat=="[35,50)"], probs = c(0.5,0.025,0.975)) # quartile
t3 <- quantile(master.c.inc.c$rinc.ART100pct[master.c.inc.c$end.noART==y3.end & master.c.inc.c$Gender==2 & master.c.inc.c$Agecat=="[35,50)"], probs = c(0.5,0.025,0.975)) # quartile
t1.m <- quantile(master.c.inc.c$rinc.ART100pct[master.c.inc.c$end.noART==y1.end & master.c.inc.c$Gender==0 & master.c.inc.c$Agecat=="[35,50)"], probs = c(0.5,0.025,0.975)) # quartile
t2.m <- quantile(master.c.inc.c$rinc.ART100pct[master.c.inc.c$end.noART==y2.end & master.c.inc.c$Gender==0 & master.c.inc.c$Agecat=="[35,50)"], probs = c(0.5,0.025,0.975)) # quartile
t3.m <- quantile(master.c.inc.c$rinc.ART100pct[master.c.inc.c$end.noART==y3.end & master.c.inc.c$Gender==0 & master.c.inc.c$Agecat=="[35,50)"], probs = c(0.5,0.025,0.975)) # quartile
t1.f <- quantile(master.c.inc.c$rinc.ART100pct[master.c.inc.c$end.noART==y1.end & master.c.inc.c$Gender==1 & master.c.inc.c$Agecat=="[35,50)"], probs = c(0.5,0.025,0.975)) # quartile
t2.f <- quantile(master.c.inc.c$rinc.ART100pct[master.c.inc.c$end.noART==y2.end & master.c.inc.c$Gender==1 & master.c.inc.c$Agecat=="[35,50)"], probs = c(0.5,0.025,0.975)) # quartile
t3.f <- quantile(master.c.inc.c$rinc.ART100pct[master.c.inc.c$end.noART==y3.end & master.c.inc.c$Gender==1 & master.c.inc.c$Agecat=="[35,50)"], probs = c(0.5,0.025,0.975)) # quartile
ART100pct <- as.data.frame(rbind(t1,t2,t3,t1.f,t2.f,t3.f,t1.m,t2.m,t3.m))
ART100pct$startyear <- c(y1.start,y2.start,y3.start,y1.start,y2.start,y3.start,y1.start,y2.start,y3.start)
ART100pct$endyear <- c(y1.end,y2.end,y3.end,y1.end,y2.end,y3.end,y1.end,y2.end,y3.end)
ART100pct$scenario <- "ART100pct"
ART100pct$Agecat <- "35-49"

results.final.incidence.3549 <- rbind(age909090 ,ART100pct)
results.final.incidence.3549$Gender <- c("All","All","All","Female","Female","Female","Male","Male","Male","All","All","All","Female","Female","Female","Male","Male","Male")

results.final.incidence.age.gender <- rbind(results.final.incidence.1524, results.final.incidence.2534, results.final.incidence.3549)
names(results.final.incidence.age.gender)
results.final.incidence.age.gender <- results.final.incidence.age.gender[order(results.final.incidence.age.gender$scenario, results.final.incidence.age.gender$endyear),]

write.csv(results.final.incidence.age.gender, "results.final.incidence.age.gender.csv")

#Plot the table
head(results.final.incidence.age.gender)
table(results.final.incidence.age.gender$endyear)
results.final.incidence.age.gender$horizon <- paste0(results.final.incidence.age.gender$startyear,"-",
                                                     results.final.incidence.age.gender$endyear)
table(results.final.incidence.age.gender$horizon)


labs.scenario <- c("909090"="Age-targeted 909090 (Scenario 4)", "ART100pct"="100% ART (Scenario 5)")

ggplot(data=subset(results.final.incidence.age.gender, horizon!="2010-2016")) +
  geom_point(aes(x=Gender, y=`50%`, color=Agecat),position = position_dodge(width = 0.3))+
  geom_errorbar(aes(x=Gender, ymin=`2.5%`, ymax=`97.5%`,color=Agecat,width=0.2),position = position_dodge(width = 0.3)) +
  facet_grid(horizon~scenario,labeller=labeller(scenario =labs.scenario)) +
  xlab("")+
  ylab("Reduction (%) in cumulative incidence relative to status-quo ART (scenario 1)")+
  geom_hline(yintercept=0, linetype=2) +
  theme_bw(base_size=16) +
  guides(color = guide_legend(keywidth = 2, keyheight = 1, title="Age group")) +
  scale_y_continuous(breaks = seq(-5,50,10),limits=c(-12,55),expand = c(0,0)) +
  theme(legend.position="bottom") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  theme(strip.background = element_rect(colour="black", fill="white")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"))

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\BMGF Lancet Paper\\Submission\\Lancet HIV\\Resubmission2\\Figures")
ggsave("cumulativeindcidencebyageandgender.pdf", height=10, width=8)

##########################################################################
#reduction in cumulative incidence
##########################################################################
reporthivbyageandgender.master$Uninfected.Population = reporthivbyageandgender.master$Population-reporthivbyageandgender.master$Infected

#set years for time horizon over which cumulative incidence will be calculated
y1.start = 2010
y2.start = 2016
y3.start = 2016
y1.end = 2016
y2.end = 2030
y3.end = 2050

new.infections.1 <- aggregate(Newly.Infected ~ sim.id+scenario, subset(reporthivbyageandgender.master, Year>y1.start & Year<y1.end+1.5), FUN=sum)
new.infections.2 <- aggregate(Newly.Infected ~ sim.id+scenario, subset(reporthivbyageandgender.master, Year>y2.start & Year<y2.end+1.5), FUN=sum)
new.infections.3 <- aggregate(Newly.Infected ~ sim.id+scenario, subset(reporthivbyageandgender.master, Year>y3.start & Year<y3.end+1.5), FUN=sum)
new.infections.1.g <- aggregate(Newly.Infected ~ Gender+sim.id+scenario, subset(reporthivbyageandgender.master, Year>y1.start & Year<y1.end+1.5), FUN=sum)
new.infections.2.g <- aggregate(Newly.Infected ~ Gender+sim.id+scenario, subset(reporthivbyageandgender.master, Year>y2.start & Year<y2.end+1.5), FUN=sum)
new.infections.3.g <- aggregate(Newly.Infected ~ Gender+sim.id+scenario, subset(reporthivbyageandgender.master, Year>y3.start & Year<y3.end+1.5), FUN=sum)

pop.uninfected.1 <- aggregate(Uninfected.Population ~ sim.id+scenario, subset(reporthivbyageandgender.master,Year>y1.start-0.5 & Year<y1.end+1), FUN=sum) 
pop.uninfected.2 <- aggregate(Uninfected.Population ~ sim.id+scenario, subset(reporthivbyageandgender.master,Year>y2.start-0.5 & Year<y2.end+1), FUN=sum) 
pop.uninfected.3 <- aggregate(Uninfected.Population ~ sim.id+scenario, subset(reporthivbyageandgender.master,Year>y3.start-0.5 & Year<y3.end+1), FUN=sum) 
pop.uninfected.1.g <- aggregate(Uninfected.Population ~ Gender+sim.id+scenario, subset(reporthivbyageandgender.master,Year>y1.start-0.5 & Year<y1.end+1), FUN=sum) 
pop.uninfected.2.g <- aggregate(Uninfected.Population ~ Gender+sim.id+scenario, subset(reporthivbyageandgender.master,Year>y2.start-0.5 & Year<y2.end+1), FUN=sum) 
pop.uninfected.3.g <- aggregate(Uninfected.Population ~ Gender+sim.id+scenario, subset(reporthivbyageandgender.master,Year>y3.start-0.5 & Year<y3.end+1), FUN=sum) 

cumulative.incidence.1 <- merge(new.infections.1, pop.uninfected.1, by=c("sim.id","scenario"))
cumulative.incidence.2 <- merge(new.infections.2, pop.uninfected.2, by=c("sim.id","scenario"))
cumulative.incidence.3 <- merge(new.infections.3, pop.uninfected.3, by=c("sim.id","scenario"))
cumulative.incidence.1.g <- merge(new.infections.1.g, pop.uninfected.1.g, by=c("Gender","sim.id","scenario"))
cumulative.incidence.2.g <- merge(new.infections.2.g, pop.uninfected.2.g, by=c("Gender","sim.id","scenario"))
cumulative.incidence.3.g <- merge(new.infections.3.g, pop.uninfected.3.g, by=c("Gender","sim.id","scenario"))

cumulative.incidence.1$c.inc <- cumulative.incidence.1$Newly.Infected/cumulative.incidence.1$Uninfected.Population
cumulative.incidence.2$c.inc <- cumulative.incidence.2$Newly.Infected/cumulative.incidence.2$Uninfected.Population
cumulative.incidence.3$c.inc <- cumulative.incidence.3$Newly.Infected/cumulative.incidence.3$Uninfected.Population
cumulative.incidence.1.g$c.inc <- cumulative.incidence.1.g$Newly.Infected/cumulative.incidence.1.g$Uninfected.Population
cumulative.incidence.2.g$c.inc <- cumulative.incidence.2.g$Newly.Infected/cumulative.incidence.2.g$Uninfected.Population
cumulative.incidence.3.g$c.inc <- cumulative.incidence.3.g$Newly.Infected/cumulative.incidence.3.g$Uninfected.Population
cumulative.incidence.1$start <- y1.start
cumulative.incidence.1$end <- y1.end
cumulative.incidence.2$start <- y2.start
cumulative.incidence.2$end <- y2.end
cumulative.incidence.3$start <- y3.start
cumulative.incidence.3$end <- y3.end
cumulative.incidence.1.g$start <- y1.start
cumulative.incidence.1.g$end <- y1.end
cumulative.incidence.2.g$start <- y2.start
cumulative.incidence.2.g$end <- y2.end
cumulative.incidence.3.g$start <- y3.start
cumulative.incidence.3.g$end <- y3.end

w.1 <- reshape(cumulative.incidence.1, timevar = "scenario", idvar = "sim.id", direction = "wide")
w.2 <- reshape(cumulative.incidence.2, timevar = "scenario", idvar = "sim.id", direction = "wide")
w.3 <- reshape(cumulative.incidence.3, timevar = "scenario", idvar = "sim.id", direction = "wide")
w.1.g <- reshape(cumulative.incidence.1.g, timevar = "scenario", idvar = c("sim.id","Gender") , direction = "wide")
w.2.g <- reshape(cumulative.incidence.2.g, timevar = "scenario", idvar = c("sim.id","Gender") , direction = "wide")
w.3.g <- reshape(cumulative.incidence.3.g, timevar = "scenario", idvar = c("sim.id","Gender") , direction = "wide")

master.c.inc <- rbind(w.1,w.2,w.3)
master.c.inc$Gender <- 2
master.c.inc.g <- rbind(w.1.g,w.2.g,w.3.g)
master.c.inc.c <- rbind(master.c.inc, master.c.inc.g)
names(master.c.inc.c)

#calculate ratio (AR%)
master.c.inc.c$rinc.noartnovmmc <- (1 - (master.c.inc.c$c.inc.baseline / master.c.inc.c$c.inc.noART))*100
master.c.inc.c$rinc.noart <- (1 - (master.c.inc.c$c.inc.noART / master.c.inc.c$c.inc.noARTnoVMMC))*100
master.c.inc.c$rinc.909090 <- (1 - (master.c.inc.c$c.inc.909090 / master.c.inc.c$c.inc.noART))*100
master.c.inc.c$rinc.ART100pct <- (1 - (master.c.inc.c$c.inc.ART100pct / master.c.inc.c$c.inc.noART))*100
head(master.c.inc.c)

t1 <- quantile(master.c.inc.c$rinc.noartnovmmc[master.c.inc.c$end.noART==y1.end & master.c.inc.c$Gender==2], probs = c(0.5,0.025,0.975)) # quartile
t2 <- quantile(master.c.inc.c$rinc.noartnovmmc[master.c.inc.c$end.noART==y2.end & master.c.inc.c$Gender==2], probs = c(0.5,0.025,0.975)) # quartile
t3 <- quantile(master.c.inc.c$rinc.noartnovmmc[master.c.inc.c$end.noART==y3.end & master.c.inc.c$Gender==2], probs = c(0.5,0.025,0.975)) # quartile
t1.m <- quantile(master.c.inc.c$rinc.noartnovmmc[master.c.inc.c$end.noART==y1.end & master.c.inc.c$Gender==0], probs = c(0.5,0.025,0.975)) # quartile
t2.m <- quantile(master.c.inc.c$rinc.noartnovmmc[master.c.inc.c$end.noART==y2.end & master.c.inc.c$Gender==0], probs = c(0.5,0.025,0.975)) # quartile
t3.m <- quantile(master.c.inc.c$rinc.noartnovmmc[master.c.inc.c$end.noART==y3.end & master.c.inc.c$Gender==0], probs = c(0.5,0.025,0.975)) # quartile
t1.f <- quantile(master.c.inc.c$rinc.noartnovmmc[master.c.inc.c$end.noART==y1.end & master.c.inc.c$Gender==1], probs = c(0.5,0.025,0.975)) # quartile
t2.f <- quantile(master.c.inc.c$rinc.noartnovmmc[master.c.inc.c$end.noART==y2.end & master.c.inc.c$Gender==1], probs = c(0.5,0.025,0.975)) # quartile
t3.f <- quantile(master.c.inc.c$rinc.noartnovmmc[master.c.inc.c$end.noART==y3.end & master.c.inc.c$Gender==1], probs = c(0.5,0.025,0.975)) # quartile
novmmcnoart <- as.data.frame(rbind(t1,t2,t3,t1.f,t2.f,t3.f,t1.m,t2.m,t3.m))
novmmcnoart$startyear <- c(y1.start,y2.start,y3.start,y1.start,y2.start,y3.start,y1.start,y2.start,y3.start)
novmmcnoart$endyear <- c(y1.end,y2.end,y3.end,y1.end,y2.end,y3.end,y1.end,y2.end,y3.end)
novmmcnoart$scenario <- "ART+VMMC"
t1 <- quantile(master.c.inc.c$rinc.noart[master.c.inc.c$end.noART==y1.end & master.c.inc.c$Gender==2], probs = c(0.5,0.025,0.975)) # quartile
t2 <- quantile(master.c.inc.c$rinc.noart[master.c.inc.c$end.noART==y2.end & master.c.inc.c$Gender==2], probs = c(0.5,0.025,0.975)) # quartile
t3 <- quantile(master.c.inc.c$rinc.noart[master.c.inc.c$end.noART==y3.end & master.c.inc.c$Gender==2], probs = c(0.5,0.025,0.975)) # quartile
t1.m <- quantile(master.c.inc.c$rinc.noart[master.c.inc.c$end.noART==y1.end & master.c.inc.c$Gender==0], probs = c(0.5,0.025,0.975)) # quartile
t2.m <- quantile(master.c.inc.c$rinc.noart[master.c.inc.c$end.noART==y2.end & master.c.inc.c$Gender==0], probs = c(0.5,0.025,0.975)) # quartile
t3.m <- quantile(master.c.inc.c$rinc.noart[master.c.inc.c$end.noART==y3.end & master.c.inc.c$Gender==0], probs = c(0.5,0.025,0.975)) # quartile
t1.f <- quantile(master.c.inc.c$rinc.noart[master.c.inc.c$end.noART==y1.end & master.c.inc.c$Gender==1], probs = c(0.5,0.025,0.975)) # quartile
t2.f <- quantile(master.c.inc.c$rinc.noart[master.c.inc.c$end.noART==y2.end & master.c.inc.c$Gender==1], probs = c(0.5,0.025,0.975)) # quartile
t3.f <- quantile(master.c.inc.c$rinc.noart[master.c.inc.c$end.noART==y3.end & master.c.inc.c$Gender==1], probs = c(0.5,0.025,0.975)) # quartile
noart <- as.data.frame(rbind(t1,t2,t3,t1.f,t2.f,t3.f,t1.m,t2.m,t3.m))
noart$startyear <- c(y1.start,y2.start,y3.start,y1.start,y2.start,y3.start,y1.start,y2.start,y3.start)
noart$endyear <- c(y1.end,y2.end,y3.end,y1.end,y2.end,y3.end,y1.end,y2.end,y3.end)
noart$scenario <- "VMMConly"
t1 <- quantile(master.c.inc.c$rinc.909090[master.c.inc.c$end.noART==y1.end & master.c.inc.c$Gender==2], probs = c(0.5,0.025,0.975)) # quartile
t2 <- quantile(master.c.inc.c$rinc.909090[master.c.inc.c$end.noART==y2.end & master.c.inc.c$Gender==2], probs = c(0.5,0.025,0.975)) # quartile
t3 <- quantile(master.c.inc.c$rinc.909090[master.c.inc.c$end.noART==y3.end & master.c.inc.c$Gender==2], probs = c(0.5,0.025,0.975)) # quartile
t1.m <- quantile(master.c.inc.c$rinc.909090[master.c.inc.c$end.noART==y1.end & master.c.inc.c$Gender==0], probs = c(0.5,0.025,0.975)) # quartile
t2.m <- quantile(master.c.inc.c$rinc.909090[master.c.inc.c$end.noART==y2.end & master.c.inc.c$Gender==0], probs = c(0.5,0.025,0.975)) # quartile
t3.m <- quantile(master.c.inc.c$rinc.909090[master.c.inc.c$end.noART==y3.end & master.c.inc.c$Gender==0], probs = c(0.5,0.025,0.975)) # quartile
t1.f <- quantile(master.c.inc.c$rinc.909090[master.c.inc.c$end.noART==y1.end & master.c.inc.c$Gender==1], probs = c(0.5,0.025,0.975)) # quartile
t2.f <- quantile(master.c.inc.c$rinc.909090[master.c.inc.c$end.noART==y2.end & master.c.inc.c$Gender==1], probs = c(0.5,0.025,0.975)) # quartile
t3.f <- quantile(master.c.inc.c$rinc.909090[master.c.inc.c$end.noART==y3.end & master.c.inc.c$Gender==1], probs = c(0.5,0.025,0.975)) # quartile
age909090 <- as.data.frame(rbind(t1,t2,t3,t1.f,t2.f,t3.f,t1.m,t2.m,t3.m))
age909090$startyear <- c(y1.start,y2.start,y3.start,y1.start,y2.start,y3.start,y1.start,y2.start,y3.start)
age909090$endyear <- c(y1.end,y2.end,y3.end,y1.end,y2.end,y3.end,y1.end,y2.end,y3.end)
age909090$scenario <- "909090"
t1 <- quantile(master.c.inc.c$rinc.ART100pct[master.c.inc.c$end.noART==y1.end & master.c.inc.c$Gender==2], probs = c(0.5,0.025,0.975)) # quartile
t2 <- quantile(master.c.inc.c$rinc.ART100pct[master.c.inc.c$end.noART==y2.end & master.c.inc.c$Gender==2], probs = c(0.5,0.025,0.975)) # quartile
t3 <- quantile(master.c.inc.c$rinc.ART100pct[master.c.inc.c$end.noART==y3.end & master.c.inc.c$Gender==2], probs = c(0.5,0.025,0.975)) # quartile
t1.m <- quantile(master.c.inc.c$rinc.ART100pct[master.c.inc.c$end.noART==y1.end & master.c.inc.c$Gender==0], probs = c(0.5,0.025,0.975)) # quartile
t2.m <- quantile(master.c.inc.c$rinc.ART100pct[master.c.inc.c$end.noART==y2.end & master.c.inc.c$Gender==0], probs = c(0.5,0.025,0.975)) # quartile
t3.m <- quantile(master.c.inc.c$rinc.ART100pct[master.c.inc.c$end.noART==y3.end & master.c.inc.c$Gender==0], probs = c(0.5,0.025,0.975)) # quartile
t1.f <- quantile(master.c.inc.c$rinc.ART100pct[master.c.inc.c$end.noART==y1.end & master.c.inc.c$Gender==1], probs = c(0.5,0.025,0.975)) # quartile
t2.f <- quantile(master.c.inc.c$rinc.ART100pct[master.c.inc.c$end.noART==y2.end & master.c.inc.c$Gender==1], probs = c(0.5,0.025,0.975)) # quartile
t3.f <- quantile(master.c.inc.c$rinc.ART100pct[master.c.inc.c$end.noART==y3.end & master.c.inc.c$Gender==1], probs = c(0.5,0.025,0.975)) # quartile
ART100pct <- as.data.frame(rbind(t1,t2,t3,t1.f,t2.f,t3.f,t1.m,t2.m,t3.m))
ART100pct$startyear <- c(y1.start,y2.start,y3.start,y1.start,y2.start,y3.start,y1.start,y2.start,y3.start)
ART100pct$endyear <- c(y1.end,y2.end,y3.end,y1.end,y2.end,y3.end,y1.end,y2.end,y3.end)
ART100pct$scenario <- "ART100pct"

results.final.cumincidence <- rbind(novmmcnoart, noart, age909090 ,ART100pct)
results.final.cumincidence$Gender <- c("All","All","All","Female","Female","Female","Male","Male","Male","All","All","All","Female","Female","Female","Male","Male","Male","All","All","All","Female","Female","Female","Male","Male","Male","All","All","All","Female","Female","Female","Male","Male","Male")

##########################################################################
#Number of deaths averted - calculates deaths over different time horizons
#and divides by number of person-years (of entire population) over that time
##########################################################################
reporthivbyageandgender.master$Uninfected.Population = reporthivbyageandgender.master$Population-reporthivbyageandgender.master$Infected

#set years for time horizon over which cumulative incidence will be calculated
y1.start = 2010
y2.start = 2016
y3.start = 2016
y1.end = 2016
y2.end = 2030
y3.end = 2050

head(reporthivbyageandgender.master)
new.infections.1 <- aggregate(Died_from_HIV ~ sim.id+scenario, subset(reporthivbyageandgender.master, Year>y1.start & Year<y1.end+1.5), FUN=sum)
new.infections.2 <- aggregate(Died_from_HIV ~ sim.id+scenario, subset(reporthivbyageandgender.master, Year>y2.start & Year<y2.end+1.5), FUN=sum)
new.infections.3 <- aggregate(Died_from_HIV ~ sim.id+scenario, subset(reporthivbyageandgender.master, Year>y3.start & Year<y3.end+1.5), FUN=sum)
new.infections.1.g <- aggregate(Died_from_HIV ~ Gender+sim.id+scenario, subset(reporthivbyageandgender.master, Year>y1.start & Year<y1.end+1.5), FUN=sum)
new.infections.2.g <- aggregate(Died_from_HIV ~ Gender+sim.id+scenario, subset(reporthivbyageandgender.master, Year>y2.start & Year<y2.end+1.5), FUN=sum)
new.infections.3.g <- aggregate(Died_from_HIV ~ Gender+sim.id+scenario, subset(reporthivbyageandgender.master, Year>y3.start & Year<y3.end+1.5), FUN=sum)

pop.uninfected.1 <- aggregate(Population ~ sim.id+scenario, subset(reporthivbyageandgender.master,Year>y1.start-0.5 & Year<y1.end+1), FUN=sum) 
pop.uninfected.2 <- aggregate(Population ~ sim.id+scenario, subset(reporthivbyageandgender.master,Year>y2.start-0.5 & Year<y2.end+1), FUN=sum) 
pop.uninfected.3 <- aggregate(Population ~ sim.id+scenario, subset(reporthivbyageandgender.master,Year>y3.start-0.5 & Year<y3.end+1), FUN=sum) 
pop.uninfected.1.g <- aggregate(Population ~ Gender+sim.id+scenario, subset(reporthivbyageandgender.master,Year>y1.start-0.5 & Year<y1.end+1), FUN=sum) 
pop.uninfected.2.g <- aggregate(Population ~ Gender+sim.id+scenario, subset(reporthivbyageandgender.master,Year>y2.start-0.5 & Year<y2.end+1), FUN=sum) 
pop.uninfected.3.g <- aggregate(Population ~ Gender+sim.id+scenario, subset(reporthivbyageandgender.master,Year>y3.start-0.5 & Year<y3.end+1), FUN=sum) 

cumulative.incidence.1 <- merge(new.infections.1, pop.uninfected.1, by=c("sim.id","scenario"))
cumulative.incidence.2 <- merge(new.infections.2, pop.uninfected.2, by=c("sim.id","scenario"))
cumulative.incidence.3 <- merge(new.infections.3, pop.uninfected.3, by=c("sim.id","scenario"))
cumulative.incidence.1.g <- merge(new.infections.1.g, pop.uninfected.1.g, by=c("Gender","sim.id","scenario"))
cumulative.incidence.2.g <- merge(new.infections.2.g, pop.uninfected.2.g, by=c("Gender","sim.id","scenario"))
cumulative.incidence.3.g <- merge(new.infections.3.g, pop.uninfected.3.g, by=c("Gender","sim.id","scenario"))

head(cumulative.incidence.3.g)

cumulative.incidence.1$c.inc <- cumulative.incidence.1$Died_from_HIV/cumulative.incidence.1$Population
cumulative.incidence.2$c.inc <- cumulative.incidence.2$Died_from_HIV/cumulative.incidence.2$Population
cumulative.incidence.3$c.inc <- cumulative.incidence.3$Died_from_HIV/cumulative.incidence.3$Population
cumulative.incidence.1.g$c.inc <- cumulative.incidence.1.g$Died_from_HIV/cumulative.incidence.1.g$Population
cumulative.incidence.2.g$c.inc <- cumulative.incidence.2.g$Died_from_HIV/cumulative.incidence.2.g$Population
cumulative.incidence.3.g$c.inc <- cumulative.incidence.3.g$Died_from_HIV/cumulative.incidence.3.g$Population
cumulative.incidence.1$start <- y1.start
cumulative.incidence.1$end <- y1.end
cumulative.incidence.2$start <- y2.start
cumulative.incidence.2$end <- y2.end
cumulative.incidence.3$start <- y3.start
cumulative.incidence.3$end <- y3.end
cumulative.incidence.1.g$start <- y1.start
cumulative.incidence.1.g$end <- y1.end
cumulative.incidence.2.g$start <- y2.start
cumulative.incidence.2.g$end <- y2.end
cumulative.incidence.3.g$start <- y3.start
cumulative.incidence.3.g$end <- y3.end

w.1 <- reshape(cumulative.incidence.1, timevar = "scenario", idvar = "sim.id", direction = "wide")
w.2 <- reshape(cumulative.incidence.2, timevar = "scenario", idvar = "sim.id", direction = "wide")
w.3 <- reshape(cumulative.incidence.3, timevar = "scenario", idvar = "sim.id", direction = "wide")
w.1.g <- reshape(cumulative.incidence.1.g, timevar = "scenario", idvar = c("sim.id","Gender") , direction = "wide")
w.2.g <- reshape(cumulative.incidence.2.g, timevar = "scenario", idvar = c("sim.id","Gender") , direction = "wide")
w.3.g <- reshape(cumulative.incidence.3.g, timevar = "scenario", idvar = c("sim.id","Gender") , direction = "wide")

master.c.inc <- rbind(w.1,w.2,w.3)
master.c.inc$Gender <- 2
master.c.inc.g <- rbind(w.1.g,w.2.g,w.3.g)
master.c.inc.c <- rbind(master.c.inc, master.c.inc.g)
names(master.c.inc.c)

#calculate ratio (AR%)
master.c.inc.c$rinc.noartnovmmc <- (1 - (master.c.inc.c$c.inc.baseline / master.c.inc.c$c.inc.noART))*100
master.c.inc.c$rinc.noart <- (1 - (master.c.inc.c$c.inc.noART / master.c.inc.c$c.inc.noARTnoVMMC))*100
master.c.inc.c$rinc.909090 <- (1 - (master.c.inc.c$c.inc.909090 / master.c.inc.c$c.inc.noART))*100
master.c.inc.c$rinc.ART100pct <- (1 - (master.c.inc.c$c.inc.ART100pct / master.c.inc.c$c.inc.noART))*100
head(master.c.inc.c)

t1 <- quantile(master.c.inc.c$rinc.noartnovmmc[master.c.inc.c$end.noART==y1.end & master.c.inc.c$Gender==2], probs = c(0.5,0.025,0.975)) # quartile
t2 <- quantile(master.c.inc.c$rinc.noartnovmmc[master.c.inc.c$end.noART==y2.end & master.c.inc.c$Gender==2], probs = c(0.5,0.025,0.975)) # quartile
t3 <- quantile(master.c.inc.c$rinc.noartnovmmc[master.c.inc.c$end.noART==y3.end & master.c.inc.c$Gender==2], probs = c(0.5,0.025,0.975)) # quartile
t1.m <- quantile(master.c.inc.c$rinc.noartnovmmc[master.c.inc.c$end.noART==y1.end & master.c.inc.c$Gender==0], probs = c(0.5,0.025,0.975)) # quartile
t2.m <- quantile(master.c.inc.c$rinc.noartnovmmc[master.c.inc.c$end.noART==y2.end & master.c.inc.c$Gender==0], probs = c(0.5,0.025,0.975)) # quartile
t3.m <- quantile(master.c.inc.c$rinc.noartnovmmc[master.c.inc.c$end.noART==y3.end & master.c.inc.c$Gender==0], probs = c(0.5,0.025,0.975)) # quartile
t1.f <- quantile(master.c.inc.c$rinc.noartnovmmc[master.c.inc.c$end.noART==y1.end & master.c.inc.c$Gender==1], probs = c(0.5,0.025,0.975)) # quartile
t2.f <- quantile(master.c.inc.c$rinc.noartnovmmc[master.c.inc.c$end.noART==y2.end & master.c.inc.c$Gender==1], probs = c(0.5,0.025,0.975)) # quartile
t3.f <- quantile(master.c.inc.c$rinc.noartnovmmc[master.c.inc.c$end.noART==y3.end & master.c.inc.c$Gender==1], probs = c(0.5,0.025,0.975)) # quartile
novmmcnoart <- as.data.frame(rbind(t1,t2,t3,t1.f,t2.f,t3.f,t1.m,t2.m,t3.m))
novmmcnoart$startyear <- c(y1.start,y2.start,y3.start,y1.start,y2.start,y3.start,y1.start,y2.start,y3.start)
novmmcnoart$endyear <- c(y1.end,y2.end,y3.end,y1.end,y2.end,y3.end,y1.end,y2.end,y3.end)
novmmcnoart$scenario <- "ART+VMMC"
t1 <- quantile(master.c.inc.c$rinc.noart[master.c.inc.c$end.noART==y1.end & master.c.inc.c$Gender==2], probs = c(0.5,0.025,0.975)) # quartile
t2 <- quantile(master.c.inc.c$rinc.noart[master.c.inc.c$end.noART==y2.end & master.c.inc.c$Gender==2], probs = c(0.5,0.025,0.975)) # quartile
t3 <- quantile(master.c.inc.c$rinc.noart[master.c.inc.c$end.noART==y3.end & master.c.inc.c$Gender==2], probs = c(0.5,0.025,0.975)) # quartile
t1.m <- quantile(master.c.inc.c$rinc.noart[master.c.inc.c$end.noART==y1.end & master.c.inc.c$Gender==0], probs = c(0.5,0.025,0.975)) # quartile
t2.m <- quantile(master.c.inc.c$rinc.noart[master.c.inc.c$end.noART==y2.end & master.c.inc.c$Gender==0], probs = c(0.5,0.025,0.975)) # quartile
t3.m <- quantile(master.c.inc.c$rinc.noart[master.c.inc.c$end.noART==y3.end & master.c.inc.c$Gender==0], probs = c(0.5,0.025,0.975)) # quartile
t1.f <- quantile(master.c.inc.c$rinc.noart[master.c.inc.c$end.noART==y1.end & master.c.inc.c$Gender==1], probs = c(0.5,0.025,0.975)) # quartile
t2.f <- quantile(master.c.inc.c$rinc.noart[master.c.inc.c$end.noART==y2.end & master.c.inc.c$Gender==1], probs = c(0.5,0.025,0.975)) # quartile
t3.f <- quantile(master.c.inc.c$rinc.noart[master.c.inc.c$end.noART==y3.end & master.c.inc.c$Gender==1], probs = c(0.5,0.025,0.975)) # quartile
noart <- as.data.frame(rbind(t1,t2,t3,t1.f,t2.f,t3.f,t1.m,t2.m,t3.m))
noart$startyear <- c(y1.start,y2.start,y3.start,y1.start,y2.start,y3.start,y1.start,y2.start,y3.start)
noart$endyear <- c(y1.end,y2.end,y3.end,y1.end,y2.end,y3.end,y1.end,y2.end,y3.end)
noart$scenario <- "VMMConly"
t1 <- quantile(master.c.inc.c$rinc.909090[master.c.inc.c$end.noART==y1.end & master.c.inc.c$Gender==2], probs = c(0.5,0.025,0.975)) # quartile
t2 <- quantile(master.c.inc.c$rinc.909090[master.c.inc.c$end.noART==y2.end & master.c.inc.c$Gender==2], probs = c(0.5,0.025,0.975)) # quartile
t3 <- quantile(master.c.inc.c$rinc.909090[master.c.inc.c$end.noART==y3.end & master.c.inc.c$Gender==2], probs = c(0.5,0.025,0.975)) # quartile
t1.m <- quantile(master.c.inc.c$rinc.909090[master.c.inc.c$end.noART==y1.end & master.c.inc.c$Gender==0], probs = c(0.5,0.025,0.975)) # quartile
t2.m <- quantile(master.c.inc.c$rinc.909090[master.c.inc.c$end.noART==y2.end & master.c.inc.c$Gender==0], probs = c(0.5,0.025,0.975)) # quartile
t3.m <- quantile(master.c.inc.c$rinc.909090[master.c.inc.c$end.noART==y3.end & master.c.inc.c$Gender==0], probs = c(0.5,0.025,0.975)) # quartile
t1.f <- quantile(master.c.inc.c$rinc.909090[master.c.inc.c$end.noART==y1.end & master.c.inc.c$Gender==1], probs = c(0.5,0.025,0.975)) # quartile
t2.f <- quantile(master.c.inc.c$rinc.909090[master.c.inc.c$end.noART==y2.end & master.c.inc.c$Gender==1], probs = c(0.5,0.025,0.975)) # quartile
t3.f <- quantile(master.c.inc.c$rinc.909090[master.c.inc.c$end.noART==y3.end & master.c.inc.c$Gender==1], probs = c(0.5,0.025,0.975)) # quartile
age909090 <- as.data.frame(rbind(t1,t2,t3,t1.f,t2.f,t3.f,t1.m,t2.m,t3.m))
age909090$startyear <- c(y1.start,y2.start,y3.start,y1.start,y2.start,y3.start,y1.start,y2.start,y3.start)
age909090$endyear <- c(y1.end,y2.end,y3.end,y1.end,y2.end,y3.end,y1.end,y2.end,y3.end)
age909090$scenario <- "909090"
t1 <- quantile(master.c.inc.c$rinc.ART100pct[master.c.inc.c$end.noART==y1.end & master.c.inc.c$Gender==2], probs = c(0.5,0.025,0.975)) # quartile
t2 <- quantile(master.c.inc.c$rinc.ART100pct[master.c.inc.c$end.noART==y2.end & master.c.inc.c$Gender==2], probs = c(0.5,0.025,0.975)) # quartile
t3 <- quantile(master.c.inc.c$rinc.ART100pct[master.c.inc.c$end.noART==y3.end & master.c.inc.c$Gender==2], probs = c(0.5,0.025,0.975)) # quartile
t1.m <- quantile(master.c.inc.c$rinc.ART100pct[master.c.inc.c$end.noART==y1.end & master.c.inc.c$Gender==0], probs = c(0.5,0.025,0.975)) # quartile
t2.m <- quantile(master.c.inc.c$rinc.ART100pct[master.c.inc.c$end.noART==y2.end & master.c.inc.c$Gender==0], probs = c(0.5,0.025,0.975)) # quartile
t3.m <- quantile(master.c.inc.c$rinc.ART100pct[master.c.inc.c$end.noART==y3.end & master.c.inc.c$Gender==0], probs = c(0.5,0.025,0.975)) # quartile
t1.f <- quantile(master.c.inc.c$rinc.ART100pct[master.c.inc.c$end.noART==y1.end & master.c.inc.c$Gender==1], probs = c(0.5,0.025,0.975)) # quartile
t2.f <- quantile(master.c.inc.c$rinc.ART100pct[master.c.inc.c$end.noART==y2.end & master.c.inc.c$Gender==1], probs = c(0.5,0.025,0.975)) # quartile
t3.f <- quantile(master.c.inc.c$rinc.ART100pct[master.c.inc.c$end.noART==y3.end & master.c.inc.c$Gender==1], probs = c(0.5,0.025,0.975)) # quartile
ART100pct <- as.data.frame(rbind(t1,t2,t3,t1.f,t2.f,t3.f,t1.m,t2.m,t3.m))
ART100pct$startyear <- c(y1.start,y2.start,y3.start,y1.start,y2.start,y3.start,y1.start,y2.start,y3.start)
ART100pct$endyear <- c(y1.end,y2.end,y3.end,y1.end,y2.end,y3.end,y1.end,y2.end,y3.end)
ART100pct$scenario <- "ART100pct"

results.final.mortality <- rbind(novmmcnoart, noart, age909090 ,ART100pct)
results.final.mortality$Gender <- c("All","All","All","Female","Female","Female","Male","Male","Male","All","All","All","Female","Female","Female","Male","Male","Male","All","All","All","Female","Female","Female","Male","Male","Male","All","All","All","Female","Female","Female","Male","Male","Male")

###########################################################################################################
#Table 2 for paper
###########################################################################################################
#set years for time in which to calculate incidence. ignore start years, only there for consistent headers
y1.start = 2010
y2.start = 2016
y3.start = 2016
y1.end = 2016
y2.end = 2030
y3.end = 2050

#add up new infections by age for each year + by gender
new.infections.1 <- aggregate(Newly.Infected ~ Agecat+sim.id+scenario, subset(reporthivbyageandgender.master, Year==y1.end+0.5 | Year==y1.end+1), FUN=sum)
new.infections.2 <- aggregate(Newly.Infected ~ Agecat+sim.id+scenario, subset(reporthivbyageandgender.master, Year==y2.end+0.5 | Year==y2.end+1), FUN=sum)
new.infections.3 <- aggregate(Newly.Infected ~ Agecat+sim.id+scenario, subset(reporthivbyageandgender.master, Year==y3.end+0.5 | Year==y3.end+1), FUN=sum)
new.infections.1.g <- aggregate(Newly.Infected ~ Agecat+Gender+sim.id+scenario, subset(reporthivbyageandgender.master, Year==y1.end+0.5 | Year==y1.end+1), FUN=sum)
new.infections.2.g <- aggregate(Newly.Infected ~ Agecat+Gender+sim.id+scenario, subset(reporthivbyageandgender.master, Year==y2.end+0.5 | Year==y2.end+1), FUN=sum)
new.infections.3.g <- aggregate(Newly.Infected ~ Agecat+Gender+sim.id+scenario, subset(reporthivbyageandgender.master, Year==y3.end+0.5 | Year==y3.end+1), FUN=sum)

#add up new infections overall and by gender
pop.uninfected.1 <- aggregate(Newly.Infected ~ sim.id+scenario, subset(reporthivbyageandgender.master,Year==y1.end+0.5 | Year==y1.end+1), FUN=sum) 
pop.uninfected.2 <- aggregate(Newly.Infected ~ sim.id+scenario, subset(reporthivbyageandgender.master,Year==y2.end+0.5 | Year==y2.end+1), FUN=sum) 
pop.uninfected.3 <- aggregate(Newly.Infected ~ sim.id+scenario, subset(reporthivbyageandgender.master,Year==y3.end+0.5 | Year==y3.end+1), FUN=sum) 
pop.uninfected.1.g <- aggregate(Newly.Infected ~ Gender+sim.id+scenario, subset(reporthivbyageandgender.master,Year==y1.end+0.5 | Year==y1.end+1), FUN=sum) 
pop.uninfected.2.g <- aggregate(Newly.Infected ~ Gender+sim.id+scenario, subset(reporthivbyageandgender.master,Year==y2.end+0.5 | Year==y2.end+1), FUN=sum) 
pop.uninfected.3.g <- aggregate(Newly.Infected ~ Gender+sim.id+scenario, subset(reporthivbyageandgender.master,Year==y3.end+0.5 | Year==y3.end+1), FUN=sum) 

cumulative.incidence.1 <- merge(new.infections.1, pop.uninfected.1, by=c("sim.id","scenario"))
cumulative.incidence.2 <- merge(new.infections.2, pop.uninfected.2, by=c("sim.id","scenario"))
cumulative.incidence.3 <- merge(new.infections.3, pop.uninfected.3, by=c("sim.id","scenario"))
cumulative.incidence.1.g <- merge(new.infections.1.g, pop.uninfected.1.g, by=c("Gender","sim.id","scenario"))
cumulative.incidence.2.g <- merge(new.infections.2.g, pop.uninfected.2.g, by=c("Gender","sim.id","scenario"))
cumulative.incidence.3.g <- merge(new.infections.3.g, pop.uninfected.3.g, by=c("Gender","sim.id","scenario"))

#calculates proportion of new infections from each age group
cumulative.incidence.1$c.inc <- cumulative.incidence.1$Newly.Infected.x/cumulative.incidence.1$Newly.Infected.y
cumulative.incidence.2$c.inc <- cumulative.incidence.2$Newly.Infected.x/cumulative.incidence.2$Newly.Infected.y
cumulative.incidence.3$c.inc <- cumulative.incidence.3$Newly.Infected.x/cumulative.incidence.3$Newly.Infected.y
cumulative.incidence.1.g$c.inc <- cumulative.incidence.1.g$Newly.Infected.x/cumulative.incidence.1.g$Newly.Infected.y
cumulative.incidence.2.g$c.inc <- cumulative.incidence.2.g$Newly.Infected.x/cumulative.incidence.2.g$Newly.Infected.y
cumulative.incidence.3.g$c.inc <- cumulative.incidence.3.g$Newly.Infected.x/cumulative.incidence.3.g$Newly.Infected.y
cumulative.incidence.1$start <- y1.start
cumulative.incidence.1$end <- y1.end
cumulative.incidence.2$start <- y2.start
cumulative.incidence.2$end <- y2.end
cumulative.incidence.3$start <- y3.start
cumulative.incidence.3$end <- y3.end
cumulative.incidence.1.g$start <- y1.start
cumulative.incidence.1.g$end <- y1.end
cumulative.incidence.2.g$start <- y2.start
cumulative.incidence.2.g$end <- y2.end
cumulative.incidence.3.g$start <- y3.start
cumulative.incidence.3.g$end <- y3.end

w.1 <- reshape(cumulative.incidence.1, timevar = "scenario", idvar = c("Agecat", "sim.id"), direction = "wide")
w.2 <- reshape(cumulative.incidence.2, timevar = "scenario", idvar = c("Agecat", "sim.id"), direction = "wide")
w.3 <- reshape(cumulative.incidence.3, timevar = "scenario", idvar = c("Agecat", "sim.id"), direction = "wide")
w.1.g <- reshape(cumulative.incidence.1.g, timevar = "scenario", idvar = c("Agecat","sim.id","Gender") , direction = "wide")
w.2.g <- reshape(cumulative.incidence.2.g, timevar = "scenario", idvar = c("Agecat","sim.id","Gender") , direction = "wide")
w.3.g <- reshape(cumulative.incidence.3.g, timevar = "scenario", idvar = c("Agecat","sim.id","Gender") , direction = "wide")

master.c.inc <- rbind(w.1,w.2,w.3)
master.c.inc$Gender <- 2
master.c.inc.g <- rbind(w.1.g,w.2.g,w.3.g)
master.c.inc.c <- rbind(master.c.inc, master.c.inc.g)
names(master.c.inc.c)

#calculate proportion of new cases in each age category, here you just copy to new variable to be consistent 
master.c.inc.c$rinc.baseline <- master.c.inc.c$c.inc.baseline
master.c.inc.c$rinc.noartnovmmc <- master.c.inc.c$c.inc.noARTnoVMMC
master.c.inc.c$rinc.noart <- master.c.inc.c$c.inc.noART
master.c.inc.c$rinc.909090 <- master.c.inc.c$c.inc.909090
master.c.inc.c$rinc.ART100pct <- master.c.inc.c$c.inc.ART100pct
head(master.c.inc.c)

#calculate the quantiles among 15-24 year olds to see how the % of cases in that age category changes over time by intervention
t1 <- quantile(master.c.inc.c$rinc.baseline[master.c.inc.c$end.noART==y1.end & master.c.inc.c$Gender==2 & master.c.inc.c$Agecat=="[15,25)"], probs = c(0.5,0.025,0.975)) # quartile
t2 <- quantile(master.c.inc.c$rinc.baseline[master.c.inc.c$end.noART==y2.end & master.c.inc.c$Gender==2 & master.c.inc.c$Agecat=="[15,25)"], probs = c(0.5,0.025,0.975)) # quartile
t3 <- quantile(master.c.inc.c$rinc.baseline[master.c.inc.c$end.noART==y3.end & master.c.inc.c$Gender==2 & master.c.inc.c$Agecat=="[15,25)"], probs = c(0.5,0.025,0.975)) # quartile
t1.m <- quantile(master.c.inc.c$rinc.baseline[master.c.inc.c$end.noART==y1.end & master.c.inc.c$Gender==0 & master.c.inc.c$Agecat=="[15,25)"], probs = c(0.5,0.025,0.975)) # quartile
t2.m <- quantile(master.c.inc.c$rinc.baseline[master.c.inc.c$end.noART==y2.end & master.c.inc.c$Gender==0 & master.c.inc.c$Agecat=="[15,25)"], probs = c(0.5,0.025,0.975)) # quartile
t3.m <- quantile(master.c.inc.c$rinc.baseline[master.c.inc.c$end.noART==y3.end & master.c.inc.c$Gender==0 & master.c.inc.c$Agecat=="[15,25)"], probs = c(0.5,0.025,0.975)) # quartile
t1.f <- quantile(master.c.inc.c$rinc.baseline[master.c.inc.c$end.noART==y1.end & master.c.inc.c$Gender==1 & master.c.inc.c$Agecat=="[15,25)"], probs = c(0.5,0.025,0.975)) # quartile
t2.f <- quantile(master.c.inc.c$rinc.baseline[master.c.inc.c$end.noART==y2.end & master.c.inc.c$Gender==1 & master.c.inc.c$Agecat=="[15,25)"], probs = c(0.5,0.025,0.975)) # quartile
t3.f <- quantile(master.c.inc.c$rinc.baseline[master.c.inc.c$end.noART==y3.end & master.c.inc.c$Gender==1 & master.c.inc.c$Agecat=="[15,25)"], probs = c(0.5,0.025,0.975)) # quartile
baseline <- as.data.frame(rbind(t1,t2,t3,t1.f,t2.f,t3.f,t1.m,t2.m,t3.m))
baseline$startyear <- c(y1.start,y2.start,y3.start,y1.start,y2.start,y3.start,y1.start,y2.start,y3.start)
baseline$endyear <- c(y1.end,y2.end,y3.end,y1.end,y2.end,y3.end,y1.end,y2.end,y3.end)
baseline$scenario <- "ART+VMMC"
t1 <- quantile(master.c.inc.c$rinc.noartnovmmc[master.c.inc.c$end.noART==y1.end & master.c.inc.c$Gender==2 & master.c.inc.c$Agecat=="[15,25)"], probs = c(0.5,0.025,0.975)) # quartile
t2 <- quantile(master.c.inc.c$rinc.noartnovmmc[master.c.inc.c$end.noART==y2.end & master.c.inc.c$Gender==2 & master.c.inc.c$Agecat=="[15,25)"], probs = c(0.5,0.025,0.975)) # quartile
t3 <- quantile(master.c.inc.c$rinc.noartnovmmc[master.c.inc.c$end.noART==y3.end & master.c.inc.c$Gender==2 & master.c.inc.c$Agecat=="[15,25)"], probs = c(0.5,0.025,0.975)) # quartile
t1.m <- quantile(master.c.inc.c$rinc.noartnovmmc[master.c.inc.c$end.noART==y1.end & master.c.inc.c$Gender==0 & master.c.inc.c$Agecat=="[15,25)"], probs = c(0.5,0.025,0.975)) # quartile
t2.m <- quantile(master.c.inc.c$rinc.noartnovmmc[master.c.inc.c$end.noART==y2.end & master.c.inc.c$Gender==0 & master.c.inc.c$Agecat=="[15,25)"], probs = c(0.5,0.025,0.975)) # quartile
t3.m <- quantile(master.c.inc.c$rinc.noartnovmmc[master.c.inc.c$end.noART==y3.end & master.c.inc.c$Gender==0 & master.c.inc.c$Agecat=="[15,25)"], probs = c(0.5,0.025,0.975)) # quartile
t1.f <- quantile(master.c.inc.c$rinc.noartnovmmc[master.c.inc.c$end.noART==y1.end & master.c.inc.c$Gender==1 & master.c.inc.c$Agecat=="[15,25)"], probs = c(0.5,0.025,0.975)) # quartile
t2.f <- quantile(master.c.inc.c$rinc.noartnovmmc[master.c.inc.c$end.noART==y2.end & master.c.inc.c$Gender==1 & master.c.inc.c$Agecat=="[15,25)"], probs = c(0.5,0.025,0.975)) # quartile
t3.f <- quantile(master.c.inc.c$rinc.noartnovmmc[master.c.inc.c$end.noART==y3.end & master.c.inc.c$Gender==1 & master.c.inc.c$Agecat=="[15,25)"], probs = c(0.5,0.025,0.975)) # quartile
novmmcnoart <- as.data.frame(rbind(t1,t2,t3,t1.f,t2.f,t3.f,t1.m,t2.m,t3.m))
novmmcnoart$startyear <- c(y1.start,y2.start,y3.start,y1.start,y2.start,y3.start,y1.start,y2.start,y3.start)
novmmcnoart$endyear <- c(y1.end,y2.end,y3.end,y1.end,y2.end,y3.end,y1.end,y2.end,y3.end)
novmmcnoart$scenario <- "noVMMCnoART"
t1 <- quantile(master.c.inc.c$rinc.noart[master.c.inc.c$end.noART==y1.end & master.c.inc.c$Gender==2 & master.c.inc.c$Agecat=="[15,25)"], probs = c(0.5,0.025,0.975)) # quartile
t2 <- quantile(master.c.inc.c$rinc.noart[master.c.inc.c$end.noART==y2.end & master.c.inc.c$Gender==2 & master.c.inc.c$Agecat=="[15,25)"], probs = c(0.5,0.025,0.975)) # quartile
t3 <- quantile(master.c.inc.c$rinc.noart[master.c.inc.c$end.noART==y3.end & master.c.inc.c$Gender==2 & master.c.inc.c$Agecat=="[15,25)"], probs = c(0.5,0.025,0.975)) # quartile
t1.m <- quantile(master.c.inc.c$rinc.noart[master.c.inc.c$end.noART==y1.end & master.c.inc.c$Gender==0 & master.c.inc.c$Agecat=="[15,25)"], probs = c(0.5,0.025,0.975)) # quartile
t2.m <- quantile(master.c.inc.c$rinc.noart[master.c.inc.c$end.noART==y2.end & master.c.inc.c$Gender==0 & master.c.inc.c$Agecat=="[15,25)"], probs = c(0.5,0.025,0.975)) # quartile
t3.m <- quantile(master.c.inc.c$rinc.noart[master.c.inc.c$end.noART==y3.end & master.c.inc.c$Gender==0 & master.c.inc.c$Agecat=="[15,25)"], probs = c(0.5,0.025,0.975)) # quartile
t1.f <- quantile(master.c.inc.c$rinc.noart[master.c.inc.c$end.noART==y1.end & master.c.inc.c$Gender==1 & master.c.inc.c$Agecat=="[15,25)"], probs = c(0.5,0.025,0.975)) # quartile
t2.f <- quantile(master.c.inc.c$rinc.noart[master.c.inc.c$end.noART==y2.end & master.c.inc.c$Gender==1 & master.c.inc.c$Agecat=="[15,25)"], probs = c(0.5,0.025,0.975)) # quartile
t3.f <- quantile(master.c.inc.c$rinc.noart[master.c.inc.c$end.noART==y3.end & master.c.inc.c$Gender==1 & master.c.inc.c$Agecat=="[15,25)"], probs = c(0.5,0.025,0.975)) # quartile
noart <- as.data.frame(rbind(t1,t2,t3,t1.f,t2.f,t3.f,t1.m,t2.m,t3.m))
noart$startyear <- c(y1.start,y2.start,y3.start,y1.start,y2.start,y3.start,y1.start,y2.start,y3.start)
noart$endyear <- c(y1.end,y2.end,y3.end,y1.end,y2.end,y3.end,y1.end,y2.end,y3.end)
noart$scenario <- "VMMConly"
t1 <- quantile(master.c.inc.c$rinc.909090[master.c.inc.c$end.noART==y1.end & master.c.inc.c$Gender==2 & master.c.inc.c$Agecat=="[15,25)"], probs = c(0.5,0.025,0.975)) # quartile
t2 <- quantile(master.c.inc.c$rinc.909090[master.c.inc.c$end.noART==y2.end & master.c.inc.c$Gender==2 & master.c.inc.c$Agecat=="[15,25)"], probs = c(0.5,0.025,0.975)) # quartile
t3 <- quantile(master.c.inc.c$rinc.909090[master.c.inc.c$end.noART==y3.end & master.c.inc.c$Gender==2 & master.c.inc.c$Agecat=="[15,25)"], probs = c(0.5,0.025,0.975)) # quartile
t1.m <- quantile(master.c.inc.c$rinc.909090[master.c.inc.c$end.noART==y1.end & master.c.inc.c$Gender==0 & master.c.inc.c$Agecat=="[15,25)"], probs = c(0.5,0.025,0.975)) # quartile
t2.m <- quantile(master.c.inc.c$rinc.909090[master.c.inc.c$end.noART==y2.end & master.c.inc.c$Gender==0 & master.c.inc.c$Agecat=="[15,25)"], probs = c(0.5,0.025,0.975)) # quartile
t3.m <- quantile(master.c.inc.c$rinc.909090[master.c.inc.c$end.noART==y3.end & master.c.inc.c$Gender==0 & master.c.inc.c$Agecat=="[15,25)"], probs = c(0.5,0.025,0.975)) # quartile
t1.f <- quantile(master.c.inc.c$rinc.909090[master.c.inc.c$end.noART==y1.end & master.c.inc.c$Gender==1 & master.c.inc.c$Agecat=="[15,25)"], probs = c(0.5,0.025,0.975)) # quartile
t2.f <- quantile(master.c.inc.c$rinc.909090[master.c.inc.c$end.noART==y2.end & master.c.inc.c$Gender==1 & master.c.inc.c$Agecat=="[15,25)"], probs = c(0.5,0.025,0.975)) # quartile
t3.f <- quantile(master.c.inc.c$rinc.909090[master.c.inc.c$end.noART==y3.end & master.c.inc.c$Gender==1 & master.c.inc.c$Agecat=="[15,25)"], probs = c(0.5,0.025,0.975)) # quartile
age909090 <- as.data.frame(rbind(t1,t2,t3,t1.f,t2.f,t3.f,t1.m,t2.m,t3.m))
age909090$startyear <- c(y1.start,y2.start,y3.start,y1.start,y2.start,y3.start,y1.start,y2.start,y3.start)
age909090$endyear <- c(y1.end,y2.end,y3.end,y1.end,y2.end,y3.end,y1.end,y2.end,y3.end)
age909090$scenario <- "909090"
t1 <- quantile(master.c.inc.c$rinc.ART100pct[master.c.inc.c$end.noART==y1.end & master.c.inc.c$Gender==2 & master.c.inc.c$Agecat=="[15,25)"], probs = c(0.5,0.025,0.975)) # quartile
t2 <- quantile(master.c.inc.c$rinc.ART100pct[master.c.inc.c$end.noART==y2.end & master.c.inc.c$Gender==2 & master.c.inc.c$Agecat=="[15,25)"], probs = c(0.5,0.025,0.975)) # quartile
t3 <- quantile(master.c.inc.c$rinc.ART100pct[master.c.inc.c$end.noART==y3.end & master.c.inc.c$Gender==2 & master.c.inc.c$Agecat=="[15,25)"], probs = c(0.5,0.025,0.975)) # quartile
t1.m <- quantile(master.c.inc.c$rinc.ART100pct[master.c.inc.c$end.noART==y1.end & master.c.inc.c$Gender==0 & master.c.inc.c$Agecat=="[15,25)"], probs = c(0.5,0.025,0.975)) # quartile
t2.m <- quantile(master.c.inc.c$rinc.ART100pct[master.c.inc.c$end.noART==y2.end & master.c.inc.c$Gender==0 & master.c.inc.c$Agecat=="[15,25)"], probs = c(0.5,0.025,0.975)) # quartile
t3.m <- quantile(master.c.inc.c$rinc.ART100pct[master.c.inc.c$end.noART==y3.end & master.c.inc.c$Gender==0 & master.c.inc.c$Agecat=="[15,25)"], probs = c(0.5,0.025,0.975)) # quartile
t1.f <- quantile(master.c.inc.c$rinc.ART100pct[master.c.inc.c$end.noART==y1.end & master.c.inc.c$Gender==1 & master.c.inc.c$Agecat=="[15,25)"], probs = c(0.5,0.025,0.975)) # quartile
t2.f <- quantile(master.c.inc.c$rinc.ART100pct[master.c.inc.c$end.noART==y2.end & master.c.inc.c$Gender==1 & master.c.inc.c$Agecat=="[15,25)"], probs = c(0.5,0.025,0.975)) # quartile
t3.f <- quantile(master.c.inc.c$rinc.ART100pct[master.c.inc.c$end.noART==y3.end & master.c.inc.c$Gender==1 & master.c.inc.c$Agecat=="[15,25)"], probs = c(0.5,0.025,0.975)) # quartile
ART100pct <- as.data.frame(rbind(t1,t2,t3,t1.f,t2.f,t3.f,t1.m,t2.m,t3.m))
ART100pct$startyear <- c(y1.start,y2.start,y3.start,y1.start,y2.start,y3.start,y1.start,y2.start,y3.start)
ART100pct$endyear <- c(y1.end,y2.end,y3.end,y1.end,y2.end,y3.end,y1.end,y2.end,y3.end)
ART100pct$scenario <- "ART100pct"

results.final.agedist <- rbind(baseline,novmmcnoart, noart, age909090 ,ART100pct)
results.final.agedist$Gender <- c("All","All","All","Female","Female","Female","Male","Male","Male","All","All","All","Female","Female","Female","Male","Male","Male","All","All","All","Female","Female","Female","Male","Male","Male","All","All","All","Female","Female","Female","Male","Male","Male","All","All","All","Female","Female","Female","Male","Male","Male")

#combine results
results.final.mortality$metric <- "mortality"
results.final.cumincidence$metric <- "cumincidence"
results.final.incidence$metric <- "incidence"
results.final.agedist$metric <- "proportionofcasesfrom1524"
results.final <- rbind(results.final.incidence, results.final.cumincidence, results.final.mortality, results.final.agedist)
write.csv(results.final, "results.final.point.est.raw.csv")

head(results.final)
names(results.final)

results.final.w <- reshape(results.final[,c(1:3,5:8)], timevar = c("endyear"), idvar = c("scenario", "Gender","metric"), direction = "wide")
head(results.final.w)
results.final.w2 <- reshape(results.final.w, timevar = c("scenario"), idvar = c("Gender","metric"), direction = "wide")

write.csv(results.final.w, "results.final.point.est.raw.wide.csv")
write.csv(results.final.w2, "results.final.point.est.raw.wide2.csv")

###########################################################################################################
#Acute HIV infection
###########################################################################################################

#Baseline scenario
wd13 <- "C:\\Users\\aakullian\\Dropbox (IDM)\\HIV\\EMOD\\Swaziland_05June2018\\Scenarios_calibration_run3_v05Sep2018_newcampaign_FINAL\\ScenariosOct_FINAL\\Baseline_TransmissionReport\\TransmissionReport\\"
setwd(paste0(wd13))
files <- list.files(full.names = F)
head(files)
f <- paste0(wd13, files)
length(f)

for (i in seq(1,250,1)){
  transmissionreport <- read.csv(f[i])
  transmissionreport <- transmissionreport[,c(2,6,7,18,19,29)]
  transmissionreport$src_age_int <- floor(transmissionreport$SRC_AGE / 365.25)
  transmissionreport$dest_age_int <- floor(transmissionreport$DEST_AGE / 365.25)
  transmissionreport$year <- floor(transmissionreport$YEAR)
  transmissionreport$year.cat <- cut(transmissionreport$year, breaks=c(1980, 2005, 2016, 2030, 2050), right=F, dig.lab=4)
  transmissionreport$scenario <- "baseline"
  transmissionreport$sim.id <- paste0(files[i])
  transmissionreport$event <- 1
  counts <- aggregate(event ~ SRC_GENDER + src_age_int + dest_age_int + year.cat, transmissionreport, sum)
  InfectionsByYear <- aggregate(event ~ year.cat+SRC_GENDER, counts, sum)
  TransmissionMatrix <- merge(counts, InfectionsByYear, by=c("year.cat","SRC_GENDER"))
  TransmissionMatrix$prop <- TransmissionMatrix$event.x/TransmissionMatrix$event.y
  TransmissionMatrix$scenario <- "baseline"
  TransmissionMatrix$sim.id <- paste0(files[i])
  if (i==1){
    TransmissionMatrix.master <- TransmissionMatrix
    transmissionreport.master <- transmissionreport} 
  else{
    TransmissionMatrix.master <- rbind(TransmissionMatrix.master, TransmissionMatrix)
    transmissionreport.master <- rbind(transmissionreport.master, transmissionreport)
    print(paste0("Why hello there, I'm working on baseline Transmission report folder",i))
  }
}

TransmissionMatrix.master.baseline <- TransmissionMatrix.master
transmissionreport.master.baseline <- transmissionreport.master
head(transmissionreport.master.baseline)

#No ART
wd14 <- "C:\\Users\\aakullian\\Dropbox (IDM)\\HIV\\EMOD\\Swaziland_05June2018\\Scenarios_calibration_run3_v05Sep2018_newcampaign_FINAL\\ScenariosOct_FINAL\\noART_TransmissionReport\\TransmissionReport\\"
setwd(paste0(wd14))
files <- list.files(full.names = F)
head(files)
f <- paste0(wd14, files)
length(f)

for (i in seq(1,250,1)){
  transmissionreport <- read.csv(f[i])
  transmissionreport <- transmissionreport[,c(2,6,7,18,19,29)]
  transmissionreport$src_age_int <- floor(transmissionreport$SRC_AGE / 365.25)
  transmissionreport$dest_age_int <- floor(transmissionreport$DEST_AGE / 365.25)
  transmissionreport$year <- floor(transmissionreport$YEAR)
  transmissionreport$year.cat <- cut(transmissionreport$year, breaks=c(1980, 2005, 2016, 2030, 2050), right=F, dig.lab=4)
  transmissionreport$scenario <- "noART"
  transmissionreport$sim.id <- paste0(files[i])
  transmissionreport$event <- 1
  counts <- aggregate(event ~ SRC_GENDER + src_age_int + dest_age_int + year.cat, transmissionreport, sum)
  InfectionsByYear <- aggregate(event ~ year.cat+SRC_GENDER, counts, sum)
  TransmissionMatrix <- merge(counts, InfectionsByYear, by=c("year.cat","SRC_GENDER"))
  TransmissionMatrix$prop <- TransmissionMatrix$event.x/TransmissionMatrix$event.y
  TransmissionMatrix$scenario <- "noART"
  TransmissionMatrix$sim.id <- paste0(files[i])
  if (i==1){
    TransmissionMatrix.master <- TransmissionMatrix
    transmissionreport.master <- transmissionreport} 
  else{
    TransmissionMatrix.master <- rbind(TransmissionMatrix.master, TransmissionMatrix)
    transmissionreport.master <- rbind(transmissionreport.master, transmissionreport)
    print(paste0("Why hello there, I'm working on noART Transmission report folder",i))
  }
}

TransmissionMatrix.master.noART <- TransmissionMatrix.master
transmissionreport.master.noART <- transmissionreport.master

#No ART no VMMC
wd15 <- "C:\\Users\\aakullian\\Dropbox (IDM)\\HIV\\EMOD\\Swaziland_05June2018\\Scenarios_calibration_run3_v05Sep2018_newcampaign_FINAL\\ScenariosOct_FINAL\\noARTnoVMMC_TransmissionReport\\TransmissionReport\\"
setwd(paste0(wd15))
files <- list.files(full.names = F)
head(files)
f <- paste0(wd15, files)
length(f)

for (i in seq(1,250,1)){
  transmissionreport <- read.csv(f[i])
  transmissionreport <- transmissionreport[,c(2,6,7,18,19,29)]
  transmissionreport$src_age_int <- floor(transmissionreport$SRC_AGE / 365.25)
  transmissionreport$dest_age_int <- floor(transmissionreport$DEST_AGE / 365.25)
  transmissionreport$year <- floor(transmissionreport$YEAR)
  transmissionreport$year.cat <- cut(transmissionreport$year, breaks=c(1980, 2005, 2016, 2030, 2050), right=F, dig.lab=4)
  transmissionreport$scenario <- "noARTnoVMMC"
  transmissionreport$sim.id <- paste0(files[i])
  transmissionreport$event <- 1
  counts <- aggregate(event ~ SRC_GENDER + src_age_int + dest_age_int + year.cat, transmissionreport, sum)
  InfectionsByYear <- aggregate(event ~ year.cat+SRC_GENDER, counts, sum)
  TransmissionMatrix <- merge(counts, InfectionsByYear, by=c("year.cat","SRC_GENDER"))
  TransmissionMatrix$prop <- TransmissionMatrix$event.x/TransmissionMatrix$event.y
  TransmissionMatrix$scenario <- "noARTnoVMMC"
  TransmissionMatrix$sim.id <- paste0(files[i])
  if (i==1){
    TransmissionMatrix.master <- TransmissionMatrix
    transmissionreport.master <- transmissionreport} 
  else{
    TransmissionMatrix.master <- rbind(TransmissionMatrix.master, TransmissionMatrix)
    transmissionreport.master <- rbind(transmissionreport.master, transmissionreport)
    print(paste0("Why hello there, I'm working on noARTnoVMMC Transmission report folder",i))
  }
}

TransmissionMatrix.master.noARTnoVMMC <- TransmissionMatrix.master
transmissionreport.master.noARTnoVMMC <- transmissionreport.master

#Stich together multiple scnearios
transmissionreport.master <- transmissionreport.master.baseline

TransmissionMatrix.master <- rbind(TransmissionMatrix.master.noART,
                                   TransmissionMatrix.master.noARTnoVMMC,
                                   TransmissionMatrix.master.baseline)

table(transmissionreport.master$scenario)
table(TransmissionMatrix.master$scenario)

head(transmissionreport.master)
#Count all the rows by year and acute stage==1
transmissionreport.master$event <- 1
transmissionreport.master$acute <- 0
transmissionreport.master$acute[transmissionreport.master$SRC_STAGE==1] <- 1
total.events <- aggregate(event ~ year+SRC_GENDER+sim.id+scenario, subset(transmissionreport.master), sum)
acute.events <- aggregate(acute ~ year+SRC_GENDER+sim.id+scenario, subset(transmissionreport.master), sum)
transmissions <- merge(total.events, acute.events, by=c("year","SRC_GENDER","sim.id", "scenario"))
transmissions$ProportionAcute <- transmissions$acute / transmissions$event
max(transmissions$year)
head(transmissions,20)
names(transmissions)

year=2050
#Proprtion Acute plot
swaziland.acute <- ggplot(data=subset(transmissions)) +
  #geom_line(aes(x=year, y=ProportionAcute*100, group=sim.id, color="grey")) +
  geom_smooth(aes(x=year, y=ProportionAcute*100),method="loess", span=0.05, se = T, size=1, linetype=1) +
  xlab("Year")+
  ylab("Percent of transmissions from acute stage")+
  theme_bw(base_size=16) +
  guides(fill = guide_legend(keywidth = 2, keyheight = 1)) +
  scale_y_continuous(limits=c(0,80), expand = c(0,0)) +
  scale_x_continuous(breaks = seq(1985,year,5),limits=c(1985,year),expand = c(0,0)) +
  theme(legend.position="bottom") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(strip.background = element_rect(colour="black", fill="white")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"))
swaziland.acute

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\BMGF Lancet Paper\\Final_Plots")
ggsave("swaziland_pctAcuteTRansmissions_scenarios2.jpg", height=8, width=8)

##########################################################################
#Age distribution of cases / incidence
##########################################################################
#baseline
setwd(paste0(wd1))
files <- list.files(full.names = F)
head(files)
f <- paste0(wd1, files)
length(f)

for (i in seq(1,100,1)){
  reporthivbyageandgender <- read.csv(f[i])
  reporthivbyageandgender$scenario <- "baseline"
  reporthivbyageandgender$sim.id <- paste0(files[i])
  reporthivbyageandgender <- reporthivbyageandgender[,c(1,3,4,7,8,9,10,23,24)]
  reporthivbyageandgender <- subset(reporthivbyageandgender, Age > 10 & Age < 50)
  if (i==1){
    reporthivbyageandgender.master <- reporthivbyageandgender
  } else {
    reporthivbyageandgender.master <- rbind(reporthivbyageandgender.master, reporthivbyageandgender)
    print(paste0("Why hello there, I'm working on baseline folder",i))
  }
}

reporthivbyageandgender.master.baseline <- reporthivbyageandgender.master
head(reporthivbyageandgender.master.baseline)

#noART (only VMMC)
setwd(paste0(wd3))
files <- list.files(full.names = F)
head(files)
f <- paste0(wd3, files)
length(f)

for (i in seq(1,100,1)){
  reporthivbyageandgender <- read.csv(f[i])
  reporthivbyageandgender$scenario <- "noART"
  reporthivbyageandgender$sim.id <- paste0(files[i])
  reporthivbyageandgender <- reporthivbyageandgender[,c(1,3,4,7,8,9,10,23,24)]
  reporthivbyageandgender <- subset(reporthivbyageandgender, Age > 10 & Age < 50)
  if (i==1){
    reporthivbyageandgender.master <- reporthivbyageandgender
  } else {
    reporthivbyageandgender.master <- rbind(reporthivbyageandgender.master, reporthivbyageandgender)
    print(paste0("Why hello there, I'm working on noART folder",i))
  }
}

reporthivbyageandgender.master.noART <- reporthivbyageandgender.master
head(reporthivbyageandgender.master.noART)

#NoARTnoVMMC
setwd(paste0(wd2))
files <- list.files(full.names = F)
head(files)
f <- paste0(wd2, files)
length(f)

for (i in seq(1,100,1)){
  reporthivbyageandgender <- read.csv(f[i])
  reporthivbyageandgender$scenario <- "noARTnoVMMC"
  reporthivbyageandgender$sim.id <- paste0(files[i])
  reporthivbyageandgender <- reporthivbyageandgender[,c(1,3,4,7,8,9,10,23,24)]
  reporthivbyageandgender <- subset(reporthivbyageandgender, Age > 10 & Age < 50)
  if (i==1){
    reporthivbyageandgender.master <- reporthivbyageandgender
  } else {
    reporthivbyageandgender.master <- rbind(reporthivbyageandgender.master, reporthivbyageandgender)
    print(paste0("Why hello there, I'm working on noARTnoVMMC folder",i))
  }
}

reporthivbyageandgender.master.noartnovmmc <- reporthivbyageandgender.master
head(reporthivbyageandgender.master.noartnovmmc)

#ART100pct
setwd(paste0(wd6))
files <- list.files(full.names = F)
head(files)
f <- paste0(wd6, files)
length(f)

for (i in seq(1,100,1)){
  reporthivbyageandgender <- read.csv(f[i])
  reporthivbyageandgender$scenario <- "ART100pct"
  reporthivbyageandgender$sim.id <- paste0(files[i])
  reporthivbyageandgender <- reporthivbyageandgender[,c(1,3,4,7,8,9,10,23,24)]
  reporthivbyageandgender <- subset(reporthivbyageandgender, Age > 10 & Age < 50)
  if (i==1){
    reporthivbyageandgender.master <- reporthivbyageandgender
  } else {
    reporthivbyageandgender.master <- rbind(reporthivbyageandgender.master, reporthivbyageandgender)
    print(paste0("Why hello there, I'm working on ART100pct folder",i))
  }
}
reporthivbyageandgender.master.art100 <- reporthivbyageandgender.master
head(reporthivbyageandgender.master.art100)

reporthivbyageandgender.master <- rbind(reporthivbyageandgender.master.art100,
                                        reporthivbyageandgender.master.noartnovmmc,
                                        reporthivbyageandgender.master.noART,
                                        reporthivbyageandgender.master.baseline)
table(reporthivbyageandgender.master$scenario)

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\HIV\\EMOD\\Swaziland_05June2018\\Scenarios_calibration_run3_v05Sep2018_newcampaign_FINAL\\ScenariosOct_FINAL")
write.csv(reporthivbyageandgender.master, "reporthivbyageandgender.master.age.shifts.csv")

#Age distribution of incidence over time
reporthivbyageandgender.master$Year2 <- floor((reporthivbyageandgender.master$Year-0.5))
reporthivbyageandgender.master$Uninfected.Population <- reporthivbyageandgender.master$Population - reporthivbyageandgender.master$Infected
trajectories_IR.1a <- aggregate(Newly.Infected ~ Year2+Gender+Age+sim.id+scenario, subset(reporthivbyageandgender.master), FUN=sum) #sums number of new infections in each year
trajectories_IR.2 <- aggregate(Uninfected.Population ~ Year+Gender+Age+sim.id+scenario, subset(reporthivbyageandgender.master), FUN=sum)
trajectories_IR.2$Year2 <- floor(trajectories_IR.2$Year-0.5)
trajectories_IR.2 <- trajectories_IR.2[!duplicated(trajectories_IR.2[c("Year2","Gender","Age", "sim.id","scenario")]),] #remove second instance of duplicate rows
trajectories_IR.2 <- trajectories_IR.2[-match("Year",names(trajectories_IR.2))]
trajectories_IR <- merge(trajectories_IR.1a, trajectories_IR.2, by=c("Year2","Gender","Age","sim.id","scenario"))
trajectories_IR$incidence <- trajectories_IR$Newly.Infected / trajectories_IR$Uninfected.Population

#Proportion of cases in each age-group
head(trajectories_IR)
newlyinfected.tot <- aggregate(Newly.Infected ~ Year2+Gender+scenario+sim.id, subset(trajectories_IR), sum)
newlyinfected.merge <- merge(trajectories_IR, newlyinfected.tot, by=c("Year2","Gender","sim.id","scenario"))
newlyinfected.merge$proportion = newlyinfected.merge$Newly.Infected.x / newlyinfected.merge$Newly.Infected.y
newlyinfected.merge <- aggregate(proportion ~ Year2+Age+Gender+scenario, subset(newlyinfected.merge), median)
head(newlyinfected.merge,40)

#plot proportion of new infections by age
labs <- c("1" = "Women", "0" = "Men")
table(newlyinfected.merge$scenario)
labs.scenario <- c("baseline"="Baseline","noART"="VMMC only", "ART100pct"="100% ART","noARTnoVMMC"="No interventions" )

#baseline 1985 - 2005 (pre-intervention)
year1=1985
year2=2005
newinf.age.year <-  ggplot() +
  geom_line(data=subset(newlyinfected.merge, Year2>=year1 & Year2 <=year2 & scenario=="baseline"), aes(y=proportion*100, x=Age, group=Year2, color=Year2), size=1.2)+ 
  #scale_colour_gradientn(colours=rainbow(2), breaks = seq(year1,year2,round(((year2-year1)/2),0))) +
  scale_colour_gradient(low="pink", high="blue", breaks=c(year1,round((year2+year1)/2,0),year2)) +
  facet_grid(scenario~Gender,labeller=labeller(Gender = labs, scenario=labs.scenario)) +
  theme(legend.position="bottom") +
  xlab("Age")+
  ylab("Percent of new infections")+
  scale_x_continuous(expand=c(0,0)) +
  #scale_y_continuous(expand=c(0,0), limits=c(0,7),breaks=seq(0,7,1)) +
  theme_bw(base_size=16) +
  guides(fill = guide_legend(keywidth = 2, keyheight = 1, nrow=1)) +
  #theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(legend.position="bottom") +
  theme(strip.background = element_rect(colour="black", fill="white"))+
  theme(legend.title=element_blank()) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
newinf.age.year

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\Adam Incidence Trends\\Final_Plots")
ggsave("swaziland.modelled.pctnewinfections.byage.jpg", height=6, width=10)

year1=1985
year2=2005

newinf.age.year <-  ggplot() +
  geom_line(data=subset(newlyinfected.merge, Year2>=year1 & Year2 <=year2), aes(y=proportion*100, x=Age, group=Year2, color=Year2), size=1.2)+ 
  #scale_colour_gradientn(colours=rainbow(2), breaks = seq(year1,year2,round(((year2-year1)/2),0))) +
  scale_colour_gradient(low="pink", high="blue", breaks=c(year1,round((year2+year1)/2,0),year2)) +
  facet_grid(Gender~scenario,labeller=labeller(Gender = labs, scenario=labs.scenario)) +
  theme(legend.position="bottom") +
  xlab("Age")+
  ylab("Percent of new infections")+
  scale_x_continuous(expand=c(0,0),limits=c(15,49),breaks=seq(15,49,5)) +
  #scale_y_continuous(expand=c(0,0), limits=c(0,7),breaks=seq(0,7,1)) +
  theme_bw(base_size=16) +
  guides(fill = guide_legend(keywidth = 2, keyheight = 1, nrow=1)) +
  #theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(legend.position="bottom") +
  theme(strip.background = element_rect(colour="black", fill="white"))+
  theme(legend.title=element_blank()) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
newinf.age.year

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\Adam Incidence Trends\\Final_Plots")
ggsave("swaziland.modelled.pctnewinfections.byage.scenario.19902005.jpg", height=6, width=10)

year1=1985
year2=2030

newinf.age.year <-  ggplot() +
  geom_line(data=subset(newlyinfected.merge, Year2>=year1 & Year2 <=year2), aes(y=proportion*100, x=Age, group=Year2, color=Year2), size=1.2)+ 
  #scale_colour_gradientn(colours=rainbow(2), breaks = seq(year1,year2,round(((year2-year1)/2),0))) +
  scale_colour_gradient(low="pink", high="blue", breaks=c(year1,round((year2+year1)/2,0),year2)) +
  facet_grid(Gender~scenario,labeller=labeller(Gender = labs, scenario=labs.scenario)) +
  theme(legend.position="bottom") +
  xlab("Age")+
  ylab("Percent of new infections")+
  scale_x_continuous(expand=c(0,0),limits=c(15,49),breaks=seq(15,49,5)) +
  #scale_y_continuous(expand=c(0,0), limits=c(0,7),breaks=seq(0,7,1)) +
  theme_bw(base_size=16) +
  guides(fill = guide_legend(keywidth = 2, keyheight = 1, nrow=1)) +
  #theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(legend.position="bottom") +
  theme(strip.background = element_rect(colour="black", fill="white"))+
  theme(legend.title=element_blank()) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
newinf.age.year

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\Adam Incidence Trends\\Final_Plots")
ggsave("swaziland.modelled.pctnewinfections.byage.scenario.20052030.jpg", height=6, width=10)

#age-specific incidence over time
incidence.by.age <- aggregate(incidence ~ Year2+Gender+Age+scenario, subset(trajectories_IR), median)

year1=1985
year2=2030

inc.age.year <-  ggplot() +
  geom_line(data=subset(incidence.by.age, Year2>=year1 & Year2 <=year2), aes(y=incidence*100, x=Age, group=Year2, color=Year2), size=1.2)+ 
  #scale_colour_gradientn(colours=rainbow(2), breaks = seq(year1,year2,round(((year2-year1)/2),0))) +
  scale_colour_gradient(low="pink", high="blue", breaks=c(year1,round((year2+year1)/2,0),year2)) +
  facet_grid(Gender~scenario,labeller=labeller(Gender = labs, scenario=labs.scenario)) +
  theme(legend.position="bottom") +
  xlab("Age")+
  ylab("Incidence (per 100 py)")+
  scale_x_continuous(expand=c(0,0)) +
  #scale_y_continuous(expand=c(0,0), limits=c(0,7),breaks=seq(0,7,1)) +
  theme_bw(base_size=16) +
  guides(fill = guide_legend(keywidth = 2, keyheight = 1, nrow=1)) +
  #theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(legend.position="bottom") +
  theme(strip.background = element_rect(colour="black", fill="white"))+
  theme(legend.title=element_blank()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
#theme(legend.position="none")
inc.age.year

#################################################################################################
#Bring in Transmission data
#################################################################################################
#baseline transmission report
wd13 <- "C:\\Users\\aakullian\\Dropbox (IDM)\\HIV\\EMOD\\Swaziland_05June2018\\Scenarios_calibration_run3_v05Sep2018_newcampaign_FINAL\\ScenariosOct_FINAL\\Baseline_TransmissionReport\\TransmissionReport\\"
setwd(paste0(wd13))
files <- list.files(full.names = F)
head(files)
f <- paste0(wd13, files)
length(f)
i = 2

transmissionreport.master <- read.csv(f[i])
transmissionreport.master$scenario <- "baseline"
transmissionreport.master$sim.id <- paste0(files[i])
names(transmissionreport.master)
table(transmissionreport.master$SRC_STAGE)
table(transmissionreport.master$SRC_VIRAL_LOAD)

for (i in seq(1,250,1)){
  transmissionreport <- read.csv(f[i])
  transmissionreport <- transmissionreport[,c(1,2,4,6,7,16,18,19,29)]
  transmissionreport$scenario <- "baseline"
  transmissionreport$sim.id <- paste0(files[i])
  if (i==1){
    transmissionreport.master <- transmissionreport
  } else {
    transmissionreport.master <- rbind(transmissionreport.master, transmissionreport)
    print(paste0("Why hello there, I'm working on baseline transmission report folder",i))
  }
}

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\HIV\\EMOD\\Swaziland_05June2018\\Scenarios_calibration_run3_v05Sep2018_newcampaign_FINAL\\ScenariosOct_FINAL\\Baseline_TransmissionReport\\TransmissionReport\\")
saveRDS(transmissionreport.master, "transmissionreport.master.rds") 

#bring in saved dataset
setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\HIV\\EMOD\\Swaziland_05June2018\\Scenarios_calibration_run3_v05Sep2018_newcampaign_FINAL\\ScenariosOct_FINAL\\Baseline_TransmissionReport\\TransmissionReport\\")
transmissionreport.master <- readRDS("transmissionreport.master.rds")

#create variables
head(transmissionreport.master)
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
#transmissionreport <- transmissionreport[,c(3,4,6,7,9,11,12,13,14,15)]

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

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\BMGF Lancet Paper\\Final_Plots")
ggsave("proportionoftransmissionsacuteandearly.jpg", height=8, width=10)

###########################################################################################################
#Transmission pair formation
###########################################################################################################
head(transmissionreport,5)

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

year1=2005
year2=2050
labs.gender <- c("0" = "Male","1" = "Female")
p.age.dist.transmissions <- ggplot(subset(transmissionreport.n.median, year.fl >= year1 & year.fl <= year2)) +
  geom_smooth(aes(x=src_age_int, y=prop*100, color=year.fl, group=year.fl),method="loess", span=0.8, se = F, size=1, linetype=1) +
  facet_grid(~SRC_GENDER, labeller=labeller(SRC_GENDER = labs.gender)) +
  scale_colour_gradientn(colours=rainbow(2),guide = guide_colourbar(barwidth=30, barheight=1), breaks=seq(2010, year2, 10)) +
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

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\BMGF Lancet Paper\\Final_Plots")
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

##########################################
#Changes in HIV incidence over time by age and gender
##########################################
# wd14 <- "C:\\Users\\aakullian\\Dropbox (IDM)\\HIV\\EMOD\\Swaziland_05June2018\\Scenarios_calibration_run3_v05Sep2018_newcampaign_FINAL\\ScenariosOct_FINAL\\baseline\\ReportHIVByAgeAndGender\\"
# setwd(paste0(wd14))
# files <- list.files(full.names = F)
# f <- paste0(wd14, files)
# reporthivbyageandgender <- read.csv(f[1])
# names(reporthivbyageandgender)
# reporthivbyageandgender <- reporthivbyageandgender[,c(1,3,7,8,9,10,11)]

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\HIV\\EMOD\\Swaziland_05June2018\\Scenarios_calibration_run3_v05Sep2018_newcampaign_FINAL\\ScenariosOct_FINAL")
reporthivbyageandgender.master.final <- read.csv("reporthivbyageandgender.master.final.csv")

reporthivbyageandgender.master2 <- subset(reporthivbyageandgender.master.final, scenario=="baseline")
table(reporthivbyageandgender.master2$scenario)

reporthivbyageandgender.master2$Year.cat <- cut(reporthivbyageandgender.master2$Year, seq(1980, 2050, 5), right=F)
reporthivbyageandgender.master2$age.cat <- cut(reporthivbyageandgender.master2$Age, c(15,25,35,100), right=F)
table(reporthivbyageandgender.master2$Year.cat)
table(reporthivbyageandgender.master2$age.cat)
head(reporthivbyageandgender.master2)

reporthivbyageandgender.master2$Year2 <- floor((reporthivbyageandgender.master2$Year-0.5))
reporthivbyageandgender.master2$Uninfected.Population <- reporthivbyageandgender.master2$Population - reporthivbyageandgender.master2$Infected
trajectories_IR.1a <- aggregate(Newly.Infected ~ Year2+Gender+age.cat+sim.id, subset(reporthivbyageandgender.master2, Age>10 & Age < 60), FUN=sum) #sums number of new infections in each year
trajectories_IR.2 <- aggregate(Uninfected.Population ~ Year+Gender+age.cat+sim.id, subset(reporthivbyageandgender.master2, Age>10 & Age < 60), FUN=sum)
trajectories_IR.2$Year2 <- floor(trajectories_IR.2$Year-0.5)
trajectories_IR.2 <- trajectories_IR.2[!duplicated(trajectories_IR.2[c("Year2","Gender","age.cat","sim.id")]),] #remove second instance of duplicate rows
trajectories_IR.2 <- trajectories_IR.2[-match("Year",names(trajectories_IR.2))]
trajectories_IR <- merge(trajectories_IR.1a, trajectories_IR.2, by=c("Year2","Gender","age.cat","sim.id"))
trajectories_IR$incidence <- trajectories_IR$Newly.Infected / trajectories_IR$Uninfected.Population
head(trajectories_IR)

year1=2000
year2=2016
#labs <- c("1" = "Women", "0" = "Men")
labs <- c("[15,25)" = "15-24","[25,35)" = "25-34", "[35,100)" = "35+")
head(trajectories_IR)
trajectories_IR_plot <- subset(trajectories_IR, Year2 >= year1 & Year2 <= year2)
p.age.dist.acquisitions <- ggplot() +
  geom_smooth(data=trajectories_IR_plot,aes(x=Year2, y=incidence*100, color=factor(Gender), group=Gender),method="loess", span=0.4, se = F, size=2, linetype=1) +
  geom_line(data=trajectories_IR_plot,aes(x=Year2, y=incidence*100, group=interaction(Gender, sim.id), color=factor(Gender)), size = 1, alpha=0.02) +
  facet_grid(~age.cat, labeller=labeller(age.cat = labs)) +
  geom_vline(xintercept=2004, linetype=2) +
  #geom_vline(xintercept=2011, linetype=1) +
  scale_color_manual(values = c("red","blue"), labels=c("Men","Women")) +
  xlab("Year") +
  ylab("Incidence per 100 py") +
  scale_x_continuous(limits = c(year1, year2+0.2), breaks=seq(year1,year2,5),expand=c(0,0)) +
  scale_y_continuous(limits = c(0, 7), breaks=seq(0,7,1),expand=c(0,0)) +
  theme_bw(base_size=20) +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1)) +
  guides(fill = guide_legend(keywidth = 2, keyheight = 1)) +
  theme(legend.title=element_blank())+
  theme(legend.position="bottom") +
  theme(strip.background = element_rect(colour="black", fill="white")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p.age.dist.acquisitions

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\BMGF Lancet Paper\\Final_Plots")
ggsave("HIVincidenceovertimebyageandgender.jpg", height=7, width=10)

#Rearrange age and gender so panels are by gender
year1=2000
year2=2016
labs <- c("1" = "Women", "0" = "Men")
#labs <- c("[15,25)" = "15-24","[25,35)" = "25-34", "[35,100)" = "35+")
head(trajectories_IR)
trajectories_IR_plot <- subset(trajectories_IR, Year2 >= year1 & Year2 <= year2)
p.age.dist.acquisitions <- ggplot() +
  geom_smooth(data=trajectories_IR_plot,aes(x=Year2, y=incidence*100, color=factor(age.cat), group=age.cat),method="loess", span=0.4, se = F, size=2, linetype=1) +
  geom_line(data=trajectories_IR_plot,aes(x=Year2, y=incidence*100, group=interaction(age.cat, sim.id), color=factor(age.cat)), size = 1, alpha=0.02) +
  facet_grid(~Gender, labeller=labeller(Gender = labs)) +
  geom_vline(xintercept=2004, linetype=2) +
  #geom_vline(xintercept=2011, linetype=1) +
  scale_color_manual(values = c("red","purple","blue"), labels=c("15-24","25-34","35+")) +
  xlab("Year") +
  ylab("Incidence per 100 py") +
  scale_x_continuous(limits = c(year1, year2+0.2), breaks=seq(year1,year2,1),expand=c(0,0)) +
  scale_y_continuous(limits = c(0, 7), breaks=seq(0,7,1),expand=c(0,0)) +
  theme_bw(base_size=20) +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1)) +
  guides(fill = guide_legend(keywidth = 2, keyheight = 1)) +
  theme(legend.title=element_blank())+
  theme(legend.position="bottom") +
  theme(strip.background = element_rect(colour="black", fill="white")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p.age.dist.acquisitions

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\BMGF Lancet Paper\\Final_Plots")
ggsave("HIVincidenceovertimebyageandgender_bygender.jpg", height=7, width=10)

############################################
#Age distribution of incidence over time
############################################
reporthivbyageandgender.master2 <- subset(reporthivbyageandgender.master.final, scenario=="baseline")
table(reporthivbyageandgender.master2$scenario)

reporthivbyageandgender.master2$Year2 <- floor((reporthivbyageandgender.master2$Year-0.5))
reporthivbyageandgender.master2$Uninfected.Population <- reporthivbyageandgender.master2$Population - reporthivbyageandgender.master2$Infected
trajectories_IR.1a <- aggregate(Newly.Infected ~ Year2+Gender+Age+sim.id, subset(reporthivbyageandgender.master2, Age>10 & Age < 60), FUN=sum) #sums number of new infections in each year
trajectories_IR.2 <- aggregate(Uninfected.Population ~ Year+Gender+Age+sim.id, subset(reporthivbyageandgender.master2, Age>10 & Age < 60), FUN=sum)
trajectories_IR.2$Year2 <- floor(trajectories_IR.2$Year-0.5)
trajectories_IR.2 <- trajectories_IR.2[!duplicated(trajectories_IR.2[c("Year2","Gender","Age","sim.id")]),] #remove second instance of duplicate rows
trajectories_IR.2 <- trajectories_IR.2[-match("Year",names(trajectories_IR.2))]
trajectories_IR <- merge(trajectories_IR.1a, trajectories_IR.2, by=c("Year2","Gender","Age","sim.id"))
trajectories_IR$incidence <- trajectories_IR$Newly.Infected / trajectories_IR$Uninfected.Population
# trajectories_IR$year.cat <- cut(trajectories_IR$Year2, seq(1990,2050,1), right=F, dig.lab=4)
# trajectories_IR.3 <- aggregate(cbind(Newly.Infected,Uninfected.Population) ~ year.cat+Gender+Age+sim.id, subset(trajectories_IR, Age>10 & Age < 60), FUN=sum)
# trajectories_IR.3$incidence <- trajectories_IR.3$Newly.Infected / trajectories_IR.3$Uninfected.Population
# trajectories_IR.3$year <- as.numeric(substr(trajectories_IR.3$year.cat, 2,5))
# table(trajectories_IR.3$year)
# head(trajectories_IR.3)

year1=2005
year2=2050
lateyear=2050
labs <- c("1" = "Women", "0" = "Men")
p.age.dist.acquisitions <- ggplot(subset(trajectories_IR, Year2>=year1 & Year2<=year2)) +
  #geom_line(aes(x=Age, y=incidence*100, color=year.cat, group=interaction(year.cat,sim.id), alpha=0.001)) +
  geom_smooth(aes(x=Age, y=incidence*100, color=Year2, group=Year2),method="loess", span=0.7, se = F, size=1, linetype=1) +
  facet_grid(~Gender, labeller=labeller(Gender = labs)) +
  scale_colour_gradientn(colours=rainbow(2),guide = guide_colourbar(barwidth=30, barheight=1), breaks=seq(year1, year2, 10)) +
  xlab("Age") +
  ylab("Incidence per 100 py") +
  scale_x_continuous(limits = c(15, 55), breaks=seq(15,50,5),expand=c(0,0)) +
  scale_y_continuous(limits = c(0, 6), breaks=seq(0,6,1),expand=c(0,0)) +
  theme_bw(base_size=24) +
  guides(fill = guide_legend(keywidth = 2, keyheight = 1)) +
  theme(legend.title=element_blank())+
  theme(legend.position="bottom") +
  theme(strip.background = element_rect(colour="black", fill="white")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p.age.dist.acquisitions

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\BMGF Lancet Paper\\Final_Plots")
ggsave("HIVincidenceovertimebyageandgender_agedist20052050.jpg", height=7, width=10)

reporthivbyageandgender.master2 <- subset(reporthivbyageandgender.master.final, scenario=="baseline")
table(reporthivbyageandgender.master2$scenario)

reporthivbyageandgender.master2$Year2 <- floor((reporthivbyageandgender.master2$Year-0.5))
reporthivbyageandgender.master2$Uninfected.Population <- reporthivbyageandgender.master2$Population - reporthivbyageandgender.master2$Infected
trajectories_IR.1a <- aggregate(Newly.Infected ~ Year2+Gender+Age+sim.id, subset(reporthivbyageandgender.master2, Age>10 & Age < 60), FUN=sum) #sums number of new infections in each year
trajectories_IR.2 <- aggregate(Uninfected.Population ~ Year+Gender+Age+sim.id, subset(reporthivbyageandgender.master2, Age>10 & Age < 60), FUN=sum)
trajectories_IR.2$Year2 <- floor(trajectories_IR.2$Year-0.5)
trajectories_IR.2 <- trajectories_IR.2[!duplicated(trajectories_IR.2[c("Year2","Gender","Age","sim.id")]),] #remove second instance of duplicate rows
trajectories_IR.2 <- trajectories_IR.2[-match("Year",names(trajectories_IR.2))]
trajectories_IR <- merge(trajectories_IR.1a, trajectories_IR.2, by=c("Year2","Gender","Age","sim.id"))
trajectories_IR$incidence <- trajectories_IR$Newly.Infected / trajectories_IR$Uninfected.Population

year1=1980
year2=2004
head(trajectories_IR)
labs <- c("1" = "Women", "0" = "Men")
p.age.dist.acquisitions <- ggplot(subset(trajectories_IR, Year2 >= year1 & Year2 <=year2 & Age < 25 | (Age==25 & Gender==0 & Year2 >= year1 & Year2 <=year2))) +
  #geom_line(aes(x=Age, y=incidence*100, color=year.cat, group=interaction(year.cat,sim.id), alpha=0.001)) +
  geom_smooth(aes(x=Year2, y=incidence*100, color=Age, group=Age),method="loess", span=0.5, se = F, size=2, linetype=1) +
  facet_grid(~Gender, labeller=labeller(Gender = labs)) +
  scale_colour_gradientn(colours=rainbow(2), limits=c(15,55)) +
  xlab("Year") +
  ylab("Incidence per 100 py") +
  scale_x_continuous(limits = c(1980, 2016), breaks=seq(1980,2016,5),expand=c(0,0)) +
  scale_y_continuous(limits = c(0, 6), breaks=seq(0,6,1),expand=c(0,0)) +
  theme_bw(base_size=20) +
  guides(fill = guide_legend(keywidth = 2, keyheight = 1)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(legend.title=element_blank())+
  theme(legend.position="bottom") +
  theme(strip.background = element_rect(colour="black", fill="white")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p.age.dist.acquisitions

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\BMGF Lancet Paper\\Final_Plots")
ggsave("HIVincidenceovertimebyageandgender_overtime_19802004_young.jpg", height=8, width=12)

year1=1980
year2=2030
head(trajectories_IR)
labs <- c("1" = "Women", "0" = "Men")
p.age.dist.acquisitions <- ggplot(subset(trajectories_IR, Year2 >= year1 & Year2 <=year2 & Age < 25 | (Age==25 & Gender==0 & Year2 >= year1 & Year2 <=year2))) +
  #geom_line(aes(x=Age, y=incidence*100, color=year.cat, group=interaction(year.cat,sim.id), alpha=0.001)) +
  geom_smooth(aes(x=Year2, y=incidence*100, color=Age, group=Age),method="loess", span=0.5, se = F, size=2, linetype=1) +
  facet_grid(~Gender, labeller=labeller(Gender = labs)) +
  scale_colour_gradientn(colours=rainbow(2), limits=c(15,55)) +
  xlab("Year") +
  ylab("Incidence per 100 py") +
  scale_x_continuous(limits = c(1980, 2016), breaks=seq(1980,2016,5),expand=c(0,0)) +
  scale_y_continuous(limits = c(0, 6), breaks=seq(0,6,1),expand=c(0,0)) +
  theme_bw(base_size=20) +
  guides(fill = guide_legend(keywidth = 2, keyheight = 1)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(legend.title=element_blank())+
  theme(legend.position="bottom") +
  theme(strip.background = element_rect(colour="black", fill="white")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p.age.dist.acquisitions

year1=1980
year2=2030
head(trajectories_IR)
labs <- c("1" = "Women", "0" = "Men")
p.age.dist.acquisitions <- ggplot(subset(trajectories_IR, Year2 >= year1 & Year2 <=year2)) +
  #geom_line(aes(x=Age, y=incidence*100, color=year.cat, group=interaction(year.cat,sim.id), alpha=0.001)) +
  geom_smooth(aes(x=Year2, y=incidence*100, color=Age, group=Age),method="loess", span=0.5, se = F, size=2, linetype=1) +
  facet_grid(~Gender, labeller=labeller(Gender = labs)) +
  scale_colour_gradientn(colours=rainbow(2), limits=c(15,55)) +
  xlab("Year") +
  ylab("Incidence per 100 py") +
  scale_x_continuous(limits = c(1980, 2016), breaks=seq(1980,2016,5),expand=c(0,0)) +
  scale_y_continuous(limits = c(0, 6), breaks=seq(0,6,1),expand=c(0,0)) +
  theme_bw(base_size=20) +
  guides(fill = guide_legend(keywidth = 2, keyheight = 1)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(legend.title=element_blank())+
  theme(legend.position="bottom") +
  theme(strip.background = element_rect(colour="black", fill="white")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p.age.dist.acquisitions

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\BMGF Lancet Paper\\Final_Plots")
ggsave("HIVincidenceovertimebyageandgender_overtime_19802016_young.jpg", height=8, width=12)

#Change in incidence by broad age group (over and under 25)
reporthivbyageandgender.master2 <- subset(reporthivbyageandgender.master.final, scenario=="baseline")
table(reporthivbyageandgender.master2$scenario)

reporthivbyageandgender.master2$under25 <- NA
reporthivbyageandgender.master2$under25[reporthivbyageandgender.master2$Age>10 & reporthivbyageandgender.master2$Age<25] <- 1
reporthivbyageandgender.master2$under25[reporthivbyageandgender.master2$Age>20] <- 0
head(reporthivbyageandgender.master2)

reporthivbyageandgender.master2$Year2 <- floor((reporthivbyageandgender.master2$Year-0.5))
reporthivbyageandgender.master2$Uninfected.Population <- reporthivbyageandgender.master2$Population - reporthivbyageandgender.master2$Infected
trajectories_IR.1a <- aggregate(Newly.Infected ~ Year2+under25+Gender+sim.id, subset(reporthivbyageandgender.master2, Age>10 & Age < 100), FUN=sum) #sums number of new infections in each year
trajectories_IR.2 <- aggregate(Uninfected.Population ~ Year+under25+Gender+sim.id, subset(reporthivbyageandgender.master2, Age>10 & Age < 100), FUN=sum)
trajectories_IR.2$Year2 <- floor(trajectories_IR.2$Year-0.5)
trajectories_IR.2 <- trajectories_IR.2[!duplicated(trajectories_IR.2[c("Year2","under25","Gender", "sim.id")]),] #remove second instance of duplicate rows
trajectories_IR.2 <- trajectories_IR.2[-match("Year",names(trajectories_IR.2))]
trajectories_IR <- merge(trajectories_IR.1a, trajectories_IR.2, by=c("Year2","under25","Gender","sim.id"))
trajectories_IR$incidence <- trajectories_IR$Newly.Infected / trajectories_IR$Uninfected.Population
summary(trajectories_IR$incidence)
head(trajectories_IR)
trajectories_IR.w <- reshape(trajectories_IR, idvar = c("Year2","Gender","sim.id"), timevar = "under25", direction = "wide")
trajectories_IR.w$IRR.0.1 <- trajectories_IR.w$incidence.0 / trajectories_IR.w$incidence.1
head(trajectories_IR.w)

#take median of under / over 25 prevalence values of 250 sims
p05 <- function(x) quantile(x, prob=.05, na.rm=TRUE) ## Newly defined functions can be used
p95 <- function(x) quantile(x, prob=.95, na.rm=TRUE)
smooth_vals <- aggregate(IRR.0.1 ~ Year2+Gender, subset(trajectories_IR.w), FUN=function(x) c(median=median(x), p05 = p05(x), p95= p95(x)))
head(smooth_vals)

#Age of peak incidence
reporthivbyageandgender.master2$Year2 <- floor((reporthivbyageandgender.master2$Year-0.5))
reporthivbyageandgender.master2$Uninfected.Population <- reporthivbyageandgender.master2$Population - reporthivbyageandgender.master2$Infected
trajectories_IR.1a <- aggregate(Newly.Infected ~ Year2+Age+Gender+sim.id, subset(reporthivbyageandgender.master2, Age>10 & Age < 100), FUN=sum) #sums number of new infections in each year
trajectories_IR.2 <- aggregate(Uninfected.Population ~ Year+Age+Gender+sim.id, subset(reporthivbyageandgender.master2, Age>10 & Age < 100), FUN=sum)
trajectories_IR.2$Year2 <- floor(trajectories_IR.2$Year-0.5)
trajectories_IR.2 <- trajectories_IR.2[!duplicated(trajectories_IR.2[c("Year2","Age","Gender", "sim.id")]),] #remove second instance of duplicate rows
trajectories_IR.2 <- trajectories_IR.2[-match("Year",names(trajectories_IR.2))]
trajectories_IR <- merge(trajectories_IR.1a, trajectories_IR.2, by=c("Year2","Age","Gender","sim.id"))
trajectories_IR$incidence <- trajectories_IR$Newly.Infected / trajectories_IR$Uninfected.Population
summary(trajectories_IR$incidence)

p05 <- function(x) quantile(x, prob=.05, na.rm=TRUE) ## Newly defined functions can be used
p95 <- function(x) quantile(x, prob=.95, na.rm=TRUE)
smooth_vals <- aggregate(incidence ~ Gender+Year2+Age, subset(trajectories_IR), median)
head(smooth_vals)
smooth_vals.max <- aggregate(incidence ~ Gender+Year2, subset(smooth_vals), max)
smooth_vals.max.merge <- merge(smooth_vals.max, smooth_vals)
smooth_vals.max.merge[smooth_vals.max.merge$Gender==1,]




##########################################
#Changes in prevalence and  viral load suppression over time by age and gender
##########################################


#Swaziland prevalence over time by age and gender
labs <- c("1" = "Women", "0" = "Men")
head(prevalence)
swaziland.prev <- ggplot(subset(prevalence, Year < 2011)) +
  geom_point(size=2, aes(x=Age, y=NationalPrevalence*100, color=factor(Year))) + 
  geom_line(size=2, aes(x=Age, y=NationalPrevalence*100, color=factor(Year))) + 
  geom_errorbar(aes(x=Age, ymin=lb*100, ymax=ub*100,  color=factor(Year)), width=2, size=1) +
  scale_color_manual(values = c("blue","purple", "red")) +
  facet_grid(~Gender, labeller=labeller(Gender = labs)) +
  xlab("Year")+
  ylab("Prevalence")+
  theme_bw(base_size=20) +
  guides(fill = guide_legend(keywidth = 2, keyheight = 1)) +
  theme(legend.title=element_blank())+
  scale_y_continuous(breaks = seq(0,60,10),limits=c(0,60),expand=c(0,0)) +
  scale_x_continuous(breaks = seq(15,60,5),expand=c(0,0)) +
  theme(legend.position="bottom") +
  theme(strip.background = element_rect(colour="black", fill="white")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
swaziland.prev

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\BMGF Lancet Paper\\Final_Plots")
ggsave("SwazilandHIVprevalencebyageandgender1.jpg", height=6, width=10)

#calculate HIV prevalence
reporthivbyageandgender.master2 <- subset(reporthivbyageandgender.master.final, scenario=="baseline")
table(reporthivbyageandgender.master2$scenario)
reporthivbyageandgender.master2$under25 <- NA
reporthivbyageandgender.master2$under25[reporthivbyageandgender.master2$Age>10 & reporthivbyageandgender.master2$Age<25] <- 1
reporthivbyageandgender.master2$under25[reporthivbyageandgender.master2$Age>20] <- 0
table(reporthivbyageandgender.master2$Age, reporthivbyageandgender.master2$under25)

Population <- aggregate(Population ~ Year+Gender+Age+sim.id, subset(reporthivbyageandgender.master2, Age>10 & Age < 100), FUN=sum)
Infected <- aggregate(Infected ~ Year+Gender+Age+sim.id, subset(reporthivbyageandgender.master2, Age>10 & Age < 100), FUN=sum)
prev <- merge(Infected, Population, by=c("Year","Gender","Age","sim.id"))
prev$prevalence <- prev$Infected / prev$Population
summary(prev$prevalence)
head(prev)

#Loop through prev by age over time (ppdv ~= prevalence before ART scale-up)
labs.age <- c("15-19","20-24","25-29","30-34","35-39","40-44","45-49","50-54","55-59","60-64")
year1 = 1990
year2 = 2004
prev.plot <- ggplot(subset(prev, Year == 1990)) +
  geom_line(size=1.2, linetype=1, aes(x = Age, y=prevalence*100, color=Year, group=interaction(Year,sim.id)), alpha=0.03)+ 
  geom_smooth(aes(x=Age, y=prevalence*100, color=Year, group=Year),method="loess", span=0.5, se = F, size=2, linetype=1) +
  scale_colour_gradientn(colours=rainbow(3),limits=c(year1, year2), breaks = c(year1,round((year1+year2)/2,0),year2)) +
  facet_grid(~Gender,labeller=labeller(Gender = labs)) +
  xlab("Age") +
  ylab("HIV prevalence") +
  scale_x_continuous(limits = c(15, 60), breaks=seq(15,60,5),expand=c(0,0), labels=labs.age) +
  scale_y_continuous(limits = c(0, 60), breaks=seq(0,60,5),expand=c(0,0)) +
  theme_bw(base_size=20) +
  theme(legend.title=element_blank())+
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1)) +
  theme(legend.position="bottom") +
  theme(strip.background = element_rect(colour="black", fill="white")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
prev.plot

head(prev)

#under versus over 25 
Population <- aggregate(Population ~ Year+Gender+under25+sim.id, subset(reporthivbyageandgender.master2, Age>10 & Age < 100), FUN=sum)
Infected <- aggregate(Infected ~ Year+Gender+under25+sim.id, subset(reporthivbyageandgender.master2, Age>10 & Age < 100), FUN=sum)
prev <- merge(Infected, Population, by=c("Year","under25","sim.id"))
prev$prevalence <- prev$Infected / prev$Population
summary(prev$prevalence)
head(prev)

#take median of under / over 25 prevalence values of 250 sims
p05 <- function(x) quantile(x, prob=.05, na.rm=TRUE) ## Newly defined functions can be used
p95 <- function(x) quantile(x, prob=.95, na.rm=TRUE)
smooth_vals <- aggregate(prevalence ~ Year+under25, subset(prev), FUN=function(x) c(median=median(x), p05 = p05(x), p95= p95(x)))
smooth_vals[smooth_vals$under25==1,]

#By 5 year age group
Population <- aggregate(Population ~ Year+Gender+Age+sim.id, subset(reporthivbyageandgender.master2, Age>10 & Age < 100), FUN=sum)
Infected <- aggregate(Infected ~ Year+Gender+Age+sim.id, subset(reporthivbyageandgender.master2, Age>10 & Age < 100), FUN=sum)
prev <- merge(Infected, Population, by=c("Year","Gender","Age","sim.id"))
prev$prevalence <- prev$Infected / prev$Population
summary(prev$prevalence)
head(prev)

#Age of peak prevalence
p05 <- function(x) quantile(x, prob=.05, na.rm=TRUE) ## Newly defined functions can be used
p95 <- function(x) quantile(x, prob=.95, na.rm=TRUE)
smooth_vals <- aggregate(prevalence ~ Gender+Year+Age, subset(prev), median)
head(smooth_vals)
smooth_vals.max <- aggregate(prevalence ~ Gender+Year, subset(smooth_vals), max)
smooth_vals.max.merge <- merge(smooth_vals.max, smooth_vals)
smooth_vals.max.merge[smooth_vals.max.merge$Gender==0,]



##########################################################
#calculate population prevalence of unsuppressed viral load
##########################################################
reporthivbyageandgender.master2 <- subset(reporthivbyageandgender.master.final, scenario=="baseline")
table(reporthivbyageandgender.master2$scenario)

OnARTnumbers <- aggregate(On_ART ~ Year+Gender+Age+sim.id, subset(reporthivbyageandgender.master2, Age>10 & Age < 100), FUN=sum) #sums prev infections in each year
Population <- aggregate(Population ~ Year+Gender+Age+sim.id, subset(reporthivbyageandgender.master2, Age>10 & Age < 100), FUN=sum)
Infected <- aggregate(Infected ~ Year+Gender+Age+sim.id, subset(reporthivbyageandgender.master2, Age>10 & Age < 100), FUN=sum)

ppdv <- merge(OnARTnumbers, Population, by=c("Year","Gender","Age","sim.id"))
ppdv <- merge(ppdv, Infected, by=c("Year","Gender","Age","sim.id"))
ppdv$suppressed <- ppdv$On_ART*0.92
ppdv$unsuppr <- ppdv$Infected - ppdv$suppressed
ppdv$ppdv <- ppdv$unsuppr / ppdv$Population
ppdv$prevalence <- ppdv$Infected / ppdv$Population
summary(ppdv$ppdv)
summary(ppdv$prevalence)
head(ppdv)


#Loop through ppdv by age over time (ppdv ~= prevalence before ART scale-up)
year1 = 1990
year2 = 2004
prev.plot <- ggplot(subset(ppdv, Year == 1990)) +
  geom_line(size=1.2, linetype=1, aes(x = Age, y=ppdv*100, color=Year, group=interaction(Year,sim.id)), alpha=0.03)+ 
  geom_smooth(aes(x=Age, y=ppdv*100, color=Year, group=Year),method="loess", span=0.5, se = F, size=2, linetype=1) +
  scale_colour_gradientn(colours=rainbow(3),limits=c(year1, year2), breaks = c(year1,round((year1+year2)/2,0),year2)) +
  facet_grid(~Gender,labeller=labeller(Gender = labs)) +
  xlab("Age") +
  ylab("% with detectable viral load") +
  scale_x_continuous(limits = c(15, 60), breaks=seq(15,59,5),expand=c(0,0)) +
  scale_y_continuous(limits = c(0, 60), breaks=seq(0,60,5),expand=c(0,0)) +
  theme_bw(base_size=20) +
  theme(legend.title=element_blank())+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(legend.position="bottom") +
  theme(strip.background = element_rect(colour="black", fill="white")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
prev.plot

predict(loess(incidence~Year2,baseline.smooth, span=0.2))

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\BMGF Lancet Paper\\Final_Plots")
ggsave("pppv.plot.1990.jpg", height=8, width=12)

prev.plot <- ggplot(subset(ppdv, Year==1990 | Year==1995)) +
  geom_line(size=1.2, linetype=1, aes(x = Age, y=ppdv*100, color=Year, group=interaction(Year,sim.id)), alpha=0.03)+ 
  geom_smooth(aes(x=Age, y=ppdv*100, color=Year, group=Year),method="loess", span=0.5, se = F, size=2, linetype=1) +
  scale_colour_gradientn(colours=rainbow(3),limits=c(year1, year2), breaks = c(year1,round((year1+year2)/2,0),year2)) +
  facet_grid(~Gender,labeller=labeller(Gender = labs)) +
  xlab("Age") +
  ylab("% with detectable viral load") +
  scale_x_continuous(limits = c(15, 60), breaks=seq(15,59,5),expand=c(0,0)) +
  scale_y_continuous(limits = c(0, 60), breaks=seq(0,60,5),expand=c(0,0)) +
  theme_bw(base_size=20) +
  theme(legend.title=element_blank())+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(legend.position="bottom") +
  theme(strip.background = element_rect(colour="black", fill="white")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
prev.plot

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\BMGF Lancet Paper\\Final_Plots")
ggsave("pppv.plot.1995.jpg", height=8, width=12)

prev.plot <- ggplot(subset(ppdv, Year==1990 | Year==1995 | Year==2000)) +
  geom_line(size=1.2, linetype=1, aes(x = Age, y=ppdv*100, color=Year, group=interaction(Year,sim.id)), alpha=0.03)+ 
  geom_smooth(aes(x=Age, y=ppdv*100, color=Year, group=Year),method="loess", span=0.5, se = F, size=2, linetype=1) +
  scale_colour_gradientn(colours=rainbow(3),limits=c(year1, year2), breaks = c(year1,round((year1+year2)/2,0),year2)) +
  facet_grid(~Gender,labeller=labeller(Gender = labs)) +
  xlab("Age") +
  ylab("% with detectable viral load") +
  scale_x_continuous(limits = c(15, 60), breaks=seq(15,59,5),expand=c(0,0)) +
  scale_y_continuous(limits = c(0, 60), breaks=seq(0,60,5),expand=c(0,0)) +
  theme_bw(base_size=20) +
  theme(legend.title=element_blank())+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(legend.position="bottom") +
  theme(strip.background = element_rect(colour="black", fill="white")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
prev.plot

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\BMGF Lancet Paper\\Final_Plots")
ggsave("pppv.plot.2000.jpg", height=8, width=12)

prev.plot <- ggplot(subset(ppdv, Year==1990 | Year==1995 | Year==2000 | Year == 2005)) +
  geom_line(size=1.2, linetype=1, aes(x = Age, y=ppdv*100, color=Year, group=interaction(Year,sim.id)), alpha=0.03)+ 
  geom_smooth(aes(x=Age, y=ppdv*100, color=Year, group=Year),method="loess", span=0.5, se = F, size=2, linetype=1) +
  scale_colour_gradientn(colours=rainbow(3),limits=c(year1, year2), breaks = c(year1,round((year1+year2)/2,0),year2)) +
  facet_grid(~Gender,labeller=labeller(Gender = labs)) +
  xlab("Age") +
  ylab("% with detectable viral load") +
  scale_x_continuous(limits = c(15, 60), breaks=seq(15,59,5),expand=c(0,0)) +
  scale_y_continuous(limits = c(0, 60), breaks=seq(0,60,5),expand=c(0,0)) +
  theme_bw(base_size=20) +
  theme(legend.title=element_blank())+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(legend.position="bottom") +
  theme(strip.background = element_rect(colour="black", fill="white")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
prev.plot

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\BMGF Lancet Paper\\Final_Plots")
ggsave("pppv.plot.2005.jpg", height=8, width=12)

prev.plot <- ggplot(subset(ppdv, Year==1990 | Year==1995 | Year==2000 | Year == 2005 | Year == 2010)) +
  geom_line(size=1.2, linetype=1, aes(x = Age, y=ppdv*100, color=Year, group=interaction(Year,sim.id)), alpha=0.03)+ 
  geom_smooth(aes(x=Age, y=ppdv*100, color=Year, group=Year),method="loess", span=0.5, se = F, size=2, linetype=1) +
  scale_colour_gradientn(colours=rainbow(3),limits=c(year1, year2), breaks = c(year1,round((year1+year2)/2,0),year2)) +
  facet_grid(~Gender,labeller=labeller(Gender = labs)) +
  xlab("Age") +
  ylab("% with detectable viral load") +
  scale_x_continuous(limits = c(15, 60), breaks=seq(15,59,5),expand=c(0,0)) +
  scale_y_continuous(limits = c(0, 60), breaks=seq(0,60,5),expand=c(0,0)) +
  theme_bw(base_size=20) +
  theme(legend.title=element_blank())+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(legend.position="bottom") +
  theme(strip.background = element_rect(colour="black", fill="white")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
prev.plot

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\BMGF Lancet Paper\\Final_Plots")
ggsave("pppv.plot.2010.jpg", height=8, width=12)

prev.plot <- ggplot(subset(ppdv, Year==1990 | Year==1995 | Year==2000 | Year == 2005 | Year == 2016)) +
  geom_line(size=1.2, linetype=1, aes(x = Age, y=ppdv*100, color=Year, group=interaction(Year,sim.id)), alpha=0.03)+ 
  geom_smooth(aes(x=Age, y=ppdv*100, color=Year, group=Year),method="loess", span=0.5, se = F, size=2, linetype=1) +
  scale_colour_gradientn(colours=rainbow(3),limits=c(year1, year2), breaks = c(year1,round((year1+year2)/2,0),year2)) +
  facet_grid(~Gender,labeller=labeller(Gender = labs)) +
  xlab("Age") +
  ylab("% with detectable viral load") +
  scale_x_continuous(limits = c(15, 60), breaks=seq(15,59,5),expand=c(0,0)) +
  scale_y_continuous(limits = c(0, 60), breaks=seq(0,60,5),expand=c(0,0)) +
  theme_bw(base_size=20) +
  theme(legend.title=element_blank())+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(legend.position="bottom") +
  theme(strip.background = element_rect(colour="black", fill="white")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
prev.plot

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\BMGF Lancet Paper\\Final_Plots")
ggsave("pppv.plot.2016.jpg", height=8, width=12)

prev.plot <- ggplot(subset(ppdv, Year==1990 | Year==1995 | Year==2000 | Year == 2005 | Year == 2016 | Year == 2030)) +
  geom_line(size=1.2, linetype=1, aes(x = Age, y=ppdv*100, color=Year, group=interaction(Year,sim.id)), alpha=0.03)+ 
  geom_smooth(aes(x=Age, y=ppdv*100, color=Year, group=Year),method="loess", span=0.5, se = F, size=2, linetype=1) +
  scale_colour_gradientn(colours=rainbow(3),limits=c(year1, year2), breaks = c(year1,round((year1+year2)/2,0),year2)) +
  facet_grid(~Gender,labeller=labeller(Gender = labs)) +
  xlab("Age") +
  ylab("% with detectable viral load") +
  scale_x_continuous(limits = c(15, 60), breaks=seq(15,59,5),expand=c(0,0)) +
  scale_y_continuous(limits = c(0, 60), breaks=seq(0,60,5),expand=c(0,0)) +
  theme_bw(base_size=20) +
  theme(legend.title=element_blank())+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(legend.position="bottom") +
  theme(strip.background = element_rect(colour="black", fill="white")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
prev.plot

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\BMGF Lancet Paper\\Final_Plots")
ggsave("pppv.plot.2030.jpg", height=8, width=12)

#calculate population prevalence of unsuppressed viral load
OnARTnumbers <- aggregate(On_ART ~ Year+Gender+Age, subset(reporthivbyageandgender, Age>10 & Age < 100), FUN=sum) #sums prev infections in each year
Population <- aggregate(Population ~ Year+Gender+Age, subset(reporthivbyageandgender, Age>10 & Age < 100), FUN=sum)
Infected <- aggregate(Infected ~ Year+Gender+Age, subset(reporthivbyageandgender, Age>10 & Age < 100), FUN=sum)

ppdv <- merge(OnARTnumbers, Population, by=c("Year","Gender","Age"))
ppdv <- merge(ppdv, Infected, by=c("Year","Gender","Age"))
ppdv$suppressed <- ppdv$On_ART*0.92
ppdv$unsuppr <- ppdv$Infected - ppdv$suppressed
ppdv$ppdv <- ppdv$unsuppr / ppdv$Population
ppdv$prevalence <- ppdv$Infected / ppdv$Population
summary(ppdv$ppdv)
summary(ppdv$prevalence)
head(ppdv)
summary(ppdv$ppdv*100)

#calculate population prevalence of unsuppressed viral load by broad age.group
head(reporthivbyageandgender)
reporthivbyageandgender$age.cat <- cut(reporthivbyageandgender$Age, c(15,25,35,100), right=F)
OnARTnumbers <- aggregate(On_ART ~ Year+Gender+age.cat, subset(reporthivbyageandgender, Age>10 & Age < 100), FUN=sum) #sums prev infections in each year
Population <- aggregate(Population ~ Year+Gender+age.cat, subset(reporthivbyageandgender, Age>10 & Age < 100), FUN=sum)
Infected <- aggregate(Infected ~ Year+Gender+age.cat, subset(reporthivbyageandgender, Age>10 & Age < 100), FUN=sum)

ppdv <- merge(OnARTnumbers, Population, by=c("Year","Gender","age.cat"))
ppdv <- merge(ppdv, Infected, by=c("Year","Gender","age.cat"))
ppdv$suppressed <- ppdv$On_ART*0.92
ppdv$unsuppr <- ppdv$Infected - ppdv$suppressed
ppdv$ppdv <- ppdv$unsuppr / ppdv$Population
ppdv$prevalence <- ppdv$Infected / ppdv$Population
summary(ppdv$ppdv)
summary(ppdv$prevalence)
head(ppdv)
summary(ppdv$ppdv*100)

#over time by age
year1 = 1985
year2 = 2020
summary(ppdv$Age)
ppdv.plot <- ggplot() +
  geom_line(data = subset(ppdv, Year >= year1 & Year <=year2), size=1.2, linetype=1, aes(x = Year, y=ppdv*100, color=factor(age.cat), group=factor(age.cat)))+ 
  #scale_colour_gradientn(colours=rainbow(3)) +
  facet_grid(~Gender,labeller=labeller(Gender = labs)) +
  xlab("Year of HIV acquisition") +
  ylab("Population prevalence of unsuppressed viral load") +
  #scale_x_continuous(limits = c(15, 60), breaks=seq(15,60,5),expand=c(0,0)) +
  #scale_y_continuous(limits = c(0, 60), breaks=seq(0,60,5),expand=c(0,0)) +
  theme_bw(base_size=16) +
  # guides(fill = guide_legend(keywidth = 2, keyheight = 1)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(legend.position="bottom") +
  theme(strip.background = element_rect(colour="black", fill="white")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
ppdv.plot

year1 = 1985
year2 = 2020
limits = c(0,max(ppdv$ppdv))
breaks = seq(0, max(ppdv$ppdv), max(ppdv$ppdv)/6)
breaks=round(breaks, 3)
table(ppdv$Age)

ppdv.raster.plot <-  ggplot() +
  geom_raster(data=subset(ppdv, Year >= year1 & Year <=year2), aes(y=Age, x=Year, fill=ppdv*100))+ 
  scale_fill_gradientn(name="Prevalence",
                       colours=c("darkblue","blue","yellow","darkred"),
                       limits=limits*100,
                       breaks=breaks*100) +
  stat_contour(data=subset(ppdv, Year >= year1 & Year <=year2), aes(y=Age, x=Year, z=ppdv*100),
               color="black", size=0.5, linetype=1, binwidth=5) +
  facet_grid(~Gender,labeller=labeller(Gender = labs)) +
  theme(legend.position="bottom") +
  xlab("Year")+
  ylab("Age")+
  scale_x_continuous(expand=c(0,0), breaks=seq(year1,year2-1,5), limits=c(year1,year2)) +
  scale_y_continuous(expand=c(0,0)) +
  theme_bw(base_size=20) +
  guides(fill = guide_legend(keywidth = 2, keyheight = 1, nrow=1)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(legend.position="bottom") +
  theme(strip.background = element_rect(colour="black", fill="white"))
ppdv.raster.plot

setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\Manuscripts\\Ongoing Manuscripts Folder\\BMGF Lancet Paper\\Final_Plots")
ggsave("pppv.plot.raster.jpg", height=8, width=12)





##########################################
#Extra plots
##########################################
#How many transmissions does each newly infected person produce within the first year?
names(transmissionreport)
transmissionreport$transmit.in.1yr <- 0
transmissionreport$transmit.in.1yr[transmissionreport$yrs.to.transmit<1] <- 1
head(transmissionreport,20)
transmit.1yr <- aggregate(transmit.in.1yr ~ age.cat + src.gender + year.acq.int, transmissionreport, FUN=sum)
nrow(transmit.1yr)
denominator <- aggregate(acquisition ~ age.cat + src.gender + year.acq.int, transmissionreport, FUN=sum)
nrow(denominator)
transmit.1yr <- merge(transmit.1yr, denominator, by=c("age.cat", "src.gender","year.acq.int"))
head(transmit.1yr)
transmit.1yr$prop <- transmit.1yr$transmit.in.1yr / transmit.1yr$acquisition

transmissions.by.year.1yr <- ggplot() +
  geom_line(data = transmit.1yr, size=1.2, linetype=1, aes(x = year.acq.int, y=prop, color=age.cat, group=age.cat))+ 
  #geom_ribbon(data= subset(agg.median.age.at.hiv.m), aes(x = year, ymin=AgeAtVisit[,2], ymax=AgeAtVisit[,3], fill = sex), alpha=0.3)+ 
  #geom_point(data= subset(p.transmit.1yr), aes(x=year.acq.int, y=transmit.in.1yr*100, color=age.cat)) +
  facet_grid(~src.gender) +
  xlab("Year of HIV acquisition") +
  ylab("Number of transmissions per infected individual") +
  scale_x_continuous(limits = c(1990, 2049), breaks=seq(1990,2049,5),expand=c(0,0)) +
  #scale_y_continuous(limits = c(0, 60), breaks=seq(0,60,5),expand=c(0,0)) +
  theme_bw(base_size=16) +
  # guides(fill = guide_legend(keywidth = 2, keyheight = 1)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(legend.position="bottom") +
  theme(strip.background = element_rect(colour="black", fill="white")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
transmissions.by.year.1yr



