######################################################################################################
library(kableExtra)
library(tidyverse)
library(janitor)
library(lubridate)
library(MASS)
library(glmnet)
library(parallel)
library(data.table)
library(networkD3)
library(RColorBrewer)
library(MESS)
library(directlabels)

## Load output from EMOD
#setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\GitHub\\EMOD_eswatini\\Output")
setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\GitHub\\EMOD_eswatini\\eSwatini2\\Eswatini_UTT_n10")
load(file = "inputfiles2.rda") #brings in transmissionreport.csv and reporthivbyageandgender.csv

#create "date of acquisition" variable as date when the dest_ID acquired HIV (whether or not there was a subsequent transmission)
date.of.acquisition <- data.frame("SRC_ID"=transmissionreport$DEST_ID, "acq.year"=transmissionreport$YEAR, "sim.id"=transmissionreport$sim.id, "scenario"=transmissionreport$scenario) #source id used to merge in with source of transmitters
network.data <- left_join(transmissionreport, date.of.acquisition, by=c("SRC_ID","sim.id","scenario")) #join to include date of acquisition

network.data.m <- network.data %>%    
  mutate(src.age.cat = cut(src_age_int, c(15,20,25,30,35,40,45,50,55,60,200), right=F)) %>%
  mutate(dest.age.cat = cut(dest_age_int, c(15,20,25,30,35,40,45,50,55,60,200), right=F)) %>%
  mutate(src.age.cat3 = cut(src_age_int, c(15,25,35,200), right=F)) %>%
  mutate(dest.age.cat3 = cut(dest_age_int, c(15,25,35,200), right=F)) %>%
  mutate(acq.age=src_age_int-(YEAR-acq.year)) %>% #age of source when they acquired HIV
  mutate(TimeToTransmission = YEAR - acq.year) %>%
  mutate(MonthsToTransmission = ceiling((YEAR - acq.year)*12)) %>% #month within which transmission occurred (e.g., 1 = transmission ocrurred within < 1 month of acq)
  mutate(YearsToTransmission = ceiling((YEAR - acq.year))) %>% #year within which transmission occurred (e.g., 1 = transmission occurred within < 1 year of acq)
  mutate(era = case_when(Year2 < 2016 ~ "Pre-UTT", TRUE ~ "UTT era")) %>% 
  mutate(Year= Year2) %>%
  mutate(src.stage.f = factor(SRC_STAGE, levels = c(1,2,3,4), labels = c("Acute","Latent","AIDS","on ART"))) %>%
  dplyr::select(-c(SRC_STAGE, Year2)) 
network.data.m$scenario = factor(network.data.m$scenario, levels=c('baseline','ART100pct'))
table(network.data.m$scenario)

labs <- c("1" = "Women", "0" = "Men")
labs.age3 <- c("[15,25)"="15-24","[25,35)"="25-34","[35,200)"="35+")

var_guide_baseline <- data.frame(
  var_name = c(names(network.data.m)),
  descr = c(
    "Source transmission ID",
    "Source sex (0 = Male, 1 = Female)",
    "Destination transmission ID",
    "Destiation sex (0 = Male, 1 = Female)",
    "Age of source (years)",
    "Age of destination (years)",
    "Year in sim time-steps",
    "Simulation ID (N=250)",
    "Model scenario (baseline=status quo ART, ART100pct=All HIV+ initiate ART within ~ 6 mo, viral load suppression within 1 year of acquisition",
    "Year transmitter (dest) acquired HIV (NA if seeded in)",
    "5-yr age category of source",
    "5-yr age category of destination",
    "10-yr age category of source",
    "10-yr age category of destination",
    "Age when transmitter acquired HIV (NA if seeded in)",
    "Time (years with decimal) from HIV acquisition to transmission (from source persepctive)",
    "Time (months) from HIV acquisition to transmission (from source persepctive)",
    "Time (years) within which HIV acquisition to transmission occured (from source persepctive)",
    "era (pre/post UTT)",
    "Year transmission occured",
    "Source HIV stage (Acute = untreated acute stage (2.9 months), Latent=untreated latent HIV (>3 months), AIDS=untreated AIDS, on ART=On ART"
    )
  )

network.vars <-kable(
  var_guide_baseline,
  col.names = c("Variable name",
    "Description")
) %>% 
  kable_styling()
network.vars

names(reporthivbyageandgender)
reporthivbyageandgender <- reporthivbyageandgender %>% 
  mutate(SRC_GENDER = Gender, Year = Year2, NeverOnART = case_when(IP_Key.ARTstate == "NeverOnART" ~ 1, TRUE ~ 0),
         onART = case_when(IP_Key.ARTstate == "OnART" ~ 1, TRUE ~ 0), src.age.cat.art = cut(Age, c(15,25,35,45), right=F),
         src.age.cat3 = cut(Age, c(15,25,35,200), right=F),
         src.age.cat = cut(Age, c(15,20,25,30,35,40,45,50,55,60,200), right=F))  %>%
    dplyr::select(-c(Year2, Gender, Age)) 

var_guide_reporthiv <- data.frame(
  var_name = c(names(reporthivbyageandgender)),
  descr = c(
    "Is Circumcised (1=Yes, 0 =No)",
    "ART status among HIV+ (NeverOnART=Never initiated ART, OffART=Initiated but currently off, OnART=On ART)",
    "HIV stage",
    "Newly infected with HIV (n)",
    "Transmitters (n)",
    "Uninfected (n)",
    "Currently infected with HIV (n)",
    "Total population (N)",
    "Simulation ID (N=250)",
    "Model scenario (baseline=status quo ART, ART100pct=All HIV+ initiate ART within ~ 6 mo, viral load suppression within 1 year of acquisition",
    "Sex (0 = Male, 1 = Female)",
    "Calendar year",
    "Never initiated ART (HIV+)",
    "on ART (HIV+)",
    "Age category (ART)",
    "Age category (5-year bins)",
    "Age category (10-year bins)"
  )
)

reportHIV.vars <- kable(
  var_guide_reporthiv,
  col.names = c("Variable name",
    "Description")
) %>% 
  kable_styling()

#ART coverage
table(reporthivbyageandgender$IP_Key.ARTstate)
artcov <- reporthivbyageandgender %>%
  group_by(SRC_GENDER,src.age.cat.art, Year, scenario, sim.id) %>%
  summarise(cov=sum(Population[onART==1])/sum(Infected))

#Bring in ART coverage data
setwd("C:\\Users\\aakullian\\Dropbox (IDM)\\HIV\\EMOD\\Swaziland_05June2018\\Data\\Swaziland_point_estimates")
infile1 <- "SWAZILAND_calibration_nationalARTprevalence"
artdata <- read.csv(paste0('./', infile1, '.csv'))
artdata$SRC_GENDER <- artdata$Gender
artdata$src.age.cat.art <- artdata$AgeBin
artdata$ART_coverage <- artdata$NationalARTPrevalence
labs <- c("1" = "Women", "0" = "Men")

#ART scale-up to 2016
table(artcov$src.age.cat.art)
table(artdata$src.age.cat.art)
ggplot() +
  geom_point(data=subset(artcov, scenario=="baseline" & is.na(cov)==F), aes(x=Year, y=cov*100, color=src.age.cat.art, group=interaction(src.age.cat.art,sim.id)),alpha=0.1) +
  geom_line(data=subset(artcov, scenario=="baseline" & is.na(cov)==F), aes(x=Year, y=cov*100, color=src.age.cat.art, group=interaction(src.age.cat.art,sim.id)),alpha=0.1) +
  geom_point(data=subset(artdata), aes(x=Year, y=ART_coverage*100, color=factor(src.age.cat.art)), size=2) +
  geom_errorbar(data = subset(artdata), aes(x=Year, ymin=lb*100, ymax=ub*100, color=factor(src.age.cat.art)), width=1, size=1.2) +
  facet_grid(~SRC_GENDER) +
  theme(legend.position="bottom") +
  ylab("ART Coverage (%)")+
  theme_bw(base_size=18) +
  #scale_y_continuous(limits=c(0,100), breaks = seq(0,90,10), expand=c(0,0)) +
  scale_x_continuous(limits = c(2000,2050), breaks = seq(2000,2050,5), expand=c(0,0)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(strip.background = element_rect(colour="black", fill="white")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"))

#proportion of prevalent infections by stage
head(reporthivbyageandgender)
table(reporthivbyageandgender$src.stage.f)
reporthivbyageandgender <- reporthivbyageandgender %>% mutate(src.stage.f = factor(HIV_Stage, levels = c("ACUTE","AIDS","LATENT","NOT_INFECTED","ON_ART"), labels = c("Acute", "AIDS", "Latent","Not Infected","on ART"))) 
reporthivbyageandgender$src.stage.f <- factor(reporthivbyageandgender$src.stage.f, levels = c("Acute", "Latent", "AIDS","on ART","Not Infected"))

prev.stage <- subset(reporthivbyageandgender)  %>%
  group_by(scenario, sim.id, Year, src.stage.f) %>%
  summarize(Infected=sum(Infected), Transmitters=sum(Transmitters), Uninfected=sum(Uninfected.Population),Population=sum(Population)) %>%
  mutate(transm_per_infected=(Transmitters/Infected), na.rm=T) %>% #divide by 100 so can plot per 100 across all metrics and still have it be per 1 infection
  group_by(scenario, sim.id, Year) %>%
  mutate(prev=Infected/sum(Population), na.rm=T, inc=Transmitters/sum(Uninfected), na.rm=T, propinc=Transmitters/sum(Transmitters)) %>%
  filter(src.stage.f != "Not Infected") %>% droplevels() %>%
  group_by(scenario, Year, src.stage.f) %>%
  summarize(med_prev=median(prev, na.rm=T), lb_prev = quantile(prev, probs = 0.025, na.rm=T), ub_prev = quantile(prev, probs = 0.975, na.rm=T),
            med_inc=median(inc, na.rm=T), lb_inc = quantile(inc, probs = 0.025, na.rm=T), ub_inc = quantile(inc, probs = 0.975, na.rm=T),
            med_transmperinfected=median(transm_per_infected, na.rm=T), lb_transmperinfected = quantile(transm_per_infected, probs = 0.025, na.rm=T),
            ub_transmperinfected = quantile(transm_per_infected, probs = 0.975, na.rm=T),
            med_propinc=median(propinc, na.rm=T), lb_propinc = quantile(propinc, probs = 0.025, na.rm=T),ub_propinc = quantile(propinc, probs = 0.975, na.rm=T)) %>%
  pivot_longer(
    !scenario:src.stage.f, 
    names_to = c(".value", "metric"), 
    names_sep = "_", 
    values_drop_na = TRUE
  )

prev.stage$metric_f = factor(prev.stage$metric, levels=c('prev','transmperinfected','inc','propinc'))
labs=c('prev'="Prevalence (%)",'transmperinfected'="Transmissions (per 100 prevalent infections)",'inc'="Incidence (per 100 py)",'propinc'="Attributable fraction (% of incidence)")

ggplot(subset(prev.stage,scenario=="baseline" & !(metric=="transmperinfected" & Year==2004) & Year>1989), aes(x=Year)) + 
  geom_line(aes(y=med*100, color=src.stage.f), size=2) +
  geom_ribbon(aes(ymin=lb*100, ymax=ub*100, fill=src.stage.f), size=2, alpha=0.2) +
  scale_x_continuous(expand = c(0,0), breaks=seq(1980,2045,5)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_color_viridis_d(guide = guide_legend(title = "Scenario")) +
  scale_fill_viridis_d(guide = guide_legend(title = "Scenario")) +
  facet_wrap(~metric_f, scales="free_y", labeller=labeller(metric_f = labs)) +
  #coord_cartesian(xlim=c(2000, 2050)) + 
  ggtitle("") + ylab("")+
  theme(strip.background = element_rect(colour="white", fill="white"), axis.text.x = element_text(angle = 90),legend.title=element_blank(),axis.line = element_line(colour = "black"),panel.background = element_blank(),legend.position = "bottom") +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+ annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)

#Plots are stratified on UTT era - chose one and then plot

transmission.table.simple <- subset(network.data.m, era=="UTT era")  %>%
  group_by(SRC_GENDER, DEST_GENDER, src.age.cat3, dest.age.cat3) %>%
  summarize(transmissions=n())

##
transmission.table.simple$source.gender.age <- paste(transmission.table.simple$SRC_GENDER, transmission.table.simple$src.age.cat3)
transmission.table.simple$destination.gender.age <- paste(transmission.table.simple$DEST_GENDER, transmission.table.simple$dest.age.cat3)

flow.data <- data.frame("orig_reg" = transmission.table.simple$source.gender.age, "dest_reg" = transmission.table.simple$destination.gender.age, "Transmissions" = transmission.table.simple$transmissions)
table(flow.data$orig_reg)

flow.data$source <- as.numeric(as.factor(flow.data$orig_reg))-1 #must be zero indexed so subtract 1
flow.data$target <- as.numeric(as.factor(flow.data$dest_reg))-1 #must be zero indexed so subtract 1
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

#Time to transmission
network.data.m$TimeToTransmission <- network.data.m$YEAR - network.data.m$acq.year #time from HIV acquisition until transmission event
hist(network.data.m$TimeToTransmission[network.data.m$SRC_GENDER==0])
hist(network.data.m$TimeToTransmission[network.data.m$SRC_GENDER==1])

TimeToTransmissionIQR <- subset(network.data.m, scenario=="baseline" & is.na(TimeToTransmission)==F) %>% group_by(era, SRC_GENDER) %>%
  summarize(median = quantile(TimeToTransmission, probs = 0.5),
            lb = quantile(TimeToTransmission, probs = 0.25),
            ub = quantile(TimeToTransmission, probs = 0.75))
head(TimeToTransmissionIQR)

labs <- c("1" = "Women", "0" = "Men")
ggplot(subset(network.data.m, scenario=="baseline"), aes(x=TimeToTransmission, group=interaction(sim.id,SRC_GENDER), color=factor(SRC_GENDER))) +
  geom_density(alpha=0.2, adjust=1/10) +
  #geom_density(data=subset(network.data.m, scenario=="baseline" & TimeToTransmission<1),alpha=0.2, adjust=1/10) +
  facet_grid(era~SRC_GENDER, labeller=labeller(SRC_GENDER = labs)) +
  scale_color_manual(values = c("blue3","red2"), labels=c("Men","Women")) +
  #scale_fill_manual(values = c("blue3","red2"), labels=c("Men","Women")) +
  ggtitle("Distribution of time to transmission") +
  theme_bw(base_size=24) +  xlab("Years")+ 
  scale_x_continuous(breaks=seq(0,39,5),limits=c(0,40)) + 
  #scale_y_continuous(expand = c(0,0), breaks=seq(0,0.15,0.05)) + coord_cartesian(ylim=c(0, 0.16)) +
  geom_vline(data = subset(TimeToTransmissionIQR), aes(xintercept = median, color=factor(SRC_GENDER)),linetype=2) +
  geom_vline(data = subset(TimeToTransmissionIQR), aes(xintercept = lb, color=factor(SRC_GENDER)), alpha=0.2,linetype=2) +
  geom_vline(data = subset(TimeToTransmissionIQR), aes(xintercept = ub, color=factor(SRC_GENDER)), alpha=0.2,linetype=2) +
  theme(strip.background = element_rect(colour="black", fill="white")) +
  guides(fill = guide_legend(reverse=TRUE, nrow=1)) + labs(fill = "Years from HIV acquisition") +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), legend.title=element_blank(),
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),
        panel.background = element_blank(), legend.position = "bottom") +  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)

# Month by month time to transmission density (among all transmission)
labs <- c("1" = "Women", "0" = "Men")
ggplot(subset(network.data.m, scenario=="baseline" & MonthsToTransmission<=12), aes(x=MonthsToTransmission, group=interaction(sim.id,SRC_GENDER), color=factor(SRC_GENDER))) +
  geom_density(alpha=0.2, adjust=1/10) +
  facet_grid(era~SRC_GENDER, labeller=labeller(SRC_GENDER = labs)) +
  scale_color_manual(values = c("blue3","red2"), labels=c("Men","Women")) +
  #scale_fill_manual(values = c("blue3","red2"), labels=c("Men","Women")) +
  ggtitle("Distribution of time to transmission") +
  theme_bw(base_size=24) +  xlab("Years")+ 
  #scale_x_continuous(expand = c(0,0),breaks=seq(0,39,5),limits=c(0,40)) + 
  #scale_y_continuous(expand = c(0,0), breaks=seq(0,0.15,0.05)) 
  theme(strip.background = element_rect(colour="black", fill="white")) +
  guides(fill = guide_legend(reverse=TRUE, nrow=1)) + labs(fill = "Years from HIV acquisition") +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), legend.title=element_blank(),
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),
        panel.background = element_blank(), legend.position = "bottom") +  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)

#proportion of HIV transmissions from transmitters with < 1 year from infection by age/year continuous
network.data.m$TimeToTransmission.1yr <- 0
network.data.m$TimeToTransmission.1yr[network.data.m$TimeToTransmission<1] <- 1
network.data.m$TimeToTransmission.1yr[is.na(network.data.m$TimeToTransmission)==T] <- NA
table(network.data.m$TimeToTransmission.1yr)

head(network.data.m)
network.data.m$yearcat5 <- cut(network.data.m$Year, seq(1980,2050,5), right=F)
table(network.data.m$yearcat5, network.data.m$Year)

network.data.m$yearcat5num <- unclass(network.data.m$yearcat5)
network.data.m$src.age.catnum <- unclass(network.data.m$src.age.cat)
network.data.m$dest.age.catnum <- unclass(network.data.m$dest.age.cat)
table(network.data.m$dest.age.catnum, network.data.m$dest.age.cat)
table(network.data.m$src.age.catnum, network.data.m$src.age.cat)

TimeToTransmissionProportions.sources <- subset(network.data.m,  is.na(TimeToTransmission)==F) %>%
  group_by(scenario, sim.id, yearcat5num, src.age.catnum, SRC_GENDER) %>%
  summarise(TimeToTransmission1yr = mean(TimeToTransmission.1yr)) %>%
  group_by(scenario, yearcat5num, src.age.catnum, SRC_GENDER) %>%
  summarize(median = quantile(TimeToTransmission1yr, probs = 0.5),
            lb = quantile(TimeToTransmission1yr, probs = 0.025),
            ub = quantile(TimeToTransmission1yr, probs = 0.975)) %>%
  mutate(age=src.age.catnum, sex=SRC_GENDER, type="source")

TimeToTransmissionProportions.dest <- subset(network.data.m,  is.na(TimeToTransmission)==F) %>%
  group_by(scenario, sim.id, yearcat5num, dest.age.catnum, DEST_GENDER) %>%
  summarise(TimeToTransmission1yr = mean(TimeToTransmission.1yr)) %>%
  group_by(scenario, yearcat5num, dest.age.catnum, DEST_GENDER) %>%
  summarize(median = quantile(TimeToTransmission1yr, probs = 0.5),
            lb = quantile(TimeToTransmission1yr, probs = 0.025),
            ub = quantile(TimeToTransmission1yr, probs = 0.975)) %>%
  mutate(age=dest.age.catnum, sex=DEST_GENDER, type="dest")

TimeToTransmissionProportions.c <- rbind(TimeToTransmissionProportions.sources, TimeToTransmissionProportions.dest)

years = unique(network.data.m$yearcat5)
ages = unique(network.data.m$src.age.cat)

labs <- c("1" = "Women", "0" = "Men")
data=TimeToTransmissionProportions.c
ggplot(data) +
  geom_tile(data=data, aes(y=age, x=yearcat5num, fill=median))+ 
  scale_fill_gradient(name="Proportion of transmissions within 1 year of infection", low="light blue", high="dark red") +
  facet_grid(type~sex, labeller=labeller(sex = labs)) +
  stat_contour(data=data, aes(y=age, x=yearcat5num, z=median),
               color="black", size=1, linetype=1, breaks=c(0.25,0.5,0.75)) +
  theme(legend.position="bottom") +
  xlab("year")+ylab("age")+
  scale_x_continuous(limits=c(2,14), expand=c(0,0),breaks=seq(2,13,1)) +
  scale_y_continuous(expand=c(0,0), breaks=seq(1,10,1)) +
  guides(fill = guide_legend(keywidth = 2, keyheight = 1, nrow=1)) +
  theme_bw(base_size=24) +
  theme(legend.position="bottom") +
  theme(strip.background = element_rect(colour="black", fill="white"))+
  theme(legend.title=element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"))


labs <- c("1" = "Women", "0" = "Men")
data=TimeToTransmissionProportions.sources
ggplot() +
  geom_tile(data=data, aes(y=src.age.catnum, x=yearcat5num, fill=median))+ 
  scale_fill_gradient(name="Proportion of transmissions within 1 year of infection", low="light blue", high="dark red") +
  facet_grid(~SRC_GENDER, labeller=labeller(SRC_GENDER = labs)) +
  stat_contour(data=data, aes(y=src.age.catnum, x=yearcat5num, z=median),
               color="black", size=1, linetype=1, breaks=c(0.25,0.5,0.75)) +
  #geom_text_contour(aes(y=src.age.catnum, x=yearcat5num, z=median), size=7, rotate=F) +
  theme(legend.position="bottom") +
  xlab("year")+ylab("age")+
  #scale_x_continuous(expand=c(0,0), breaks=seq(1980,2050,5)) +
  #scale_y_continuous(expand=c(0,0), breaks=seq(15,100,5)) +
  guides(fill = guide_legend(keywidth = 2, keyheight = 1, nrow=1)) +
  theme_bw(base_size=24) +
  theme(legend.position="bottom") +
  theme(strip.background = element_rect(colour="black", fill="white"))+
  theme(legend.title=element_blank()) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"))

data=TimeToTransmissionProportions.dest
ggplot() +
  geom_tile(data=data, aes(y=dest.age.catnum, x=yearcat5num, fill=median))+ 
  scale_fill_gradient(name="Proportion of transmissions within 1 year of infection", low="light blue", high="dark red") +
  facet_grid(~DEST_GENDER, labeller=labeller(DEST_GENDER = labs)) +
  stat_contour(data=data, aes(y=dest.age.catnum, x=yearcat5num, z=median),
               color="black", size=1, linetype=1, breaks=c(0.25,0.5,0.75)) +
  theme(legend.position="bottom") +
  xlab("year")+ylab("age")+
  #scale_x_continuous(expand=c(0,0), breaks=seq(1980,2050,5)) +
  #scale_y_continuous(expand=c(0,0), breaks=seq(15,100,5)) +
  guides(fill = guide_legend(keywidth = 2, keyheight = 1, nrow=1)) +
  theme_bw(base_size=24) +
  theme(legend.position="bottom") +
  theme(strip.background = element_rect(colour="black", fill="white"))+
  theme(legend.title=element_blank()) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) 

ggplot(subset(TimeToTransmissionProportions, TimeToTransmissionCat=="[0,1)" & scenario=="baseline")) + 
  geom_point(aes(y=median, x=Year, color=factor(SRC_GENDER), group=factor(SRC_GENDER)), size=2) +
  geom_line(aes(y=median, x=Year, color=factor(SRC_GENDER)), size=2) +
  #geom_errorbar(aes(x=Year, ymin = lb, ymax = ub, color=factor(SRC_GENDER)), width = 0.2) +
  geom_ribbon(aes(ymin=lb, ymax=ub, x=Year, fill=factor(SRC_GENDER), group=factor(SRC_GENDER)),alpha=0.3) +
  facet_grid(~src.age.cat3, labeller = labeller(src.age.cat3=labs.age3))+
  scale_x_continuous(expand = c(0,0), breaks=seq(1980,2050,5)) +
  scale_y_continuous(expand = c(0,0), breaks=seq(0,1,0.05)) + 
  coord_cartesian(ylim=c(0, 1)) +  theme_bw(base_size=20) +
  xlab("Year")+ ylab("Proportion") +
  scale_color_manual(values = c("blue", "red"), labels = c("Male", "Female"))+
  scale_fill_manual(values = c("blue", "red"), labels = c("Male", "Female"))+
  ggtitle("Proportion of HIV transmissions from transmitters infected for < 1 year") +
  theme(strip.background = element_rect(colour="black", fill="white")) +
  theme(panel.border = element_blank(),
        legend.title=element_blank(),
        panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 90),
        axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        legend.position = "bottom")

#########################################################################################
t1 <- network.data.m %>% #count number of transmissions by scenario, sim.id, gender, acquisition year, and months to transmission
  filter(!is.na(acq.year)) %>% mutate(acq.year = floor(acq.year)) %>%
  group_by(scenario, sim.id, SRC_GENDER, acq.year, MonthsToTransmission) %>%
  summarize(transmissions = n()) 
t1.all <- network.data.m %>% #now do the same but across both genders
  filter(!is.na(acq.year)) %>% mutate(acq.year = floor(acq.year)) %>%
  group_by(scenario, sim.id, acq.year, MonthsToTransmission) %>%
  summarize(transmissions = n()) 
#create data frame to calculate the number of individuals who acquire HIV (index infections) by year and gender
date.of.acquisition2 <- data.frame("SRC_ID"=transmissionreport$DEST_ID, "SRC_GENDER"=transmissionreport$DEST_GENDER, "acq.year"=transmissionreport$YEAR,  "sim.id"=transmissionreport$sim.id, "scenario"=transmissionreport$scenario) #source id used to merge in with source of transmitters
#create data frame to calculate the number of individuals who acquire HIV (index infections) by year (both genders)
date.of.acquisition2.all <- data.frame("SRC_ID"=transmissionreport$DEST_ID, "acq.year"=transmissionreport$YEAR,  "sim.id"=transmissionreport$sim.id, "scenario"=transmissionreport$scenario) #source id used to merge in with source of transmitters, this is aggregated across both genders
#count number of index acquisitions by year
r1 <- date.of.acquisition2  %>% mutate(acq.year = floor(acq.year)) %>% 
  group_by(scenario, sim.id, SRC_GENDER, acq.year) %>% summarize(count=n()) 
r1.all <- date.of.acquisition2.all  %>% mutate(acq.year = floor(acq.year)) %>% #do for both genders
  group_by(scenario, sim.id, acq.year) %>% summarize(count=n()) 
#combine number of index acquisitions with number of transmissions
ttable <- left_join(t1,r1,by=c("scenario","sim.id","SRC_GENDER","acq.year"))
ttable.all <- left_join(t1.all,r1.all,by=c("scenario","sim.id","acq.year"))
ttable <- ttable %>% mutate(reff.transperperson = transmissions/count) 
ttable.all <- ttable.all %>% mutate(reff.transperperson = transmissions/count, SRC_GENDER=2)  #calculate reff as transmissions per number of index acquisitions
ttable.merge <- rbind(ttable,ttable.all)

timeframe = 30*12 #time horizon in months to calculate Reff (default is 30 years * 12 months)
#fill in missing month and add zeros if there were not transmission events during that month (this takes about 1 min per year)
ttable.selection1 <- complete(subset(ttable.merge, MonthsToTransmission<timeframe & acq.year%in%c(1990, 2000, 2010, 2020)), 
                             MonthsToTransmission = 1:timeframe, fill = list(reff.transperperson = 0))

ttable.selection <- ttable.selection1  %>%
  group_by(scenario, sim.id, SRC_GENDER, acq.year) %>% #get cumulative sum
  mutate(reffx.cum = cumsum(reff.transperperson), scenario_f = factor(scenario, levels=c('baseline',"ART100pct"))) %>%
  group_by(scenario, sim.id, SRC_GENDER, acq.year) %>%
  gather(Metric, Value, reff.transperperson:reffx.cum)

ttable.selection.med <- ttable.selection %>% group_by(scenario, SRC_GENDER, acq.year, MonthsToTransmission, Metric) %>% summarize(med.Value=mean(Value))

#Cumulative Reff over time (years)
labs <- c("2"="All", "1" = "Women", "0" = "Men")
labs.metric <- c("reff.transperperson"="Year-specific R effective", "reffx.cum" = "cumulative R-effective")
labs.scenario <- c("baseline"="90-90-90", "ART100pct" = "100% ART init.")
ggplot(subset(ttable.selection, scenario=="baseline" & acq.year==2010)) + 
  #geom_point(aes(y=Value, x=MonthsToTransmission/12, color=factor(SRC_GENDER), group=interaction(sim.id, SRC_GENDER)), size=1,  alpha=0.01) +
  geom_line(aes(y=Value, x=MonthsToTransmission/12, color=factor(SRC_GENDER), group=interaction(sim.id, SRC_GENDER)), size=1, alpha=0.01) +
  #geom_line(data=subset(ttable.selection.med, scenario=="baseline" & acq.year==2010), aes(y=med.Value, x=MonthsToTransmission/12, color=factor(SRC_GENDER), group=(SRC_GENDER)), size=1) +
  #geom_smooth(aes(x=MonthsToTransmission/12, y=Value, group=factor(SRC_GENDER), color=factor(SRC_GENDER)),method="loess", span=0.05, se = F, size=2, linetype=1) +
  facet_wrap(scenario_f~Metric, labeller=labeller(scenario_f=labs.scenario), scales="free")+
  #geom_hline(yintercept=1, linetype=2) +
  scale_x_continuous(breaks=seq(0,timeframe/12,1)) + 
  #scale_y_continuous(expand = c(0,0)) + #, breaks=seq(0,1.5,0.1)) + coord_cartesian(ylim=c(0, 1.5)) +
  xlab("Years since HIV acquisition")+ ylab("Reff") +
  scale_color_manual(values = c("blue", "red","purple"), labels = c("Male", "Female","All"))+
  #scale_fill_manual(values = c("blue", "red"), labels = c("Male", "Female"))+
  ggtitle("Number of secondary transmissions per infected individual by year from HIV acquisition") +
  theme(strip.background = element_rect(colour="black", fill="white")) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        legend.title=element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        legend.position = "bottom") +  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)

#Reff by month and cumulative Reff over time (months)
ggplot(subset(ttable.selection, MonthsToTransmission<13 & acq.year==2010)) + 
  #geom_point(aes(y=Value, x=MonthsToTransmission, color=factor(SRC_GENDER), group=interaction(sim.id, SRC_GENDER)), size=1,  alpha=0.03) +
  geom_line(aes(y=Value, x=MonthsToTransmission, color=factor(SRC_GENDER), group=interaction(sim.id, SRC_GENDER)), size=1, alpha=0.03) +
  #geom_line(data=subset(ttable.selection.med,MonthsToTransmission<25& scenario=="baseline"), aes(y=med.Value, x=MonthsToTransmission, color=factor(SRC_GENDER), group=interaction(SRC_GENDER)), size=1) +
  #geom_smooth(aes(x=MonthsToTransmission, y=Value, group=factor(SRC_GENDER), color=factor(SRC_GENDER)),method="loess", span=0.2, se = F, size=2, linetype=1) +
  facet_grid(acq.year~Metric, scales="free", labeller=labeller(Metric=labs.metric, scenario_f=labs.scenario))+
  #geom_hline(yintercept=1) +
  scale_x_continuous(expand = c(0,0), breaks=seq(0,timeframe,1)) + 
  scale_y_continuous(expand = c(0,0)) +
  xlab("Months since HIV acquisition")+ ylab("Reff") +
  scale_color_manual(values = c("blue", "red","purple"), labels = c("Male", "Female","All"))+
  #scale_fill_manual(values = c("blue", "red"), labels = c("Male", "Female"))+
  ggtitle("Number of transmissions per infected individual by month from HIV acquisition") +
  theme(strip.background = element_rect(colour="black", fill="white")) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        legend.title=element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        legend.position = "bottom") +  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)

########################################################################################################
#iterate R effective over all year/time horizon combinations
t1 <- network.data.m %>% 
  filter(!is.na(acq.year)) %>% mutate(acq.year = floor(acq.year)) %>%
  group_by(scenario, sim.id, SRC_GENDER, acq.year, YearsToTransmission) %>%
  summarize(transmissions = n()) 
t1.all <- network.data.m %>% #now aggregate across both genders
  filter(!is.na(acq.year)) %>% mutate(acq.year = floor(acq.year)) %>%
  group_by(scenario, sim.id, acq.year, YearsToTransmission) %>%
  summarize(transmissions = n()) 
date.of.acquisition2 <- data.frame("SRC_ID"=transmissionreport$DEST_ID, "SRC_GENDER"=transmissionreport$DEST_GENDER, "acq.year"=transmissionreport$YEAR,  "sim.id"=transmissionreport$sim.id, "scenario"=transmissionreport$scenario) #source id used to merge in with source of transmitters
date.of.acquisition2.all <- data.frame("SRC_ID"=transmissionreport$DEST_ID, "acq.year"=transmissionreport$YEAR,  "sim.id"=transmissionreport$sim.id, "scenario"=transmissionreport$scenario) #source id used to merge in with source of transmitters, this is aggregated across both genders
r1 <- date.of.acquisition2  %>% mutate(acq.year = floor(acq.year)) %>%
  group_by(scenario, sim.id, SRC_GENDER, acq.year) %>% summarize(count=n()) 
r1.all <- date.of.acquisition2.all  %>% mutate(acq.year = floor(acq.year)) %>% #do for both genders
  group_by(scenario, sim.id, acq.year) %>% summarize(count=n()) 
ttable <- left_join(t1,r1,by=c("scenario","sim.id","SRC_GENDER","acq.year"))
ttable.all <- left_join(t1.all,r1.all,by=c("scenario","sim.id","acq.year"))
ttable <- ttable %>% mutate(transperperson = transmissions/count) 
ttable.all <- ttable.all %>% mutate(transperperson = transmissions/count, SRC_GENDER=2)  
ttable.merge <- rbind(ttable,ttable.all)

timeframe = 30 #time horizon in years to calculate Reff
ttable.years.selection <- complete(subset(ttable.merge, YearsToTransmission<timeframe+1 & acq.year%in%seq(1980,2020,5)), YearsToTransmission = 1:timeframe, fill = list(transperperson = 0)) #fill in missing month and add zeros if there were not transmission events during that month
ttable.yearscumsum.cont <- 
  ttable.years.selection %>% 
  group_by(scenario, sim.id, SRC_GENDER, acq.year) %>% #get cumulative sum
  mutate(reff = cumsum(transperperson), scenario_f = factor(scenario, levels=c('baseline',"ART100pct"))) %>%
  group_by(scenario_f, SRC_GENDER, acq.year, YearsToTransmission) %>%    
  summarize(med.transperperson=mean(transperperson), med.reff=mean(reff))  #get mean across 250 simulations 
#categorize continuous yeartstotransmission into fewer categories for better plotting
cuts <- c(0,1,2,3,4,5,10,15,20,25,30)
ttable.yearscumsum.cont$med.transperperson.cat <- cut(ttable.yearscumsum.cont$YearsToTransmission, cuts, right=T)
ttable.yearscumsum.cont$TransmissionYear <- ttable.yearscumsum.cont$YearsToTransmission + ttable.yearscumsum.cont$acq.year
ttable.yearscumsum <-  ttable.yearscumsum.cont %>% 
  group_by(scenario_f, SRC_GENDER, acq.year, med.transperperson.cat) %>% 
  summarise(med.transperperson =sum(med.transperperson), med.reff = max(med.reff), YearsToTransmission=max(YearsToTransmission))
ttable.yearscumsum$YearsToTransmission.rev <- factor(ttable.yearscumsum$YearsToTransmission, levels = rev(levels(factor(ttable.yearscumsum$YearsToTransmission))))

#Reff by year from acquisition and acquisition year
colfunc <- colorRampPalette(c("white", "white"))
colfunc(10)
labs <- c("2"="All", "1" = "Women", "0" = "Men")
ggplot(subset(ttable.yearscumsum)) + 
  geom_area(aes(y=med.transperperson, x=acq.year, fill=(YearsToTransmission.rev)), size=1, position="stack") +
  geom_line(data=subset(ttable.yearscumsum, YearsToTransmission==30), aes(y=med.reff, x=acq.year), color="black", size=1) +
  facet_grid(scenario_f~SRC_GENDER, labeller=labeller(SRC_GENDER=labs, scenario_f=labs.scenario))+
  #scale_x_continuous(expand = c(0,0), breaks=c(1980,2020,5)) + 
  scale_y_continuous(expand = c(0,0), breaks=seq(0,2.9,0.5)) + coord_cartesian(ylim=c(0,3)) +
  xlab("Year of index HIV acquisition")+ ylab("R effective") +
  geom_hline(yintercept=1, linetype=2) +   theme(text=element_text(size=24))+
  #scale_color_manual(values = c("blue", "red"), labels = c("Male", "Female"))+
  scale_fill_manual(values=colfunc(10))+
  ggtitle("Number of HIV transmissions per infected by year of index acquisition and time from acquisition") +
  theme(strip.background = element_rect(colour="black", fill="white")) +
  guides(fill = guide_legend(reverse=TRUE, nrow=1)) + labs(fill = "Years from HIV acquisition") +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),
        panel.background = element_blank(), legend.position = "bottom") +  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)

labs <- c("2"="All", "1" = "Women", "0" = "Men")
ggplot(subset(ttable.yearscumsum)) + 
  geom_area(aes(y=med.transperperson, x=acq.year, fill=factor(YearsToTransmission.rev)), size=1, position="stack") +
  geom_line(data=subset(ttable.yearscumsum, YearsToTransmission==30), aes(y=med.reff, x=acq.year), color="black", size=1) +
  facet_grid(scenario_f~SRC_GENDER, labeller=labeller(SRC_GENDER=labs, scenario_f=labs.scenario))+
  #scale_x_continuous(expand = c(0,0), breaks=c(1980,2020,5)) + 
  scale_y_continuous(expand = c(0,0), breaks=seq(0,2.9,0.5)) + coord_cartesian(ylim=c(0,3)) +
  xlab("Year of index HIV acquisition")+ ylab("R effective") +
  geom_hline(yintercept=1, linetype=2) +   theme(text=element_text(size=24))+
  #scale_color_manual(values = c("blue", "red"), labels = c("Male", "Female"))+
  #scale_fill_manual(values = c("blue", "red"), labels = c("Male", "Female"))+
  ggtitle("Number of HIV transmissions per infected by year of index acquisition and time from acquisition") +
  theme(strip.background = element_rect(colour="black", fill="white")) +
  guides(fill = guide_legend(reverse=TRUE, nrow=1)) + labs(fill = "Years from HIV acquisition") +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),
        panel.background = element_blank(), legend.position = "bottom") +  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)

#Reff cross section
labs <- c("2"="All", "1" = "Women", "0" = "Men")
ggplot(subset(ttable.yearscumsum, SRC_GENDER==2 & acq.year%in%c(1990,2000,2020))) + 
  geom_line(aes(y=med.reff, x=YearsToTransmission, group=acq.year), size=1) +
  geom_point(aes(y=med.reff, x=YearsToTransmission, color=factor(YearsToTransmission.rev)), size=4) +
  facet_grid(scenario_f~SRC_GENDER, labeller=labeller(SRC_GENDER=labs, scenario_f=labs.scenario))+
  scale_x_continuous(breaks=cuts) + scale_y_continuous(expand = c(0,0), breaks=seq(0,2,0.5)) + coord_cartesian(ylim=c(0, 1.9)) +
  xlab("Year of HIV acquisition")+ ylab("R effective") +
  geom_hline(yintercept=1, linetype=2) +   theme(text=element_text(size=24))+
  #scale_color_manual(values = c("blue", "red"), labels = c("Male", "Female"))+
  #scale_fill_manual(values = c("blue", "red"), labels = c("Male", "Female"))+
  ggtitle("Number of transmissions per infected individual by year of HIV acquisition and time from acquisition") +
  theme(strip.background = element_rect(colour="black", fill="white")) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),
        panel.background = element_blank(), legend.position = "none") +  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)


