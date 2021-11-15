################################################################################
## Adam Akullian and Allen Roberts
## Function to generate summary plots and output from EMOD experiments
## Assumes outputs are stored in "ReportHIVByAgeAndGender"
## Time series of population, prevalence, incidence, and ART coverage with
## various stratifications by age, gender, and risk group
## Outputs are "scenario_summary_plots.pdf" and "scenario_summary.RData" in a given output directory
################################################################################

make_results_summary <- function(scenario, 
                                 suite_name, 
                                 data_dir,
                                 output_dir, 
                                 report_interval = 0.5, 
                                 make_plots = FALSE,
                                 make_vaccine_cohort_plots=FALSE) {
  
  ## Libraries
  library(ggplot2)
  library(tidyverse)
  library(lemon)
  library(parallel)
  library(foreach)
  library(doParallel)
  library(mgsub)
  
  ## Suppress unhelpful warnings
  options(dplyr.summarise.inform = FALSE) 
  
  ## Argument error handling
  if(missing(scenario)) {
    
    stop("Error: No scenario argument supplied")
    
  }
  
  if(missing(suite_name)) {
    
    stop("Error: No suite_name argument supplied")
    
  }
  
  if(missing(data_dir)) {
    
    stop("Error: No data_dir argument supplied")
  }
  
  if(missing(output_dir)) {
    
    stop("Error: No output_dir argument supplied")
  }
  
  ## Set up parallels - this might be redundant here, but it's necessary for the 
  ## functions to run
  numCores <- detectCores()
  registerDoParallel(numCores)
  
  ## Plot settings
  theme_set(theme_classic())
  
  ## Load simulation run - assumes they are stored in ReportHIVByAgeAndGender
  file_path <- file.path(data_dir, suite_name, scenario, "ReportHIVByAgeAndGender")
  file_list <- list.files(file_path, pattern = ".csv")
  
  ## Create output directories
  dir.create(file.path(output_dir, suite_name, scenario), recursive = TRUE)
  
  ## Function to combine results after parallel processing
  comb <- function(x, ...) {  
    mapply(rbind,x,...,SIMPLIFY=FALSE)
  }
  
    outputs <- foreach(i=1:length(file_list), .combine = 'comb', .multicombine = TRUE) %dopar% {
    #outputs <- foreach(i=1:20, .combine = 'comb', .multicombine = TRUE) %dopar% {
    
    ## Reload packages - sometimes necessary when parallel computing
    library(ggplot2)
    library(tidyverse)
    library(lemon)
    
    print(i)
    ff <- file_list[i]
    run_num <- as.integer(gsub("[^0-9]", "", unlist(strsplit(ff, "_"))[2]))
    df <- read.csv(file.path(file_path, ff), stringsAsFactors = FALSE)
    df$run_num <- run_num
    
    # Format variables
    df <- df %>% mutate(age_group = cut(df$Age, breaks = c(10, 14, 19, 24, 29, 34, 39, 44, 49), labels = c("10-14", "15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49")),
             gender = factor(df$Gender, levels = c(0,1), labels = c("Men", "Women")),
             year_floor = floor(Year))  %>%
      rename(year = Year)
    
    df <- df[df$Age>=10 & df$Age<=49,]
    
    # df <- df %>%
    #   mutate(age_group = cut(df$Age, breaks = c(14, 19, 24, 29, Inf), labels = c("15-19", "20-24", "25-29", "30+")),
    #          gender = factor(df$Gender, levels = c(0,1), labels = c("Men", "Women")),
    #          year_floor = floor(Year))  %>%
    #   rename(year = Year)
    # 
    
    ## Person-time at risk
    df$HasHIV <- 1
    df$HasHIV[df$HIV_Stage=="NOT_INFECTED"] <- 0
    df$person_years_at_risk <- ifelse(df$HasHIV == 0, df$Population*report_interval, df$Newly.Infected*report_interval/2)
    
    #create on ART variable
    df$On_ART <- 0
    df$On_ART[df$HIV_Stage=="ON_ART"] <- df$Population[df$HIV_Stage=="ON_ART"]
    
    ## Risk groups
    if("IP_Key.Risk" %in% colnames(df)){
      df$risk <- "LOW"
      df$risk[df$IP_Key.Risk=="MEDIUM"|df$IP_Key.Risk=="HIGH"] <- "MEDIUM/HIGH"
    }
  
    ## Total population over time
    pop <- (df %>%
              group_by(year, run_num) %>%
              summarise(pop = sum(Population))
    )
    
    ## Prevalence, ages 15-49, by gender
    prev_adults_gender <- (df %>%
                             filter(Age >= 15 & Age <= 49) %>%
                             group_by(gender, year, run_num) %>%
                             summarise(prev = 100*sum(Infected)/sum(Population))) 
    
    ## Prevalence - age/gender
    prev_age_gender <- (df %>% 
                          group_by(gender, age_group, year, run_num) %>%
                          summarise(prev = 100*sum(Infected)/sum(Population))
    )
    
    ## Incidence, adults 15-49, by gender
    inc_adults_gender <- (df %>%
                            filter(Age >= 15 & Age <= 49) %>%
                            group_by(gender, year_floor, run_num) %>%
                            summarise(new_infections = sum(Newly.Infected),
                                      person_years = sum(person_years_at_risk)) %>%
                            mutate(inc = 100*new_infections/person_years,
                                   year = year_floor) %>%
                            select(-year_floor)
                          
    )
    
    ## Incidence - age/gender
    inc_age_gender <- (df %>% 
                         group_by(gender, age_group, year_floor, run_num) %>%
                         summarise(new_infections = sum(Newly.Infected),
                                   person_years = sum(person_years_at_risk)) %>%
                         mutate(inc = 100*new_infections/person_years,
                                year = year_floor) %>%
                         select(-year_floor)
    )
    
    ## Incidence - age/gender/risk group
      inc_age_gender_risk <- (df %>% 
                           group_by(gender, age_group, risk, year_floor, run_num) %>%
                           summarise(new_infections = sum(Newly.Infected),
                                     person_years = sum(person_years_at_risk)) %>%
                           mutate(inc = 100*new_infections/person_years,
                                  year = year_floor) %>%
                           select(-year_floor)
      )
    
    ## ART coverage - adults 15-49 by gender
    art_adults_gender <- (df %>% 
                            filter(Age >= 15 & Age <= 49) %>%
                            group_by(gender, year, run_num) %>%
                            summarise(art_cov = 100 * sum(On_ART)/sum(Infected))
    )
    
    ## ART coverage - adults 15-49 by gender
    art_adults_age_gender <- (df %>% 
                            filter(Age >= 15 & Age <= 49) %>%
                            group_by(gender, Age, year, run_num) %>%
                            summarise(art_cov = 100 * sum(On_ART)/sum(Infected))
    )
    
    ## VMMC coverage by age and gender
    vmmc_age <- (subset(df, gender=="Men") %>%
                            filter(Age >= 10 & Age <= 49) %>%
                            group_by(Age, year, run_num) %>%
                            summarise(vmmc_cov = 100 * sum(Population[IsCircumcised==1])/sum(Population))
    )

    return(list("pop" = pop, 
                "prev_adults_gender" = prev_adults_gender, 
                "prev_age_gender" = prev_age_gender, 
                "inc_adults_gender" = inc_adults_gender, 
                "inc_age_gender" = inc_age_gender, 
                "art_adults_gender" = art_adults_gender,
                "art_adults_age_gender" = art_adults_age_gender,
                "vmmc_age"=vmmc_age))
  }
  
  ## Save output
  save(outputs, file = file.path(output_dir, suite_name, scenario, "scenario_summary.RData"))
  
  # ## Bring in reference data on prevalence, incidence, and ART coverage
  # prevdata <- read.csv("..\\Calibration\\Data\\SWAZILAND_nationalprevalence_bothgenders.csv")
  # head(prevdata)
  # table(prevdata$AgeBin)
  # prevdata1549 <- prevdata %>%
  #   rename(gender=Gender, year=Year) %>%
  #   mutate(gender=recode(gender, Male="Men", Female="Women")) %>%
  #   filter(AgeBin %in% c("[15:50)")) %>%
  #   group_by(year, gender) %>%
  #   summarise(alpha=sum(Prevalence*effective_count), beta=sum((1-Prevalence)*effective_count), prev=alpha/(alpha+beta), lb=qbeta(0.025, alpha, beta), ub=qbeta(0.975, alpha, beta))
  # 
  # prevdata_agebin <- prevdata %>%
  #   rename(gender=Gender, year=Year) %>%
  #   mutate(gender=recode(gender, Male="Men", Female="Women")) %>%
  #   mutate(age_group=paste(substr(AgeBin, 2,3),"-",as.numeric(substr(AgeBin, 5,6))-1, sep="")) %>%
  #   filter(!age_group %in% c("15-49","50-54","55-59")) %>%
  #   group_by(year, gender, age_group) %>%
  #   summarise(alpha=sum(Prevalence*effective_count), beta=sum((1-Prevalence)*effective_count), prev=alpha/(alpha+beta), lb=qbeta(0.025, alpha, beta), ub=qbeta(0.975, alpha, beta))
  
  if(make_plots == TRUE) {
    ## Plots
    pdf(file = file.path(output_dir, suite_name, scenario, "scenario_summary_plots.pdf"), 
        height = 5, width = 7)
    
    ## Total population over time
    print(ggplot(data = outputs$pop, aes(x = year, y = pop)) +
            geom_line(aes(group = run_num), alpha = 0.5) +
            labs(x = "Year", y = "Population size") +
            ggtitle("Total population")
    )

    ## Prevalence - adults by gender
    print(ggplot(data = outputs$prev_adults_gender, aes(x = year, y = prev, color = gender)) +
            geom_line(aes(group = interaction(run_num, gender)), alpha = 0.5) +
            #geom_point(data=prevdata1549, aes(x=year, y= prev*100),color="black") +
            #geom_errorbar(data=prevdata1549, aes(x=year, ymin= lb*100, ymax=ub*100),color="black") +
            facet_wrap(~gender) +
            labs(x = "Year", y = "Prevalence (%)", color = "Gender") +
            ggtitle("HIV prevalence, ages 15-49, by gender") +
            guides(colour = guide_legend(override.aes = list(alpha = 1)))
    )

    ## Prevalence by age and gender
    print(ggplot(data = subset(outputs$prev_age_gender), aes(x = year, y = prev)) +
            geom_line(aes(group = interaction(run_num, gender), color=gender),  alpha = 0.5) +
            #geom_point(data=subset(prevdata_agebin, gender=="Men"), aes(x=year, y= prev*100),color="black") +
            #geom_errorbar(data=subset(prevdata_agebin, gender=="Men"), aes(x=year, ymin= lb*100, ymax=ub*100),color="black") +
            labs(x = "Year", y = "Prevalence (%)", color = "Men") +
            facet_rep_wrap(~age_group, repeat.tick.labels = TRUE) +
            ggtitle("HIV prevalence by age and gender") +
            guides(colour = guide_legend(override.aes = list(alpha = 1)))
    )

    ## Incidence - adults 15-49, by gender
    print(ggplot(data = outputs$inc_adults_gender, aes(x = year, y = inc, color = gender)) +
            geom_line(aes(group = interaction(run_num, gender)), alpha = 0.5) +
            labs(x = "Year", y = "Incidence Rate (per 100 PY)", color = "Gender") +
            ggtitle("HIV incidence, ages 15-49 by gender") +
            guides(colour = guide_legend(override.aes = list(alpha = 1)))
    )


    ## Incidence - age and gender
    print(ggplot(data = outputs$inc_age_gender, aes(x = year, y = inc, color = gender)) +
            geom_line(aes(group = interaction(run_num, gender)), alpha = 0.5) +
            labs(x = "Year", y = "Incidence Rate (per 100 PY)", color = "Gender") +
            facet_rep_wrap(~age_group, repeat.tick.labels = TRUE) +
            ggtitle("HIV incidence by age and gender") +
            guides(colour = guide_legend(override.aes = list(alpha = 1)))
    )

# #    if("IP_Key.Risk" %in% colnames(df)){
#     print(ggplot(data=subset(outputs$inc_age_gender_risk), aes(x = year, y = inc, color = as.factor(age_group))) +
#             geom_point()+geom_line()+
#             #geom_ribbon(aes(ymin=IRR_lb, ymax=IRR_ub)) +
#             labs(x = "Year", y = "Incidence Rate (per 100 PY)") +
#             facet_rep_wrap(gender~risk, repeat.tick.labels = TRUE) +
#             ggtitle("HIV incidence by risk-group") +
#             guides(colour = guide_legend(override.aes = list(alpha = 1)))
#     )
# #    }

    if("HasIntervention.Placebo" %in% colnames(df)){
      print(ggplot(data=vaccine_trial_efficacy, aes(x = age_group)) +
              geom_point(aes(y = ve_mean, color = gender))+
              geom_errorbar(aes(ymin=ve_lb, ymax=ve_ub)) +
              #labs(x = "Year", y = "Vaccine efficacy") +
              facet_rep_wrap(~gender, repeat.tick.labels = TRUE) +
              ggtitle("") +
              guides(colour = guide_legend(override.aes = list(alpha = 1)))
      )
    }

    ## ART coverage among adults 15-49 by gender
    print(ggplot(data = outputs$art_adults_gender, aes(x = year, y = art_cov, color = gender)) +
            geom_line(aes(group = interaction(run_num, gender)), alpha = 0.5) +
            labs(x = "Year", y = "ART Coverage (%)", color = "Gender") +
            ggtitle("ART coverage by gender") +
            guides(colour = guide_legend(override.aes = list(alpha = 1)))
    )

    ## ART coverage among adults 15-49 by age and gender
    print(ggplot(data = subset(outputs$art_adults_age_gender), aes(x = year, y = art_cov, color = Age)) +
            geom_line(aes(group = interaction(run_num, Age)), alpha = 0.5) +
            facet_rep_wrap(~gender, repeat.tick.labels = TRUE) +
            labs(x = "Year", y = "ART Coverage (%)", color = "Gender") +
            ggtitle("ART coverage by age and gender") +
            guides(colour = guide_legend(override.aes = list(alpha = 1)))
    )
    
    ## VMMC coverage among men
    print(ggplot(data = subset(outputs$vmmc_age), aes(x = year, y = vmmc_cov)) +
            geom_line(aes(color = Age, group = interaction(run_num, Age)), alpha = 0.5) +
            geom_point(aes(color = Age, group = Age), alpha = 0.5) +
            geom_hline(yintercept=80, linetype="dashed")+
            labs(x = "Year", y = "Circumcision Coverage (%)", color = "Age") +
            scale_colour_gradientn(colours = rainbow(10), breaks=seq(10,45,5)) +
            #facet_rep_wrap(~gender, repeat.tick.labels = TRUE) +
            ggtitle("Circumcision coverage by age-group") +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
            scale_x_continuous(limits=c(1980,2050), breaks=seq(1980,2050,5),expand=c(0,0)) +
            scale_y_continuous(limits=c(0,100), breaks=seq(0,100,10),expand=c(0,0)) +
            guides(colour = guide_legend(override.aes = list(alpha = 1)))
    )
    
    dev.off()
  }  
}
