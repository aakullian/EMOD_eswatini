################################################################################
## Adam Akullian and Allen Roberts
## Compare simulation results across scenarios
## Uses summary outputs computed for each scenario across simulation runs
## Saves a dataset (scenario_comparison.RData) and plots (scenario_comparison_plots.pdf)
################################################################################

#rm(list = ls())

# Libraries
library(ggplot2)
library(tidyverse)
library(viridisLite)
library(parallel)
library(foreach)
library(doParallel)

## Plot settings
theme_set(theme_classic())

## Suppress uninformative warnings
options(dplyr.summarise.inform=F)

## Function to combine results after parallel processing
comb <- function(x, ...) {  
  mapply(rbind,x,...,SIMPLIFY=FALSE)
}

## Set up parallels
numCores <- detectCores()
registerDoParallel(numCores)

## Output directory - where scenario results are stored
output_dir <- "SummaryOutput"

## Suite name. This example assumes that all the scenario results are located within the
## same suite, but it ideally would allow arbitrary specification of 
## simulations (eg, across suites, or only certain runs within a given suite)
suite_name <- "Eswatini_Age1519_CampaignDuration_Sweep"

## Identify scenarios for a given suite
## In this test example, they are stored locally, but 
## they could be located on COMPS
scenarios <- data.frame("name" = list.dirs(file.path(output_dir, suite_name), full.names = FALSE, recursive = FALSE))

#########################################################################################################################
## For each scenario, compute summary statistics and return outputs
#########################################################################################################################
outputs <- foreach(i=1:nrow(scenarios), .combine = 'comb', .multicombine = TRUE) %dopar% {
  
  library(ggplot2)
  library(tidyverse)
  library(viridisLite)
  library(parallel)
  library(foreach)
  library(doParallel)
  
  scenario_name <- scenarios$name[i]
  
  ## Load outputs
  load(file.path(output_dir, suite_name, scenario_name, "scenario_summary.RData"))
  
  ## Prevalence, ages 15-49, by gender
  prev_adults_gender <- outputs$prev_adults_gender %>%
                           group_by(gender, year) %>%
                           summarise(prev_mean = mean(prev),
                                     prev_median = median(prev),
                                     prev_lower = quantile(prev, probs = 0.025),
                                     prev_upper = quantile(prev, probs = 0.975)) %>%
                           mutate(scenario = sub(".*__", "", scenario_name))
  
  ## Incidence, adults 15-49, by gender
  inc_adults_gender <- outputs$inc_adults_gender %>%
                          group_by(gender, year) %>%
                          summarise(inc_mean = mean(inc),
                                    inc_median = median(inc),
                                    inc_lower = quantile(inc, probs = 0.025),
                                    inc_upper = quantile(inc, probs = 0.975)) %>%
                          mutate(scenario = sub(".*__", "", scenario_name))
  
  ## Incidence, by age and gender
  inc_age_gender <- outputs$inc_age_gender %>%
    group_by(age_group, gender, year) %>%
    summarise(inc_mean = mean(inc),
              inc_median = median(inc),
              inc_lower = quantile(inc, probs = 0.025),
              inc_upper = quantile(inc, probs = 0.975)) %>%
    mutate(scenario = sub(".*__", "", scenario_name))
  
  ## Cumulative incidence, adults 15-49, by gender
  t_horiz_start <- 2020
  t_horiz_end <- seq(2020,2050,1)
  
  for(i in 1:length(t_horiz_end)){
    df <- outputs$inc_adults_gender %>%
      group_by(gender, run_num) %>%
      summarise(new_infections_sum=sum(new_infections[year>=t_horiz_start & year<=t_horiz_end[i]]), 
                py_sum=sum(person_years[year>=t_horiz_start & year<=t_horiz_end[i]]), cuminc=100*new_infections_sum/py_sum) %>%
      mutate(horizon_start = t_horiz_start, horizon_end = t_horiz_end[i],scenario = sub(".*__", "", scenario_name))
      if(i==1){cum_inc_adults_gender<-df} 
      if(i>1){cum_inc_adults_gender<-rbind(df, cum_inc_adults_gender)}
  }
  
  ## Cumulative incidence, adults by age and gender
  for(i in 1:length(t_horiz_end)){
    df <- outputs$inc_age_gender %>%
      group_by(age_group, gender, run_num) %>%
      summarise(new_infections_sum=sum(new_infections[year>=t_horiz_start & year<=t_horiz_end[i]]), 
                py_sum=sum(person_years[year>=t_horiz_start & year<=t_horiz_end[i]]), cuminc=100*new_infections_sum/py_sum) %>%
      mutate(horizon_start = t_horiz_start, horizon_end = t_horiz_end[i],scenario = sub(".*__", "", scenario_name))
    if(i==1){cum_inc_age_gender<-df} 
    if(i>1){cum_inc_age_gender<-rbind(df, cum_inc_age_gender)}
  }

  return(list("prev_adults_gender" = prev_adults_gender,
              "inc_adults_gender" = inc_adults_gender,
              "cum_inc_adults_gender" = cum_inc_adults_gender,
              "cum_inc_age_gender"=cum_inc_age_gender))
}

## Save
save(outputs, file = file.path(output_dir, suite_name, "scenario_comparison.RData"))

## Plots
pdf(file = file.path(output_dir, suite_name, "scenario_comparison_plots.pdf"), 
    height = 5, width = 7)

## Prevalence, adults 15-49, by gender across scenarios
print(outputs$prev_adults_gender %>%
        ggplot(aes(x = year, y = prev_mean)) +
        geom_line(aes(color = factor(scenario))) +
        #geom_ribbon(aes(ymin = prev_lower, ymax = prev_upper, fill = factor(scenario)), alpha = 0.3, show.legend = FALSE) +
        #scale_color_manual(guide = guide_legend(title = ""),values=c("blue","green","red"), labels = c("Status-quo", "StopVMMC2020","VMMC151980pct")) +
        #scale_fill_manual(guide = guide_legend(title = ""),values=c("blue","green","red"), labels = c("Status-quo", "StopVMMC2020","VMMC151980pct")) +
        #scale_color_manual(values = c("red", "blue")) +
        labs(x = "Year", y = "Prevalence") +
        facet_wrap(~gender) +
        ggtitle("Prevalence among adults age 15-49") +
        theme(legend.position="bottom") 
      )

## Incidence, adults 15-49, by gender across scenarios
print(subset(outputs$inc_adults_gender, year>2019) %>%
        ggplot(aes(x = year, y = inc_mean)) +
        geom_line(aes(color = factor(scenario))) +
        #geom_ribbon(aes(ymin = inc_lower, ymax = inc_upper, fill = factor(scenario)), alpha = 0.3, show.legend = FALSE) +
        #scale_color_manual(guide = guide_legend(title = ""),values=c("blue","green","red"), labels = c("Status-quo", "StopVMMC2020","VMMC151980pct")) +
        scale_y_continuous(expand=c(0,0), limits=c(0,2), breaks=seq(0,2,0.1)) + 
        labs(x = "Year", y = "Incidence (per 100 PY)") +
        facet_wrap(~gender) +
        #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
        geom_hline(yintercept = 0.1, linetype="dashed") +
        ggtitle("Incidence among adults age 15-49") +
        theme(legend.position="bottom") 
)

## Cumulative incidence and vaccine efficacy adults 15-49, by gender across scenarios
effectiveness.plot <- outputs$cum_inc_adults_gender %>% 
  pivot_wider(names_from = scenario, values_from = c(new_infections_sum, py_sum, cuminc)) %>%
  select(c(gender, run_num, horizon_start, horizon_end, contains("cuminc"))) %>%
  mutate_at(vars(contains("s801519")), funs(rel = 1-(./cuminc_noVMMC), abs = cuminc_noVMMC-.)) %>%
  select(-cuminc_noVMMC) %>%
  select(gender, run_num, horizon_start, horizon_end, contains("rel")|contains("abs")) %>%
  group_by(gender, horizon_start, horizon_end) %>%
  summarise_at(vars(contains("s801519")), funs(med=median(.), lb=quantile(., probs = 0.025), ub=quantile(., probs = 0.975))) %>%
  pivot_longer(cols = contains("s801519"), names_to = c("outcome","campaign", "program_years","metric","est"), names_sep="_", values_to = "value", values_drop_na = TRUE) %>%
  pivot_wider(names_from = c(est), values_from = c(value)) %>%
  mutate(program_years = fct_relevel(program_years, "statusquo", "5yr", "10yr", "15yr", "20yr", "30yr"))

dodge <- position_dodge(3)

print(subset(effectiveness.plot, metric=="rel" & horizon_end %in% seq(2020,2050,5)) %>%
        ggplot(aes(x = horizon_end, y = med)) +
        geom_point(aes(color=program_years), size=1, position = dodge) +
        geom_errorbar(aes(ymin = lb, ymax = ub, color=program_years), position = dodge) + 
        scale_color_viridis_d(guide = guide_legend(title = "")) +
        #geom_line(aes(color=efficacy_pct), size=1) +
        #geom_ribbon(aes(ymin = lb, ymax = ub, fill=efficacy), alpha = 0.3, show.legend = FALSE) +
        geom_hline(yintercept = 0, linetype="dashed") +
        labs(x = "Year", y = "Reduction in cumulative incidence") +
        facet_grid(metric~gender) +
        ggtitle("Reduction in adult (15-49) incidence (%) over increasing time-horizons") +
        theme(legend.position="bottom") 
)

print(subset(effectiveness.plot, metric=="abs" & horizon_end %in% seq(2020,2050,5)) %>%
        ggplot(aes(x = horizon_end, y = med)) +
        geom_point(aes(color=program_years), size=1, position = dodge) +
        geom_errorbar(aes(ymin = lb, ymax = ub, color=program_years), position = dodge) + 
        scale_color_viridis_d(guide = guide_legend(title = "")) +
        #geom_line(aes(color=efficacy_pct), size=1) +
        #geom_ribbon(aes(ymin = lb, ymax = ub, fill=efficacy), alpha = 0.3, show.legend = FALSE) +
        geom_hline(yintercept = 0, linetype="dashed") +
        labs(x = "Year", y = "Reduction in cumulative incidence") +
        facet_grid(metric~gender) +
        ggtitle("Reduction in adult (15-49) incidence (%) over increasing time-horizons") +
        theme(legend.position="bottom") 
)

# effectiveness.plot.age <- outputs$cum_inc_age_gender %>%  select(-scenario_f) %>%
#   pivot_wider(names_from = scenario, values_from = c(new_infections_sum, py_sum, cuminc)) %>%
#   select(c(age_group, gender, run_num, horizon_start, horizon_end, contains("cuminc"))) %>%
#   mutate_at(vars(contains("Vaccine")), funs(rel = 1-(./cuminc_baseline), abs = cuminc_baseline-.)) %>%
#   select(-cuminc_baseline) %>%
#   select(age_group, gender, run_num, horizon_start, horizon_end, contains("rel")|contains("abs")) %>%
#   group_by(age_group, gender, horizon_start, horizon_end) %>%
#   summarise_at(vars(contains("Vaccine")), funs(med=median(.), lb=quantile(., probs = 0.025), ub=quantile(., probs = 0.975))) %>%
#   pivot_longer(cols = contains("Vaccine"), names_to = c("incidence_name", "efficacy","metric","est"), names_sep="_", values_to = "value", values_drop_na = TRUE) %>%
#   pivot_wider(names_from = c(est), values_from = c(value)) %>%
#   mutate(efficacy_pct = gsub("[^0-9.-]", "",efficacy))
# 
# print(subset(vaccine_efficacy_age, metric=="abs" & horizon_end>2028) %>%
#         ggplot(aes(x = horizon_end, y = med)) +
#         geom_line(aes(color=efficacy_pct), size=1) +
#         #geom_ribbon(aes(ymin = lb, ymax = ub, fill=efficacy), alpha = 0.3, show.legend = FALSE) +
#         scale_color_viridis_d(guide = guide_legend(title = ""),labels = c("50%", "55%","60%","65%","70%")) +
#         geom_hline(yintercept = 0, linetype="dashed") +
#         labs(x = "Year", y = "Reduction in cumulative incidence") +
#         facet_grid(gender~age_group) +
#         #ggtitle("Population-level vaccine effectiveness in adults by age \n over varying time-horizons from 2028-2050") +
#         theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#         theme(legend.position="bottom", text=element_text(size=21)) 
# )
# ggsave(".//SummaryOutput//VaccineRCTTargetBy1yrAge//abscuminc_age.pdf", height=10, width=15)

# #########################################################################################################################
# ### Compare vaccine efficacy (within trial) across scenarios
# #########################################################################################################################
# ## For each scenario, compute summary statistics and return outputs
# outputs <- foreach(i=2:nrow(scenarios), .combine = 'comb', .multicombine = TRUE) %dopar% { #only bring together vaccine scenarios
#   
#   library(ggplot2)
#   library(tidyverse)
#   library(viridisLite)
#   library(parallel)
#   library(foreach)
#   library(doParallel)
#   library(lemon)
#   
#   scenario_name <- scenarios$name[i]
#   
#   ## Load outputs
#   load(file.path(output_dir, suite_name, scenario_name, "scenario_summary.RData"))
#   
#   if("vaccine_trial_efficacy" %in% names(outputs)){    
#     #Vaccine trial efficacy
#     vaccine_trial_efficacy <- outputs$vaccine_trial_efficacy %>%
#       group_by(gender) %>%
#       summarise(inc_mean_vaccine0 = mean(inc_arm_vaccine0),
#                 inc_lb_vaccine0 = quantile(inc_arm_vaccine0, probs = 0.025,na.rm = TRUE),
#                 inc_ub_vaccine0 = quantile(inc_arm_vaccine0, probs = 0.975,na.rm = TRUE),
#                 inc_mean_vaccine1 = mean(inc_arm_vaccine1),
#                 inc_lb_vaccine1 = quantile(inc_arm_vaccine1, probs = 0.025,na.rm = TRUE),
#                 inc_ub_vaccine1 = quantile(inc_arm_vaccine1, probs = 0.975,na.rm = TRUE),
#                 ve_mean = mean(ve),
#                 ve_lb = quantile(ve, probs = 0.025,na.rm = TRUE),
#                 ve_ub = quantile(ve, probs = 0.975,na.rm = TRUE)) %>%
#       mutate(scenario = scenario_name)
#   }
#   
#   return(list("vaccine_trial_efficacy"=vaccine_trial_efficacy))
# }
# 
# ## Save
# save(outputs, file = file.path(output_dir, suite_name, "scenario_comparison_VE.RData"))
# 
# load(file = file.path(output_dir, suite_name, "scenario_comparison_VE.RData"))
# 
# outputs$vaccine_trial_efficacy$scenario_f <- "baseline"
# outputs$vaccine_trial_efficacy$scenario_f[grep("VaccineEfficacy", outputs$vaccine_trial_efficacy$scenario)] <- str_sub(outputs$vaccine_trial_efficacy$scenario[grep("VaccineEfficacy", outputs$vaccine_trial_efficacy$scenario)], start=-17)
# outputs$vaccine_trial_efficacy$scenario <- outputs$vaccine_trial_efficacy$scenario_f
# outputs$vaccine_trial_efficacy$pca_efficacy <- as.numeric(str_sub(outputs$vaccine_trial_efficacy$scenario, -2,-1))/100
# 
# labels=c("VaccineEfficacy50"="50%","VaccineEfficacy55"= "55%","VaccineEfficacy60"="60%","VaccineEfficacy65"="65%","VaccineEfficacy70"="70%")
# print(ggplot(data=outputs$vaccine_trial_efficacy, aes(x = gender)) +
#         geom_errorbar(aes(ymin=ve_lb, ymax=ve_ub, color=gender)) +
#         geom_point(aes(y = ve_mean, color = gender))+
#         geom_hline(aes(yintercept=pca_efficacy), linetype="dashed") +
#         labs(x = "age group", y = "observed vaccine efficacy") +
#         facet_grid(~scenario, labeller=labeller(scenario = labels)) +
#         ggtitle("") +
#         theme_bw(base_size=16) +
#         theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#         scale_y_continuous(breaks=seq(-0.2,1,0.1), limits=c(0,1), expand=c(0,0)) +
#         theme(legend.title=element_blank()) +
#         guides(colour = guide_legend(override.aes = list(alpha = 1)))
# )
# 
dev.off()
