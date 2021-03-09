###############################################
##  Plotting malaria variables for figure 3  ##
###############################################
#
# Description:
#   Plotting crude incidence and malaria proportion of outpatient visits
#       from figure 3a and 3b.
#   Time series for each of the 70 districts for each malaria variable
#
#   Also making Figure 1 on health facility reports received (complete or not)
#
#
#  Sebastian Rodriguez (sebastian@rodriguez.cr)
#  Last edited Mar 09, 2021
#





rm(list = ls(all = TRUE))


library("ggplot2")
library("plyr")
library("dplyr")
library("zoo")
library("gridExtra")

library("lubridate")


cases <- read.csv("~/Box/NU-malaria-team/projects/smc_impact/data/outputs/U5_DS_cases_seasonal_smc_good_rows_kalman_imputes.csv", header  = TRUE, strip.white = TRUE, stringsAsFactors = FALSE)
cases$U5_pop <- cases$District.Pop * .18
cases$Date <- as.yearmon(cases$Date)




## Figure 3a
ggplot(cases, aes(x = Date, y =  conf_rdt_mic_u5 / (U5_pop/1000),
                  group = as.factor(District))) +
    geom_line(alpha = 0.25, size = 1, show.legend = FALSE, color = "blue") + ylab("Cases per 1000") +
    ggtitle("Confirmed cases of malaria from RDT or Microscopy among children under 5") +
    scale_x_yearmon("Date", breaks = sort(unique(cases$Date))[c(seq(1,48,6), 48)],
                     labels = sort(unique(cases$Date))[c(seq(1,48,6), 48)]) +
    theme_bw() + theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45, hjust = 1), panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))



## Figure 3b
ggplot(cases, aes(x = Date, y = conf_rdt_mic_u5 / allout_u5,
                  group = as.factor(District))) +
    geom_line(alpha = 0.25, size = 1, show.legend = FALSE, color = "blue") + ylab("") +
    ggtitle("Proportion of all-cause outpatient visits due to malaria in under-5 population") +
    scale_x_yearmon("Date", breaks = sort(unique(cases$Date))[c(seq(1,48,6), 48)],
                     labels = sort(unique(cases$Date))[c(seq(1,48,6), 48)]) +
    theme_bw() + theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45, hjust = 1), panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))



