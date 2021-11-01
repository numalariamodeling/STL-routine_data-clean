

rm(list = ls(all = TRUE))


library("ggplot2")
library("plyr")
library("dplyr")
library("zoo")
library("gridExtra")

library("lubridate")
library("maditr")





################################################################################
## FIRST U5


cases <- read.csv("~/Box/NU-malaria-team/projects/smc_impact/data/outputs/U5_DS_cases_seasonal_smc_good_rows_MA_imputes_testing_w_rep_weights.csv", header  = TRUE, strip.white = TRUE, stringsAsFactors = FALSE)
cases$U5_pop <- cases$District.Pop * .18
cases$Date <- as.yearmon(cases$Date)




cases_national <- ddply(cases, .(Date), summarize,
                        conf_rdt_mic_u5 = sum(conf_rdt_mic_u5),
                        U5_pop = sum(U5_pop))




ggplot(cases_national, aes(x = Date, y = conf_rdt_mic_u5 / (U5_pop/1000))) +
    geom_line(alpha = 0.25, size = 1, show.legend = FALSE, color = "blue") + ylab("National cases per 1000 among children under-5") +
    ggtitle("") +
    scale_x_yearmon("Date", breaks = sort(unique(cases_national$Date))[c(seq(1,48,6), 48)],
                    labels = sort(unique(cases_national$Date))[c(seq(1,48,6), 48)]) +
    theme_bw() + theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45, hjust = 1),
                       panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))



################################################################################
## NOW OV5


cases <- read.csv("~/Box/NU-malaria-team/projects/smc_impact/data/outputs/Ov5_DS_cases_seasonal_smc_good_rows_MA_imputes_w_rep_weights.csv", header  = TRUE, strip.white = TRUE, stringsAsFactors = FALSE)
cases$Ov5_pop <- cases$District.Pop * .82
cases$Date <- as.yearmon(cases$Date)




cases_national <- ddply(cases, .(Date), summarize,
                        conf_rdt_mic_ov5 = sum(conf_rdt_mic_ov5),
                        Ov5_pop = sum(Ov5_pop))




ggplot(cases_national, aes(x = Date, y = conf_rdt_mic_ov5 / (Ov5_pop/1000))) +
    geom_line(alpha = 0.25, size = 1, show.legend = FALSE, color = "blue") + ylab("National cases per 1000 among individuals over-5") +
    ggtitle("") +
    scale_x_yearmon("Date", breaks = sort(unique(cases_national$Date))[c(seq(1,48,6), 48)],
                    labels = sort(unique(cases_national$Date))[c(seq(1,48,6), 48)]) +
    theme_bw() + theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45, hjust = 1),
                       panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))



