
rm(list = ls(all = TRUE))


library("ggplot2")
library("plyr")
library("dplyr")
library("zoo")
library("gridExtra")

library("lubridate")
library("maditr")





HF_cases <- read.csv("~/Box/NU-malaria-team/projects/smc_impact/data/outputs/U5_HF_cases_smc_coords_imputed_rdts_and_allout_testing_MA.csv", stringsAsFactors = FALSE)
HF_cases$Date <- as.yearmon(HF_cases$Date)



tmp_total_HF <- ddply(HF_cases, c(.(District), .(Date), .(month)), summarise,
                      total_HF = length(unique(UID)))



average_HF_counts_per_month <- ddply(HF_cases, c(.(District), .(UID), .(month)), summarise,
                                     average_counts_HF = mean(conf_rdt_mic_u5, na.rm = T),
                                     total_counts_HF = sum(conf_rdt_mic_u5, na.rm = T),
                                     na_counts_HF = sum(is.na(conf_rdt_u5)))

average_District_counts_per_month <- ddply(HF_cases, c(.(District), .(month)), summarise,
                                           average_counts_D = mean(conf_rdt_mic_u5, na.rm = T),
                                           total_counts_D = sum(conf_rdt_mic_u5, na.rm = T),
                                           na_counts_D = sum(is.na(conf_rdt_u5)))


average_counts_per_month <- left_join(average_District_counts_per_month, average_HF_counts_per_month,
                                      by = c("District", "month"))
average_counts_per_month <- left_join(average_counts_per_month, tmp_total_HF,
                                      by = c("District", "month"))

average_counts_per_month <- average_counts_per_month[-which(average_counts_per_month$na_counts_HF == 4),]
average_counts_per_month$weights <- average_counts_per_month$total_counts_HF/average_counts_per_month$total_counts_D




tmp_total_HF_weighted <- ddply(average_counts_per_month, c(.(District), .(month), .(total_HF)), summarise,
                               total_HF_weighted = sum(weights))




