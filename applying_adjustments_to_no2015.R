

rm(list = ls(all = TRUE))


library("ggplot2")
library("plyr")
library("dplyr")
library("zoo")
library("gridExtra")

library("lubridate")
library("maditr")



cases <- read.csv("~/Box/NU-malaria-team/projects/smc_impact/data/outputs/U5_DS_cases_seasonal_smc_good_rows_MA_imputes_testing_w_rep_weights.csv", header  = TRUE, strip.white = TRUE, stringsAsFactors = FALSE)
cases$U5_pop <- cases$District.Pop * .18
cases$Date <- as.yearmon(cases$Date)

cases$rep_rate <- cases$reporting_HF / cases$total_HF


cases_no2015 <- cases[which(cases$year != 2015),]


untested_ind <- which(cases_no2015$susp_u5 >= cases_no2015$test_rdt_mic_u5)
unconfirmed_ind <- which((cases_no2015$maltreat_u5 >= cases_no2015$conf_rdt_mic_u5) &
                             (cases_no2015$susp_u5 < cases_no2015$test_rdt_mic_u5))

cases_no2015$presumed_u5 <- NA
cases_no2015[untested_ind, "presumed_u5"] <- cases_no2015[untested_ind, "susp_u5"] - cases_no2015[untested_ind, "test_rdt_mic_u5"]
cases_no2015[unconfirmed_ind, "presumed_u5"] <- cases_no2015[unconfirmed_ind, "maltreat_u5"] - cases_no2015[unconfirmed_ind, "conf_rdt_mic_u5"]


sum(is.na(cases_no2015$presumed_u5))




cases_no2015$N1_adj <- cases_no2015$conf_rdt_mic_u5 + (cases_no2015$presumed_u5 * (cases_no2015$conf_rdt_mic_u5 / cases_no2015$test_rdt_mic_u5))



# cases_no2015$N2_adj <- cases_no2015$N1_adj + (cases_no2015$N1_adj * (1 - cases_no2015$rep_rate))
cases_no2015$N2_adj <- cases_no2015$N1_adj / cases_no2015$rep_rate




ggplot(cases_no2015, aes(x = Date, y = N2_adj / (U5_pop/1000),
                         group = as.factor(District))) +
    geom_line(alpha = 0.25, size = 1, show.legend = FALSE, color = "blue") + ylab("Cases per 1000") +
    ggtitle("Adjusted malaria cases among u5") +
    scale_x_yearmon("Date", breaks = sort(unique(cases$Date))[c(seq(1,48,6), 48)],
                    labels = sort(unique(cases$Date))[c(seq(1,48,6), 48)]) +
    theme_bw() + theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45, hjust = 1),
                       panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))




cases_no2015$N2_weighted_adj <- cases_no2015$N1_adj / cases_no2015$weighted_rep_rate


ggplot(cases_no2015, aes(x = Date, y = N2_weighted_adj / (U5_pop/1000),
                         group = as.factor(District))) +
    geom_line(alpha = 0.25, size = 1, show.legend = FALSE, color = "blue") + ylab("Cases per 1000") +
    ggtitle("Adjusted malaria cases among u5 (weighted rep rate)") +
    scale_x_yearmon("Date", breaks = sort(unique(cases$Date))[c(seq(1,48,6), 48)],
                    labels = sort(unique(cases$Date))[c(seq(1,48,6), 48)]) +
    theme_bw() + theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45, hjust = 1),
                       panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))






#############################################################################################################







# cases$cases_rep_adj <- cases$conf_rdt_mic_u5 + (cases$conf_rdt_mic_u5 * (1 - cases$rep_rate))
cases$cases_rep_adj <- cases$conf_rdt_mic_u5 / cases$rep_rate



ggplot(cases, aes(x = Date, y = cases_rep_adj / (U5_pop/1000),
                         group = as.factor(District))) +
    geom_line(alpha = 0.25, size = 1, show.legend = FALSE, color = "blue") + ylab("Cases per 1000") +
    ggtitle("Reporting adjusted malaria cases among u5") +
    scale_x_yearmon("Date", breaks = sort(unique(cases$Date))[c(seq(1,48,6), 48)],
                    labels = sort(unique(cases$Date))[c(seq(1,48,6), 48)]) +
    theme_bw() + theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45, hjust = 1),
                       panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))



ggplot(cases, aes(x = Date, y = rep_rate,
                  group = as.factor(District))) +
    geom_line(alpha = 0.25, size = 1, show.legend = FALSE, color = "blue") + ylab("") +
    ggtitle("Health facility reporting rate by District") +
    scale_x_yearmon("Date", breaks = sort(unique(cases$Date))[c(seq(1,48,6), 48)],
                    labels = sort(unique(cases$Date))[c(seq(1,48,6), 48)]) +
    theme_bw() + theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45, hjust = 1),
                       panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))






cases$cases_rep_weighted_adj <- cases$conf_rdt_mic_u5 / cases$weighted_rep_rate



ggplot(cases, aes(x = Date, y = cases_rep_weighted_adj / (U5_pop/1000),
                  group = as.factor(District))) +
    geom_line(alpha = 0.25, size = 1, show.legend = FALSE, color = "blue") + ylab("Cases per 1000") +
    ggtitle("Reporting adjusted malaria cases among u5 (w/ rep weights)") +
    scale_x_yearmon("Date", breaks = sort(unique(cases$Date))[c(seq(1,48,6), 48)],
                    labels = sort(unique(cases$Date))[c(seq(1,48,6), 48)]) +
    theme_bw() + theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45, hjust = 1),
                       panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))



ggplot(cases, aes(x = Date, y = weighted_rep_rate,
                  group = as.factor(District))) +
    geom_line(alpha = 0.25, size = 1, show.legend = FALSE, color = "blue") + ylab("") +
    ggtitle("Health facility weighted reporting rate by District") +
    scale_x_yearmon("Date", breaks = sort(unique(cases$Date))[c(seq(1,48,6), 48)],
                    labels = sort(unique(cases$Date))[c(seq(1,48,6), 48)]) +
    theme_bw() + theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45, hjust = 1),
                       panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))















ggplot(cases, aes(x = Date, y = conf_rdt_mic_u5 / (U5_pop/1000),
                  group = as.factor(District))) +
    geom_line(alpha = 0.25, size = 1, show.legend = FALSE, color = "blue") + ylab("Cases per 1000") +
    ggtitle("Crude malaria incidence among u5") +
    scale_x_yearmon("Date", breaks = sort(unique(cases$Date))[c(seq(1,48,6), 48)],
                    labels = sort(unique(cases$Date))[c(seq(1,48,6), 48)]) +
    theme_bw() + theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45, hjust = 1),
                       panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))















#############################################################################################################





# med_2014 <- read.csv("~/Box/NU-malaria-team/projects/hbhi_burkina/DS DHS estimates/medfever/med_BF14DS.csv", stringsAsFactors = FALSE)
# med_2017 <- read.csv("~/Box/NU-malaria-team/projects/hbhi_burkina/DS DHS estimates/medfever/med_BF17DS.csv", stringsAsFactors = FALSE)


med_2014 <- read.csv("~/OneDrive/Desktop/medfever_region_2014.csv", stringsAsFactors = FALSE)
med_2017 <- read.csv("~/OneDrive/Desktop/medfever_region_2017.csv", stringsAsFactors = FALSE)



med_2014 <- med_2014[,-4]
names(med_2014)[3] <- "medfever_2014"
med_2017 <- med_2017[,-4]
names(med_2017)[3] <- "medfever_2017"

medfever_DHS <- inner_join(med_2014, med_2017, by = c("NOMDEP", "NOMREGION"))
medfever_DHS <- cbind(medfever_DHS[,1:3], medfever_DHS[,3], medfever_DHS[,4], medfever_DHS[,4])
names(medfever_DHS) <- c("District", "Region", "2015", "2016", "2017", "2018")

medfever_DHS <- melt(medfever_DHS, id.vars = c("District", "Region"),
                     variable.name = "year", value.name = "medfever_regional")



district_map <- cbind(sort(unique(cases$District)),
                      sort(unique(medfever_DHS$District))[c(1:40,42,41,43:70)])

for (i in 1:nrow(district_map))
{
    DS <- district_map[i, 2]
    new_name_DS <- district_map[i, 1]
    
    medfever_DHS[which(medfever_DHS$District == DS), "District"] <- new_name_DS
}

medfever_DHS$year <- as.numeric(as.character(medfever_DHS$year))

# medfever_DHS <- left_join(medfever_DHS[,-2], unique(cases[,c("Region", "District", "year", "U5_pop")]), by = c("District", "year"))


# medfever_region <- ddply(medfever_DHS[!is.na(medfever_DHS$medfever),], c(.(Region), .(year)),
#                          summarise, medfever = sum(medfever * U5_pop)/ sum(U5_pop))
# 
# 
# medfever_DHS_tmp_1 <- medfever_DHS[!is.na(medfever_DHS$medfever),]
# medfever_DHS_tmp_2 <- medfever_DHS[is.na(medfever_DHS$medfever),]
# 
# medfever_DHS_tmp_2 <- left_join(medfever_DHS_tmp_2[,-3], medfever_region, by = c("year", "Region"))
# medfever_DHS_tmp_2 <- medfever_DHS_tmp_2[, c(1:2, 5, 3:4)]
# 
# medfever_DHS <- rbind(medfever_DHS_tmp_1, medfever_DHS_tmp_2)

# names(medfever_region)[3] <- "medfever_regional"
# medfever_DHS <- left_join(medfever_DHS, medfever_region, by = c("year", "Region"))


cases <- left_join(cases, medfever_DHS[,-2], by = c("District", "year"))

for (D in unique(cases$District))
{
    cases[which(cases$District == D & cases$Date %in% as.yearmon(seq(as.Date("2016-05-01"), as.Date("2016-12-01"), by="months"))), "medfever_regional"] <- cases[which(cases$District == D & cases$Date == "Jan 2017"), "medfever_regional"]
}



# cases$cases_trtseeking_adj <- cases$conf_rdt_mic_u5 + (cases$conf_rdt_mic_u5 * (1 - cases$medfever))
cases$cases_trtseeking_adj <- cases$conf_rdt_mic_u5 + (cases$conf_rdt_mic_u5 * (1 - cases$medfever_regional))



ggplot(cases, aes(x = Date, y = cases_trtseeking_adj / (U5_pop/1000),
                  group = as.factor(District))) +
    geom_line(alpha = 0.25, size = 1, show.legend = FALSE, color = "blue") + ylab("Cases per 1000") +
    ggtitle("Treatment seeking adjusted malaria cases among u5") +
    scale_x_yearmon("Date", breaks = sort(unique(cases$Date))[c(seq(1,48,6), 48)],
                    labels = sort(unique(cases$Date))[c(seq(1,48,6), 48)]) +
    theme_bw() + theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45, hjust = 1),
                       panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))









#############################################################################################################






cases$cases_reporting_trtseeking_adj <- cases$cases_rep_adj + (cases$cases_rep_adj * (1 - cases$medfever_regional))




ggplot(cases, aes(x = Date, y = cases_reporting_trtseeking_adj / (U5_pop/1000),
                  group = as.factor(District))) +
    geom_line(alpha = 0.25, size = 1, show.legend = FALSE, color = "blue") + ylab("Cases per 1000") +
    ggtitle("Malaria cases among u5 adjusted for reporting and treatment seeking") +
    scale_x_yearmon("Date", breaks = sort(unique(cases$Date))[c(seq(1,48,6), 48)],
                    labels = sort(unique(cases$Date))[c(seq(1,48,6), 48)]) +
    theme_bw() + theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45, hjust = 1),
                       panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))







cases_no2015 <- left_join(cases_no2015, medfever_DHS[,c(1:3, 6)], by = c("District", "year"))
cases_no2015$N3_adj <- cases_no2015$N2_adj + (cases_no2015$N2_adj * (1 - cases_no2015$medfever_regional))


ggplot(cases_no2015, aes(x = Date, y = N3_adj / (U5_pop/1000),
                         group = as.factor(District))) +
    geom_line(alpha = 0.25, size = 1, show.legend = FALSE, color = "blue") + ylab("Cases per 1000") +
    ggtitle("Adjusted malaria cases among u5") +
    scale_x_yearmon("Date", breaks = sort(unique(cases$Date))[c(seq(1,48,6), 48)],
                    labels = sort(unique(cases$Date))[c(seq(1,48,6), 48)]) +
    theme_bw() + theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45, hjust = 1),
                       panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))














cases$cases_weighted_reporting_trtseeking_adj <- cases$cases_rep_weighted_adj + (cases$cases_rep_weighted_adj * (1 - cases$medfever_regional))




ggplot(cases, aes(x = Date, y = cases_weighted_reporting_trtseeking_adj / (U5_pop/1000),
                  group = as.factor(District))) +
    geom_line(alpha = 0.25, size = 1, show.legend = FALSE, color = "blue") + ylab("Cases per 1000") +
    ggtitle("Malaria cases among u5 adjusted for weighted reporting and treatment seeking") +
    scale_x_yearmon("Date", breaks = sort(unique(cases$Date))[c(seq(1,48,6), 48)],
                    labels = sort(unique(cases$Date))[c(seq(1,48,6), 48)]) +
    theme_bw() + theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45, hjust = 1),
                       panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))







cases_no2015 <- left_join(cases_no2015, medfever_DHS[,c(1:3, 6)], by = c("District", "year"))
cases_no2015$N3_adj <- cases_no2015$N2_adj + (cases_no2015$N2_adj * (1 - cases_no2015$medfever_regional))


ggplot(cases_no2015, aes(x = Date, y = N3_adj / (U5_pop/1000),
                         group = as.factor(District))) +
    geom_line(alpha = 0.25, size = 1, show.legend = FALSE, color = "blue") + ylab("Cases per 1000") +
    ggtitle("Adjusted malaria cases among u5 (w/ rep weights)") +
    scale_x_yearmon("Date", breaks = sort(unique(cases$Date))[c(seq(1,48,6), 48)],
                    labels = sort(unique(cases$Date))[c(seq(1,48,6), 48)]) +
    theme_bw() + theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45, hjust = 1),
                       panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))






