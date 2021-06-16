##############################################
##  Seasonal Trend Decomposition SI figures ##
##############################################
#
# Description:
#
#
#
#
#
#  Sebastian Rodriguez (sebastian@rodriguez.cr)
#  Last edited Dec 13, 2020
#
#



rm(list = ls(all = TRUE))


require("ggplot2")
require("plyr")
require("dplyr")
require("zoo")
require("gridExtra")
require("lubridate")


require("stlplus")
require("reshape2")
require("trend")


require("RColorBrewer")
require("spdep")
require("rgdal")
require("sf")
require("rgeos")



# Function for standardizing malaria variables

getNormalized <- function(vec)
{
    norm_vec <- (vec - mean(vec))/sd(vec)
    return(norm_vec)
}






# Loading health district dataset
# Creating under-5 population column and fixing date column

# cases <- read.csv("~/Box/NU-malaria-team/projects/smc_impact/data/outputs/U5_DS_cases_seasonal_smc_good_rows_MA_imputes.csv", header  = TRUE, strip.white = TRUE, stringsAsFactors = FALSE)
cases <- read.csv("~/Box/NU-malaria-team/projects/smc_impact/data/outputs/U5_DS_cases_seasonal_smc_good_rows_MA_imputes_testing_w_rep_weights.csv", header  = TRUE, strip.white = TRUE, stringsAsFactors = FALSE)
cases$U5_pop <- cases$District.Pop * .18
cases$Date <- as.yearmon(cases$Date)





# Finding what date SMC was first implemented in each district
SMC.init_tmp <- cases[which(cases$SMC.received == 1),
                      c("District", "Date")]
SMC.init_tmp_ordered <- SMC.init_tmp[order(SMC.init_tmp$District, SMC.init_tmp$Date),]
SMC.init <- SMC.init_tmp_ordered[!duplicated(SMC.init_tmp_ordered$District),]

# Fixing so we get 2014 start dates where appropriate
DS_SMC_2014 <- c("bogande", "bousse", "boussouma", "garango", "kaya",
                 "sebba", "seguenega", "tougan")
SMC.init[which(SMC.init$District %in% DS_SMC_2014), "Date"] <- as.yearmon("Aug 2014")

# Same for 2019
SMC.init_2019 <- data.frame(District = c("baskuy", "bogodogo", "boulmiougou",
                                         "nongr-massom", "sig-noghin"),
                            Date = as.yearmon(rep("Jul 2019", 5)))
SMC.init <- rbind(SMC.init, SMC.init_2019)
names(SMC.init)[2] <- "SMC.init_date"

# Joining SMC initialization date information
cases <- inner_join(cases, SMC.init, by = "District")




################################################################################




cases <- cases[order(cases$District, cases$Date),]


cases$mal_cases <- cases$conf_rdt_mic_u5 / (cases$U5_pop/1000)
cases$all_cases <- cases$allout_u5 / (cases$U5_pop/1000)
cases$all_nonMal_cases <- (cases$allout_u5 - cases$conf_rdt_mic_u5) / (cases$U5_pop/1000)
cases$mal_ratio <- cases$conf_rdt_mic_u5 / cases$allout_u5





################################################################################


med_2014 <- read.csv("~/OneDrive/Desktop/medfever_region_2014.csv", stringsAsFactors = FALSE)
med_2017 <- read.csv("~/OneDrive/Desktop/medfever_region_2017.csv", stringsAsFactors = FALSE)



med_2014 <- med_2014[,-4]
names(med_2014)[3] <- "medfever_2014"
med_2017 <- med_2017[,-4]
names(med_2017)[3] <- "medfever_2017"

medfever_DHS <- inner_join(med_2014, med_2017, by = c("NOMDEP", "NOMREGION"))
medfever_DHS <- cbind(medfever_DHS[,1:3], medfever_DHS[,4])
names(medfever_DHS) <- c("District", "Region", "Aug 2014", "Dec 2017")

medfever_DHS <- melt(medfever_DHS, id.vars = c("District", "Region"),
                     variable.name = "Date", value.name = "medfever_regional")



district_map <- cbind(sort(unique(cases$District)),
                      sort(unique(medfever_DHS$District))[c(1:40,42,41,43:70)])

for (i in 1:nrow(district_map))
{
    DS <- district_map[i, 2]
    new_name_DS <- district_map[i, 1]
    
    medfever_DHS[which(medfever_DHS$District == DS), "District"] <- new_name_DS
}

medfever_DHS$Date <- as.yearmon(medfever_DHS$Date)
medfever_DHS <- medfever_DHS[order(medfever_DHS$District, medfever_DHS$Date),]



medfever_DHS <- as.data.frame(medfever_DHS)

dates <- seq(as.Date("2014-08-01"), as.Date("2018-12-01"), by="months")
medfever_DHS_fitted <- data.frame()

for (D in unique(cases$District))
{
    D_data <- medfever_DHS[medfever_DHS$District == D,]
    D_data$ind <- c(1, 41)
    
    D_lm <- lm(medfever_regional ~ ind, data = D_data)
    D_fitted <- predict(D_lm, data.frame("ind" = c(1:53)))
    
    medfever_DHS_fitted_D <- data.frame("District" = rep(D, 53),
                                        "Date" = dates,
                                        "fitted_regional_medfever" = D_fitted)
    medfever_DHS_fitted <- rbind(medfever_DHS_fitted, medfever_DHS_fitted_D)
}
medfever_DHS_fitted$Date <- as.yearmon(medfever_DHS_fitted$Date)



cases <- left_join(cases, medfever_DHS_fitted, by = c("District", "Date"))

cases$cases_linear_trtseeking_adj <- (cases$mal_cases / cases$fitted_regional_medfever) / cases$weighted_rep_rate




################################################################################



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



cases <- left_join(cases, medfever_DHS[,-2], by = c("District", "year"))

for (D in unique(cases$District))
{
    cases[which(cases$District == D & cases$Date %in% as.yearmon(seq(as.Date("2016-05-01"), as.Date("2016-12-01"), by="months"))), "medfever_regional"] <- cases[which(cases$District == D & cases$Date == "Jan 2017"), "medfever_regional"]
}


cases$cases_trtseeking_adj <- cases$mal_cases / cases$medfever_regional




cases$rep_rate <- cases$reporting_HF / cases$total_HF

# cases$cases_rep_adj <- cases$conf_rdt_mic_u5  + (cases$conf_rdt_mic_u5 * (1 - cases$rep_rate))
cases$cases_rep_adj <- cases$mal_cases  / cases$rep_rate

cases$cases_rep_weighted_adj <- cases$mal_cases / cases$weighted_rep_rate



cases$cases_both_adj <- cases$cases_rep_weighted_adj / cases$medfever_regional




cases$all_cases_rep_weighted_adj <- cases$all_cases / cases$weighted_rep_rate
cases$all_nonMal_cases_rep_weighted_adj <- cases$all_nonMal_cases / cases$weighted_rep_rate











## Fixing pouytenga district for NA value in Nov 2016
# pouytenga has NA cases in Nov 2016
# Imputing with mean of previous and subsequent months

pouytenga_ind_Nov_2017 <- which(cases$Date == "Nov 2017" & cases$District == "pouytenga")
pouytenga_ind_Jan_2018 <- which(cases$Date == "Jan 2018" & cases$District == "pouytenga")


mal_cases_pouytenga_Dec_2017 <- mean(c(cases$mal_cases[pouytenga_ind_Nov_2017],
                                       cases$mal_cases[pouytenga_ind_Jan_2018]))
all_cases_pouytenga_Dec_2017 <- mean(c(cases$all_cases[pouytenga_ind_Nov_2017],
                                       cases$all_cases[pouytenga_ind_Jan_2018]))
all_nonMal_cases_pouytenga_Dec_2017 <- mean(c(cases$all_nonMal_cases[pouytenga_ind_Nov_2017],
                                              cases$all_nonMal_cases[pouytenga_ind_Jan_2018]))
mal_ratio_pouytenga_Dec_2017 <- mean(c(cases$mal_ratio[pouytenga_ind_Nov_2017],
                                       cases$mal_ratio[pouytenga_ind_Jan_2018]))

trt_adj_pouytenga_Dec_2017 <- mean(c(cases$cases_trtseeking_adj[pouytenga_ind_Nov_2017],
                                     cases$cases_trtseeking_adj[pouytenga_ind_Jan_2018]))
rep_adj_pouytenga_Dec_2017 <- mean(c(cases$cases_rep_adj[pouytenga_ind_Nov_2017],
                                     cases$cases_rep_adj[pouytenga_ind_Jan_2018]))
weight_rep_adj_pouytenga_Dec_2017 <- mean(c(cases$cases_rep_weighted_adj[pouytenga_ind_Nov_2017],
                                            cases$cases_rep_weighted_adj[pouytenga_ind_Jan_2018]))
both_adj_pouytenga_Dec_2017 <- mean(c(cases$cases_both_adj[pouytenga_ind_Nov_2017],
                                      cases$cases_both_adj[pouytenga_ind_Jan_2018]))

all_cases_weight_rep_adj_pouytenga_Dec_2017 <- mean(c(cases$all_cases_rep_weighted_adj[pouytenga_ind_Nov_2017],
                                                      cases$all_cases_rep_weighted_adj[pouytenga_ind_Jan_2018]))
all_nonMal_weight_rep_adj_pouytenga_Dec_2017 <- mean(c(cases$all_nonMal_cases_rep_weighted_adj[pouytenga_ind_Nov_2017],
                                                       cases$all_nonMal_cases_rep_weighted_adj[pouytenga_ind_Jan_2018]))



pouytenga_Dec_2017_row_tmp <- cases[pouytenga_ind_Nov_2017,]
pouytenga_Dec_2017_row_tmp$month <- 12
pouytenga_Dec_2017_row_tmp$Date <- as.yearmon("Dec 2017")
pouytenga_Dec_2017_row_tmp$mal_cases <- mal_cases_pouytenga_Dec_2017
pouytenga_Dec_2017_row_tmp$all_cases <- all_cases_pouytenga_Dec_2017
pouytenga_Dec_2017_row_tmp$all_nonMal_cases <- all_nonMal_cases_pouytenga_Dec_2017
pouytenga_Dec_2017_row_tmp$mal_ratio <- mal_ratio_pouytenga_Dec_2017
pouytenga_Dec_2017_row_tmp$cases_trtseeking_adj <- trt_adj_pouytenga_Dec_2017
pouytenga_Dec_2017_row_tmp$cases_rep_adj <- rep_adj_pouytenga_Dec_2017
pouytenga_Dec_2017_row_tmp$cases_rep_weighted_adj <- weight_rep_adj_pouytenga_Dec_2017
pouytenga_Dec_2017_row_tmp$cases_both_adj <- both_adj_pouytenga_Dec_2017
pouytenga_Dec_2017_row_tmp$all_cases_rep_weighted_adj <- all_cases_weight_rep_adj_pouytenga_Dec_2017
pouytenga_Dec_2017_row_tmp$all_nonMal_cases_rep_weighted_adj <- all_nonMal_weight_rep_adj_pouytenga_Dec_2017


cases <- rbind(cases[1:pouytenga_ind_Nov_2017,],
               pouytenga_Dec_2017_row_tmp,
               cases[pouytenga_ind_Jan_2018:nrow(cases),])





# Pre-allocating data.frame for STL results

STL_result_DF_norm <- data.frame()


# Decompose each district

for (DS in sort(unique(cases$District)))
{
    # Selecting current district + ordering by date
    
    cases_dist <- cases[which(cases$District == DS),]
    cases_dist <- cases_dist[order(cases_dist$Date),]
    
    
    
    # Create crude incidence variable
    # Standardize it, make it a timeseries object, and decompose it
    mal_cases_norm <- getNormalized(cases_dist$mal_cases)
    mal_cases_norm_ts <- ts(mal_cases_norm, start = c(2015, 1), deltat = 1/12)
    cases_stl <- stlplus(mal_cases_norm_ts, s.window="periodic")
    
    
    all_cases_norm <- getNormalized(cases_dist$all_cases)
    all_cases_norm_ts <- ts(all_cases_norm, start = c(2015, 1), deltat = 1/12)
    all_cases_stl <- stlplus(all_cases_norm_ts, s.window="periodic")
    
    
    all_nonMal_cases_norm <- getNormalized(cases_dist$all_nonMal_cases)
    all_nonMal_cases_norm_ts <- ts(all_nonMal_cases_norm, start = c(2015, 1), deltat = 1/12)
    all_nonMal_cases_stl <- stlplus(all_nonMal_cases_norm_ts, s.window="periodic")
    
    
    mal_ratio_norm <- getNormalized(cases_dist$mal_ratio)
    mal_ratio_norm_ts <- ts(mal_ratio_norm, start = c(2015, 1), deltat = 1/12)
    ratio_stl <- stlplus(mal_ratio_norm_ts, s.window="periodic")
    
    
    trt_adj_norm <- getNormalized(cases_dist$cases_trtseeking_adj)
    trt_adj_norm_ts <- ts(trt_adj_norm, start = c(2015, 1), deltat = 1/12)
    trt_adj_stl <- stlplus(trt_adj_norm_ts, s.window="periodic")
    
    
    rep_adj_norm <- getNormalized(cases_dist$cases_rep_adj)
    rep_adj_norm_ts <- ts(rep_adj_norm, start = c(2015, 1), deltat = 1/12)
    rep_adj_stl <- stlplus(rep_adj_norm_ts, s.window="periodic")
    
    
    rep_weighted_adj_norm <- getNormalized(cases_dist$cases_rep_weighted_adj)
    rep_weighted_adj_norm_ts <- ts(rep_weighted_adj_norm, start = c(2015, 1), deltat = 1/12)
    rep_weighted_adj_stl <- stlplus(rep_weighted_adj_norm_ts, s.window="periodic")
    
    
    both_adj_norm <- getNormalized(cases_dist$cases_both_adj)
    both_adj_norm_ts <- ts(both_adj_norm, start = c(2015, 1), deltat = 1/12)
    both_adj_stl <- stlplus(both_adj_norm_ts, s.window="periodic")
    
    
    all_cases_rep_weighted_adj_norm <- getNormalized(cases_dist$all_cases_rep_weighted_adj)
    all_cases_rep_weighted_adj_norm_ts <- ts(all_cases_rep_weighted_adj_norm, start = c(2015, 1), deltat = 1/12)
    all_cases_rep_weighted_adj_stl <- stlplus(all_cases_rep_weighted_adj_norm_ts, s.window="periodic")
    
    all_nonMal_rep_weighted_adj_norm <- getNormalized(cases_dist$all_nonMal_cases_rep_weighted_adj)
    all_nonMal_rep_weighted_adj_norm_ts <- ts(all_nonMal_rep_weighted_adj_norm, start = c(2015, 1), deltat = 1/12)
    all_nonMal_rep_weighted_adj_stl <- stlplus(all_nonMal_rep_weighted_adj_norm_ts, s.window="periodic")
    
    
    
    
    # Create date object to create master data-frame for decomposition outputs
    
    dates <- seq(as.Date("2015-01-01"), as.Date("2018-12-01"), by="months")
    
    cases_stl_ts <- as.data.frame(cases_stl$data[,1:4]) 
    cases_stl_ts$type <- rep("cases", nrow(cases_stl_ts))
    cases_stl_ts$lty <- rep("solid", nrow(cases_stl_ts))
    cases_stl_ts$dates <- dates
    cases_stl_ts$MK_p <- smk.test(mal_cases_norm_ts)$p.value
    cases_stl_ts$sens_slope <- sea.sens.slope(mal_cases_norm_ts)
    
    
    all_cases_stl_ts <- as.data.frame(all_cases_stl$data[,1:4]) 
    all_cases_stl_ts$type <- rep("allout", nrow(all_cases_stl_ts))
    all_cases_stl_ts$lty <- rep("solid", nrow(all_cases_stl_ts))
    all_cases_stl_ts$dates <- dates
    all_cases_stl_ts$MK_p <- smk.test(all_cases_norm_ts)$p.value
    all_cases_stl_ts$sens_slope <- sea.sens.slope(all_cases_norm_ts)
    
    
    all_nonMal_cases_stl_ts <- as.data.frame(all_nonMal_cases_stl$data[,1:4]) 
    all_nonMal_cases_stl_ts$type <- rep("all_nonMal", nrow(all_nonMal_cases_stl_ts))
    all_nonMal_cases_stl_ts$lty <- rep("solid", nrow(all_nonMal_cases_stl_ts))
    all_nonMal_cases_stl_ts$dates <- dates
    all_nonMal_cases_stl_ts$MK_p <- smk.test(all_nonMal_cases_norm_ts)$p.value
    all_nonMal_cases_stl_ts$sens_slope <- sea.sens.slope(all_nonMal_cases_norm_ts)
    
    
    ratio_stl_ts <- as.data.frame(ratio_stl$data[,1:4]) 
    ratio_stl_ts$type <- rep("ratio", nrow(ratio_stl_ts))
    ratio_stl_ts$lty <- rep("longdash", nrow(ratio_stl_ts))
    ratio_stl_ts$dates <- dates
    ratio_stl_ts$MK_p <- smk.test(mal_ratio_norm_ts)$p.value
    ratio_stl_ts$sens_slope <- sea.sens.slope(mal_ratio_norm_ts)
    
    
    trt_adj_stl_ts <- as.data.frame(trt_adj_stl$data[,1:4]) 
    trt_adj_stl_ts$type <- rep("trt_adj", nrow(trt_adj_stl_ts))
    trt_adj_stl_ts$lty <- rep("longdash", nrow(trt_adj_stl_ts))
    trt_adj_stl_ts$dates <- dates
    trt_adj_stl_ts$MK_p <- smk.test(trt_adj_norm_ts)$p.value
    trt_adj_stl_ts$sens_slope <- sea.sens.slope(trt_adj_norm_ts)
    
    
    rep_adj_stl_ts <- as.data.frame(rep_adj_stl$data[,1:4]) 
    rep_adj_stl_ts$type <- rep("rep_adj", nrow(rep_adj_stl_ts))
    rep_adj_stl_ts$lty <- rep("longdash", nrow(rep_adj_stl_ts))
    rep_adj_stl_ts$dates <- dates
    rep_adj_stl_ts$MK_p <- smk.test(rep_adj_norm_ts)$p.value
    rep_adj_stl_ts$sens_slope <- sea.sens.slope(rep_adj_norm_ts)
    
    
    rep_weighted_adj_stl_ts <- as.data.frame(rep_weighted_adj_stl$data[,1:4]) 
    rep_weighted_adj_stl_ts$type <- rep("rep_weighted_adj", nrow(rep_weighted_adj_stl_ts))
    rep_weighted_adj_stl_ts$lty <- rep("longdash", nrow(rep_weighted_adj_stl_ts))
    rep_weighted_adj_stl_ts$dates <- dates
    rep_weighted_adj_stl_ts$MK_p <- smk.test(rep_weighted_adj_norm_ts)$p.value
    rep_weighted_adj_stl_ts$sens_slope <- sea.sens.slope(rep_weighted_adj_norm_ts)
    
    
    both_adj_stl_ts <- as.data.frame(both_adj_stl$data[,1:4]) 
    both_adj_stl_ts$type <- rep("both_adj", nrow(both_adj_stl_ts))
    both_adj_stl_ts$lty <- rep("longdash", nrow(both_adj_stl_ts))
    both_adj_stl_ts$dates <- dates
    both_adj_stl_ts$MK_p <- smk.test(both_adj_norm_ts)$p.value
    both_adj_stl_ts$sens_slope <- sea.sens.slope(both_adj_norm_ts)
    
    
    all_cases_rep_weighted_adj_stl_ts <- as.data.frame(all_cases_rep_weighted_adj_stl$data[,1:4]) 
    all_cases_rep_weighted_adj_stl_ts$type <- rep("allout_rep_weighted_adj", nrow(all_cases_rep_weighted_adj_stl_ts))
    all_cases_rep_weighted_adj_stl_ts$lty <- rep("longdash", nrow(all_cases_rep_weighted_adj_stl_ts))
    all_cases_rep_weighted_adj_stl_ts$dates <- dates
    all_cases_rep_weighted_adj_stl_ts$MK_p <- smk.test(all_cases_rep_weighted_adj_norm_ts)$p.value
    all_cases_rep_weighted_adj_stl_ts$sens_slope <- sea.sens.slope(all_cases_rep_weighted_adj_norm_ts)
    
    
    all_nonMal_rep_weighted_adj_stl_ts <- as.data.frame(all_nonMal_rep_weighted_adj_stl$data[,1:4]) 
    all_nonMal_rep_weighted_adj_stl_ts$type <- rep("all_nonMal_rep_weighted_adj", nrow(all_nonMal_rep_weighted_adj_stl_ts))
    all_nonMal_rep_weighted_adj_stl_ts$lty <- rep("longdash", nrow(all_nonMal_rep_weighted_adj_stl_ts))
    all_nonMal_rep_weighted_adj_stl_ts$dates <- dates
    all_nonMal_rep_weighted_adj_stl_ts$MK_p <- smk.test(all_nonMal_rep_weighted_adj_norm_ts)$p.value
    all_nonMal_rep_weighted_adj_stl_ts$sens_slope <- sea.sens.slope(all_nonMal_rep_weighted_adj_norm_ts)
    
    
    
    
    # Create output data-frame from stl outputs of both malaria variables
    
    STL_result_DF_DS_norm <- rbind(cases_stl_ts,
                                   all_cases_stl_ts,
                                   all_nonMal_cases_stl_ts,
                                   ratio_stl_ts,
                                   trt_adj_stl_ts,
                                   rep_adj_stl_ts,
                                   rep_weighted_adj_stl_ts,
                                   both_adj_stl_ts,
                                   all_cases_rep_weighted_adj_stl_ts,
                                   all_nonMal_rep_weighted_adj_stl_ts)
    STL_result_DF_DS_norm$District <- DS
    
    
    # Adding dates of SMC rounds for plotting
    
    SMC_DF_DS <- cases_dist[, c("Date", "SMC.received")]
    SMC_DF_DS$Date <- as.Date(SMC_DF_DS$Date)
    STL_result_DF_DS_norm <- left_join(STL_result_DF_DS_norm, SMC_DF_DS,
                                       by = c("dates" = "Date"))
    
    # Appending to District data-frame to master data-frame and continuing
    
    STL_result_DF_norm <- rbind(STL_result_DF_norm, STL_result_DF_DS_norm)
}

# Adding SMC rollout group info for each district
STL_result_DF_norm <- left_join(STL_result_DF_norm, unique(cases[,c("District", "Region")]),
                                by = "District")

# Ordering by rollout group and then by district
STL_result_DF_norm <- STL_result_DF_norm[order(STL_result_DF_norm$Region, STL_result_DF_norm$District),]


# Setting line-types and district factors for plotting
STL_result_DF_norm$lty <- factor(STL_result_DF_norm$lty, levels = c("solid", "longdash"))
STL_result_DF_norm$District <- factor(STL_result_DF_norm$District, levels = sort(unique(cases$District)))






################################################################################

# STL_result_DF_norm_cases <- unique(STL_result_DF_norm[which(STL_result_DF_norm$type == "cases"), c("District", "Region", "MK_p", "sens_slope")])
# names(STL_result_DF_norm_cases)[3:4] <- c("MK_p_cases", "sens_slope_cases")
# 
# sum(STL_result_DF_norm_cases$MK_p_cases > .05)
# 
# 
# 
# STL_result_DF_norm_ratio <- unique(STL_result_DF_norm[which(STL_result_DF_norm$type == "ratio"), c("District", "Region", "MK_p", "sens_slope")])
# names(STL_result_DF_norm_ratio)[3:4] <- c("MK_p_ratio", "sens_slope_ratio")
# 
# sum(STL_result_DF_norm_ratio$MK_p_ratio > .05)
# 
# STL_result_DF_slopes <- full_join(STL_result_DF_norm_cases, STL_result_DF_norm_ratio, by = c("District", "Region"))



STL_result_DF_slopes <- unique(STL_result_DF_norm[,c(5,8:10,12)])







burkina_shape <- readOGR("~/Box/NU-malaria-team/data/burkina_shapefiles/Burkina Faso Health Districts SHP (130715)/BFA.shp")
burkina_shape_DF <- fortify(burkina_shape, region = "District")


district_list <- cbind(sort(levels(STL_result_DF_slopes$District)), sort(unique(burkina_shape_DF$id)))
setDistrict <- function(Row)
{
    district <- district_list[which(district_list[,1] == Row[4]), 2]
    
    return(district)
}
STL_result_DF_slopes$id <- apply(STL_result_DF_slopes, 1, setDistrict)


STL_result_DF_slopes$plotting_sens_slope <- STL_result_DF_slopes$sens_slope
STL_result_DF_slopes[which(STL_result_DF_slopes$MK_p > .05), "plotting_sens_slope"] <- NA


burkina_shape_DF_sens <- inner_join(burkina_shape_DF, STL_result_DF_slopes, by = "id")



colr <- brewer.pal(9, "RdYlBu")

facet_reporting_p <- ggplot(data = burkina_shape_DF_sens, aes(x = long, y = lat, group = group)) +
    geom_polygon(aes(fill = plotting_sens_slope), color = "white") + 
    coord_equal() + scale_fill_gradientn("Sen's slope coefficient", colors = colr) + theme_void() +
    theme(plot.title = element_text(hjust = 0.5)) + facet_wrap(~type, ncol = 3)










# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# STL_result_DF_slopes$sign <- sign(STL_result_DF_slopes$sens_slope_cases) * sign(STL_result_DF_slopes$sens_slope_ratio)
# STL_result_DF_slopes$X <- ifelse(STL_result_DF_slopes$sens_slope_cases > 0 & STL_result_DF_slopes$sens_slope_ratio < 0 &
#                                      STL_result_DF_slopes$MK_p_cases <= .05 & STL_result_DF_slopes$MK_p_ratio <= .05, -1, 0)
# 
# 
# 
# sum(STL_result_DF_slopes$sign == -1)
# 
# sum(STL_result_DF_slopes$X == -1)
# 
# 
# 
# 
# 
# 
# 
