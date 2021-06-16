#############################################
##  Seasonal Trend Decomposition Figure 7  ##
#############################################
#
# Description:
#   Creating the output for figure 7 in the main text
#   Seasonal trend decomposition for all districts grouped by SMC rollout group
#   Decomposing crude incidence and malaria share of outpatient visits variables
#   Only trend component is plotted, for individual districts
#       and aggregate rollout group
#   All data was standardized for comparison between groups and
#       for comparison to figure 6
#       
#   Figure was further cleaned and assembled in Adobe Illustrator
#
#
#
#
#  Sebastian Rodriguez (sebastian@rodriguez.cr)
#  Last edited Mar 09, 2021
#



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
require("maditr")




# Function for standardizing malaria variables

getNormalized <- function(vec)
{
    norm_vec <- (vec - mean(vec))/sd(vec)
    return(norm_vec)
}






# Loading health district dataset
# Creating under-5 population column and fixing date column

# cases <- read.csv("~/Box/NU-malaria-team/projects/smc_impact/data/outputs/U5_DS_cases_seasonal_smc_good_rows_MA_imputes_testing.csv", header  = TRUE, strip.white = TRUE, stringsAsFactors = FALSE)
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


#########################################################################3

cases <- cases[order(cases$District, cases$Date),]

cases$mal_cases <- cases$conf_rdt_mic_u5 / (cases$U5_pop/1000)
cases$all_cases <- cases$allout_u5 / (cases$U5_pop/1000)
cases$all_nonMal_cases <- (cases$allout_u5 - cases$conf_rdt_mic_u5) / (cases$U5_pop/1000)
cases$mal_ratio <- cases$conf_rdt_mic_u5 / cases$allout_u5




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


cases$cases_trtseeking_adj_raw <- cases$conf_rdt_mic_u5 / cases$medfever_regional
cases$cases_trtseeking_adj <- cases$cases_trtseeking_adj_raw / (cases$U5_pop / 1000)



cases$rep_rate <- cases$reporting_HF / cases$total_HF

# cases$cases_rep_adj <- cases$conf_rdt_mic_u5  + (cases$conf_rdt_mic_u5 * (1 - cases$rep_rate))
cases$cases_rep_adj_raw <- cases$conf_rdt_mic_u5  / cases$rep_rate
cases$cases_rep_adj <- cases$cases_rep_adj_raw / (cases$U5_pop / 1000)


cases$cases_rep_weighted_adj_raw <- cases$conf_rdt_mic_u5 / cases$weighted_rep_rate
cases$cases_rep_weighted_adj <- cases$cases_rep_weighted_adj_raw / (cases$U5_pop / 1000)


cases$cases_both_adj_raw <- cases$cases_rep_weighted_adj / cases$medfever_regional
cases$cases_both_adj <- cases$cases_both_adj_raw / (cases$U5_pop / 1000)



#########################################################################3




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


cases <- rbind(cases[1:pouytenga_ind_Nov_2017,],
               pouytenga_Dec_2017_row_tmp,
               cases[pouytenga_ind_Jan_2018:nrow(cases),])






############################################################################



cases$mal_cases_norm <- getNormalized(cases$mal_cases)
cases$all_cases_norm <- getNormalized(cases$all_cases)
cases$all_nonMal_cases_norm <- getNormalized(cases$all_nonMal_cases)
cases$mal_ratio_norm <- getNormalized(cases$mal_ratio)
cases$cases_rep_weighted_adj_norm <- getNormalized(cases$cases_rep_weighted_adj)
cases$both_adj_norm <- getNormalized(cases$cases_both_adj)




############################################################################





# Performing STL on each individual district
# For both malaria variables

STL_result_DF <- data.frame()

for (DS in sort(unique(cases$District)))
{
    cases_dist <- cases[which(cases$District == DS),]
    cases_dist <- cases_dist[order(cases_dist$Date),]
    
 
    
    mal_cases_ts <- ts(cases_dist$mal_cases_norm, start = c(2015, 1), deltat = 1/12)
    cases_stl <- stlplus(mal_cases_ts, s.window = "periodic")
    
    
    all_cases_ts <- ts(cases_dist$all_cases_norm, start = c(2015, 1), deltat = 1/12)
    all_cases_stl <- stlplus(all_cases_ts, s.window="periodic")
    
    
    all_nonMal_cases_ts <- ts(cases_dist$all_nonMal_cases_norm, start = c(2015, 1), deltat = 1/12)
    all_nonMal_cases_stl <- stlplus(all_nonMal_cases_ts, s.window="periodic")
    
    
    mal_ratio_ts <- ts(cases_dist$mal_ratio_norm, start = c(2015, 1), deltat = 1/12)
    ratio_stl <- stlplus(mal_ratio_ts, s.window = "periodic")
    
    
    cases_rep_weighted_adj_ts <- ts(cases_dist$cases_rep_weighted_adj_norm, start = c(2015, 1), deltat = 1/12)
    cases_rep_weighted_adj_stl <- stlplus(cases_rep_weighted_adj_ts, s.window = "periodic")
    
    
    both_adj_ts <- ts(cases_dist$both_adj_norm, start = c(2015, 1), deltat = 1/12)
    both_adj_stl <- stlplus(both_adj_ts, s.window="periodic")
    
    
    
    # Assembling data frame for decompositon of all districts
    
    dates <- seq(as.Date("2015-01-01"), as.Date("2018-12-01"), by = "months")
    
    
    cases_stl_ts <- as.data.frame(cases_stl$data[,1:4])
    cases_stl_ts$type <- rep("cases", nrow(cases_stl_ts))
    cases_stl_ts$dates <- dates
    
    
    all_cases_stl_ts <- as.data.frame(all_cases_stl$data[,1:4])
    all_cases_stl_ts$type <- rep("all_cases", nrow(all_cases_stl_ts))
    all_cases_stl_ts$dates <- dates
    
    
    all_nonMal_cases_stl_ts <- as.data.frame(all_nonMal_cases_stl$data[,1:4])
    all_nonMal_cases_stl_ts$type <- rep("all_nonMal_cases", nrow(all_nonMal_cases_stl_ts))
    all_nonMal_cases_stl_ts$dates <- dates
    
    
    ratio_stl_ts <- as.data.frame(ratio_stl$data[,1:4])
    ratio_stl_ts$type <- rep("ratio", nrow(ratio_stl_ts))
    ratio_stl_ts$dates <- dates
    
    
    cases_rep_weighted_adj_stl_ts <- as.data.frame(cases_rep_weighted_adj_stl$data[,1:4])
    cases_rep_weighted_adj_stl_ts$type <- rep("rep_weighted_adj", nrow(cases_rep_weighted_adj_stl_ts))
    cases_rep_weighted_adj_stl_ts$dates <- dates
    
    
    both_adj_stl_ts <- as.data.frame(both_adj_stl$data[,1:4])
    both_adj_stl_ts$type <- rep("both_adj_cases", nrow(both_adj_stl_ts))
    both_adj_stl_ts$dates <- dates

    
    
    
    STL_result_DF_DS <- rbind(cases_stl_ts,
                              all_cases_stl_ts,
                              all_nonMal_cases_stl_ts,
                              ratio_stl_ts,
                              cases_rep_weighted_adj_stl_ts,
                              both_adj_stl_ts)
    STL_result_DF_DS$District <- DS
    STL_result_DF_DS$Region <- unique(cases_dist$Region)
    STL_result_DF <- rbind(STL_result_DF,
                           STL_result_DF_DS)
    
    
}






####################################################################




# Creating dataset for aggregated districts by Region

cases_Region_agg <- ddply(cases, c(.(Region), .(Date)), summarise,
                          conf_rdt_mic_u5 = sum(conf_rdt_mic_u5),
                          allout_u5 = sum(allout_u5),
                          cases_rep_weighted_adj_raw = sum(cases_rep_weighted_adj_raw),
                          cases_both_adj_raw = sum(cases_both_adj_raw),
                          U5_pop = sum(U5_pop))

cases_Region_agg <- cases_Region_agg[order(cases_Region_agg$Region, cases_Region_agg$Date),]


cases_Region_agg$mal_cases <- cases_Region_agg$conf_rdt_mic_u5 / (cases_Region_agg$U5_pop/1000)
cases_Region_agg$all_cases <- cases_Region_agg$allout_u5 / (cases_Region_agg$U5_pop/1000)
cases_Region_agg$all_nonMal_cases <- (cases_Region_agg$allout_u5 - cases_Region_agg$conf_rdt_mic_u5) / (cases_Region_agg$U5_pop/1000)
cases_Region_agg$mal_ratio <- cases_Region_agg$conf_rdt_mic_u5 / cases_Region_agg$allout_u5


cases_Region_agg$cases_rep_weighted_adj = cases_Region_agg$cases_rep_weighted_adj_raw / (cases_Region_agg$U5_pop/1000)
cases_Region_agg$cases_both_adj = cases_Region_agg$cases_both_adj_raw / (cases_Region_agg$U5_pop/1000)




## Standardizing both malaria variables across all districts

cases_Region_agg$mal_cases_norm <- getNormalized(cases_Region_agg$mal_cases)
cases_Region_agg$all_cases_norm <- getNormalized(cases_Region_agg$all_cases)
cases_Region_agg$all_nonMal_cases_norm <- getNormalized(cases_Region_agg$all_nonMal_cases)
cases_Region_agg$mal_ratio_norm <- getNormalized(cases_Region_agg$mal_ratio)
cases_Region_agg$cases_rep_weighted_adj_norm <- getNormalized(cases_Region_agg$cases_rep_weighted_adj)
cases_Region_agg$both_adj_norm <- getNormalized(cases_Region_agg$cases_both_adj)




# Performing STL on each aggregated SMC rollout group
# i.e. aggregate of all districts belonging to the same rollout group
# For both malaria variables 

STL_result_DF_Region <- data.frame()

for (R in sort(unique(cases_Region_agg$Region)))
{
    cases_R <- cases_Region_agg[which(cases_Region_agg$Region == R),]
    cases_R <- cases_R[order(cases_R$Date),]
    
    
    mal_cases_ts <- ts(cases_R$mal_cases_norm, start = c(2015, 1), deltat = 1/12)
    cases_stl <- stlplus(mal_cases_ts, s.window = "periodic")
    
    
    all_cases_ts <- ts(cases_R$all_cases_norm, start = c(2015, 1), deltat = 1/12)
    all_cases_stl <- stlplus(all_cases_ts, s.window="periodic")
    
    
    all_nonMal_cases_ts <- ts(cases_R$all_nonMal_cases_norm, start = c(2015, 1), deltat = 1/12)
    all_nonMal_cases_stl <- stlplus(all_nonMal_cases_ts, s.window="periodic")
    
    
    mal_ratio_ts <- ts(cases_R$mal_ratio_norm, start = c(2015, 1), deltat = 1/12)
    ratio_stl <- stlplus(mal_ratio_ts, s.window = "periodic")
    
    
    cases_rep_weighted_adj_ts <- ts(cases_R$cases_rep_weighted_adj_norm, start = c(2015, 1), deltat = 1/12)
    cases_rep_weighted_adj_stl <- stlplus(cases_rep_weighted_adj_ts, s.window = "periodic")
    
    
    both_adj_ts <- ts(cases_R$both_adj_norm, start = c(2015, 1), deltat = 1/12)
    both_adj_stl <- stlplus(both_adj_ts, s.window="periodic")
    
    
    
    
    
    # Assembling data frame for decompositon of all SMC rollout groups
    
    dates <- seq(as.Date("2015-01-01"), as.Date("2018-12-01"), by = "months")
    
    
    cases_stl_ts <- as.data.frame(cases_stl$data[,1:4])
    cases_stl_ts$type <- rep("cases", nrow(cases_stl_ts))
    cases_stl_ts$dates <- dates
    
    
    all_cases_stl_ts <- as.data.frame(all_cases_stl$data[,1:4])
    all_cases_stl_ts$type <- rep("all_cases", nrow(all_cases_stl_ts))
    all_cases_stl_ts$dates <- dates
    
    
    all_nonMal_cases_stl_ts <- as.data.frame(all_nonMal_cases_stl$data[,1:4])
    all_nonMal_cases_stl_ts$type <- rep("all_nonMal_cases", nrow(all_nonMal_cases_stl_ts))
    all_nonMal_cases_stl_ts$dates <- dates
    
    
    ratio_stl_ts <- as.data.frame(ratio_stl$data[,1:4])
    ratio_stl_ts$type <- rep("ratio", nrow(ratio_stl_ts))
    ratio_stl_ts$dates <- dates
    
    
    cases_rep_weighted_adj_stl_ts <- as.data.frame(cases_rep_weighted_adj_stl$data[,1:4])
    cases_rep_weighted_adj_stl_ts$type <- rep("rep_weighted_adj", nrow(cases_rep_weighted_adj_stl_ts))
    cases_rep_weighted_adj_stl_ts$dates <- dates
    
    
    both_adj_stl_ts <- as.data.frame(both_adj_stl$data[,1:4])
    both_adj_stl_ts$type <- rep("both_adj_cases", nrow(both_adj_stl_ts))
    both_adj_stl_ts$dates <- dates
    
    
    STL_result_DF_Region_R <- rbind(cases_stl_ts,
                                    all_cases_stl_ts,
                                    all_nonMal_cases_stl_ts,
                                    ratio_stl_ts,
                                    cases_rep_weighted_adj_stl_ts,
                                    both_adj_stl_ts)
    STL_result_DF_Region_R$District <- paste("master", R, sep = " ")
    STL_result_DF_Region_R$Region <- R
    STL_result_DF_Region <- rbind(STL_result_DF_Region,
                                  STL_result_DF_Region_R)

}





#########################################################################


# Joining both dataframes for plotting

STL_result_DF_master <- rbind(STL_result_DF, STL_result_DF_Region)





for (t in unique(STL_result_DF$type))
{
    STL_result_DF_type <- STL_result_DF[which(STL_result_DF$type == t),]
    
    
    
    plot_type_trend_for_grid <- ggplot(data = STL_result_DF_type,
                                        aes(x = dates, y = trend)) +
        geom_line(aes(group = District,
                      color = as.factor(District)),
                  show.legend = FALSE) +
        xlab("") + ylab("") + facet_wrap(~Region, scales = "free") +
        theme_classic()
    
    
    
    pdf(paste("~/OneDrive/Desktop/clean_trends_20210608/", t, "_trends_per_region.pdf", sep = ""), width = 11, height = 7)
    print(plot_type_trend_for_grid)
    dev.off()

}




















