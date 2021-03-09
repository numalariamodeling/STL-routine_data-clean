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



# Function for standardizing malaria variables

getNormalized <- function(vec)
{
    norm_vec <- (vec - mean(vec))/sd(vec)
    return(norm_vec)
}






# Loading health district dataset
# Creating under-5 population column and fixing date column

cases <- read.csv("/Users/sebas/Box/NU-malaria-team/projects/smc_impact/data/outputs/U5_DS_cases_seasonal_smc_good_rows_kalman_imputes.csv", header  = TRUE, strip.white = TRUE, stringsAsFactors = FALSE)
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
cases$mal_ratio <- cases$conf_rdt_mic_u5 / cases$allout_u5



## Fixing pouytenga district for NA value in Nov 2016
# pouytenga has NA cases in Nov 2016
# Imputing with mean of previous and subsequent months

pouytenga_ind_Nov_2017 <- which(cases$Date == "Nov 2017" & cases$District == "pouytenga")
pouytenga_ind_Jan_2018 <- which(cases$Date == "Jan 2018" & cases$District == "pouytenga")

mal_cases_pouytenga_Dec_2017 <- mean(c(cases$mal_cases[pouytenga_ind_Nov_2017],
                                       cases$mal_cases[pouytenga_ind_Jan_2018]))
mal_ratio_pouytenga_Dec_2017 <- mean(c(cases$mal_ratio[pouytenga_ind_Nov_2017],
                                       cases$mal_ratio[pouytenga_ind_Jan_2018]))

pouytenga_Dec_2017_row_tmp <- cases[pouytenga_ind_Nov_2017,]
pouytenga_Dec_2017_row_tmp$month <- 12
pouytenga_Dec_2017_row_tmp$Date <- as.yearmon("Dec 2017")
pouytenga_Dec_2017_row_tmp$mal_cases <- mal_cases_pouytenga_Dec_2017
pouytenga_Dec_2017_row_tmp$mal_ratio <- mal_ratio_pouytenga_Dec_2017


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
    
    
    
    
    # Do the same for the malaria share variable
    mal_ratio_norm <- getNormalized(cases_dist$mal_ratio)
    mal_ratio_norm_ts <- ts(mal_ratio_norm, start = c(2015, 1), deltat = 1/12)
    ratio_stl <- stlplus(mal_ratio_norm_ts, s.window="periodic")
    
    
    
    
    # Create date object to create master data-frame for decomposition outputs
    
    dates <- seq(as.Date("2015-01-01"), as.Date("2018-12-01"), by="months")
    
    cases_stl_ts <- as.data.frame(cases_stl$data[,1:4]) 
    cases_stl_ts$type <- rep("cases", nrow(cases_stl_ts))
    cases_stl_ts$lty <- rep("solid", nrow(cases_stl_ts))
    cases_stl_ts$dates <- dates
    
    ratio_stl_ts <- as.data.frame(ratio_stl$data[,1:4]) 
    ratio_stl_ts$type <- rep("ratio", nrow(ratio_stl_ts))
    ratio_stl_ts$lty <- rep("longdash", nrow(ratio_stl_ts))
    ratio_stl_ts$dates <- dates
    
    
    
    # Create output data-frame from stl outputs of both malaria variables
    
    STL_result_DF_DS_norm <- rbind(cases_stl_ts,
                                   ratio_stl_ts)
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
STL_result_DF_norm <- left_join(STL_result_DF_norm, unique(cases[,c("District", "SMC.init_date")]),
                                by = "District")
STL_result_DF_norm$SMC_rollout_group <- year(STL_result_DF_norm$SMC.init_date)

# Ordering by rollout group and then by district
STL_result_DF_norm <- STL_result_DF_norm[order(STL_result_DF_norm$SMC_rollout_group, STL_result_DF_norm$District),]


# Setting line-types and district factors for plotting
STL_result_DF_norm$lty <- factor(STL_result_DF_norm$lty, levels = c("solid", "longdash"))
STL_result_DF_norm$District <- factor(STL_result_DF_norm$District, levels = sort(unique(cases$District)))






DS_raw_coefs_DF <- data.frame(District = unique(STL_result_DF_norm$District),
                          cases_coef = rep(0, 70),
                          cases_p = rep(0, 70),
                          ratio_coef = rep(0, 70),
                          ratio_p = rep(0, 70))

for(DS in unique(STL_result_DF_norm$District))
{
    cases_raw <- STL_result_DF_norm[which(STL_result_DF_norm$District == DS &
                                                STL_result_DF_norm$type == "cases"), "raw"]
    
    ratio_raw <- STL_result_DF_norm[which(STL_result_DF_norm$District == DS &
                                                STL_result_DF_norm$type == "ratio"), "raw"]
    
    
    cases.fit <- lm(cases_raw ~ c(1:48))
    ratio.fit <- lm(ratio_raw ~ c(1:48))
    
    DS_raw_coefs_DF[which(DS_raw_coefs_DF$District == DS), c(2,4)] <- c(coef(cases.fit)[2],
                                                                        coef(ratio.fit)[2])
    DS_raw_coefs_DF[which(DS_raw_coefs_DF$District == DS), c(3,5)] <- c(anova(cases.fit)$`Pr(>F)`[1],
                                                                        anova(ratio.fit)$`Pr(>F)`[1])
    # cases_plot <- ggplot() + geom_line(aes(x = c(1:48), y = cases_raw), col = "blue") +
    #     geom_abline(slope = coef(cases.fit)[2], intercept = coef(cases.fit)[1]) + xlab("")
    # ratio_plot <- ggplot() + geom_line(aes(x = c(1:48), y = ratio_raw), col = "red") +
    #     geom_abline(slope = coef(ratio.fit)[2], intercept = coef(ratio.fit)[1]) + xlab("")
    # 
    # grid.arrange(cases_plot, ratio_plot, ncol = 2)
}

DS_raw_coefs_DF$sign <- sign(DS_raw_coefs_DF$cases_coef) * sign(DS_raw_coefs_DF$ratio_coef)
DS_raw_coefs_DF$X <- ifelse(DS_raw_coefs_DF$cases_coef > 0 & DS_raw_coefs_DF$ratio_coef < 0 &
                                DS_raw_coefs_DF$cases_p <= .05 & DS_raw_coefs_DF$ratio_p <= .05, -1, 0)

DS_raw_coefs_DF <- left_join(DS_raw_coefs_DF, unique(cases[,c("District", "SMC.init_date")]),
                         by = "District")






DS_trend_coefs_DF <- data.frame(District = unique(STL_result_DF_norm$District),
                                cases_int = rep(0, 70),
                                cases_coef = rep(0, 70),
                                cases_p = rep(0, 70),
                                ratio_int = rep(0, 70),
                                ratio_coef = rep(0, 70),
                                ratio_p = rep(0, 70))

for(DS in unique(STL_result_DF_norm$District))
{
    cases_trend <- STL_result_DF_norm[which(STL_result_DF_norm$District == DS &
                                                STL_result_DF_norm$type == "cases"), "trend"]
    
    ratio_trend <- STL_result_DF_norm[which(STL_result_DF_norm$District == DS &
                                                STL_result_DF_norm$type == "ratio"), "trend"]
    
    
    cases.fit <- lm(cases_trend ~ c(1:48))
    ratio.fit <- lm(ratio_trend ~ c(1:48))
    
    
    DS_trend_coefs_DF[which(DS_trend_coefs_DF$District == DS), c(2,5)] <- c(coef(cases.fit)[1],
                                                                            coef(ratio.fit)[1])
    DS_trend_coefs_DF[which(DS_trend_coefs_DF$District == DS), c(3,6)] <- c(coef(cases.fit)[2],
                                                                            coef(ratio.fit)[2])
    DS_trend_coefs_DF[which(DS_trend_coefs_DF$District == DS), c(4,7)] <- c(anova(cases.fit)$`Pr(>F)`[1],
                                                                            anova(ratio.fit)$`Pr(>F)`[1])
    
    
    
    # cases_plot <- ggplot() + geom_line(aes(x = c(1:48), y = cases_trend), col = "blue") +
    #     geom_abline(slope = coef(cases.fit)[2], intercept = coef(cases.fit)[1]) + xlab("")
    # ratio_plot <- ggplot() + geom_line(aes(x = c(1:48), y = ratio_trend), col = "red") +
    #     geom_abline(slope = coef(ratio.fit)[2], intercept = coef(ratio.fit)[1]) + xlab("")
    # 
    # grid.arrange(cases_plot, ratio_plot, ncol = 2)
}

DS_trend_coefs_DF$sign <- sign(DS_trend_coefs_DF$cases_coef) * sign(DS_trend_coefs_DF$ratio_coef)
DS_trend_coefs_DF$X <- ifelse(DS_trend_coefs_DF$cases_coef > 0 & DS_trend_coefs_DF$ratio_coef < 0 &
                                  DS_trend_coefs_DF$cases_p <= .05 & DS_trend_coefs_DF$ratio_p <= .05, -1, 0)

# DS_trend_coefs_DF <- left_join(DS_trend_coefs_DF, unique(cases[,c("District", "SMC.init_date")]),
#                                by = "District")



DS_trend_coefs_plotting <- STL_result_DF_norm[,c("District", "trend", "type", "lty",
                                                 "dates", "SMC.received", "SMC_rollout_group")]

DS_trend_coefs_plotting_cases <- left_join(DS_trend_coefs_plotting[which(DS_trend_coefs_plotting$type == "cases"),],
                                           DS_trend_coefs_DF[,c(1:3)], by = "District")
names(DS_trend_coefs_plotting_cases)[8:9] <- c("int", "slope")

DS_trend_coefs_plotting_ratio <- left_join(DS_trend_coefs_plotting[which(DS_trend_coefs_plotting$type == "ratio"),],
                                           DS_trend_coefs_DF[,c(1,5:6)], by = "District")
names(DS_trend_coefs_plotting_ratio)[8:9] <- c("int", "slope")

DS_trend_coefs_plotting <- rbind(DS_trend_coefs_plotting_cases, DS_trend_coefs_plotting_ratio)



for (i in seq(1,70,20))
{
    # Select 4 district for plotting
    # Dont want alphabetical but my SMC rollout group THEN alphabetical from there
    
    DS_list <- unique(DS_trend_coefs_plotting$District)[i:(i+19)]
    
    plotting_DF <- DS_trend_coefs_plotting[which(DS_trend_coefs_plotting$District %in% DS_list),]
    plotting_DF$facet_name <- paste(plotting_DF$District,
                                    plotting_DF$SMC_rollout_group, sep = " ")
    plotting_DF$facet_name <- factor(plotting_DF$facet_name, levels = unique(plotting_DF$facet_name))
    
    
    plot_trends <- ggplot(plotting_DF, aes(x = dates, y = trend)) +
        geom_line(aes(group = type, col = District, linetype = lty), show.legend = FALSE) +
        geom_line(inherit.aes = FALSE, aes(x = dates, y = (int + c(1:48)*slope),
                                           group = type, linetype = lty),
                  col = "black", alpha = .6, show.legend = FALSE) +
        scale_color_discrete("") + xlab("") + facet_wrap(~facet_name, ncol = 4) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"),
              panel.spacing.x = unit(6, "mm"))
    
    
    
    
    # Saving SI plots
    pdf(paste("/Users/sebas/Desktop/SI_figures_cleaned/X/SI_X_", i, ".pdf", sep = ""))
    print(plot_trends)
    dev.off()
    
}



