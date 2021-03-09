#############################################
##  Seasonal Trend Decomposition Figure 6  ##
#############################################
#
# Description:
#   Creating the output for figure 6 in the main text
#   Seasonal trend decomposition for 4 representative districts for SMC implementation
#   Decomposing crude incidence and malaria share of outpatient visits variables
#   after standardization
#
#
#
#  Sebastian Rodriguez (sebastian@rodriguez.cr)
#  Last edited Mar 09, 2021
#
#



rm(list = ls(all = TRUE))


require("ggplot2")
require("plyr")
require("dplyr")
require("zoo")
require("gridExtra")

require("stlplus")


# Loading health district dataset
# Creating under-5 population column and fixing date column

cases <- read.csv("~/Box/NU-malaria-team/projects/smc_impact/data/outputs/U5_DS_cases_seasonal_smc_good_rows_kalman_imputes.csv", header  = TRUE, strip.white = TRUE, stringsAsFactors = FALSE)
cases$U5_pop <- cases$District.Pop * .18
cases$Date <- as.yearmon(cases$Date)





# Function for standardizing malaria variables

getNormalized <- function(vec)
{
    norm_vec <- (vec - mean(vec))/sd(vec)
    return(norm_vec)
}




# Setting the 4 representative districts of SMC rollout groups 2015-2018

DS_list <- c("tougouri", "koudougou", "lena", "gaoua")




# Pre-allocating data.frame for STL results
STL_result_DF_norm <- data.frame()

# Decompose each district in list
for (DS in DS_list)
{
    # Selecting current district + ordering by date
    cases_dist <- cases[which(cases$District == DS),]
    cases_dist <- cases_dist[order(cases_dist$Date),]
    
    
    
    
    # Create crude incidence variable
    # Standardize it, make it a timeseries object, and decompose it
    mal_cases <- cases_dist$conf_rdt_mic_u5 / (cases_dist$U5_pop/1000)
    mal_cases_norm <- getNormalized(mal_cases)
    mal_cases_norm_ts <- ts(mal_cases_norm, start = c(2015, 1), deltat = 1/12)
    cases_stl <- stlplus(mal_cases_norm_ts, s.window="periodic")
    
    
    # Do the same for the malaria share variable
    mal_ratio <- cases_dist$conf_rdt_mic_u5 / cases_dist$allout_u5
    mal_ratio_norm <- getNormalized(mal_ratio)
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
    # And precipitation
    SMC_DF_DS <- cases_dist[, c("Date", "SMC.received", "precip_era5")]
    SMC_DF_DS$precip_norm <- getNormalized(SMC_DF_DS$precip_era5)
    SMC_DF_DS$Date <- as.Date(SMC_DF_DS$Date)
    STL_result_DF_DS_norm <- left_join(STL_result_DF_DS_norm, SMC_DF_DS,
                                       by = c("dates" = "Date"))
    
    
    
    
    
    
    # Appending to District data-frame to master data-frame and continuing
    STL_result_DF_norm <- rbind(STL_result_DF_norm, STL_result_DF_DS_norm)
}

# Setting line-types and district factors for plotting
STL_result_DF_norm$lty <- factor(STL_result_DF_norm$lty, levels = c("solid", "longdash"))
STL_result_DF_norm$District <- factor(STL_result_DF_norm$District, levels = DS_list)



color_list <- c("#E41A1C", "#4DAF4A", "#377EB8", "#D95F02")



# Create plots for original data and the 3 extracted components from the decomposition
# We plot both malaria variables on the same plot for ease of visualization

plot_data <- ggplot(STL_result_DF_norm, aes(x = dates, y = raw)) +
    geom_line(aes(group = type, col = District, linetype = lty), size = 1.3, show.legend = FALSE) + scale_color_discrete("") +
    geom_vline(data = STL_result_DF_norm[which(STL_result_DF_norm$SMC.received == 1),],
               aes(xintercept = dates), col = "black", linetype = "solid", alpha = .25, show.legend = FALSE) +
    scale_color_manual(values = color_list) +
    xlab("") + facet_wrap(~District, ncol = 4) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          panel.spacing.x = unit(4, "mm"))


plot_seasonal <- ggplot(STL_result_DF_norm, aes(x = dates, y = seasonal)) +
    geom_line(inherit.aes = FALSE, aes(x = dates, y = precip_norm), col = "black", show.legend = FALSE) +
    geom_line(aes(group = type, col = District, linetype = as.factor(lty)), size = 1.3, show.legend = FALSE) +
    scale_color_manual(values = color_list) +
    geom_vline(data = STL_result_DF_norm[which(STL_result_DF_norm$SMC.received == 1),],
               aes(xintercept = dates), col = "black", linetype = "solid", alpha = .25, show.legend = FALSE) +
    xlab("") + facet_wrap(~District, ncol = 4) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          panel.spacing.x = unit(4, "mm"))


plot_trend <- ggplot(STL_result_DF_norm, aes(x = dates, y = trend)) +
    geom_line(aes(group = type, col = District, linetype = lty), size = 1.3, show.legend = FALSE) + scale_color_discrete("") +
    geom_vline(data = STL_result_DF_norm[which(STL_result_DF_norm$SMC.received == 1),],
               aes(xintercept = dates), col = "black", linetype = "solid", alpha = .25, show.legend = FALSE) +
    scale_color_manual(values = color_list) +
    xlab("") + facet_wrap(~District, ncol = 4) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          panel.spacing.x = unit(4, "mm"))


plot_rem <- ggplot(STL_result_DF_norm, aes(x = dates, y = remainder)) +
    geom_line(aes(group = type, col = District, linetype = lty), size = 1.3, show.legend = FALSE) + scale_color_discrete("") +
    geom_vline(data = STL_result_DF_norm[which(STL_result_DF_norm$SMC.received == 1),],
               aes(xintercept = dates), col = "black", linetype = "solid", alpha = .25, show.legend = FALSE) +
    scale_color_manual(values = color_list) +
    xlab("") + facet_wrap(~District, ncol = 4) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          panel.spacing.x = unit(4, "mm"))



grid.arrange(plot_data, plot_seasonal, plot_trend, plot_rem, ncol = 1)








