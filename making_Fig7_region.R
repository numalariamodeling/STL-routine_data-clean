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



rm(list = ls(all = TRUE))


require("ggplot2")
require("plyr")
require("dplyr")
require("stringr")
require("zoo")
require("gridExtra")
require("lubridate")


require("RColorBrewer")
require("colorRamps")

require("stlplus")






# Function for standardizing malaria variables

getNormalized <- function(vec)
{
    norm_vec <- (vec - mean(vec))/sd(vec)
    return(norm_vec)
}





# Loading health district dataset
# Creating under-5 population column and fixing date column

cases <- read.csv("~/Box/NU-malaria-team/projects/smc_impact/data/outputs/U5_DS_cases_seasonal_smc_good_rows_MA_imputes.csv",
                  header  = TRUE, strip.white = TRUE, stringsAsFactors = FALSE)
cases$U5_pop <- cases$District.Pop * .18
cases$Date <- as.yearmon(cases$Date)




# Finding what date SMC was first implemented in each district

SMC.init_tmp <- cases[which(cases$SMC.received == 1),
                      c("District", "Date")]
SMC.init_tmp_ordered <- SMC.init_tmp[order(SMC.init_tmp$District, SMC.init_tmp$Date),]
SMC.init <- SMC.init_tmp_ordered[!duplicated(SMC.init_tmp_ordered$District),]
SMC.init$Date <- year(SMC.init$Date)




# Fixing so we get 2014 start dates where appropriate

DS_SMC_2014 <- c("bogande", "bousse", "boussouma", "garango", "kaya",
                 "sebba", "seguenega", "tougan")
SMC.init[which(SMC.init$District %in% DS_SMC_2014), "Date"] <- 2014



# Same for 2019

SMC.init_2019 <- data.frame(District = c("baskuy", "bogodogo", "boulmiougou",
                                         "nongr-massom", "sig-noghin"),
                            Date = rep(2019, 5))
SMC.init <- rbind(SMC.init, SMC.init_2019)
names(SMC.init)[2] <- "SMC.init_date"

cases <- inner_join(cases, SMC.init, by = "District")
cases <- cases[order(cases$District, cases$Date),]




############################################################################






## Making malaria variables

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
pouytenga_Dec_2017_row_tmp$mal_cases <- mal_cases_pouytenga_Dec_2017
pouytenga_Dec_2017_row_tmp$mal_ratio <- mal_ratio_pouytenga_Dec_2017


cases <- rbind(cases[1:pouytenga_ind_Nov_2017,],
               pouytenga_Dec_2017_row_tmp,
               cases[pouytenga_ind_Jan_2018:nrow(cases),])






## Standardizing both malaria variables across all districts

cases$mal_cases_norm <- getNormalized(cases$mal_cases)
cases$mal_ratio_norm <- getNormalized(cases$mal_ratio)



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
    
    
    
    mal_ratio_ts <- ts(cases_dist$mal_ratio_norm, start = c(2015, 1), deltat = 1/12)
    ratio_stl <- stlplus(mal_ratio_ts, s.window = "periodic")
    
    
    
    
    # Assembling data frame for decompositon of all districts
    
    dates <- seq(as.Date("2015-01-01"), as.Date("2018-12-01"), by = "months")
    
    
    cases_stl_ts <- as.data.frame(cases_stl$data[,1:4])
    cases_stl_ts$type <- rep("cases", nrow(cases_stl_ts))
    cases_stl_ts$dates <- dates
    
    ratio_stl_ts <- as.data.frame(ratio_stl$data[,1:4])
    ratio_stl_ts$type <- rep("ratio", nrow(ratio_stl_ts))
    ratio_stl_ts$dates <- dates
    
    
    STL_result_DF_DS <- rbind(cases_stl_ts,
                              ratio_stl_ts)
    STL_result_DF_DS$District <- DS
    STL_result_DF_DS$Region <- unique(cases_dist$Region)
    STL_result_DF <- rbind(STL_result_DF,
                           STL_result_DF_DS)
    
    
}


# Setting figure descriptors for ease of plotting with ggplot2

STL_result_DF$alpha <- 0.25
STL_result_DF$lwd <- 1








####################################################################




# Creating dataset for aggregated districts by Region

cases_Region_agg <- ddply(cases, c(.(Region), .(Date)), summarise,
                          conf_rdt_mic_u5 = sum(conf_rdt_mic_u5),
                          allout_u5 = sum(allout_u5),
                          U5_pop = sum(U5_pop),
                          SMC.received = mean(SMC.received))

# If a district only receives 3 rounds of SMC (not 4)

cases_Region_agg[cases_Region_agg$SMC.received > 0, "SMC.received"] <- 1



## Making malaria variables

cases_Region_agg$mal_cases <- cases_Region_agg$conf_rdt_mic_u5 / (cases_Region_agg$U5_pop/1000)
cases_Region_agg$mal_ratio <- cases_Region_agg$conf_rdt_mic_u5 / cases_Region_agg$allout_u5



## Standardizing both malaria variables across all districts

cases_Region_agg$mal_cases_norm <- getNormalized(cases_Region_agg$mal_cases)
cases_Region_agg$mal_ratio_norm <- getNormalized(cases_Region_agg$mal_ratio)






# Performing STL on each aggregated SMC rollout group
# i.e. aggregate of all districts belonging to the same rollout group
# For both malaria variables 

STL_result_DF_Region <- data.frame()

for (R in sort(unique(cases_Region_agg$Region)))
{
    cases_region <- cases_Region_agg[which(cases_Region_agg$Region == R),]
    cases_region <- cases_region[order(cases_region$Date),]
    
    
    mal_cases_ts <- ts(cases_region$mal_cases_norm, start = c(2015, 1), deltat = 1/12)
    cases_stl <- stlplus(mal_cases_ts, s.window = "periodic")
    
    
    mal_ratio_ts <- ts(cases_region$mal_ratio_norm, start = c(2015, 1), deltat = 1/12)
    ratio_stl <- stlplus(mal_ratio_ts, s.window = "periodic")
    
    
    
    # Assembling data frame for decompositon of all SMC rollout groups
    
    dates <- seq(as.Date("2015-01-01"), as.Date("2018-12-01"), by = "months")
    
    
    cases_stl_ts <- as.data.frame(cases_stl$data[,1:4])
    cases_stl_ts$type <- rep("cases", nrow(cases_stl_ts))
    cases_stl_ts$dates <- dates
    
    ratio_stl_ts <- as.data.frame(ratio_stl$data[,1:4])
    ratio_stl_ts$type <- rep("ratio", nrow(ratio_stl_ts))
    ratio_stl_ts$dates <- dates
    
    
    STL_result_DF_Region_R <- rbind(cases_stl_ts,
                                    ratio_stl_ts)
    STL_result_DF_Region_R$District <- paste("master", R, sep = " ")
    STL_result_DF_Region_R$Region <- R
    STL_result_DF_Region <- rbind(STL_result_DF_Region,
                                  STL_result_DF_Region_R)
    
}


# Setting figure descriptors for ease of plotting with ggplot2

STL_result_DF_Region$alpha <- 1
STL_result_DF_Region$lwd <- 1.5





#########################################################################


# Joining both dataframes for plotting

STL_result_DF_master <- rbind(STL_result_DF, STL_result_DF_Region)
STL_result_DF_master$lty <- factor(STL_result_DF_master$lty, levels = c("solid", "longdash"))


# color_list <- c("#984EA3", "#E41A1C", "#4DAF4A", "#377EB8", "#D95F02", "#E6AB02")


# Figure 4
# Trend component of each district, grouped by Region
# Along with overall trend of aggregate Region

regions1 <- unique(STL_result_DF_master$Region)[1:6]
regions2 <- unique(STL_result_DF_master$Region)[7:13]



plot_cases_trend_for_grid <- ggplot(data = STL_result_DF_master[which(STL_result_DF_master$type == "cases" &
                                                                          STL_result_DF_master$Region %in% regions1),],
                                    aes(x = dates, y = trend)) +
    geom_line(aes(group = District,
                  color = as.factor(Region),
                  alpha = factor(alpha),
                  size = factor(lwd)),
              show.legend = FALSE) +
    scale_alpha_manual(values = unique(STL_result_DF_master$alpha)) +
    scale_size_manual(values = unique(STL_result_DF_master$lwd)) +
    # scale_color_manual(values = color_list) +
    xlab("") + facet_wrap(~Region, ncol = 1, scales = "free") +
    theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                       axis.text.y = element_text(colour = "black"),
                       axis.text.x = element_text(colour = "black"),
                       panel.border = element_rect(colour = "black", fill = NA),
                       panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plot_ratio_trend_for_grid <- ggplot(data = STL_result_DF_master[which(STL_result_DF_master$type == "ratio" &
                                                                          STL_result_DF_master$Region %in% regions1),],
                                    aes(x = dates, y = trend)) +
    geom_line(aes(group = District,
                  color = as.factor(Region),
                  alpha = factor(alpha),
                  size = factor(lwd)),
              show.legend = FALSE) +
    scale_alpha_manual(values = unique(STL_result_DF_master$alpha)) +
    scale_size_manual(values = unique(STL_result_DF_master$lwd)) +
    # scale_color_manual(values = color_list) +
    xlab("") + ylab("") + facet_wrap(~Region, ncol = 1, scales = "free") +
    theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                       axis.text.y = element_text(colour = "black"),
                       axis.text.x = element_text(colour = "black"),
                       panel.border = element_rect(colour = "black", fill = NA),
                       panel.grid.major = element_blank(), panel.grid.minor = element_blank())

grid.arrange(plot_cases_trend_for_grid, plot_ratio_trend_for_grid, ncol = 2)



plot_cases_trend_for_grid <- ggplot(data = STL_result_DF_master[which(STL_result_DF_master$type == "cases" &
                                                                          STL_result_DF_master$Region %in% regions2),],
                                    aes(x = dates, y = trend)) +
    geom_line(aes(group = District,
                  color = as.factor(Region),
                  alpha = factor(alpha),
                  size = factor(lwd)),
              show.legend = FALSE) +
    scale_alpha_manual(values = unique(STL_result_DF_master$alpha)) +
    scale_size_manual(values = unique(STL_result_DF_master$lwd)) +
    # scale_color_manual(values = color_list) +
    xlab("") + facet_wrap(~Region, ncol = 1, scales = "free") +
    theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                       axis.text.y = element_text(colour = "black"),
                       axis.text.x = element_text(colour = "black"),
                       panel.border = element_rect(colour = "black", fill = NA),
                       panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plot_ratio_trend_for_grid <- ggplot(data = STL_result_DF_master[which(STL_result_DF_master$type == "ratio" &
                                                                          STL_result_DF_master$Region %in% regions2),],
                                    aes(x = dates, y = trend)) +
    geom_line(aes(group = District,
                  color = as.factor(Region),
                  alpha = factor(alpha),
                  size = factor(lwd)),
              show.legend = FALSE) +
    scale_alpha_manual(values = unique(STL_result_DF_master$alpha)) +
    scale_size_manual(values = unique(STL_result_DF_master$lwd)) +
    # scale_color_manual(values = color_list) +
    xlab("") + ylab("") + facet_wrap(~Region, ncol = 1, scales = "free") +
    theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                       axis.text.y = element_text(colour = "black"),
                       axis.text.x = element_text(colour = "black"),
                       panel.border = element_rect(colour = "black", fill = NA),
                       panel.grid.major = element_blank(), panel.grid.minor = element_blank())

grid.arrange(plot_cases_trend_for_grid, plot_ratio_trend_for_grid, ncol = 2)










# Trend of Regions together (bottom row of Figure 4)

ggplot(data = STL_result_DF_Region,
       aes(x = dates, y = trend)) +
    geom_line(aes(group = District,
                  color = as.factor(Region),
                  size = factor(lwd))) +
    scale_alpha_manual(values = unique(STL_result_DF_master$alpha)) +
    scale_size_manual(values = unique(STL_result_DF_master$lwd)) +
    # scale_color_manual(values = color_list) +
    xlab("") + facet_wrap(~type, scales = "free") +
    theme_bw() + theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_text(colour = "black"),
                       axis.text.x = element_text(colour = "black", angle = 45, hjust = 1),
                       panel.border = element_rect(colour = "black", fill = NA),
                       panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                       axis.line = element_line(colour = "black"))



## Figures were assembled and cleaned further in Adobe Illustrator





