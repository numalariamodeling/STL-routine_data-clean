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
# cases$all_cases <- cases$allout_u5 / (cases$U5_pop/1000)
# cases$all_nonMal_cases <- (cases$allout_u5 - cases$conf_rdt_mic_u5) / (cases$U5_pop/1000)
# cases$mal_ratio <- cases$conf_rdt_mic_u5 / cases$allout_u5




# med_2014 <- read.csv("~/Box/NU-malaria-team/projects/hbhi_burkina/DS DHS estimates/medfever/med_BF14DS.csv", stringsAsFactors = FALSE)
# med_2017 <- read.csv("~/Box/NU-malaria-team/projects/hbhi_burkina/DS DHS estimates/medfever/med_BF17DS.csv", stringsAsFactors = FALSE)
# 
# 
# 
# med_2014 <- med_2014[,-c(1,4)]
# names(med_2014)[2] <- "medfever_2014"
# med_2017 <- med_2017[,-c(1,4)]
# names(med_2017)[2] <- "medfever_2017"
# 
# medfever_DHS <- inner_join(med_2014, med_2017, by = "NOMDEP")
# medfever_DHS <- cbind(medfever_DHS[,1:2], medfever_DHS[,2], medfever_DHS[,3], medfever_DHS[,3])
# names(medfever_DHS) <- c("District", "2015", "2016", "2017", "2018")
# 
# medfever_DHS <- melt(medfever_DHS, id.vars = "District", variable.name = "year", value.name = "medfever")
# 
# 
# 
# district_map <- cbind(sort(unique(cases$District)),
#                       sort(unique(medfever_DHS$District))[c(1:40,42,41,43:70)])
# 
# for (i in 1:nrow(district_map))
# {
#     DS <- district_map[i, 2]
#     new_name_DS <- district_map[i, 1]
#     
#     medfever_DHS[which(medfever_DHS$District == DS), "District"] <- new_name_DS
# }
# 
# medfever_DHS$year <- as.numeric(as.character(medfever_DHS$year))
# 
# medfever_DHS <- left_join(medfever_DHS, unique(cases[,c("Region", "District", "year", "U5_pop")]), by = c("District", "year"))
# 
# 
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
# 
# names(medfever_region)[3] <- "medfever_regional"
# medfever_DHS <- left_join(medfever_DHS, medfever_region, by = c("year", "Region"))
# 
# 
# cases <- left_join(cases, medfever_DHS[,c(1:3, 6)], by = c("District", "year"))
# 
# for (D in unique(cases$District))
# {
#     cases[which(cases$District == D & cases$Date %in% as.yearmon(seq(as.Date("2016-05-01"), as.Date("2016-12-01"), by="months"))), "medfever"] <- cases[which(cases$District == D & cases$Date == "Jan 2017"), "medfever"]
#     cases[which(cases$District == D & cases$Date %in% as.yearmon(seq(as.Date("2016-05-01"), as.Date("2016-12-01"), by="months"))), "medfever_regional"] <- cases[which(cases$District == D & cases$Date == "Jan 2017"), "medfever_regional"]
# }





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

cases$medfever_regional_step1 <- cases$medfever_regional
for (D in unique(cases$District))
{
    cases[which(cases$District == D & cases$Date %in% as.yearmon(seq(as.Date("2016-05-01"), as.Date("2016-12-01"), by="months"))), "medfever_regional_step1"] <- cases[which(cases$District == D & cases$Date == "Jan 2017"), "medfever_regional_step1"]
}


cases$cases_step1_trtseeking_adj <- (cases$mal_cases / cases$medfever_regional_step1) / cases$weighted_rep_rate

cases$cases_step2_trtseeking_adj <- (cases$mal_cases / cases$medfever_regional) / cases$weighted_rep_rate





# cases$rep_rate <- cases$reporting_HF / cases$total_HF
# 
# cases$cases_rep_adj <- cases$conf_rdt_mic_u5  + (cases$conf_rdt_mic_u5 * (1 - cases$rep_rate))
# cases$cases_rep_adj <- cases$mal_cases  / cases$rep_rate
# 
# cases$cases_rep_weighted_adj <- cases$mal_cases / cases$weighted_rep_rate
# 
# 
# 
# cases$cases_both_adj <- cases$cases_rep_weighted_adj / cases$fitted_regional_medfever




# cases$all_cases_rep_weighted_adj <- cases$all_cases / cases$weighted_rep_rate
# cases$all_nonMal_cases_rep_weighted_adj <- cases$all_nonMal_cases / cases$weighted_rep_rate





################################################################################




## Fixing pouytenga district for NA value in Nov 2016
# pouytenga has NA cases in Nov 2016
# Imputing with mean of previous and subsequent months

pouytenga_ind_Nov_2017 <- which(cases$Date == "Nov 2017" & cases$District == "pouytenga")
pouytenga_ind_Jan_2018 <- which(cases$Date == "Jan 2018" & cases$District == "pouytenga")


mal_cases_pouytenga_Dec_2017 <- mean(c(cases$mal_cases[pouytenga_ind_Nov_2017],
                                       cases$mal_cases[pouytenga_ind_Jan_2018]))
# all_cases_pouytenga_Dec_2017 <- mean(c(cases$all_cases[pouytenga_ind_Nov_2017],
#                                        cases$all_cases[pouytenga_ind_Jan_2018]))
# all_nonMal_cases_pouytenga_Dec_2017 <- mean(c(cases$all_nonMal_cases[pouytenga_ind_Nov_2017],
#                                               cases$all_nonMal_cases[pouytenga_ind_Jan_2018]))
# mal_ratio_pouytenga_Dec_2017 <- mean(c(cases$mal_ratio[pouytenga_ind_Nov_2017],
#                                        cases$mal_ratio[pouytenga_ind_Jan_2018]))

linear_trt_adj_pouytenga_Dec_2017 <- mean(c(cases$cases_linear_trtseeking_adj[pouytenga_ind_Nov_2017],
                                            cases$cases_linear_trtseeking_adj[pouytenga_ind_Jan_2018]))
step1_trt_adj_pouytenga_Dec_2017 <- mean(c(cases$cases_step1_trtseeking_adj[pouytenga_ind_Nov_2017],
                                           cases$cases_step1_trtseeking_adj[pouytenga_ind_Jan_2018]))
step2_trt_adj_pouytenga_Dec_2017 <- mean(c(cases$cases_step2_trtseeking_adj[pouytenga_ind_Nov_2017],
                                           cases$cases_step2_trtseeking_adj[pouytenga_ind_Jan_2018]))
# rep_adj_pouytenga_Dec_2017 <- mean(c(cases$cases_rep_adj[pouytenga_ind_Nov_2017],
#                                      cases$cases_rep_adj[pouytenga_ind_Jan_2018]))
# weight_rep_adj_pouytenga_Dec_2017 <- mean(c(cases$cases_rep_weighted_adj[pouytenga_ind_Nov_2017],
#                                             cases$cases_rep_weighted_adj[pouytenga_ind_Jan_2018]))
# both_adj_pouytenga_Dec_2017 <- mean(c(cases$cases_both_adj[pouytenga_ind_Nov_2017],
#                                       cases$cases_both_adj[pouytenga_ind_Jan_2018]))




pouytenga_Dec_2017_row_tmp <- cases[pouytenga_ind_Nov_2017,]
pouytenga_Dec_2017_row_tmp$month <- 12
pouytenga_Dec_2017_row_tmp$Date <- as.yearmon("Dec 2017")
pouytenga_Dec_2017_row_tmp$mal_cases <- mal_cases_pouytenga_Dec_2017
pouytenga_Dec_2017_row_tmp$cases_linear_trtseeking_adj <- linear_trt_adj_pouytenga_Dec_2017
pouytenga_Dec_2017_row_tmp$cases_step1_trtseeking_adj <- step1_trt_adj_pouytenga_Dec_2017
pouytenga_Dec_2017_row_tmp$cases_step2_trtseeking_adj <- step2_trt_adj_pouytenga_Dec_2017
# pouytenga_Dec_2017_row_tmp$cases_both_adj <- both_adj_pouytenga_Dec_2017


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
    
    
    linear_trt_adj_norm <- getNormalized(cases_dist$cases_linear_trtseeking_adj)
    linear_trt_adj_norm_ts <- ts(linear_trt_adj_norm, start = c(2015, 1), deltat = 1/12)
    linear_trt_adj_stl <- stlplus(linear_trt_adj_norm_ts, s.window="periodic")
    

    step1_trt_adj_norm <- getNormalized(cases_dist$cases_step1_trtseeking_adj)
    step1_trt_adj_norm_ts <- ts(step1_trt_adj_norm, start = c(2015, 1), deltat = 1/12)
    step1_trt_adj_stl <- stlplus(step1_trt_adj_norm_ts, s.window="periodic")
    
    
    step2_trt_adj_norm <- getNormalized(cases_dist$cases_step2_trtseeking_adj)
    step2_trt_adj_norm_ts <- ts(step2_trt_adj_norm, start = c(2015, 1), deltat = 1/12)
    step2_trt_adj_stl <- stlplus(step2_trt_adj_norm_ts, s.window="periodic")

    
    
    # Create date object to create master data-frame for decomposition outputs
    
    dates <- seq(as.Date("2015-01-01"), as.Date("2018-12-01"), by="months")
    
    cases_stl_ts <- as.data.frame(cases_stl$data[,1:4]) 
    cases_stl_ts$type <- rep("cases", nrow(cases_stl_ts))
    cases_stl_ts$dates <- dates
    

    linear_trt_adj_stl_ts <- as.data.frame(linear_trt_adj_stl$data[,1:4])
    linear_trt_adj_stl_ts$type <- rep("linear_trt_adj", nrow(linear_trt_adj_stl_ts))
    linear_trt_adj_stl_ts$dates <- dates
    

    step1_trt_adj_stl_ts <- as.data.frame(step1_trt_adj_stl$data[,1:4])
    step1_trt_adj_stl_ts$type <- rep("step1_trt_adj", nrow(step1_trt_adj_stl_ts))
    step1_trt_adj_stl_ts$dates <- dates
    
    
    step2_trt_adj_stl_ts <- as.data.frame(step2_trt_adj_stl$data[,1:4])
    step2_trt_adj_stl_ts$type <- rep("step2_trt_adj", nrow(step2_trt_adj_stl_ts))
    step2_trt_adj_stl_ts$dates <- dates
    
    
    
    
    # Create output data-frame from stl outputs of both malaria variables
    
    STL_result_DF_DS_norm <- rbind(cases_stl_ts,
                                   linear_trt_adj_stl_ts,
                                   step1_trt_adj_stl_ts,
                                   step2_trt_adj_stl_ts)
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
STL_result_DF_norm <- left_join(STL_result_DF_norm, unique(cases[,c("District", "Region", "SMC.init_date")]), by = "District")
STL_result_DF_norm$SMC_rollout_group <- year(STL_result_DF_norm$SMC.init_date)

# Ordering by rollout group and then by district
STL_result_DF_norm <- STL_result_DF_norm[order(STL_result_DF_norm$Region, STL_result_DF_norm$District),]


# Setting line-types and district factors for plotting
STL_result_DF_norm$District <- factor(STL_result_DF_norm$District, levels = sort(unique(cases$District)))



#######################################################################################


# Create plots for original data and the 3 extracted components from the decomposition
# We plot both malaria variables on the same plot for ease of visualization
# Make 4x4 grid (4 districts) spanning all 70 districts



for (i in seq(1,70,8))
{
    # Select 4 district for plotting
    # Dont want alphabetical but my SMC rollout group THEN alphabetical from there
    
    
    DS_list <- unique(STL_result_DF_norm$District)[i:(i+7)]
    
    plotting_DF <- STL_result_DF_norm[which(STL_result_DF_norm$District %in% DS_list),]
    # plotting_DF$facet_name <- paste(plotting_DF$District, plotting_DF$SMC_rollout_group, sep = " ")
    # plotting_DF$facet_name <- factor(plotting_DF$facet_name, levels = unique(plotting_DF$facet_name))
    
    
    plot_trends <- ggplot(plotting_DF[which(plotting_DF$type %in% c("cases", "linear_trt_adj", "step1_trt_adj", "step2_trt_adj")),],
                           aes(x = dates, y = trend)) +
        geom_line(aes(group = type, linetype = "solid",
                      color = factor(type, levels = c("cases", "linear_trt_adj", "step1_trt_adj", "step2_trt_adj")))) +
        scale_color_manual("", values = c("#913058", "#F6851F", "#00A08A", "#8971B3")) + scale_linetype_identity("") + 
        xlab("") + facet_wrap(~District, ncol = 4) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"),
              panel.spacing.x = unit(6, "mm"))

    
    pdf(paste("~/OneDrive/Desktop/SI_figures_cleaned_for_Bea/trt_seeking_trends/trt_seeking_trends_", i, ".pdf", sep = ""),
        width = 11, height = 7)
    print(plot_trends)
    dev.off()
    
}





# plot1 <- ggplot(STL_result_DF_norm[which(STL_result_DF_norm$District %in% c("tougouri", "koudougou", "lena", "gaoua") &
#                                              STL_result_DF_norm$type %in% c("cases", "trt_adj", "rep_adj", "rep_weighted_adj", "both_adj")),],
#                 aes(x = dates, y = trend)) +
#     geom_line(aes(group = type, linetype = "solid",
#                   color = factor(type, levels = c("cases", "trt_adj", "rep_adj", "rep_weighted_adj", "both_adj")))) +
#     scale_color_manual("", values = c("#913058", "#F6851F", "#00A08A", "#8971B3", "#D61B5A")) +
#     scale_linetype_identity("") +
#     xlab("") + facet_wrap(~District, ncol = 4) +
#     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#           panel.background = element_blank(), axis.line = element_line(colour = "black"),
#           panel.spacing.x = unit(6, "mm"))




