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

cases <- read.csv("~/Box/NU-malaria-team/projects/smc_impact/data/outputs/U5_DS_cases_seasonal_smc_good_rows_MA_imputes_testing.csv", header  = TRUE, strip.white = TRUE, stringsAsFactors = FALSE)
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

# cases$mal_cases <- cases$conf_rdt_mic_u5 / (cases$U5_pop/1000)
# cases$all_cases <- cases$allout_u5 / (cases$U5_pop/1000)
# cases$mal_ratio <- cases$conf_rdt_mic_u5 / cases$allout_u5




med_2014 <- read.csv("~/Box/NU-malaria-team/projects/hbhi_burkina/DS DHS estimates/medfever/med_BF14DS.csv", stringsAsFactors = FALSE)
med_2017 <- read.csv("~/Box/NU-malaria-team/projects/hbhi_burkina/DS DHS estimates/medfever/med_BF17DS.csv", stringsAsFactors = FALSE)



med_2014 <- med_2014[,-c(1,4)]
names(med_2014)[2] <- "medfever_2014"
med_2017 <- med_2017[,-c(1,4)]
names(med_2017)[2] <- "medfever_2017"

medfever_DHS <- inner_join(med_2014, med_2017, by = "NOMDEP")
medfever_DHS <- cbind(medfever_DHS[,1:2], medfever_DHS[,2], medfever_DHS[,3], medfever_DHS[,3])
names(medfever_DHS) <- c("District", "2015", "2016", "2017", "2018")

medfever_DHS <- melt(medfever_DHS, id.vars = "District", variable.name = "year", value.name = "medfever")



district_map <- cbind(sort(unique(cases$District)),
                      sort(unique(medfever_DHS$District))[c(1:40,42,41,43:70)])

for (i in 1:nrow(district_map))
{
    DS <- district_map[i, 2]
    new_name_DS <- district_map[i, 1]
    
    medfever_DHS[which(medfever_DHS$District == DS), "District"] <- new_name_DS
}

medfever_DHS$year <- as.numeric(as.character(medfever_DHS$year))

medfever_DHS <- left_join(medfever_DHS, unique(cases[,c("Region", "District", "year", "U5_pop")]), by = c("District", "year"))


medfever_region <- ddply(medfever_DHS[!is.na(medfever_DHS$medfever),], c(.(Region), .(year)),
                         summarise, medfever = sum(medfever * U5_pop)/ sum(U5_pop))


medfever_DHS_tmp_1 <- medfever_DHS[!is.na(medfever_DHS$medfever),]
medfever_DHS_tmp_2 <- medfever_DHS[is.na(medfever_DHS$medfever),]

medfever_DHS_tmp_2 <- left_join(medfever_DHS_tmp_2[,-3], medfever_region, by = c("year", "Region"))
medfever_DHS_tmp_2 <- medfever_DHS_tmp_2[, c(1:2, 5, 3:4)]

medfever_DHS <- rbind(medfever_DHS_tmp_1, medfever_DHS_tmp_2)


cases <- left_join(cases, medfever_DHS[,1:3], by = c("District", "year"))




cases$cases_trtseeking_adj <- cases$conf_rdt_mic_u5 + (cases$conf_rdt_mic_u5 * (1 - cases$medfever))




cases$rep_rate <- cases$reporting_HF / cases$total_HF

# cases$cases_rep_adj <- cases$conf_rdt_mic_u5 + (cases$conf_rdt_mic_u5 * (1 - cases$rep_rate))
cases$cases_rep_adj <- cases$conf_rdt_mic_u5 / cases$rep_rate


cases$cases_both_adj <- cases$cases_rep_adj + (cases$cases_rep_adj * (1 - cases$medfever))
















#########################################################################3




## Fixing pouytenga district for NA value in Nov 2016
# pouytenga has NA cases in Nov 2016
# Imputing with mean of previous and subsequent months

pouytenga_ind_Nov_2017 <- which(cases$Date == "Nov 2017" & cases$District == "pouytenga")
pouytenga_ind_Jan_2018 <- which(cases$Date == "Jan 2018" & cases$District == "pouytenga")

trt_adj_pouytenga_Dec_2017 <- mean(c(cases$cases_trtseeking_adj[pouytenga_ind_Nov_2017],
                                     cases$cases_trtseeking_adj[pouytenga_ind_Jan_2018]))
rep_adj_pouytenga_Dec_2017 <- mean(c(cases$cases_rep_adj[pouytenga_ind_Nov_2017],
                                     cases$cases_rep_adj[pouytenga_ind_Jan_2018]))
both_adj_pouytenga_Dec_2017 <- mean(c(cases$cases_both_adj[pouytenga_ind_Nov_2017],
                                      cases$cases_both_adj[pouytenga_ind_Jan_2018]))

pouytenga_Dec_2017_row_tmp <- cases[pouytenga_ind_Nov_2017,]
pouytenga_Dec_2017_row_tmp$month <- 12
pouytenga_Dec_2017_row_tmp$Date <- as.yearmon("Dec 2017")
pouytenga_Dec_2017_row_tmp$cases_trtseeking_adj <- trt_adj_pouytenga_Dec_2017
pouytenga_Dec_2017_row_tmp$cases_rep_adj <- rep_adj_pouytenga_Dec_2017
pouytenga_Dec_2017_row_tmp$cases_both_adj <- both_adj_pouytenga_Dec_2017


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
    trt_adj_norm <- getNormalized(cases_dist$cases_trtseeking_adj)
    trt_adj_norm_ts <- ts(trt_adj_norm, start = c(2015, 1), deltat = 1/12)
    trt_adj_stl <- stlplus(trt_adj_norm_ts, s.window="periodic")
    
    
    rep_adj_norm <- getNormalized(cases_dist$cases_rep_adj)
    rep_adj_norm_ts <- ts(rep_adj_norm, start = c(2015, 1), deltat = 1/12)
    rep_adj_stl <- stlplus(rep_adj_norm_ts, s.window="periodic")
    
    
    
    # Do the same for the malaria share variable
    both_adj_norm <- getNormalized(cases_dist$cases_both_adj)
    both_adj_norm_ts <- ts(both_adj_norm, start = c(2015, 1), deltat = 1/12)
    both_adj_stl <- stlplus(both_adj_norm_ts, s.window="periodic")
    
    
    
    
    # Create date object to create master data-frame for decomposition outputs
    
    dates <- seq(as.Date("2015-01-01"), as.Date("2018-12-01"), by="months")
    
    trt_adj_stl_ts <- as.data.frame(trt_adj_stl$data[,1:4]) 
    trt_adj_stl_ts$type <- rep("trt_adj", nrow(trt_adj_stl_ts))
    trt_adj_stl_ts$lty <- rep("solid", nrow(trt_adj_stl_ts))
    trt_adj_stl_ts$dates <- dates
    
    
    rep_adj_stl_ts <- as.data.frame(rep_adj_stl$data[,1:4]) 
    rep_adj_stl_ts$type <- rep("rep_adj", nrow(rep_adj_stl_ts))
    rep_adj_stl_ts$lty <- rep("longdash", nrow(rep_adj_stl_ts))
    rep_adj_stl_ts$dates <- dates
    
    
    both_adj_stl_ts <- as.data.frame(both_adj_stl$data[,1:4]) 
    both_adj_stl_ts$type <- rep("both_adj", nrow(both_adj_stl_ts))
    both_adj_stl_ts$lty <- rep("dotted", nrow(both_adj_stl_ts))
    both_adj_stl_ts$dates <- dates
    
    
    
    # Create output data-frame from stl outputs of both malaria variables
    
    STL_result_DF_DS_norm <- rbind(trt_adj_stl_ts,
                                   rep_adj_stl_ts,
                                   both_adj_stl_ts)
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
STL_result_DF_norm <- left_join(STL_result_DF_norm, unique(cases[,c("District", "SMC.init_date")]), by = "District")
STL_result_DF_norm$SMC_rollout_group <- year(STL_result_DF_norm$SMC.init_date)

# Ordering by rollout group and then by district
STL_result_DF_norm <- STL_result_DF_norm[order(STL_result_DF_norm$SMC_rollout_group, STL_result_DF_norm$District),]


# Setting line-types and district factors for plotting
STL_result_DF_norm$lty <- factor(STL_result_DF_norm$lty, levels = c("solid", "longdash", "dotted"))
STL_result_DF_norm$District <- factor(STL_result_DF_norm$District, levels = sort(unique(cases$District)))



#######################################################################################


# Create plots for original data and the 3 extracted components from the decomposition
# We plot both malaria variables on the same plot for ease of visualization
# Make 4x4 grid (4 districts) spanning all 70 districts



for (i in seq(1,70,16))
{
    # Select 4 district for plotting
    # Dont want alphabetical but my SMC rollout group THEN alphabetical from there
    
    ## DS_list <- sort(unique(STL_result_DF_norm$District))[i:(i+3)]
    DS_list <- unique(STL_result_DF_norm$District)[i:(i+15)]
    
    plotting_DF <- STL_result_DF_norm[which(STL_result_DF_norm$District %in% DS_list),]
    plotting_DF$facet_name <- paste(plotting_DF$District, plotting_DF$SMC_rollout_group, sep = " ")
    plotting_DF$facet_name <- factor(plotting_DF$facet_name, levels = unique(plotting_DF$facet_name))
    
    # 
    # plot_data <- ggplot(plotting_DF, aes(x = dates, y = raw)) +
    #     geom_vline(data = plotting_DF[which(plotting_DF$SMC.received == 1),],
    #                aes(xintercept = dates), col = "black", size = .3, alpha = .5, linetype = "solid", show.legend = FALSE) +
    #     geom_line(aes(group = type, col = District, linetype = lty), show.legend = FALSE) + scale_color_discrete("") +
    #     xlab("") + facet_wrap(~facet_name, ncol = 4) + ylab("data") +
    #     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    #           panel.background = element_blank(), axis.line = element_line(colour = "black"),
    #           panel.spacing.x = unit(6, "mm"))
    # 
    # 
    # plot_seasonal <- ggplot(plotting_DF, aes(x = dates, y = seasonal)) +
    #     geom_vline(data = plotting_DF[which(plotting_DF$SMC.received == 1),],
    #                aes(xintercept = dates), col = "black", size = .3, alpha = .5, linetype = "solid", show.legend = FALSE) +
    #     geom_line(aes(group = type, col = District, linetype = lty), show.legend = FALSE) + scale_color_discrete("") +
    #     xlab("") + facet_wrap(~facet_name, ncol = 4) +
    #     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    #           panel.background = element_blank(), axis.line = element_line(colour = "black"),
    #           panel.spacing.x = unit(6, "mm"))
    # 
    
    plot_trend <- ggplot(plotting_DF, aes(x = dates, y = trend)) +
        geom_vline(data = plotting_DF[which(plotting_DF$SMC.received == 1),],
                   aes(xintercept = dates), col = "black", size = .3, alpha = .5, linetype = "solid", show.legend = FALSE) +
        geom_line(aes(group = type, col = District, linetype = lty), show.legend = FALSE) + scale_color_discrete("") +
        scale_linetype_identity() +
        xlab("") + facet_wrap(~facet_name, ncol = 4) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"),
              panel.spacing.x = unit(6, "mm"))
    
    # 
    # plot_rem <- ggplot(plotting_DF, aes(x = dates, y = remainder)) +
    #     geom_vline(data = plotting_DF[which(plotting_DF$SMC.received == 1),],
    #                aes(xintercept = dates), col = "black", size = .3, alpha = .5, linetype = "solid", show.legend = FALSE) +
    #     geom_line(aes(group = type, col = District, linetype = lty), show.legend = FALSE) + scale_color_discrete("") +
    #     xlab("") + facet_wrap(~facet_name, ncol = 4) +
    #     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    #           panel.background = element_blank(), axis.line = element_line(colour = "black"),
    #           panel.spacing.x = unit(6, "mm"))
    # 
    
    
    
    # Saving SI plots
    pdf(paste("~/OneDrive/Desktop/SI_figures_cleaned_for_Bea/adjusted/SI_decomp_", i, ".pdf", sep = ""))
    print(plot_trend)
    # grid.arrange(plot_data, plot_seasonal, plot_trend, plot_rem, ncol = 1)
    dev.off()
    
}


ggplot(STL_result_DF_norm[which(STL_result_DF_norm$District == "baskuy"),],
       aes(x = dates, y = trend)) +
    geom_vline(data = plotting_DF[which(plotting_DF$SMC.received == 1),],
               aes(xintercept = dates), col = "black", size = .3, alpha = .5, linetype = "solid", show.legend = FALSE) +
    geom_line(aes(group = type, color = factor(type, levels = c("trt_adj", "rep_adj", "both_adj")),
                  linetype = factor(type, levels = c("trt_adj", "rep_adj", "both_adj")))) +
    scale_color_discrete("") + scale_linetype("") + xlab("") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          panel.spacing.x = unit(6, "mm"))









