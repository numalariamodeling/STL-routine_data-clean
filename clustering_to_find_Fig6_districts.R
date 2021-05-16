#######################################################
##  Plotting health facility reporting for figure 1  ##
#######################################################
#
# Description:
#   Finding representative districts for all SMC rollout groups based on
#       lowest Euclidian distance to other time series from districts in same SMC rollout group
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


library("TSclust")
library("maditr")




# Function for standardizing malaria variables

getNormalized <- function(vec)
{
    norm_vec <- (vec - mean(vec))/sd(vec)
    return(norm_vec)
}





# Loading health district dataset
# Creating under-5 population column and fixing date column

cases <- read.csv("~/Box/NU-malaria-team/projects/smc_impact/data/outputs/U5_DS_cases_seasonal_smc_good_rows_MA_imputes.csv", header  = TRUE, strip.white = TRUE, stringsAsFactors = FALSE)
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


## For 2014 SMC group

# First on normalized cases

cases_2014 <- cases[which(cases$SMC.init_date == 2014), c("District", "Date", "mal_cases_norm")]
cases_2014_cast <- dcast(cases_2014, formula = Date~District, value.var = "mal_cases_norm")
cases_2014_cast <- cases_2014_cast[,-1] # removing dates since sorted

tsdist <- diss(cases_2014_cast, "EUCL")
cases_2014_mean_distances <- colMeans(as.matrix(tsdist))




# Then on normalized proportion

prop_2014 <- cases[which(cases$SMC.init_date == 2014), c("District", "Date", "mal_ratio_norm")]
prop_2014_cast <- dcast(prop_2014, formula = Date~District, value.var = "mal_ratio_norm")
prop_2014_cast <- prop_2014_cast[,-1] # removing dates since sorted

tsdist <- diss(prop_2014_cast, "EUCL")
prop_2014_mean_distances <- colMeans(as.matrix(tsdist))

mean_2014_distances <- rbind(cases_2014_mean_distances, prop_2014_mean_distances)
row.names(mean_2014_distances) <- c("cases", "proportion")





## For 2015 SMC group

cases_2015 <- cases[which(cases$SMC.init_date == 2015), c("District", "Date", "mal_cases_norm")]
cases_2015_cast <- dcast(cases_2015, formula = Date~District, value.var = "mal_cases_norm")
cases_2015_cast <- cases_2015_cast[,-1] # removing dates since sorted

tsdist <- diss(cases_2015_cast, "EUCL")
cases_2015_mean_distances <- colMeans(as.matrix(tsdist))



prop_2015 <- cases[which(cases$SMC.init_date == 2015), c("District", "Date", "mal_ratio_norm")]
prop_2015_cast <- dcast(prop_2015, formula = Date~District, value.var = "mal_ratio_norm")
prop_2015_cast <- prop_2015_cast[,-1] # removing dates since sorted

tsdist <- diss(prop_2015_cast, "EUCL")
prop_2015_mean_distances <- colMeans(as.matrix(tsdist))

mean_2015_distances <- rbind(cases_2015_mean_distances, prop_2015_mean_distances)
row.names(mean_2015_distances) <- c("cases", "proportion")





## For 2016 SMC group

cases_2016 <- cases[which(cases$SMC.init_date == 2016), c("District", "Date", "mal_cases_norm")]
cases_2016_cast <- dcast(cases_2016, formula = Date~District,
                         value.var = "mal_cases_norm", fun.aggregate = sum)
cases_2016_cast <- cases_2016_cast[,-1] # removing dates since sorted

tsdist <- diss(cases_2016_cast, "EUCL")
cases_2016_mean_distances <- colMeans(as.matrix(tsdist))



prop_2016 <- cases[which(cases$SMC.init_date == 2016), c("District", "Date", "mal_ratio_norm")]
prop_2016_cast <- dcast(prop_2016, formula = Date~District,
                        value.var = "mal_ratio_norm", fun.aggregate = sum)
prop_2016_cast <- prop_2016_cast[,-1] # removing dates since sorted

tsdist <- diss(prop_2016_cast, "EUCL")
prop_2016_mean_distances <- colMeans(as.matrix(tsdist))

mean_2016_distances <- rbind(cases_2016_mean_distances, prop_2016_mean_distances)
row.names(mean_2016_distances) <- c("cases", "proportion")






## For 2017 SMC group

cases_2017 <- cases[which(cases$SMC.init_date == 2017), c("District", "Date", "mal_cases_norm")]
cases_2017_cast <- dcast(cases_2017, formula = Date~District, value.var = "mal_cases_norm")
cases_2017_cast <- cases_2017_cast[,-1] # removing dates since sorted

tsdist <- diss(cases_2017_cast, "EUCL")
cases_2017_mean_distances <- colMeans(as.matrix(tsdist))



prop_2017 <- cases[which(cases$SMC.init_date == 2017), c("District", "Date", "mal_ratio_norm")]
prop_2017_cast <- dcast(prop_2017, formula = Date~District, value.var = "mal_ratio_norm")
prop_2017_cast <- prop_2017_cast[,-1] # removing dates since sorted

tsdist <- diss(prop_2017_cast, "EUCL")
prop_2017_mean_distances <- colMeans(as.matrix(tsdist))

mean_2017_distances <- rbind(cases_2017_mean_distances, prop_2017_mean_distances)
row.names(mean_2017_distances) <- c("cases", "proportion")






## For 2018 SMC group

cases_2018 <- cases[which(cases$SMC.init_date == 2018), c("District", "Date", "mal_cases_norm")]
cases_2018_cast <- dcast(cases_2018, formula = Date~District, value.var = "mal_cases_norm")
cases_2018_cast <- cases_2018_cast[,-1] # removing dates since sorted

tsdist <- diss(cases_2018_cast, "EUCL")
cases_2018_mean_distances <- colMeans(as.matrix(tsdist))



prop_2018 <- cases[which(cases$SMC.init_date == 2018), c("District", "Date", "mal_ratio_norm")]
prop_2018_cast <- dcast(prop_2018, formula = Date~District, value.var = "mal_ratio_norm")
prop_2018_cast <- prop_2018_cast[,-1] # removing dates since sorted

tsdist <- diss(prop_2018_cast, "EUCL")
prop_2018_mean_distances <- colMeans(as.matrix(tsdist))

mean_2018_distances <- rbind(cases_2018_mean_distances, prop_2018_mean_distances)
row.names(mean_2018_distances) <- c("cases", "proportion")




## For 2019 SMC group

cases_2019 <- cases[which(cases$SMC.init_date == 2019), c("District", "Date", "mal_cases_norm")]
cases_2019_cast <- dcast(cases_2019, formula = Date~District, value.var = "mal_cases_norm")
cases_2019_cast <- cases_2019_cast[,-1] # removing dates since sorted

tsdist <- diss(cases_2019_cast, "EUCL")
cases_2019_mean_distances <- colMeans(as.matrix(tsdist))



prop_2019 <- cases[which(cases$SMC.init_date == 2019), c("District", "Date", "mal_ratio_norm")]
prop_2019_cast <- dcast(prop_2019, formula = Date~District, value.var = "mal_ratio_norm")
prop_2019_cast <- prop_2019_cast[,-1] # removing dates since sorted

tsdist <- diss(prop_2019_cast, "EUCL")
prop_2019_mean_distances <- colMeans(as.matrix(tsdist))

mean_2019_distances <- rbind(cases_2019_mean_distances, prop_2019_mean_distances)
row.names(mean_2019_distances) <- c("cases", "proportion")



most_rep <- rbind(colnames(mean_2015_distances)[apply(mean_2015_distances, 1, which.min)],
                  colnames(mean_2016_distances)[apply(mean_2016_distances, 1, which.min)],
                  colnames(mean_2017_distances)[apply(mean_2017_distances, 1, which.min)],
                  colnames(mean_2018_distances)[apply(mean_2018_distances, 1, which.min)])
row.names(most_rep) <- c("2015", "2016", "2017", "2018")
colnames(most_rep) <- c("cases", "proportion")
most_rep_EUCL <- as.data.frame(most_rep)





