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

pouytenga_Dec_2017_row_tmp <- cases[pouytenga_ind_Nov_2017,]
pouytenga_Dec_2017_row_tmp$month <- 12
pouytenga_Dec_2017_row_tmp$Date <- as.yearmon("Dec 2017")
pouytenga_Dec_2017_row_tmp$mal_cases <- mal_cases_pouytenga_Dec_2017
pouytenga_Dec_2017_row_tmp$all_cases <- all_cases_pouytenga_Dec_2017
pouytenga_Dec_2017_row_tmp$all_nonMal_cases <- all_nonMal_cases_pouytenga_Dec_2017
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
  
  
  all_cases_norm <- getNormalized(cases_dist$all_cases)
  all_cases_norm_ts <- ts(all_cases_norm, start = c(2015, 1), deltat = 1/12)
  all_cases_stl <- stlplus(all_cases_norm_ts, s.window="periodic")
  
  
  all_nonMal_cases_norm <- getNormalized(cases_dist$all_nonMal_cases)
  all_nonMal_cases_norm_ts <- ts(all_nonMal_cases_norm, start = c(2015, 1), deltat = 1/12)
  all_nonMal_cases_stl <- stlplus(all_nonMal_cases_norm_ts, s.window="periodic")
  
  
  mal_ratio_norm <- getNormalized(cases_dist$mal_ratio)
  mal_ratio_norm_ts <- ts(mal_ratio_norm, start = c(2015, 1), deltat = 1/12)
  ratio_stl <- stlplus(mal_ratio_norm_ts, s.window="periodic")

  
  
  # Create date object to create master data-frame for decomposition outputs
  
  dates <- seq(as.Date("2015-01-01"), as.Date("2018-12-01"), by="months")
  
  cases_stl_ts <- as.data.frame(cases_stl$data[,1:4]) 
  cases_stl_ts$type <- rep("cases", nrow(cases_stl_ts))
  cases_stl_ts$lty <- rep("solid", nrow(cases_stl_ts))
  cases_stl_ts$dates <- dates
  
  
  all_cases_stl_ts <- as.data.frame(all_cases_stl$data[,1:4]) 
  all_cases_stl_ts$type <- rep("allout", nrow(all_cases_stl_ts))
  all_cases_stl_ts$lty <- rep("solid", nrow(all_cases_stl_ts))
  all_cases_stl_ts$dates <- dates
  
  
  all_nonMal_cases_stl_ts <- as.data.frame(all_nonMal_cases_stl$data[,1:4]) 
  all_nonMal_cases_stl_ts$type <- rep("all_nonMal", nrow(all_nonMal_cases_stl_ts))
  all_nonMal_cases_stl_ts$lty <- rep("solid", nrow(all_nonMal_cases_stl_ts))
  all_nonMal_cases_stl_ts$dates <- dates
  
  
  ratio_stl_ts <- as.data.frame(ratio_stl$data[,1:4]) 
  ratio_stl_ts$type <- rep("ratio", nrow(ratio_stl_ts))
  ratio_stl_ts$lty <- rep("longdash", nrow(ratio_stl_ts))
  ratio_stl_ts$dates <- dates
  
  
  
  # Create output data-frame from stl outputs of both malaria variables
  
  STL_result_DF_DS_norm <- rbind(cases_stl_ts,
                                 all_cases_stl_ts,
                                 all_nonMal_cases_stl_ts,
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
STL_result_DF_norm <- left_join(STL_result_DF_norm, unique(cases[,c("District", "SMC.init_date")]), by = "District")
STL_result_DF_norm$SMC_rollout_group <- year(STL_result_DF_norm$SMC.init_date)

# Ordering by rollout group and then by district
STL_result_DF_norm <- STL_result_DF_norm[order(STL_result_DF_norm$SMC_rollout_group, STL_result_DF_norm$District),]


# Setting line-types and district factors for plotting
STL_result_DF_norm$lty <- factor(STL_result_DF_norm$lty, levels = c("solid", "longdash"))
STL_result_DF_norm$District <- factor(STL_result_DF_norm$District, levels = sort(unique(cases$District)))

STL_result_DF_norm$Date <- as.yearmon(STL_result_DF_norm$dates)



#######################################################################################



med_2014 <- read.csv("~/OneDrive/Desktop/medfever_region_2014.csv", stringsAsFactors = FALSE)
med_2017 <- read.csv("~/OneDrive/Desktop/medfever_region_2017.csv", stringsAsFactors = FALSE)



med_2014 <- med_2014[,-4]
names(med_2014)[3] <- "medfever_2014"
med_2017 <- med_2017[,-4]
names(med_2017)[3] <- "medfever_2017"

medfever_DHS <- inner_join(med_2014, med_2017, by = c("NOMDEP", "NOMREGION"))
medfever_DHS <- cbind(medfever_DHS[,1:3], medfever_DHS[,4])
names(medfever_DHS) <- c("District", "Region", "2014", "2017")

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

medfever_DHS <- medfever_DHS[order(medfever_DHS$District, medfever_DHS$year),]




######################################################################################




tmp_allout_percent_inc <- sapply(sort(unique(STL_result_DF_norm$District)), function(D) {
  DS_allout_trend <- STL_result_DF_norm[which(STL_result_DF_norm$District == D & STL_result_DF_norm$type == "allout"), "trend"]
  DS_allout_trend_diff <- diff(DS_allout_trend)
  DS_allout_percent_inc <- DS_allout_trend_diff / abs(DS_allout_trend[-length(DS_allout_trend)])

  return(c(NA, DS_allout_percent_inc))
})

cases$rate_allout_u5 <- cases$allout_u5 / cases$U5_pop
cases_tmp <- cbind(cases[,c("Region", "District", "Date", "allout_u5", "U5_pop", "rate_allout_u5")], c(tmp_allout_percent_inc))
names(cases_tmp)[7] <- "allout_percent_inc"
cases_tmp$medfever_estim <- NA

cases_tmp <- left_join(cases_tmp, STL_result_DF_norm[which(STL_result_DF_norm$type == "allout"),
                                                     c("trend", "District", "Date")], by = c("District", "Date"))


cases_tmp[is.na(cases_tmp$allout_percent_inc), "medfever_estim"] <- medfever_DHS[which(medfever_DHS$year == 2014), "medfever_regional"]
cases_tmp[which(cases_tmp$Date == "Dec 2017"), "medfever_estim"] <- medfever_DHS[which(medfever_DHS$year == 2017), "medfever_regional"]


for (DS in unique(cases_tmp$District))
{
  cases_tmp_DS <- cases_tmp[which(cases_tmp$District == DS),]
  for (i in 2:48)
  {
    cases_tmp_DS[i, "medfever_estim"] <- cases_tmp_DS[i - 1, "medfever_estim"] + (cases_tmp_DS[i - 1, "medfever_estim"] * cases_tmp_DS[i, "allout_percent_inc"])
  }
  cases_tmp[which(cases_tmp$District == DS), "medfever_estim"] <- cases_tmp_DS[, "medfever_estim"]
}



#########################################################################################################

## Some tests for Banfora


cases_tmp_DS$inds <- 7:54
fit <- lm(trend ~ bs(ind, knots = seq(13,54,6)), data = cases_tmp_DS)

plot(1:54, c(rep(NA,6),cases_tmp_DS$trend),col="grey")
points(c(1:54), predict(fit, newdata = list(ind = c(1:54))),col="darkgreen",lwd=2,type="l")



cases_tmp_DS$inds <- 1:48
fit <- lm(trend ~ bs(inds, knots = seq(7,48,6)), data = cases_tmp_DS)

plot(1:48, cases_tmp_DS$trend,col="grey")

points(c(1:48), predict(fit, newdata = list(trend = c(-0.7071068, rep(NA, 34), 0.7071068, rep(NA, 12)),
                                            ind = c(1:48))),col="darkgreen",lwd=2,type="l")





tmp_medfever_banfora <- data.frame(inds = c(1,36), medfever = c(0.6583144, 0.8433474))
fit <- lm(medfever ~ inds, data = tmp_medfever_banfora)

plot(tmp_medfever_banfora$inds, tmp_medfever_banfora$medfever, col="grey", ylim = c(.6, .9))
points(c(1:48), predict(fit, newdata = list(inds = c(1:48))),col="darkgreen",lwd=2,type="l")





tmp_medfever_banfora <- data.frame(inds = c(1,36), medfever = c(0.6583144, 0.8433474))
fit <- lm(trend ~ bs(inds, knots = seq(7,48,6)), data = cases_tmp_DS)
plot(tmp_medfever_banfora$inds, tmp_medfever_banfora$medfever, col="grey", ylim = c(.6, .9))
points(c(1:48), ((predict(fit, newdata = list(inds = c(1:48)))*sd(tmp_medfever_banfora$medfever)) + mean(tmp_medfever_banfora$medfever)),col="darkgreen",lwd=2,type="l")

C <- 0.025
points(c(1:48), ((predict(fit, newdata = list(inds = c(1:48)))*sd(tmp_medfever_banfora$medfever)) +
                   mean(tmp_medfever_banfora$medfever) + C),col="red",lwd=2,type="l")






#########################################################################################################





cases_tmp_2 <- cases[,c("Region", "District", "Date", "allout_u5", "U5_pop", "conf_rdt_mic_u5", "weighted_rep_rate")]
cases_tmp_2$medfever <- NA
cases_tmp_2$medfever_estim <- NA

cases_tmp_2 <- left_join(cases_tmp_2, STL_result_DF_norm[which(STL_result_DF_norm$type == "allout"),
                                                         c("trend", "District", "Date")], by = c("District", "Date"))


cases_tmp_2[which(cases_tmp_2$Date == "Jan 2015"), "medfever"] <- medfever_DHS[which(medfever_DHS$year == 2014), "medfever_regional"]
cases_tmp_2[which(cases_tmp_2$Date == "Dec 2017"), "medfever"] <- medfever_DHS[which(medfever_DHS$year == 2017), "medfever_regional"]




for (DS in unique(cases_tmp_2$District))
{
  cases_tmp_DS <- cases_tmp_2[which(cases_tmp_2$District == DS),]
  
  cases_tmp_DS$inds <- 1:48
  fit <- lm(trend ~ bs(inds, knots = seq(7,48,6)), data = cases_tmp_DS)
  
  estim_medfever <- (predict(fit, newdata = list(inds = c(1:48))) * sd(cases_tmp_DS$medfever, na.rm = T)) + mean(cases_tmp_DS$medfever, na.rm = T)

  
  y_range <- c(min(range(estim_medfever)[1], cases_tmp_DS[which(cases_tmp_DS$Date == "Jan 2015"), "medfever"]) - .05,
               max(range(estim_medfever)[2], cases_tmp_DS[which(cases_tmp_DS$Date == "Dec 2017"), "medfever"]) + .05)
  
  
  pdf(paste("~/OneDrive/Desktop/estim_medfever/", DS, ".pdf", sep = ""))
  
  plot(cases_tmp_DS$inds, cases_tmp_DS$medfever, col = "black", ylim = y_range, pch = 17, main = paste(DS, unique(cases_tmp_DS$Region), sep = ", "))
  points(c(1:48), (cases_tmp_DS$trend * sd(cases_tmp_DS$medfever, na.rm = T)) + mean(cases_tmp_DS$medfever, na.rm = T), col = "red")
  points(c(1:48), estim_medfever,
         col = "darkgreen", lwd = 2, type = "l")
  
  dev.off()
  
  
  
  cases_tmp_2[which(cases_tmp_2$District == DS), "medfever_estim"] <- estim_medfever
}









ggplot() + geom_line(data = cases_tmp_DS, aes(x = Date, y = conf_rdt_mic_u5 / (estim_medfever * weighted_rep_rate)), col = "blue") +
  geom_line(data = cases_tmp_DS, aes(x = Date, y = conf_rdt_mic_u5), col = "red")


ggplot() + geom_line(data = cases_tmp_DS, aes(x = Date, y = conf_rdt_mic_u5 / allout_u5), col = "blue")



tmp_1 <- cases_tmp_DS$conf_rdt_mic_u5 / (estim_medfever * cases_tmp_DS$weighted_rep_rate)
tmp_1_norm <- getNormalized(tmp_1)
tmp_1_norm_stl <- stlplus(ts(tmp_1_norm, start = c(2015, 1), deltat = 1/12), s.window="periodic")
tmp_1_stl_ts <- as.data.frame(tmp_1_norm_stl$data[,1:4]) 



tmp_2 <- cases_tmp_DS$conf_rdt_mic_u5 / cases_tmp_DS$allout_u5
tmp_2_norm <- getNormalized(tmp_2)
tmp_2_norm_stl <- stlplus(ts(tmp_2_norm, start = c(2015, 1), deltat = 1/12), s.window="periodic")
tmp_2_stl_ts <- as.data.frame(tmp_2_norm_stl$data[,1:4]) 





ggplot() + geom_line(data = tmp_1_stl_ts, aes(x = 1:48, y = trend), col = "blue") +
  geom_line(data = tmp_2_stl_ts, aes(x = 1:48, y = trend), col = "red")






