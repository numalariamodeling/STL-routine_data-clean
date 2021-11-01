######################################################
##  Imputing observations at health facility level  ##
######################################################
#
# Description:
#   Imputing health facility observations for confirmed rdt cases and all-cause outpatient cases
#       with Kilman smoothing. Only impute from HFs missing less than or equal to 5 obs and no more than 2
#       in a row
#
#
#  Sebastian Rodriguez (sebastian@rodriguez.cr)
#  Last edited Mar 09, 2021
#



rm(list = ls(all = TRUE))

require("plyr")
require("dplyr")
require("zoo")

library("imputeTS")
library("maditr")



# Loading health facility dataset

HF_cases <- read.csv("~/Box/NU-malaria-team/projects/smc_impact/data/outputs/U5_HF_cases_smc_coords.csv")
HF_cases$Date <- as.Date(as.yearmon(HF_cases$Date))
HF_cases <- HF_cases[order(HF_cases$UID, HF_cases$Date),]


## Fixing issue with "boulmiougou boulmiougou csps maiga"
tmp_HF_cases_csps_maiga <- HF_cases[which(HF_cases$UID == "boulmiougou boulmiougou csps maiga"),]
tmp_HF_cases_csps_maiga_unique <- tmp_HF_cases_csps_maiga[duplicated(tmp_HF_cases_csps_maiga[,c("Date","precip_era5")]),]

HF_cases <- HF_cases[-which(HF_cases$UID == "boulmiougou boulmiougou csps maiga"),]
HF_cases <- rbind(HF_cases, tmp_HF_cases_csps_maiga_unique)
HF_cases <- HF_cases[order(HF_cases$UID, HF_cases$Date),]


## Fixing some ridiculous HF outliers by changing to NA
HF_cases[which(HF_cases$UID == "zorgho zorgho csps tuire" & HF_cases$Date == "2018-09-01"), "allout_u5"] <- NA
HF_cases[which(HF_cases$UID == "nanoro kindi csps zerkoum" & HF_cases$Date == "2016-09-01"), "allout_u5"] <- NA
HF_cases[which(HF_cases$UID == "tenkodogo tenkodogo csps urbain ii" & HF_cases$Date == "2018-11-01"), "allout_u5"] <- NA
HF_cases[which(HF_cases$UID == "toma toma cma toma" & HF_cases$Date == "2018-08-01"), "allout_u5"] <- NA
HF_cases[which(HF_cases$UID == "sapone ipelce csps ipelce" & HF_cases$Date == "2016-09-01"), "allout_u5"] <- NA


HF_cases[which(HF_cases$UID == "kampti kampti csps latara" & HF_cases$Date == "2017-03-01"), "conf_rdt_ov5"] <- NA
HF_cases[which(HF_cases$UID == "nouna djibasso csps ba" & HF_cases$Date == "2016-05-01"), "conf_rdt_ov5"] <- NA





## Making NA Oct, Nov, and Dec 2016 for districts in which the prop of conf/allout > 1 for those months
Oct_districts <- c("ouargaye", "boussouma", "leo", "reo", "lena", "toma", "tougan", "garango",
                   "pouytenga", "nanoro", "sabou", "gourcy", "dano", "diebougou", "dedougou", "gorom",
                   "nouna", "sindou", "bitou", "kongoussi", "tougouri", "manga", "kombissiri", "dori",
                   "po", "sapone", "diapaga", "fada", "seguenega", "thiou", "yako", "ziniare", "zorgho")
Nov_districts <- c("ouargaye", "sapouy", "banfora", "mangodara", "koupela", "tenkodogo", "kongoussi", "tougouri",
                   "leo", "nanoro", "reo", "manga", "kombissiri", "po", "bogande", "fada", "ziniare", "sindou")
Sep_districts <- c("fada", "gorom")
Dec_districts <- c("solenzo", "tougan", "banfora", "reo", "tenado", "pama", "bousse", "ziniare", "kampti")
Jul_districts <- "po"


HF_cases[which(HF_cases$District %in% Oct_districts & HF_cases$Date == "2016-10-01" &
                   (HF_cases$conf_rdt_mic_u5 / HF_cases$allout_u5) > 1), c("conf_rdt_u5", "conf_mic_u5")] <- NA

HF_cases[which(HF_cases$District %in% Nov_districts & HF_cases$Date == "2016-11-01" &
                   (HF_cases$conf_rdt_mic_u5 / HF_cases$allout_u5) > 1), c("conf_rdt_u5", "conf_mic_u5")] <- NA

HF_cases[which(HF_cases$District %in% Sep_districts & HF_cases$Date == "2016-09-01" &
                   (HF_cases$conf_rdt_mic_u5 / HF_cases$allout_u5) > 1), c("conf_rdt_u5", "conf_mic_u5")] <- NA

HF_cases[which(HF_cases$District %in% Dec_districts & HF_cases$Date == "2016-12-01" &
                   (HF_cases$conf_rdt_mic_u5 / HF_cases$allout_u5) > 1), c("conf_rdt_u5", "conf_mic_u5")] <- NA

HF_cases[which(HF_cases$District %in% Jul_districts & HF_cases$Date == "2016-07-01" &
                   (HF_cases$conf_rdt_mic_u5 / HF_cases$allout_u5) > 1), c("conf_rdt_u5", "conf_mic_u5")] <- NA








## rounding some columns to not have FP issues later
HF_cases$precip_era5 <- round(HF_cases$precip_era5, digits = 3)
HF_cases$air_temp_era5 <- round(HF_cases$air_temp_era5, digits = 3)
HF_cases$SMC.coverage <- round(HF_cases$SMC.coverage, digits = 3)





apply(HF_cases[,c("allout_u5", "conf_rdt_u5", "conf_mic_u5")], 2, function(y) sum(!is.na(y)))



####################################################################


rle.try_allout <- ddply(HF_cases, .(UID), summarize,
                        is_NA = rle(is.na(allout_u5))[2],
                        consec = rle(is.na(allout_u5))[1])


####################################################################

consec_NA_list_allout_u5 <- split(rle.try_allout[,2:3], rle.try_allout$UID)
consec_NA_list_allout_u5_reshaped <- lapply(consec_NA_list_allout_u5, function(x) {
    x <- do.call(cbind.data.frame, x);
    names(x) <- c("is.NA", "consec vals");
    return(x) })


## Finding HFs that have leq 5 NAs and no more than 2 in a row
imputing_HFs_allout_u5_list <- sapply(consec_NA_list_allout_u5_reshaped, function(x) {
    na_rows <- which(x$is.NA == TRUE);
    num_NAs <- sum(x[na_rows, "consec vals"]);
    
    if (num_NAs > 0 & num_NAs <= 6 & !(0 %in% as.numeric(x[na_rows, "consec vals"] <= 2)))
    {
        return(TRUE)
    } else {
        return(FALSE)
    }
})

imputing_HFs_allout_u5_data <- HF_cases[which(HF_cases$UID %in% names(imputing_HFs_allout_u5_list[imputing_HFs_allout_u5_list == TRUE])),]



imputing_HFs_allout_u5_data_list <- split(imputing_HFs_allout_u5_data$allout_u5, as.character(imputing_HFs_allout_u5_data$UID))

new_imputed_HFs_allout_u5_data_list <- lapply(imputing_HFs_allout_u5_data_list, function(x){

    cases_ts <- ts(x, start = c(2015, 1), deltat = 1/12);
    imp_cases <- na_seadec(cases_ts, algorithm = "ma");
    

    return(as.numeric(round(imp_cases))); # we round to preserve whole numbers
})


result_imputed_allout_u5_DF <- data.frame("UID" = sort(rep(names(new_imputed_HFs_allout_u5_data_list), 48)),
                                          "Date" = rep(unique(HF_cases$Date), length(new_imputed_HFs_allout_u5_data_list)),
                                          "imputed_allout_u5" = unlist(new_imputed_HFs_allout_u5_data_list, use.names = F))

for (U in unique(result_imputed_allout_u5_DF$UID))
    HF_cases[which(HF_cases$UID == U), "allout_u5"] <- result_imputed_allout_u5_DF[which(result_imputed_allout_u5_DF$UID == U), "imputed_allout_u5"]




####################################################################

rle.try_conf_u5 <- ddply(HF_cases, .(UID), summarize,
                         is_NA = rle(is.na(conf_rdt_u5))[2],
                         consec_NAs = rle(is.na(conf_rdt_u5))[1])

####################################################################


consec_NA_list_conf_u5 <- split(rle.try_conf_u5[,2:3], rle.try_conf_u5$UID)
consec_NA_list_conf_u5_reshaped <- lapply(consec_NA_list_conf_u5, function(x) {
    x <- do.call(cbind.data.frame, x);
    names(x) <- c("is.NA", "consec vals");
    return(x) })


## Finding HFs that have leq 5 NAs and no more than 2 in a row
imputing_HFs_conf_u5_list <- sapply(consec_NA_list_conf_u5_reshaped, function(x) {
    na_rows <- which(x$is.NA == TRUE);
    num_NAs <- sum(x[na_rows, "consec vals"]);
    
    if (num_NAs > 0 & num_NAs <= 6 & !(0 %in% as.numeric(x[na_rows, "consec vals"] <= 3)))
    {
        return(TRUE)
    } else {
        return(FALSE)
    }
})

imputing_HFs_conf_u5_data <- HF_cases[which(HF_cases$UID %in% names(imputing_HFs_conf_u5_list[imputing_HFs_conf_u5_list == TRUE])),]



imputing_HFs_conf_u5_data_list <- split(imputing_HFs_conf_u5_data$conf_rdt_u5, as.character(imputing_HFs_conf_u5_data$UID))

new_imputed_HFs_conf_u5_data_list <- lapply(imputing_HFs_conf_u5_data_list, function(x) {
    
    cases_ts <- ts(x, start = c(2015, 1), deltat = 1/12);
    imp_cases <- na_seadec(cases_ts, algorithm = "ma");
    
    return(as.numeric(round(imp_cases))); # we round to preserve whole numbers
})


result_imputed_conf_u5_DF <- data.frame("UID" = sort(rep(names(new_imputed_HFs_conf_u5_data_list), 48)),
                                        "Date" = rep(unique(HF_cases$Date), length(new_imputed_HFs_conf_u5_data_list)),
                                        "imputed_conf_rdt_u5" = unlist(new_imputed_HFs_conf_u5_data_list, use.names = F))

for (U in unique(result_imputed_conf_u5_DF$UID))
    HF_cases[which(HF_cases$UID == U), "conf_rdt_u5"] <- result_imputed_conf_u5_DF[which(result_imputed_conf_u5_DF$UID == U), "imputed_conf_rdt_u5"]





####################################################################

rle.try_conf_mic_u5 <- ddply(HF_cases, .(UID), summarize,
                             is_NA = rle(is.na(conf_mic_u5))[2],
                             consec_NAs = rle(is.na(conf_mic_u5))[1])

####################################################################


consec_NA_list_conf_mic_u5 <- split(rle.try_conf_mic_u5[,2:3], rle.try_conf_mic_u5$UID)
consec_NA_list_conf_mic_u5_reshaped <- lapply(consec_NA_list_conf_mic_u5, function(x) {
    x <- do.call(cbind.data.frame, x);
    names(x) <- c("is.NA", "consec vals");
    return(x) })


## Finding HFs that have leq 5 NAs and no more than 2 in a row
imputing_HFs_conf_mic_u5_list <- sapply(consec_NA_list_conf_mic_u5_reshaped, function(x) {
    na_rows <- which(x$is.NA == TRUE);
    num_NAs <- sum(x[na_rows, "consec vals"]);
    
    if (num_NAs > 0 & num_NAs <= 6 & !(0 %in% as.numeric(x[na_rows, "consec vals"] <= 2)))
    {
        return(TRUE)
    } else {
        return(FALSE)
    }
})

imputing_HFs_conf_mic_u5_data <- HF_cases[which(HF_cases$UID %in% names(imputing_HFs_conf_mic_u5_list[imputing_HFs_conf_mic_u5_list == TRUE])),]


imputing_HFs_conf_mic_u5_data_list <- split(imputing_HFs_conf_mic_u5_data$conf_mic_u5, as.character(imputing_HFs_conf_mic_u5_data$UID))

new_imputed_HFs_conf_mic_u5_data_list <- lapply(imputing_HFs_conf_mic_u5_data_list, function(x) {
    
    cases_ts <- ts(x, start = c(2015, 1), deltat = 1/12);
    imp_cases <- na_seadec(cases_ts, algorithm = "ma");
    
    
    return(as.numeric(round(imp_cases))); # we round to preserve whole numbers
})


result_imputed_conf_mic_u5_DF <- data.frame("UID" = sort(rep(names(new_imputed_HFs_conf_mic_u5_data_list), 48)),
                                            "Date" = rep(unique(HF_cases$Date), length(new_imputed_HFs_conf_mic_u5_data_list)),
                                            "imputed_conf_mic_u5" = unlist(new_imputed_HFs_conf_mic_u5_data_list, use.names = F))

for (U in unique(result_imputed_conf_mic_u5_DF$UID))
    HF_cases[which(HF_cases$UID == U), "conf_mic_u5"] <- result_imputed_conf_mic_u5_DF[which(result_imputed_conf_mic_u5_DF$UID == U), "imputed_conf_mic_u5"]










####################################################################


HF_cases[which(HF_cases$allout_u5 < 0), "allout_u5"] <- 0
HF_cases[which(HF_cases$conf_rdt_u5 < 0), "conf_rdt_u5"] <- 0
HF_cases[which(HF_cases$conf_mic_u5 < 0), "conf_mic_u5"] <- 0


####################################################################


HF_cases$conf_rdt_mic_u5 <- rowSums(HF_cases[,c("conf_rdt_u5", "conf_mic_u5")], na.rm = T)



apply(HF_cases[,c("allout_u5", "conf_rdt_u5", "conf_mic_u5")], 2, function(y) sum(is.na(y)))





write.csv(HF_cases, "~/Box/NU-malaria-team/projects/smc_impact/data/outputs/U5_HF_cases_smc_coords_imputed_rdts_and_allout_MA.csv", row.names = FALSE)


