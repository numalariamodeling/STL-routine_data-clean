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
    
    if (num_NAs > 0 & num_NAs <= 5 & !(0 %in% as.numeric(x[na_rows, "consec vals"] <= 2)))
    {
        return(TRUE)
    } else {
        return(FALSE)
    }
})

imputing_HFs_allout_u5_data <- HF_cases[which(HF_cases$UID %in% names(imputing_HFs_allout_u5_list[imputing_HFs_allout_u5_list == TRUE])),]



imputing_HFs_allout_u5_data_list <- split(imputing_HFs_allout_u5_data$allout_u5, as.character(imputing_HFs_allout_u5_data$UID))

# new_imputed_HFs_allout_u5_data_list <- lapply(imputing_HFs_allout_u5_data_list, function(x) {
#     
#     # Normalizing ts first
#     na_inds <- which(is.na(x));
#     tmp_x <- x[-na_inds];
#     
#     x_mean <- mean(tmp_x);
#     x_sd <- sd(tmp_x);
#     x_norm <- (tmp_x - x_mean) / x_sd;
#     
#     x[which(!is.na(x))] <- x_norm;
#     
#     
#     cases_ts <- ts(x, start = c(2015, 1), deltat = 1/12);
#     
#     if (sum(tmp_x, na.rm = T) > 0)
#     {
#         imp_cases <- tryCatch({
#             na_seadec(cases_ts, algorithm = "kalman");
#         },
#         warning = function(cond){
#             print("ummm");
#             na_seadec(cases_ts, algorithm = "ma");
#         })
#     } else {
#         imp_cases <- na_seadec(cases_ts, algorithm = "ma");
#     }
#     
#     # un-normalizing
#     imp_cases <- (imp_cases * x_sd) + x_mean
#     
#     
#     return(as.numeric(round(imp_cases))); # we round to preserve whole numbers
# })



new_imputed_HFs_allout_u5_data_list <- list()

for (i in 1:length(imputing_HFs_allout_u5_data_list))
{
    x <- imputing_HFs_allout_u5_data_list[[i]]
    
    # Normalizing ts first
    na_inds <- which(is.na(x));
    tmp_x <- x[-na_inds];
    
    if (sum(tmp_x) != 0)
    {
        x_mean <- mean(tmp_x);
        x_sd <- sd(tmp_x);
        x_norm <- (tmp_x - x_mean) / x_sd;
        
        x[which(!is.na(x))] <- x_norm;
    }
    
    cases_ts <- ts(x, start = c(2015, 1), deltat = 1/12);
    
    if (sum(tmp_x, na.rm = T) > 0)
    {
        imp_cases <- tryCatch({
            na_seadec(cases_ts, algorithm = "kalman");
        },
        warning = function(cond){
            print("ummm");
            na_seadec(cases_ts, algorithm = "ma");
        })
    } else {
        imp_cases <- na_seadec(cases_ts, algorithm = "ma");
    }
    
    # un-normalizing
    if (sum(tmp_x) != 0)
    {
        imp_cases <- (imp_cases * x_sd) + x_mean;
    }
    
    new_imputed_HFs_allout_u5_data_list[[i]] <- imp_cases
}


names(new_imputed_HFs_allout_u5_data_list) <- names(imputing_HFs_allout_u5_data_list)


















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
    
    if (num_NAs > 0 & num_NAs <= 5 & !(0 %in% as.numeric(x[na_rows, "consec vals"] <= 2)))
    {
        return(TRUE)
    } else {
        return(FALSE)
    }
})

imputing_HFs_conf_u5_data <- HF_cases[which(HF_cases$UID %in% names(imputing_HFs_conf_u5_list[imputing_HFs_conf_u5_list == TRUE])),]



imputing_HFs_conf_u5_data_list <- split(imputing_HFs_conf_u5_data$conf_rdt_u5, as.character(imputing_HFs_conf_u5_data$UID))

new_imputed_HFs_conf_u5_data_list <- lapply(imputing_HFs_conf_u5_data_list, function(x) {
    
    # Normalizing ts first
    na_inds <- which(is.na(x));
    tmp_x <- x[-na_inds];
    
    x_mean <- mean(tmp_x);
    x_sd <- sd(tmp_x);
    x_norm <- (tmp_x - x_mean) / x_sd;
    
    x[which(!is.na(x))] <- x_norm;
    
    
    cases_ts <- ts(x, start = c(2015, 1), deltat = 1/12);

    if (sum(tmp_x, na.rm = T) > 0)
    {
        imp_cases <- tryCatch({
            na_seadec(cases_ts, algorithm = "kalman");
        },
        warning = function(cond){
            print("ummm")
            na_seadec(cases_ts, algorithm = "ma");
        })
    } else {
        imp_cases <- na_seadec(cases_ts, algorithm = "ma");
    }
    
    # un-normalizing
    imp_cases <- (imp_cases * x_sd) + x_mean;
    
    
    return(as.numeric(round(imp_cases))); # we round to preserve whole numbers
})


result_imputed_conf_u5_DF <- data.frame("UID" = sort(rep(names(new_imputed_HFs_conf_u5_data_list), 48)),
                                     "Date" = rep(unique(HF_cases$Date), length(new_imputed_HFs_conf_u5_data_list)),
                                     "imputed_conf_rdt_u5" = unlist(new_imputed_HFs_conf_u5_data_list, use.names = F))

for (U in unique(result_imputed_conf_u5_DF$UID))
    HF_cases[which(HF_cases$UID == U), "conf_rdt_u5"] <- result_imputed_conf_u5_DF[which(result_imputed_conf_u5_DF$UID == U), "imputed_conf_rdt_u5"]






####################################################################


HF_cases[which(HF_cases$allout_u5 < 0), "allout_u5"] <- 0
HF_cases[which(HF_cases$conf_rdt_u5 < 0), "conf_rdt_u5"] <- 0



####################################################################


HF_cases$conf_rdt_mic_u5 <- rowSums(HF_cases[,c("conf_rdt_u5", "conf_mic_u5")], na.rm = T)



# write.csv(HF_cases, "~/Box/NU-malaria-team/projects/smc_impact/data/outputs/U5_HF_cases_smc_coords_imputed_rdts_and_allout_kalman.csv", row.names = FALSE)


