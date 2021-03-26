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


rle.try_allout_u5 <- ddply(HF_cases, .(UID), summarize,
                           is_NA = rle(is.na(allout_u5))[2],
                           consec = rle(is.na(allout_u5))[1])


####################################################################

consec_NA_list_allout_u5 <- split(rle.try_allout_u5[,2:3], rle.try_allout_u5$UID)
consec_NA_list_allout_u5_reshaped <- lapply(consec_NA_list_allout_u5, function(x) {
    x <- do.call(cbind.data.frame, x);
    names(x) <- c("is.NA", "consec vals");
    return(x) })


## Finding HFs that have leq 5 NAs and no more than 2 in a row
imputing_HFs_allout_u5_list <- sapply(consec_NA_list_allout_u5_reshaped, function(x) {
    na_rows <- which(x$is.NA == TRUE);
    num_NAs <- sum(x[na_rows, "consec vals"]);
    
    if (num_NAs <= 5 & !(0 %in% as.numeric(x[na_rows, "consec vals"] <= 2)))
    {
        return(TRUE)
    } else {
        return(FALSE)
    }
})

imputing_HFs_allout_u5_data <- HF_cases[which(HF_cases$UID %in% names(imputing_HFs_allout_u5_list[imputing_HFs_allout_u5_list == TRUE])),]



imputing_HFs_allout_u5_data_list <- split(imputing_HFs_allout_u5_data$allout_u5, as.character(imputing_HFs_allout_u5_data$UID))

new_imputed_HFs_allout_u5_data_list <- lapply(imputing_HFs_allout_u5_data_list, function(x) {
    cases_ts <- ts(x, start = c(2015, 1), deltat = 1/12);
    
    ma_flag <- 0;
    if (sum(cases_ts, na.rm = T) > 0)
    {
        imp_cases_and_flag <- tryCatch({
            imp_cases <- na_seadec(cases_ts, algorithm = "kalman");
            
            c(ma_flag, imp_cases);
        },
        warning = function(cond){
            imp_cases <- na_seadec(cases_ts, algorithm = "ma");
            ma_flag <- 2;

            imp_cases_and_flag <- c(ma_flag, imp_cases);   
            return(imp_cases_and_flag);
        })
        imp_cases <- imp_cases_and_flag[-1];
        ma_flag <- imp_cases_and_flag[1];
    } else {
        imp_cases <- na_seadec(cases_ts, algorithm = "ma");
        ma_flag <- 1;
    }
    
    ret_list <- list(cases_ts, round(imp_cases), ma_flag);
    
    return(ret_list); # we round to preserve whole numbers
})

sum(sapply(new_imputed_HFs_allout_u5_data_list, function(x) (x[[3]]) > 0))
for (i in 1:length(new_imputed_HFs_allout_u5_data_list))
{
    if (new_imputed_HFs_allout_u5_data_list[[i]][[3]] > 0)
        print(paste(i, new_imputed_HFs_allout_u5_data_list[[i]][[3]], sep = " "))
}


ggplot_na_distribution(new_imputed_HFs_allout_u5_data_list[[1449]][[1]],
                       title = names(new_imputed_HFs_allout_u5_data_list)[1449])
ggplot_na_distribution(new_imputed_HFs_allout_u5_data_list[[1492]][[1]],
                       title = names(new_imputed_HFs_allout_u5_data_list)[1492])


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
    
    if (num_NAs <= 5 & !(0 %in% as.numeric(x[na_rows, "consec vals"] <= 2)))
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
    
    ma_flag <- 0;
    if (sum(cases_ts, na.rm = T) > 0)
    {
        imp_cases_and_flag <- tryCatch({
            imp_cases <- na_seadec(cases_ts, algorithm = "kalman");
            
            c(ma_flag, imp_cases); 
        },
        warning = function(cond){
            imp_cases <- na_seadec(cases_ts, algorithm = "ma");
            ma_flag <- 2;
            
            imp_cases_and_flag <- c(ma_flag, imp_cases);   
            return(imp_cases_and_flag);
        })
        imp_cases <- imp_cases_and_flag[-1];
        ma_flag <- imp_cases_and_flag[1];
    } else {
        imp_cases <- na_seadec(cases_ts, algorithm = "ma");
        ma_flag <- 1;
    }

    ret_list <- list(cases_ts, round(imp_cases), ma_flag);
    
    return(ret_list); # we round to preserve whole numbers
})

sum(sapply(new_imputed_HFs_conf_u5_data_list, function(x) (x[[3]]) > 0))
for (i in 1:length(new_imputed_HFs_conf_u5_data_list))
{
    if (new_imputed_HFs_conf_u5_data_list[[i]][[3]] > 0)
        print(paste(i, new_imputed_HFs_conf_u5_data_list[[i]][[3]], sep = " "))
}


ggplot_na_distribution(new_imputed_HFs_conf_u5_data_list[[1252]][[1]],
                       title = names(new_imputed_HFs_conf_u5_data_list)[1252])
ggplot_na_distribution(new_imputed_HFs_conf_u5_data_list[[1535]][[1]],
                       title = names(new_imputed_HFs_conf_u5_data_list)[1535])
ggplot_na_distribution(new_imputed_HFs_conf_u5_data_list[[1551]][[1]],
                       title = names(new_imputed_HFs_conf_u5_data_list)[1551])
ggplot_na_distribution(new_imputed_HFs_conf_u5_data_list[[1579]][[1]],
                       title = names(new_imputed_HFs_conf_u5_data_list)[1579])


####################################################################


rle.try_allout_ov5 <- ddply(HF_cases, .(UID), summarize,
                            is_NA = rle(is.na(allout_ov5))[2],
                            consec = rle(is.na(allout_ov5))[1])


####################################################################

consec_NA_list_allout_ov5 <- split(rle.try_allout_ov5[,2:3], rle.try_allout_ov5$UID)
consec_NA_list_allout_ov5_reshaped <- lapply(consec_NA_list_allout_ov5, function(x) {
    x <- do.call(cbind.data.frame, x);
    names(x) <- c("is.NA", "consec vals");
    return(x) })


## Finding HFs that have leq 5 NAs and no more than 2 in a row
imputing_HFs_allout_ov5_list <- sapply(consec_NA_list_allout_ov5_reshaped, function(x) {
    na_rows <- which(x$is.NA == TRUE);
    num_NAs <- sum(x[na_rows, "consec vals"]);
    
    if (num_NAs <= 5 & !(0 %in% as.numeric(x[na_rows, "consec vals"] <= 2)))
    {
        return(TRUE)
    } else {
        return(FALSE)
    }
})

imputing_HFs_allout_ov5_data <- HF_cases[which(HF_cases$UID %in% names(imputing_HFs_allout_ov5_list[imputing_HFs_allout_ov5_list == TRUE])),]



imputing_HFs_allout_ov5_data_list <- split(imputing_HFs_allout_ov5_data$allout_ov5, as.character(imputing_HFs_allout_ov5_data$UID))

new_imputed_HFs_allout_ov5_data_list <- lapply(imputing_HFs_allout_ov5_data_list, function(x) {
    cases_ts <- ts(x, start = c(2015, 1), deltat = 1/12);
    
    ma_flag <- 0;
    if (sum(cases_ts, na.rm = T) > 0)
    {
        imp_cases_and_flag <- tryCatch({
            imp_cases <- na_seadec(cases_ts, algorithm = "kalman");
            
            c(ma_flag, imp_cases); 
        },
        warning = function(cond){
            imp_cases <- na_seadec(cases_ts, algorithm = "ma");
            ma_flag <- 2;
            
            imp_cases_and_flag <- c(ma_flag, imp_cases);   
            return(imp_cases_and_flag);
        })
        imp_cases <- imp_cases_and_flag[-1];
        ma_flag <- imp_cases_and_flag[1];
    } else {
        imp_cases <- na_seadec(cases_ts, algorithm = "ma");
        ma_flag <- 1;
    }
    
    ret_list <- list(cases_ts, round(imp_cases), ma_flag);
    
    return(ret_list); # we round to preserve whole numbers
})

sum(sapply(new_imputed_HFs_allout_ov5_data_list, function(x) (x[[3]]) > 0))
for (i in 1:length(new_imputed_HFs_allout_ov5_data_list))
{
    if (new_imputed_HFs_allout_ov5_data_list[[i]][[3]] > 0)
        print(paste(i, new_imputed_HFs_allout_ov5_data_list[[i]][[3]], sep = " "))
}


ggplot_na_distribution(new_imputed_HFs_allout_ov5_data_list[[66]][[1]],
                       title = names(new_imputed_HFs_allout_ov5_data_list)[66])
ggplot_na_distribution(new_imputed_HFs_allout_ov5_data_list[[476]][[1]],
                       title = names(new_imputed_HFs_allout_ov5_data_list)[476])


####################################################################

rle.try_conf_ov5 <- ddply(HF_cases, .(UID), summarize,
                          is_NA = rle(is.na(conf_rdt_ov5))[2],
                          consec_NAs = rle(is.na(conf_rdt_ov5))[1])

####################################################################


consec_NA_list_conf_ov5 <- split(rle.try_conf_ov5[,2:3], rle.try_conf_ov5$UID)
consec_NA_list_conf_ov5_reshaped <- lapply(consec_NA_list_conf_ov5, function(x) {
    x <- do.call(cbind.data.frame, x);
    names(x) <- c("is.NA", "consec vals");
    return(x) })


## Finding HFs that have leq 5 NAs and no more than 2 in a row
imputing_HFs_conf_ov5_list <- sapply(consec_NA_list_conf_ov5_reshaped, function(x) {
    na_rows <- which(x$is.NA == TRUE);
    num_NAs <- sum(x[na_rows, "consec vals"]);
    
    if (num_NAs <= 5 & !(0 %in% as.numeric(x[na_rows, "consec vals"] <= 2)))
    {
        return(TRUE)
    } else {
        return(FALSE)
    }
})

imputing_HFs_conf_ov5_data <- HF_cases[which(HF_cases$UID %in% names(imputing_HFs_conf_ov5_list[imputing_HFs_conf_ov5_list == TRUE])),]



imputing_HFs_conf_ov5_data_list <- split(imputing_HFs_conf_ov5_data$conf_rdt_ov5, as.character(imputing_HFs_conf_ov5_data$UID))

new_imputed_HFs_conf_ov5_data_list <- lapply(imputing_HFs_conf_ov5_data_list, function(x) {
    cases_ts <- ts(x, start = c(2015, 1), deltat = 1/12);
    
    ma_flag <- 0;
    if (sum(cases_ts, na.rm = T) > 0)
    {
        imp_cases_and_flag <- tryCatch({
            imp_cases <- na_seadec(cases_ts, algorithm = "kalman");
            
            c(ma_flag, imp_cases); 
        },
        warning = function(cond){
            imp_cases <- na_seadec(cases_ts, algorithm = "ma");
            ma_flag <- 2;
            
            imp_cases_and_flag <- c(ma_flag, imp_cases);   
            return(imp_cases_and_flag);
        })
        imp_cases <- imp_cases_and_flag[-1];
        ma_flag <- imp_cases_and_flag[1];
    } else {
        imp_cases <- na_seadec(cases_ts, algorithm = "ma");
        ma_flag <- 1;
    }
    
    ret_list <- list(cases_ts, round(imp_cases), ma_flag);
    
    return(ret_list); # we round to preserve whole numbers
})

sum(sapply(new_imputed_HFs_conf_ov5_data_list, function(x) (x[[3]]) > 0))
for (i in 1:length(new_imputed_HFs_conf_ov5_data_list))
{
    if (new_imputed_HFs_conf_ov5_data_list[[i]][[3]] > 0)
        print(paste(i, new_imputed_HFs_conf_ov5_data_list[[i]][[3]], sep = " "))
}


ggplot_na_distribution(new_imputed_HFs_conf_ov5_data_list[[250]][[1]],
                       title = names(new_imputed_HFs_conf_ov5_data_list)[250])
ggplot_na_distribution(new_imputed_HFs_conf_ov5_data_list[[751]][[1]],
                       title = names(new_imputed_HFs_conf_ov5_data_list)[751])
ggplot_na_distribution(new_imputed_HFs_conf_ov5_data_list[[1494]][[1]],
                       title = names(new_imputed_HFs_conf_ov5_data_list)[1494])
ggplot_na_distribution(new_imputed_HFs_conf_ov5_data_list[[1554]][[1]],
                       title = names(new_imputed_HFs_conf_ov5_data_list)[1554])


####################################################################

rle.try_maladm <- ddply(HF_cases, .(UID), summarize,
                        is_NA = rle(is.na(maladm_u5))[2],
                        consec_NAs = rle(is.na(maladm_u5))[1])

####################################################################


consec_NA_list_maladm <- split(rle.try_maladm[,2:3], rle.try_maladm$UID)
consec_NA_list_maladm_reshaped <- lapply(consec_NA_list_maladm, function(x) {
    x <- do.call(cbind.data.frame, x);
    names(x) <- c("is.NA", "consec vals");
    return(x) })


## Finding HFs that have leq 5 NAs and no more than 2 in a row
imputing_HFs_maladm_list <- sapply(consec_NA_list_maladm_reshaped, function(x) {
    na_rows <- which(x$is.NA == TRUE);
    num_NAs <- sum(x[na_rows, "consec vals"]);
    
    if (num_NAs <= 5 & !(0 %in% as.numeric(x[na_rows, "consec vals"] <= 2)))
    {
        return(TRUE)
    } else {
        return(FALSE)
    }
})

imputing_HFs_maladm_data <- HF_cases[which(HF_cases$UID %in% names(imputing_HFs_maladm_list[imputing_HFs_maladm_list == TRUE])),]



imputing_HFs_maladm_data_list <- split(imputing_HFs_maladm_data$maladm_u5, as.character(imputing_HFs_maladm_data$UID))

new_imputed_HFs_maladm_data_list <- lapply(imputing_HFs_maladm_data_list, function(x) {
    cases_ts <- ts(x, start = c(2015, 1), deltat = 1/12);
    
    ma_flag <- 0;
    if (sum(cases_ts, na.rm = T) > 0)
    {
        imp_cases_and_flag <- tryCatch({
            imp_cases <- na_seadec(cases_ts, algorithm = "kalman");
            
            c(ma_flag, imp_cases); 
        },
        warning = function(cond){
            imp_cases <- na_seadec(cases_ts, algorithm = "ma");
            ma_flag <- 2;
            
            imp_cases_and_flag <- c(ma_flag, imp_cases);   
            return(imp_cases_and_flag);
        })
        imp_cases <- imp_cases_and_flag[-1];
        ma_flag <- imp_cases_and_flag[1];
    } else {
        imp_cases <- na_seadec(cases_ts, algorithm = "ma");
        ma_flag <- 1;
    }
    
    ret_list <- list(cases_ts, round(imp_cases), ma_flag);
    
    return(ret_list); # we round to preserve whole numbers
})

sum(sapply(new_imputed_HFs_maladm_data_list, function(x) (x[[3]]) > 0))
for (i in 1:length(new_imputed_HFs_maladm_data_list))
{
    if (new_imputed_HFs_maladm_data_list[[i]][[3]] > 0)
        print(paste(i, new_imputed_HFs_maladm_data_list[[i]][[3]], sep = " "))
}


ggplot_na_distribution(new_imputed_HFs_maladm_data_list[[73]][[1]],
                       title = names(new_imputed_HFs_maladm_data_list)[73])
ggplot_na_distribution(new_imputed_HFs_maladm_data_list[[120]][[1]],
                       title = names(new_imputed_HFs_maladm_data_list)[120])
ggplot_na_distribution(new_imputed_HFs_maladm_data_list[[584]][[1]],
                       title = names(new_imputed_HFs_maladm_data_list)[584])
ggplot_na_distribution(new_imputed_HFs_maladm_data_list[[715]][[1]],
                       title = names(new_imputed_HFs_maladm_data_list)[715])


####################################################################


HF_cases[which(HF_cases$allout_u5 < 0), "allout_u5"] <- 0
HF_cases[which(HF_cases$conf_rdt_u5 < 0), "conf_rdt_u5"] <- 0

HF_cases[which(HF_cases$allout_ov5 < 0), "allout_ov5"] <- 0
HF_cases[which(HF_cases$conf_rdt_ov5 < 0), "conf_rdt_ov5"] <- 0

HF_cases[which(HF_cases$maladm_u5 < 0), "maladm_u5"] <- 0



####################################################################

HF_cases$conf_rdt_mic_u5 <- rowSums(HF_cases[,c("conf_rdt_u5", "conf_mic_u5")], na.rm = T)

HF_cases$conf_rdt_mic_ov5 <- rowSums(HF_cases[,c("conf_rdt_ov5", "conf_mic_ov5")], na.rm = T)




