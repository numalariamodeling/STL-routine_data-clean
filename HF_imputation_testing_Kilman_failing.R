######################################################
##  Testing when Kalman smoothing imputation fails  ##
######################################################
#
#
#
#
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
    
    if (num_NAs > 0 & num_NAs <= 5 & !(0 %in% as.numeric(x[na_rows, "consec vals"] <= 2)))
    {
        return(TRUE)
    } else {
        return(FALSE)
    }
})

imputing_HFs_allout_u5_data <- HF_cases[which(HF_cases$UID %in% names(imputing_HFs_allout_u5_list[imputing_HFs_allout_u5_list == TRUE])),]



imputing_HFs_allout_u5_data_list <- split(imputing_HFs_allout_u5_data$allout_u5, as.character(imputing_HFs_allout_u5_data$UID))

new_imputed_HFs_allout_u5_data_list <- lapply(imputing_HFs_allout_u5_data_list, function(x) {
    
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
    
    ma_flag <- 0;
    if (sum(tmp_x, na.rm = T) > 0)
    {
        imp_cases_and_flag <- tryCatch({
            imp_cases <- na_kalman(cases_ts);
            
            c(ma_flag, imp_cases);
        },
        warning = function(cond){
            imp_cases <- na_ma(cases_ts);
            ma_flag <- 2;

            imp_cases_and_flag <- c(ma_flag, imp_cases);   
            return(imp_cases_and_flag);
        })
        imp_cases <- imp_cases_and_flag[-1];
        ma_flag <- imp_cases_and_flag[1];
    } else {
        imp_cases <- na_ma(cases_ts);
        ma_flag <- 1;
    }
    
    
    # un-normalizing
    if (sum(tmp_x) != 0)
    {
        imp_cases <- (imp_cases * x_sd) + x_mean;
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


ggplot_na_distribution(new_imputed_HFs_allout_u5_data_list[[160]][[1]],
                       title = paste(names(new_imputed_HFs_allout_u5_data_list)[160], "allout_u5", sep = " "))
ggplot_na_distribution(new_imputed_HFs_allout_u5_data_list[[281]][[1]],
                       title = paste(names(new_imputed_HFs_allout_u5_data_list)[281], "allout_u5", sep = " "))
ggplot_na_distribution(new_imputed_HFs_allout_u5_data_list[[282]][[1]],
                       title = paste(names(new_imputed_HFs_allout_u5_data_list)[282], "allout_u5", sep = " "))









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
    
    if (sum(tmp_x) != 0)
    {
        x_mean <- mean(tmp_x);
        x_sd <- sd(tmp_x);
        x_norm <- (tmp_x - x_mean) / x_sd;
        
        x[which(!is.na(x))] <- x_norm;
    }
    
    
    cases_ts <- ts(x, start = c(2015, 1), deltat = 1/12);
    
    ma_flag <- 0;
    if (sum(tmp_x, na.rm = T) > 0)
    {
        imp_cases_and_flag <- tryCatch({
            imp_cases <- na_kalman(cases_ts);
            
            c(ma_flag, imp_cases); 
        },
        warning = function(cond){
            imp_cases <- na_ma(cases_ts);
            ma_flag <- 2;
            
            imp_cases_and_flag <- c(ma_flag, imp_cases);   
            return(imp_cases_and_flag);
        })
        imp_cases <- imp_cases_and_flag[-1];
        ma_flag <- imp_cases_and_flag[1];
    } else {
        imp_cases <- na_ma(cases_ts);
        ma_flag <- 1;
    }
    
    
    # un-normalizing
    if (sum(tmp_x) != 0)
    {
        imp_cases <- (imp_cases * x_sd) + x_mean;
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


# ggplot_na_distribution(new_imputed_HFs_conf_u5_data_list[[1252]][[1]],
#                        title = paste(names(new_imputed_HFs_conf_u5_data_list)[1252], "conf_u5", sep = " "))
# ggplot_na_distribution(new_imputed_HFs_conf_u5_data_list[[1535]][[1]],
#                        title = paste(names(new_imputed_HFs_conf_u5_data_list)[1535], "conf_u5", sep = " "))
# ggplot_na_distribution(new_imputed_HFs_conf_u5_data_list[[1551]][[1]],
#                        title = paste(names(new_imputed_HFs_conf_u5_data_list)[1551], "conf_u5", sep = " "))
# ggplot_na_distribution(new_imputed_HFs_conf_u5_data_list[[1579]][[1]],
#                        title = paste(names(new_imputed_HFs_conf_u5_data_list)[1579], "conf_u5", sep = " "))




ggplot_na_distribution(new_imputed_HFs_conf_u5_data_list[[401]][[1]],
                       title = paste(names(new_imputed_HFs_conf_u5_data_list)[401], "conf_u5", sep = " "))
ggplot_na_distribution(new_imputed_HFs_conf_u5_data_list[[487]][[1]],
                       title = paste(names(new_imputed_HFs_conf_u5_data_list)[487], "conf_u5", sep = " "))
ggplot_na_distribution(new_imputed_HFs_conf_u5_data_list[[705]][[1]],
                       title = paste(names(new_imputed_HFs_conf_u5_data_list)[705], "conf_u5", sep = " "))
ggplot_na_distribution(new_imputed_HFs_conf_u5_data_list[[839]][[1]],
                       title = paste(names(new_imputed_HFs_conf_u5_data_list)[839], "conf_u5", sep = " "))


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
                       title = paste(names(new_imputed_HFs_allout_ov5_data_list)[66], "allout_ov5", sep = " "))
ggplot_na_distribution(new_imputed_HFs_allout_ov5_data_list[[476]][[1]],
                       title = paste(names(new_imputed_HFs_allout_ov5_data_list)[476], "allout_ov5", sep = " "))


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
                       title = paste(names(new_imputed_HFs_conf_ov5_data_list)[250], "conf_ov5", sep = " "))
ggplot_na_distribution(new_imputed_HFs_conf_ov5_data_list[[751]][[1]],
                       title = paste(names(new_imputed_HFs_conf_ov5_data_list)[751], "conf_ov5", sep = " "))
ggplot_na_distribution(new_imputed_HFs_conf_ov5_data_list[[1494]][[1]],
                       title = paste(names(new_imputed_HFs_conf_ov5_data_list)[1494], "conf_ov5", sep = " "))
ggplot_na_distribution(new_imputed_HFs_conf_ov5_data_list[[1554]][[1]],
                       title = paste(names(new_imputed_HFs_conf_ov5_data_list)[1554], "conf_ov5", sep = " "))


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
                       title = paste(names(new_imputed_HFs_maladm_data_list)[73], "maladm_u5", sep = " "))
ggplot_na_distribution(new_imputed_HFs_maladm_data_list[[120]][[1]],
                       title = paste(names(new_imputed_HFs_maladm_data_list)[120], "maladm_u5", sep = " "))
ggplot_na_distribution(new_imputed_HFs_maladm_data_list[[584]][[1]],
                       title = paste(names(new_imputed_HFs_maladm_data_list)[584], "maladm_u5", sep = " "))
ggplot_na_distribution(new_imputed_HFs_maladm_data_list[[715]][[1]],
                       title = paste(names(new_imputed_HFs_maladm_data_list)[715], "maladm_u5", sep = " "))


####################################################################

rle.try_maladm <- ddply(HF_cases, .(UID), summarize,
                        is_NA = rle(is.na(maladm_ov5))[2],
                        consec_NAs = rle(is.na(maladm_ov5))[1])

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



imputing_HFs_maladm_data_list <- split(imputing_HFs_maladm_data$maladm_ov5, as.character(imputing_HFs_maladm_data$UID))

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


ggplot_na_distribution(new_imputed_HFs_maladm_data_list[[761]][[1]],
                       title = paste(names(new_imputed_HFs_maladm_data_list)[761], "maladm_ov5", sep = " "))
ggplot_na_distribution(new_imputed_HFs_maladm_data_list[[860]][[1]],
                       title = paste(names(new_imputed_HFs_maladm_data_list)[860], "maladm_ov5", sep = " "))
ggplot_na_distribution(new_imputed_HFs_maladm_data_list[[936]][[1]],
                       title = paste(names(new_imputed_HFs_maladm_data_list)[936], "maladm_ov5", sep = " "))


####################################################################


HF_cases[which(HF_cases$allout_u5 < 0), "allout_u5"] <- 0
HF_cases[which(HF_cases$conf_rdt_u5 < 0), "conf_rdt_u5"] <- 0

HF_cases[which(HF_cases$allout_ov5 < 0), "allout_ov5"] <- 0
HF_cases[which(HF_cases$conf_rdt_ov5 < 0), "conf_rdt_ov5"] <- 0

HF_cases[which(HF_cases$maladm_u5 < 0), "maladm_u5"] <- 0
HF_cases[which(HF_cases$maladm_ov5 < 0), "maladm_ov5"] <- 0


####################################################################

HF_cases$conf_rdt_mic_u5 <- rowSums(HF_cases[,c("conf_rdt_u5", "conf_mic_u5")], na.rm = T)

HF_cases$conf_rdt_mic_ov5 <- rowSums(HF_cases[,c("conf_rdt_ov5", "conf_mic_ov5")], na.rm = T)









###########
## Testing with arima TS struct

na_seadec(new_imputed_HFs_maladm_data_list[[761]][[1]], algorithm = "kalman")


normalized_TS <- (new_imputed_HFs_maladm_data_list[[761]][[1]][-39] - mean(new_imputed_HFs_maladm_data_list[[761]][[1]][-39]))/sd(new_imputed_HFs_maladm_data_list[[761]][[1]][-39])
normalized_TS <- c(normalized_TS[1:38], NA, normalized_TS[39:47])
normalized_TS <- ts(normalized_TS, start = c(2015, 1), deltat = 1/12);
na_seadec(normalized_TS, algorithm = "kalman")

plot_norm <- ggplot_na_imputations(normalized_TS, na_seadec(normalized_TS, algorithm = "kalman"))
plot_orig <- ggplot_na_distribution(new_imputed_HFs_maladm_data_list[[761]][[1]])
grid.arrange(plot_norm, plot_orig, ncol=1)



na_seadec(new_imputed_HFs_maladm_data_list[[761]][[1]], algorithm = "ma")


na_seadec(new_imputed_HFs_maladm_data_list[[761]][[1]], algorithm = "kalman",
          model = "auto.arima", optim.method = "L-BFGS-B") # Same as typing nothing or NULL

na_seadec(new_imputed_HFs_maladm_data_list[[761]][[1]], algorithm = "kalman",
          model = "auto.arima", optim.method = "Nelder-Mead") # Default, same as "L-BFGS-B"

na_seadec(new_imputed_HFs_maladm_data_list[[761]][[1]], algorithm = "kalman",
          model = "auto.arima", optim.method = "Brent") # Same as "L-BFGS-B"

na_seadec(new_imputed_HFs_maladm_data_list[[761]][[1]], algorithm = "kalman",
          model = "auto.arima", optim.method = "SANN") # Same as "L-BFGS-B"



