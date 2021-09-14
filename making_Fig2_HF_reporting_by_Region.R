
rm(list = ls(all = TRUE))

require("plyr")
require("dplyr")
require("zoo")
require("reshape2")
require("ggplot2")


# Loading health facility dataset

# HF_cases <- read.csv("~/Box/NU-malaria-team/projects/smc_impact/data/outputs/U5_HF_cases_smc_coords_imputed_rdts_and_allout_MA.csv", stringsAsFactors = FALSE)
HF_cases <- read.csv("~/Box/NU-malaria-team/projects/smc_impact/data/outputs/U5_HF_cases_smc_coords_imputed_rdts_and_allout_testing_pres_MA.csv", stringsAsFactors = FALSE)


HF_cases$Date <- as.Date(as.yearmon(HF_cases$Date))
HF_cases <- HF_cases[order(HF_cases$UID, HF_cases$Date),]




#########################################################################################



bad_rows_rdt_age1 <- which(is.na(HF_cases$conf_rdt_age1) & !is.na(HF_cases$conf_rdt_age2) &
                               !is.na(HF_cases$conf_rdt_u5) &
                               HF_cases$conf_rdt_age2 == HF_cases$conf_rdt_u5)
bad_rows_rdt_age2 <- which(is.na(HF_cases$conf_rdt_age2) & !is.na(HF_cases$conf_rdt_age1) &
                               !is.na(HF_cases$conf_rdt_u5) &
                               HF_cases$conf_rdt_age1 == HF_cases$conf_rdt_u5)
bad_rows_rdt_u5 <- which(is.na(HF_cases$conf_rdt_u5))

bad_rdt_rows <- unique(c(bad_rows_rdt_age1, bad_rows_rdt_age2, bad_rows_rdt_u5))


bad_rows_allout <- which(is.na(HF_cases$allout_u5) | HF_cases$allout_u5 == 0)


bad_rows_inconsistent <- which(HF_cases$conf_rdt_mic_u5 > HF_cases$allout_u5)


bad_rows_missing <- unique(c(bad_rdt_rows, bad_rows_allout))

bad_rows_all_2 <- unique(c(bad_rows_missing, bad_rows_inconsistent))

bad_rows_missing_allout_only <- bad_rows_allout[which(!(bad_rows_allout %in% bad_rows_inconsistent) &
                                                          !(bad_rows_allout %in% bad_rdt_rows))]
bad_rows_missing_NOT_allout_only <- bad_rows_missing[which(!(bad_rows_missing %in% bad_rows_missing_allout_only))]
bad_rows_inconsistent_only <- bad_rows_inconsistent[which(!(bad_rows_inconsistent %in% bad_rows_missing))]


HF_cases_area <- HF_cases[,c("Region", "UID", "Date")]

HF_cases_area$bad_rows <- "good report"
HF_cases_area[bad_rows_missing_allout_only, "bad_rows"] <- "missing all-cause outpatient visits cases"
HF_cases_area[bad_rows_missing_NOT_allout_only, "bad_rows"] <- "missing confirmed malaria cases or missing both values"
HF_cases_area[bad_rows_inconsistent_only, "bad_rows"] <- "inconsistent report"


HF_cases_area_casted <- dcast(HF_cases_area, Region + Date ~ bad_rows)
HF_cases_area_casted$total_HFs <- rowSums(HF_cases_area_casted[,3:6])
HF_cases_area_by_Region <- melt(HF_cases_area_casted, id.vars = c("Region", "Date", "total_HFs"),
                                variable.name = "Type", value.name = "reports")
HF_cases_area_by_Region$Type <- factor(HF_cases_area_by_Region$Type,
                                       levels = c("inconsistent report",
                                                  "missing confirmed malaria cases or missing both values",
                                                  "missing all-cause outpatient visits cases",
                                                  "good report"))



ggplot(HF_cases_area_by_Region, aes(fill = Type, y = reports/total_HFs, x = Date))  + 
    geom_area() + 
    scale_x_date(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    facet_wrap(~Region) +
    theme(panel.spacing.x = unit(2, "mm"))




