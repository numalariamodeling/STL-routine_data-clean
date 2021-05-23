##############################################
##  Clean HF data before aggregating to DS  ##
##############################################
#
# Description:
#   Cleaning health facility data for exclusion criteria prior to aggregation to health district
#
#
#  Sebastian Rodriguez (sebastian@rodriguez.cr)
#  Last edited Mar 09, 2021
#




rm(list = ls(all = TRUE))

require("plyr")
require("dplyr")
require("zoo")


# Loading health facility dataset

HF_cases <- read.csv("~/Box/NU-malaria-team/projects/smc_impact/data/outputs/U5_HF_cases_smc_coords_imputed_rdts_and_allout_MA.csv", stringsAsFactors = FALSE)





tmp_total_HF <- ddply(HF_cases, c(.(District), .(Date)), summarise,
                      total_HF = length(unique(UID)))




#################################################################################################################################


good_rows <- which(!is.na(HF_cases$conf_rdt_u5) &
                       !is.na(HF_cases$allout_u5) &
                       HF_cases$allout_u5 != 0 &
                       HF_cases$conf_rdt_mic_u5 <= HF_cases$allout_u5)

HF_cases_good <- HF_cases[good_rows,]

bad_rows_age1 <- which(is.na(HF_cases_good$conf_rdt_age1) & !is.na(HF_cases_good$conf_rdt_age2) &
                           HF_cases_good$conf_rdt_age2 == HF_cases_good$conf_rdt_u5)
bad_rows_age2 <- which(is.na(HF_cases_good$conf_rdt_age2) & !is.na(HF_cases_good$conf_rdt_age1) &
                           HF_cases_good$conf_rdt_age1 == HF_cases_good$conf_rdt_u5)
HF_cases_good <- HF_cases_good[-c(bad_rows_age1, bad_rows_age2),]

# 46655 removed



tmp_reporting_HF <- ddply(HF_cases_good, c(.(District), .(Date)), summarise,
                          reporting_HF = length(unique(UID)))


#################################################################################################################################




D_cases <- ddply(HF_cases_good[,-c(11, 164:170)], 
                 c(.(Region), .(District), .(periodname), .(month), .(year), .(Date)),
                 numcolwise(sum, na.rm = TRUE))


# Get unique data for precip, air.temp, Pop, SMC_rec, num_children_smc we took off above
# merge with health district data
HF_cases_good$Number.of.children.treated.with.SMC <- round(HF_cases_good$Number.of.children.treated.with.SMC)
unique_rows <- unique(HF_cases_good[,c(1,3,7:10,11,164:168)])
D_cases <- left_join(D_cases, unique_rows,
                     by = c("Region", "District", "periodname", "month", "year", "Date"))




D_cases <- left_join(D_cases, tmp_total_HF, by = c("District", "Date"))
D_cases <- left_join(D_cases, tmp_reporting_HF, by = c("District", "Date"))



# Saving

write.csv(D_cases, "~/Box/NU-malaria-team/projects/smc_impact/data/outputs/U5_DS_cases_seasonal_smc_good_rows_MA_imputes.csv", row.names = FALSE)



