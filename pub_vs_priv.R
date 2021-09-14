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



#################################################################################################################################




D_cases <- ddply(HF_cases[,-c(8, 11, 164:170)], 
                 c(.(Region), .(District), .(year)),
                 numcolwise(sum, na.rm = TRUE))


# Get unique data for precip, air.temp, Pop, SMC_rec, num_children_smc we took off above
# merge with health district data
HF_cases$Number.of.children.treated.with.SMC <- round(HF_cases$Number.of.children.treated.with.SMC)
unique_rows <- unique(HF_cases[,c(1,3,9,11)])
D_cases <- left_join(D_cases, unique_rows,
                     by = c("Region", "District", "year"))




# Saving

write.csv(D_cases, "~/Box/NU-malaria-team/projects/smc_impact/data/outputs/supersupertmp.csv", row.names = FALSE)





PP_cases <- read.csv("~/OneDrive/Documents/public vs private checking.csv", stringsAsFactors = FALSE)[-c(281:307),-c(4:5)]


district_map <- cbind(sort(unique(D_cases$District)),
                      sort(unique(PP_cases$X.2))[c(1:40,42,41,43:70)])

for (i in 1:nrow(district_map))
{
  DS <- district_map[i, 2]
  new_name_DS <- district_map[i, 1]
  
  PP_cases[which(PP_cases$X.2 == DS), "District"] <- new_name_DS
}



library("ggplot2")


plotting_data <- left_join(D_cases, PP_cases, by = c("year" = "X", "District"))


ggplot(data = plotting_data, aes(x = conf_rdt_mic_u5,
                                 y = Malaria.cases.confirmed..RDT.Microscopy....5.years.)) +
  geom_point() + geom_abline(slope = 1, intercept = 0) +
  ylab("Conf malaria cases among U5 from national MPR and Stratification data") +
  xlab("Conf malaria cases among u5 from HMIS data")



