##################################################
##  Adding columns of SMC coverage and climate  ##
##################################################
#
# Description:
#   Joining columns for SMC coverage and climate data
#
#
#  Sebastian Rodriguez (sebastian@rodriguez.cr)
#  Last edited Mar 09, 2021
#



rm(list = ls(all = TRUE))


library("zoo")
library("ggplot2")
library("plyr")
library("dplyr")
library("rgdal")
library("reshape2")



cases_cleaned <- read.csv("~/Box/NU-malaria-team/projects/smc_impact/data/outputs/BFA_Routine_case_data_HF_aggregate_Seb_sum_columns.csv", header = TRUE, strip.white = TRUE, stringsAsFactors = FALSE)

# changing name of district column
names(cases_cleaned)[c(1,3,5)] <- c("Region", "District", "HF")

# fixing DS names to merge with SMC file
cases_cleaned[,3] <- gsub("DS ", "", cases_cleaned[,3])




SMC <- read.csv("~/Box/NU-malaria-team/projects/hbhi_burkina/simulation_inputs/SMC_with_coverage_v3.csv", header = TRUE, stringsAsFactors = FALSE)
DS_order_in_SMC <- c(1:16, 18, 17, 19:51, 53:56, 52, 57:70)

SMC_DS_list <- cbind(sort(unique(cases_cleaned$District))[DS_order_in_SMC], sort(unique(SMC$DS_Name)))
setDistrict <- function(Row)
{
    district <- SMC_DS_list[which(SMC_DS_list[,1] == Row[3]), 2]
    
    return(district)
}
cases_cleaned$District <- apply(cases_cleaned, 1, setDistrict)


# adding SMC schedule
schedule <- list(c(8:11),
                 c(7:10),    # maybe should be same as 2015?
                 c(7:10),
                 c(7:10))
Years <- c(2015, 2016, 2017, 2018)

cases_new <- cbind(cases_cleaned, "SMC received" = rep(0,nrow(cases_cleaned)),
                   "Number of children treated with SMC" = rep(0, nrow(cases_cleaned)),
                   "SMC coverage" = rep(0, nrow(cases_cleaned)))


for (i in 1:4)
{
    s <- schedule[[i]]
    Y <- Years[i]
    
    for (R in unique(cases_new$District))
    {
        #SMC[which(SMC$DS_Name == R & SMC$year == Y), "coverage_high_access_U5"]
        Rows <- SMC[which(SMC$DS_Name == R & SMC$year == Y),]
        smc_RY <- as.numeric(Rows[,"children.treated"])
        smc_cov_RY <- Rows[,"coverage_high_access_U5"]*Rows["high_access_U5"] +
                                     Rows[,"coverage_low_access_U5"]*(1 - Rows["high_access_U5"])
        
        if (sum(smc_RY))
        {
            for (m in 1:length(smc_RY))
            {
                M <- s[m]
                
                cases_new[which(cases_new$District == R &
                                    cases_new$year == Y &
                                    cases_new$month == M), 
                          "SMC received"] <- 1
                cases_new[which(cases_new$District == R &
                                    cases_new$year == Y &
                                    cases_new$month == M),
                          "Number of children treated with SMC"] <- smc_RY[m]
                cases_new[which(cases_new$District == R &
                                    cases_new$year == Y &
                                    cases_new$month == M),
                          "SMC coverage"] <- as.numeric(smc_cov_RY[m, 1])
            }
        }
    }
}



cases_new$District <- iconv(cases_new$District, to="ASCII//TRANSLIT")




# Adding climate
climate <- read.csv("~/Box/NU-malaria-team/data/africa_health_district_climate/climate/Burkina 2013-18/Burkina_long_format.csv", header = TRUE, strip.white = TRUE, stringsAsFactors = FALSE)
names(climate)[1] <- "District"
climate_DS_order <- c(1:40, 42, 41, 43:70)
climate_DS_list <- cbind(sort(unique(cases_new$District)), sort(unique(climate$District))[climate_DS_order])
setDistrict <- function(Row)
{
    district <- climate_DS_list[which(climate_DS_list[,2] == Row[1]), 1]
    
    return(district)
}
climate$District <- apply(climate, 1, setDistrict)

cases_new <- inner_join(cases_new, climate[,-2], by=c("District", "year", "month"))



# fixing names for future so we can map onto shapefile
burkina_shape <- readOGR("~/Box/NU-malaria-team/data/burkina_shapefiles/Burkina Faso Health Districts SHP (130715)/BFA.shp")
DS_order_in_shape <- c(1:40, 42, 41, 43:70)
shape_DS_list <- cbind(sort(unique(cases_new$District))[DS_order_in_shape], levels(burkina_shape$District))
setDistrict <- function(Row)
{
    district <- shape_DS_list[which(shape_DS_list[,1] == Row[3]), 2]
    
    return(district)
}
cases_new$District <- apply(cases_new, 1, setDistrict)


# fixing orgunit names
Encoding(cases_new$orgunitlevel3) <- "UTF-8"
cases_new$orgunitlevel3 <- gsub("  ", " ", cases_new$orgunitlevel3)
cases_new$orgunitlevel3 <- iconv(cases_new$orgunitlevel3, from="UTF-8", to="ASCII//TRANSLIT")
cases_new$orgunitlevel3 <- tolower(cases_new$orgunitlevel3)

Encoding(cases_new$orgunitlevel5) <- "UTF-8"
cases_new$orgunitlevel5 <- gsub("  ", " ", cases_new$orgunitlevel5)
cases_new$orgunitlevel5 <- iconv(cases_new$orgunitlevel5, from="UTF-8", to="ASCII//TRANSLIT")
cases_new$orgunitlevel5 <- tolower(cases_new$orgunitlevel5)


# fixing HF names
Encoding(cases_new$HF) <- "UTF-8"
cases_new$HF <- gsub("  ", " ", cases_new$HF)
cases_new$HF <- iconv(cases_new$HF, from="UTF-8", to="ASCII//TRANSLIT")
cases_new$HF <- tolower(cases_new$HF)





cases_new$Date <- as.yearmon(paste(cases_new$year, cases_new$month, sep="-"))

# Adding population
pop_data <- read.csv("~/Box/NU-malaria-team/projects/smc_impact/data/outputs/cases_seasonal_smc.csv", header = TRUE, strip.white = TRUE, stringsAsFactors = FALSE)
# This dataset has a few different names for Districts
pop_data[which(pop_data$Admin2 == "Bittou"), "Admin2"] <- "Bitou"
pop_data[which(pop_data$Admin2 == "Gorom-Gorom"), "Admin2"] <- "Gorom"
pop_data[which(pop_data$Admin2 == "Vigue"), "Admin2"] <- "Karangasso Vigue"
pop_data[which(pop_data$Admin2 == "N'Dorola"), "Admin2"] <- "Ndorola"
pop_data[which(pop_data$Admin2 == "Sig-Noghin"), "Admin2"] <- "Sig-noghin"
names(pop_data)[1:4] <- c("Region", "District", "year", "month")

cases_new <- inner_join(cases_new, pop_data[1:5], by=c("Region", "District", "year", "month"))
cases_new <- cases_new[,c(1:8, 166:167, 9:165)]
names(cases_new)[which(names(cases_new) == "Population")] <- "District Pop"


cases_new$air_temp_era5 <- as.numeric(cases_new$air_temp_era5)
cases_new$precip_era5 <- as.numeric(cases_new$precip_era5)

# saving
write.csv(cases_new, "~/Users/sebas/Box/NU-malaria-team/projects/smc_impact/data/outputs/U5_HF_cases_seasonal_smc.csv", row.names = FALSE)








