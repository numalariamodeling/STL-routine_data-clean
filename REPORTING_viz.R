

rm(list = ls(all = TRUE))

require("plyr")
require("dplyr")
require("zoo")
require("reshape2")
require("ggplot2")
require("lubridate")

library("RColorBrewer")
library("spdep")
library("rgdal")
library("sf")
library("rgeos")


# Loading health facility dataset

HF_cases <- read.csv("~/Box/NU-malaria-team/projects/smc_impact/data/outputs/U5_HF_cases_smc_coords_imputed_rdts_and_allout_MA.csv", stringsAsFactors = FALSE)

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
bad_rows_missing_rdt_only <- bad_rdt_rows[which(!(bad_rdt_rows %in% bad_rows_inconsistent) &
                                                    !(bad_rdt_rows %in% bad_rows_allout))]
bad_rows_missing_both <- bad_rows_missing[which(!(bad_rows_missing %in% bad_rows_missing_allout_only) &
                                                    !(bad_rows_missing %in% bad_rows_missing_rdt_only))]
bad_rows_inconsistent_only <- bad_rows_inconsistent[which(!(bad_rows_inconsistent %in% bad_rows_missing))]


HF_cases_area <- HF_cases[,c("District", "UID", "Date")]

HF_cases_area$bad_rows <- "good_report"
HF_cases_area[bad_rows_missing_allout_only, "bad_rows"] <- "missing all-cause outpatient visits"
HF_cases_area[bad_rows_missing_rdt_only, "bad_rows"] <- "missing confirmed malaria cases"
HF_cases_area[bad_rows_missing_both, "bad_rows"] <- "missing both values"
HF_cases_area[bad_rows_inconsistent_only, "bad_rows"] <- "inconsistent report"


HF_cases_area_casted <- dcast(HF_cases_area, District + Date ~ bad_rows)
HF_cases_area_casted$total_HFs <- rowSums(HF_cases_area_casted[,3:7])
# 
# HF_cases_area_by_District <- melt(HF_cases_area_casted, id.vars = c("District", "Date", "total_HFs"),
#                                   variable.name = "Type", value.name = "reports")
# HF_cases_area_by_District$Type <- factor(HF_cases_area_by_District$Type,
#                                          levels = c("inconsistent report",
#                                                     "missing both values",
#                                                     "missing confirmed malaria cases",
#                                                     "missing all-cause outpatient visits",
#                                                     "good report"))




#########################################################################################






burkina_shape <- readOGR("~/Box/NU-malaria-team/data/burkina_shapefiles/Burkina Faso Health Districts SHP (130715)/BFA.shp")
burkina_shape_DF <- fortify(burkina_shape, region = "District")




district_list <- cbind(sort(unique(HF_cases_area_casted$District)), sort(unique(burkina_shape_DF$id)))
setDistrict <- function(Row)
{
    district <- district_list[which(district_list[,1] == Row[1]), 2]
    
    return(district)
}
HF_cases_area_casted$id <- apply(HF_cases_area_casted, 1, setDistrict)




HF_cases_reporting_yearly <- ddply(HF_cases_area_casted, c(.(id), .(year(Date))), summarize,
                                   good_report = sum(good_report),
                                   total_HFs = sum(total_HFs))
names(HF_cases_reporting_yearly)[2] <- "year"
HF_cases_reporting_yearly$reporting_percent <- HF_cases_reporting_yearly$good_report / HF_cases_reporting_yearly$total_HFs



burkina_shape_DF_yearly <- inner_join(burkina_shape_DF, HF_cases_reporting_yearly, by = "id")

colr <- brewer.pal(9, "RdYlBu")
# colr <- brewer.pal(9, "BuGn")
# colr <- brewer.pal(9, "Greens")


facet_reporting_p <- ggplot(data = burkina_shape_DF_yearly, aes(x = long, y = lat, group = group)) +
    ggtitle("Yearly percent reporting") +
    geom_polygon(aes(fill = reporting_percent), color = "white") + 
    coord_equal() + scale_fill_gradientn("Percent reporting", colors = colr) + theme_void() +
    theme(plot.title = element_text(hjust = 0.5)) + facet_wrap(~year)










# seasons <- data.frame(Date = sort(unique(HF_cases_area_casted$Date)),
#                       season = c(rep("Dry cool 2015", 1),
#                                  rep("Dry hot 2015", 4),
#                                  rep("Rainy 2015", 5),
#                                  rep("Dry cool 2015", 2),
#                                  rep("Dry cool 2016", 1),
#                                  rep("Dry hot 2016", 4),
#                                  rep("Rainy 2016", 5),
#                                  rep("Dry cool 2016", 2),
#                                  rep("Dry cool 2017", 1),
#                                  rep("Dry hot 2017", 4),
#                                  rep("Rainy 2017", 5),
#                                  rep("Dry cool 2017", 2),
#                                  rep("Dry cool 2018", 1),
#                                  rep("Dry hot 2018", 4),
#                                  rep("Rainy 2018", 5),
#                                  rep("Dry cool 2018", 2)))
# HF_cases_area_casted <- left_join(HF_cases_area_casted, seasons, by = "Date")

seasons <- data.frame(Date = sort(unique(HF_cases_area_casted$Date))[-1],
                      season = c(rep("Dry hot 2015", 4),
                                 rep("Rainy 2015", 5),
                                 rep("Dry cool 2015", 3),
                                 rep("Dry hot 2016", 4),
                                 rep("Rainy 2016", 5),
                                 rep("Dry cool 2016", 3),
                                 rep("Dry hot 2017", 4),
                                 rep("Rainy 2017", 5),
                                 rep("Dry cool 2017", 3),
                                 rep("Dry hot 2018", 4),
                                 rep("Rainy 2018", 5),
                                 rep("Dry cool 2018", 2)))
HF_cases_area_casted <- inner_join(HF_cases_area_casted, seasons, by = "Date")
HF_cases_area_casted$season <- factor(as.factor(HF_cases_area_casted$season), levels = unique(seasons$season))



HF_cases_reporting_seasonally <- ddply(HF_cases_area_casted, c(.(id), .(season)), summarize,
                                   good_report = sum(good_report),
                                   total_HFs = sum(total_HFs))
names(HF_cases_reporting_seasonally)[2] <- "season"
HF_cases_reporting_seasonally$reporting_percent <- HF_cases_reporting_seasonally$good_report / HF_cases_reporting_seasonally$total_HFs



burkina_shape_DF_seasonal <- inner_join(burkina_shape_DF, HF_cases_reporting_seasonally, by = "id")



facet_reporting_p <- ggplot(data = burkina_shape_DF_seasonal, aes(x = long, y = lat, group = group)) +
    ggtitle("Seasonal percent reporting") +
    geom_polygon(aes(fill = reporting_percent), color = "white") + 
    coord_equal() + scale_fill_gradientn("Percent reporting", colors = colr) + theme_void() +
    theme(plot.title = element_text(hjust = 0.5)) + facet_wrap(~season, ncol = 3)




