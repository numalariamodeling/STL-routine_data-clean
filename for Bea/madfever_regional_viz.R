
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


cases <- read.csv("~/Box/NU-malaria-team/projects/smc_impact/data/outputs/U5_DS_cases_seasonal_smc_good_rows_MA_imputes_testing.csv", header  = TRUE, strip.white = TRUE, stringsAsFactors = FALSE)
cases$U5_pop <- cases$District.Pop * .18
cases$Date <- as.yearmon(cases$Date)



med_2014 <- read.csv("~/Box/NU-malaria-team/projects/hbhi_burkina/DS DHS estimates/medfever/med_BF14DS.csv", stringsAsFactors = FALSE)
med_2017 <- read.csv("~/Box/NU-malaria-team/projects/hbhi_burkina/DS DHS estimates/medfever/med_BF17DS.csv", stringsAsFactors = FALSE)



med_2014 <- med_2014[,-c(1,4)]
names(med_2014)[2] <- "medfever_2014"
med_2017 <- med_2017[,-c(1,4)]
names(med_2017)[2] <- "medfever_2017"

medfever_DHS <- inner_join(med_2014, med_2017, by = "NOMDEP")
names(medfever_DHS) <- c("District", "2015", "2017")

medfever_DHS <- melt(medfever_DHS, id.vars = "District", variable.name = "year", value.name = "medfever")


district_map <- cbind(sort(unique(cases$District)),
                      sort(unique(medfever_DHS$District))[c(1:40,42,41,43:70)])

for (i in 1:nrow(district_map))
{
    DS <- district_map[i, 2]
    new_name_DS <- district_map[i, 1]
    
    medfever_DHS[which(medfever_DHS$District == DS), "District"] <- new_name_DS
}



medfever_DHS$year <- as.numeric(as.character(medfever_DHS$year))

medfever_DHS <- left_join(medfever_DHS, unique(cases[,c("Region", "District", "year", "U5_pop")]), by = c("District", "year"))


medfever_region <- ddply(medfever_DHS[!is.na(medfever_DHS$medfever),], c(.(Region), .(year)),
                         summarise, medfever = sum(medfever * U5_pop)/ sum(U5_pop))


medfever_DHS_tmp_1 <- medfever_DHS[!is.na(medfever_DHS$medfever),]
medfever_DHS_tmp_2 <- medfever_DHS[is.na(medfever_DHS$medfever),]

medfever_DHS_tmp_2 <- left_join(medfever_DHS_tmp_2[,-3], medfever_region, by = c("year", "Region"))
medfever_DHS_tmp_2 <- medfever_DHS_tmp_2[, c(1:2, 5, 3:4)]

medfever_DHS <- rbind(medfever_DHS_tmp_1, medfever_DHS_tmp_2)

names(medfever_region)[3] <- "medfever_regional"
medfever_DHS <- left_join(medfever_DHS, medfever_region, by = c("year", "Region"))



medfever_DHS[which(medfever_DHS$year == 2015), "year"] <- 2014



medfever_DHS <- medfever_DHS[order(medfever_DHS$District, medfever_DHS$year),]

diff_medfever <- medfever_DHS[which(medfever_DHS$year == 2017), "medfever"] - medfever_DHS[which(medfever_DHS$year == 2014), "medfever"]
diff_medfever_regional <- medfever_DHS[which(medfever_DHS$year == 2017), "medfever_regional"] - medfever_DHS[which(medfever_DHS$year == 2014), "medfever_regional"]

medfever_diff <- medfever_DHS[which(medfever_DHS$year == 2017), ]
medfever_diff$year <- "diff"
medfever_diff$medfever <- diff_medfever
medfever_diff$medfever_regional <- diff_medfever_regional


# medfever_DHS <- rbind(medfever_DHS, medfever_diff)





#########################################################################################









burkina_shape <- readOGR("~/Box/NU-malaria-team/data/burkina_shapefiles/Burkina Faso Health Districts SHP (130715)/BFA.shp")
burkina_shape_DF <- fortify(burkina_shape, region = "District")




district_list <- cbind(sort(unique(medfever_DHS$District)), sort(unique(burkina_shape_DF$id)))
setDistrict <- function(Row)
{
    district <- district_list[which(district_list[,1] == Row[1]), 2]
    
    return(district)
}
medfever_DHS$id <- apply(medfever_DHS, 1, setDistrict)




burkina_shape_DF_medfever <- inner_join(burkina_shape_DF, medfever_DHS, by = "id")

colr <- brewer.pal(9, "RdYlBu")
# colr <- brewer.pal(9, "BuGn")
# colr <- brewer.pal(9, "Greens")


facet_medfever_p <- ggplot(data = burkina_shape_DF_medfever, aes(x = long, y = lat, group = group)) +
    ggtitle("Regional febrile treatment seeking rate") +
    geom_polygon(aes(fill = medfever_regional), color = "white") + 
    coord_equal() + scale_fill_gradientn("Treatment seeking rate", colors = colr) + theme_void() +
    theme(plot.title = element_text(hjust = 0.5)) + facet_wrap(~year)








district_list <- cbind(sort(unique(medfever_diff$District)), sort(unique(burkina_shape_DF$id)))
setDistrict <- function(Row)
{
    district <- district_list[which(district_list[,1] == Row[1]), 2]
    
    return(district)
}
medfever_diff$id <- apply(medfever_diff, 1, setDistrict)




burkina_shape_DF_medfever <- inner_join(burkina_shape_DF, medfever_diff, by = "id")


medfever_diff_p <- ggplot(data = burkina_shape_DF_medfever, aes(x = long, y = lat, group = group)) +
    ggtitle("Difference in regional febrile treatment seeking rate") +
    geom_polygon(aes(fill = medfever_regional), color = "white") + 
    coord_equal() + scale_fill_gradientn("Difference", colors = colr) + theme_void() +
    theme(plot.title = element_text(hjust = 0.5))


