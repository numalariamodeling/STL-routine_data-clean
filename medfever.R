

rm(list = ls(all = TRUE))


require("ggplot2")
require("plyr")
require("dplyr")




med_2014 <- read.csv("~/Box/NU-malaria-team/projects/hbhi_burkina/DS DHS estimates/medfever/med_BF14DS.csv", stringsAsFactors = FALSE)
med_2017 <- read.csv("~/Box/NU-malaria-team/projects/hbhi_burkina/DS DHS estimates/medfever/med_BF17DS.csv", stringsAsFactors = FALSE)

med_2014 <- med_2014[,-c(1,4)]
names(med_2014)[2] <- "medfever_2014"
med_2017 <- med_2017[,-c(1,4)]
names(med_2017)[2] <- "medfever_2017"

medfever_DHS <- inner_join(med_2014, med_2017, by = "NOMDEP")

lm.fit <- lm(medfever_2014 ~ medfever_2017, data = medfever_DHS)




medfever_DHS <- medfever_DHS[which(!is.na(medfever_DHS$medfever_2014) &
                                       !is.na(medfever_DHS$medfever_2017)),]


mean((medfever_DHS$medfever_2017 - medfever_DHS$medfever_2014) / medfever_DHS$medfever_2014)
