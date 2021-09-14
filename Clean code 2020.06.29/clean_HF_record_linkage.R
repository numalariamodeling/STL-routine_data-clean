##############################
##  Record Linkage attempt  ##
##############################
#
#
#   After finding and fixing typos on communes
#   Clean HFs and merge
#
#

rm(list = ls(all = TRUE))


library("plyr")
library("dplyr")
library("RecordLinkage")




cases <- read.csv("~/Box/NU-malaria-team/projects/smc_impact/data/outputs/U5_HF_cases_seasonal_smc.csv", header = TRUE, strip.white = TRUE, stringsAsFactors = FALSE)
names(cases)[c(4,5)] <- c("Commune", "HF.name")
cases$UID <- paste(cases$District, cases$Commune, cases$HF.name)
cases <- cases[,c(1:5,168,6:167)]
cases[,1:6] <- apply(cases[,1:6], 2, tolower)


# Easier to merge
HF_coords <- read.csv("~/Box/NU-malaria-team/projects/smc_impact/data/outputs/BF_HF_coords_v2.csv", header = TRUE, strip.white = TRUE, stringsAsFactors = FALSE)
names(HF_coords)[c(5,7)] <- c("Commune", "HF.name") 
HF_coords$UID <- paste(HF_coords$District, HF_coords$Commune, HF_coords$HF.name)
HF_coords[,1:7] <- apply(HF_coords[,1:7], 2, tolower)






cases_names <- as.data.frame(unique(cases[,c(1,3:6)]))

# Takes care of "csi wend-na-songdo" dupe with "csi wend na songdo" (but this is all NA)
cases_names$HF.name <- gsub("-", " ", cases_names$HF.name)
cases_names$HF.name <- gsub("'", "", cases_names$HF.name)

cases_names$Commune <- gsub("-", " ", cases_names$Commune)
cases_names$Commune <- gsub("'", "", cases_names$Commune)
cases_names$Commune <- gsub("\\(yat\\)", "", cases_names$Commune)
cases_names$Commune <- gsub("\\(\\)", "", cases_names$Commune)
cases_names$Commune <- gsub(" \\(commune\\)", "", cases_names$Commune)
cases_names$Commune <- gsub("mani", "manni", cases_names$Commune)
cases_names$Commune <- gsub("pobemengao", "pobe mengao", cases_names$Commune)
# Mostly sure about these edit
cases_names$Commune <- gsub("arrondissement de ", "", cases_names$Commune)
cases_names$Commune <- gsub("aribinda", "arbinda", cases_names$Commune)
cases_names$Commune <- gsub("kassou", "cassou", cases_names$Commune)
cases_names$Commune <- gsub("falagountou", "falangountou", cases_names$Commune)
cases_names$Commune <- gsub(" c", "", cases_names$Commune)




HF_coords_names <- as.data.frame(HF_coords[,c(1, 3, 5, 7, 10)])

HF_coords_names$Region <- gsub("-", " ", HF_coords_names$Region)
HF_coords_names$HF.name <- gsub("-", " ", HF_coords_names$HF.name)
HF_coords_names$HF.name <- gsub("'", "", HF_coords_names$HF.name)

HF_coords_names$Commune <- gsub("-", " ", HF_coords_names$Commune)
HF_coords_names$Commune <- gsub("_", " ", HF_coords_names$Commune)
HF_coords_names$Commune <- gsub("'", "", HF_coords_names$Commune)
HF_coords_names$Commune <- gsub("signoghin", "sig noghin", HF_coords_names$Commune)
HF_coords_names$Commune <- gsub("rurae de pabre", "rurale de pabre", HF_coords_names$Commune)
HF_coords_names$Commune <- gsub("rural ", "", HF_coords_names$Commune)
# not sure about these one (but mostly)
HF_coords_names$Commune <- gsub("barsalgho", "barsalogho", HF_coords_names$Commune)
HF_coords_names$Commune <- gsub("bongande", "bogande", HF_coords_names$Commune)
HF_coords_names$Commune <- gsub("bouroun", "bouroum", HF_coords_names$Commune)
HF_coords_names$Commune <- gsub("boundry", "boudry", HF_coords_names$Commune)
HF_coords_names$Commune <- gsub("boudokuy", "bondokuy", HF_coords_names$Commune)
HF_coords_names$Commune <- gsub("bondoukuy", "bondokuy", HF_coords_names$Commune)
HF_coords_names$Commune <- gsub("boondokuy", "bondokuy", HF_coords_names$Commune)
HF_coords_names$Commune <- gsub("bomborekui", "bomborokuy", HF_coords_names$Commune)
HF_coords_names$Commune <- gsub("cheriba", "tcheriba", HF_coords_names$Commune)
HF_coords_names$Commune <- gsub("dapelgo", "dapelogo", HF_coords_names$Commune)
HF_coords_names$Commune <- gsub("doubala", "doumbala", HF_coords_names$Commune)
HF_coords_names$Commune <- gsub("fada", "fada ngourma", HF_coords_names$Commune)
HF_coords_names$Commune <- gsub("goromg gorom", "gorom gorom", HF_coords_names$Commune)
HF_coords_names$Commune <- gsub("gombousgou", "gomboussougou", HF_coords_names$Commune)
HF_coords_names$Commune <- gsub("gorgadii", "gorgadji", HF_coords_names$Commune)
# HF_coords_names$Commune <- gsub("goughin", "gounghin", HF_coords_names$Commune)
HF_coords_names$Commune <- gsub("ioloniono", "iolonioro", HF_coords_names$Commune)
HF_coords_names$Commune <- gsub("karangasso", "karankasso", HF_coords_names$Commune)
HF_coords_names$Commune <- gsub("kanchari", "kantchari", HF_coords_names$Commune)
HF_coords_names$Commune <- gsub("kassoum", "cassoum", HF_coords_names$Commune)
HF_coords_names$Commune <- gsub("kossoum", "cassoum", HF_coords_names$Commune)
HF_coords_names$Commune <- gsub("kokologo", "kokologho", HF_coords_names$Commune)
HF_coords_names$Commune <- gsub("koundoungou", "koundougou", HF_coords_names$Commune)
HF_coords_names$Commune <- gsub("lorepeni", "loropeni", HF_coords_names$Commune)
HF_coords_names$Commune <- gsub("manssila", "mansila", HF_coords_names$Commune)
HF_coords_names$Commune <- gsub("namissignima", "namissiguima", HF_coords_names$Commune)
HF_coords_names$Commune <- gsub("niangologo", "niangoloko", HF_coords_names$Commune)
HF_coords_names$Commune <- gsub("niankoudougou", "niankorodougou", HF_coords_names$Commune)
HF_coords_names$Commune <- gsub("nongr massom", "nongremassoum", HF_coords_names$Commune)
HF_coords_names$Commune <- gsub("pompo¤", "pompoi", HF_coords_names$Commune)
HF_coords_names$Commune <- gsub("pouni sud", "pouni", HF_coords_names$Commune)
HF_coords_names$Commune <- gsub("sanga", "sangha", HF_coords_names$Commune)
HF_coords_names$Commune <- gsub("ttcheriba", "tcheriba", HF_coords_names$Commune)
HF_coords_names$Commune <- gsub("wolonkoto", "wolokoto", HF_coords_names$Commune)
HF_coords_names$Commune <- gsub("zan", "zam", HF_coords_names$Commune)
# not sure about these one (really)
HF_coords_names$Commune <- gsub("rurale de ", "", HF_coords_names$Commune)
HF_coords_names$Commune <- gsub("rurale ", "", HF_coords_names$Commune)
HF_coords_names$Commune <- gsub("arrondissement de ", "", HF_coords_names$Commune)
HF_coords_names$Commune <- gsub("de ", "", HF_coords_names$Commune)
HF_coords_names$Commune <- gsub("kayaho", "kayao", HF_coords_names$Commune)

# Strange corner case, projection grabbed Bogande district but matches HF when changed to Koupela
HF_coords_names[which(HF_coords_names$District == "bogande" & HF_coords_names$Commune == "goughin"), "District"] <- "koupela"
HF_coords_names$Commune <- gsub("goughin", "gounghin", HF_coords_names$Commune)





# Manually changing one link but how to merge them?
cases_names$HF.name <- gsub("cabinet de soins infirmiers wend lad yolsda", 
                            "cabinet de soins infirmiers wend la yolsda", cases_names$HF.name)

cases_names <- as.data.frame(unique(cases_names))
HF_coords_names <- as.data.frame(unique(HF_coords_names))





###############
##  Merging  ##
###############
# w/o blocking
# pairs <- compare.linkage(cases_names, HF_coords_names,
#                          exclude = "UID",
#                          strcmp = TRUE)
# 
# em.pairs <- emWeights(pairs)
# summary(em.pairs)



# w/ blocking
pairs <- compare.linkage(cases_names, HF_coords_names,
                         blockfld = c(1, 2, 3), exclude = "UID",
                         strcmp = TRUE)

em.pairs <- emWeights(pairs)
summary(em.pairs)

class.pairs <- emClassify(em.pairs, threshold.upper = 13.62, threshold.lower = 10)
# summary(class.pairs)
# links <- getPairs(class.pairs, show="links")
# View(links)
# possibles <- getPairs(class.pairs, show="possible")
# View(possibles)
# possibles_clean <- possibles[c(82:87, 199:251), -c(1,6,7)]

# save list of possibles
# write.csv(possibles_clean, "/Users/sebas/Documents/GitHub/smc-impact/SMC U5 code/merging and deduping/health_facility_possible_links.csv", row.names = FALSE)




############################
##  Actually merging now  ##
############################
links <- getPairs(class.pairs, show="links", single.row = TRUE)

names(HF_coords)[10] <- "UID.2"
HF_links <- left_join(links, HF_coords, by = "UID.2")


cases$Long <- NA
cases$Lat <- NA


# Doing slow for-loop
for (U in unique(cases$UID))
{
    ind <- which(HF_links$UID.1 == U)
    if(length(ind) != 0) {
        long_lat <- HF_links[ind, c("POINT_X", "POINT_Y")]
    } else {
        long_lat <- c(NA, NA)
    }
    
    # print(U)
    cases[which(cases$UID == U), c("Long", "Lat")] <- long_lat
}


# saving
write.csv(cases, "~/Box/NU-malaria-team/projects/smc_impact/data/outputs/U5_HF_cases_smc_coords.csv", row.names = FALSE)


