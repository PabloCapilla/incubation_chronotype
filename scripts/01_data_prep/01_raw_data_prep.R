###
###
#' 
#' Script for:
#' Reproductive fitness is associated with female chronotype in a songbird
#' Womack, et al. 
#' Preprint: 
#' 
#' Latest update: 2022-06-28
#' 
###
###

# Clear memory to make sure there are not files loaded that could cause problems
rm(list=ls())

##
##
##### Script aim: #####
#' This script prepares the dataset for analysis
#' 
##
##
##

##
##
##### libraries #####
##
##
pacman::p_load(openxlsx, 
               lubridate, dplyr, tidyr, suncalc,
               lme4, 
               ggplot2, extrafont)
loadfonts()

##
##
##### additional functions #####
##
##
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

##
##
##### data #####
##
##
data <- read.csv("./data/01_raw_data/Full_incubation_reproduction_dataset.csv")
data <- data %>% 
  filter(day_before_hatch < 0) # remove the day of hatching (2 observations)
head(data)

##
##
##### sunrise calculation on the inc start date #####
##
##

## site coordinates
# table with coordinates for the different sites
coor_df <- data.frame(site = c("City", "Forest"),
                      lat = c(55.868230, 56.130595),
                      lon = c(-4.282496, -4.614817))


# new columns for incubation start times
data$sunrise_inc_start <- ymd_hms(NA)
data$sunset_inc_start <- ymd_hms(NA)
data$dawn_inc_start <- ymd_hms(NA)
data$dusk_inc_start <- ymd_hms(NA)


for(i in 1:nrow(data)){
  #add column with dawn and dusk
  suntimes  <- getSunlightTimes(date = ymd(paste0(data$year[i], "-03-31"))+data$inc_start_aprildays[i], 
                                lat = ifelse(data$area[i] == "City", coor_df$lat[1], coor_df$lat[2]),
                                lon = ifelse(data$area[i] == "City", coor_df$lon[1], coor_df$lon[2]),
                                keep = c("dawn", "dusk", "sunrise", "sunset"), 
                                tz = "Europe/London")
  data$sunrise_inc_start[i] <- ymd_hms(suntimes$sunrise)
  data$sunset_inc_start[i] <- ymd_hms(suntimes$sunset)
  data$dawn_inc_start[i] <- ymd_hms(suntimes$dawn)
  data$dusk_inc_start[i] <- ymd_hms(suntimes$dusk)
}
data$sunrise_inc_start_dec <- hour(ymd_hms(data$sunrise_inc_start)) + (minute(ymd_hms(data$sunrise_inc_start))/60)
data$sunset_inc_start_dec <- hour(ymd_hms(data$sunset_inc_start)) + (minute(ymd_hms(data$sunset_inc_start))/60)
  

##
##
##### Save data #####
##
##

# data for incubation analyses
saveRDS(object = data, file = "./data/data_incubation.RDS")
