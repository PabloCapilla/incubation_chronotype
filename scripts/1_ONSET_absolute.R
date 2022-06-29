##
##### libraries #####
##
pacman::p_load(lubridate, lme4, MuMIn, dplyr, ggplot2, ggridges, extrafont, rptR)
loadfonts()


##
##### functions #####
##
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}



##
##
##### Data #####
##
##
data <- read.csv("./data/Full_incubation_reproduction_data_ROBYN_&_CIARA.csv")
data$box <- substr(as.character(data$box),1,nchar(as.character(data$box))-5) # correct format for box ID
head(data)
data <- left_join(x = data, 
                  y = data %>% 
                    group_by(year, box) %>% 
                    summarise(mean_days_april = mean(date_aprildays, na.rm = T)), 
                  by = c("year", "box")) %>% 
  filter(day_before_hatch < 0) # very little data on the day of hatching

data$inc_stage <- ifelse(data$day_before_hatch <= -10, "1_early", 
                         ifelse(data$day_before_hatch <= -5, "2_mid", "3_late"))

head(data)
table(data$area)


##
##
##### sunrise calculation on the inc start date #####
##
##

library(suncalc)
## site coordinates
# table with coordinates for the different sites
coor_df <- data.frame(site = c("City", "Forest"),
                      lat = c(55.868230, 56.130595),
                      lon = c(-4.282496, -4.614817))


# new columns
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
plot(sunrise_inc_start_dec ~ inc_start_aprildays, data = data)


##
##
##### Samnple sizes #####
##
##
head(data)
nrow(data) # total data points
length(unique(data$box)) # number of nest-boxes

# number of clutches
data %>% 
  group_by(year, box) %>% 
  filter(row_number() == 1) %>% 
  summarise(n_obs = n()) %>% 
  nrow()

# days of incubation per clutch
data %>% 
  group_by(box,year) %>% 
  summarise(days_box = n()) %>%
  ungroup(box) %>% 
  summarise(min_range = min(days_box),
            max_range = max(days_box),
            med = median(days_box))


# number of boxes per site (not clutches)
data %>% 
  group_by(site, box) %>% 
  filter(row_number() == 1) %>% 
  group_by(site) %>% 
  summarise(n_boxes = n())


summary(data$diff_sunrise)
mean(data$diff_sunrise)
sd(data$diff_sunrise)/sqrt(nrow(data))

# start of incubation
mean(data$inc_start_aprildays[data$area == "City"])
sd(data$inc_start_aprildays[data$area == "City"])/sqrt(nrow(data[data$area == "City",]))
mean(data$inc_start_aprildays[data$area == "Forest"])
sd(data$inc_start_aprildays[data$area == "Forest"])/sqrt(nrow(data[data$area == "Forest",]))

data %>% 
  group_by(year, area, box) %>% 
  summarise(mean_box = mean(inc_start_aprildays)) %>% 
  group_by(area) %>% 
  summarise(mean_date = mean(mean_box),
            n_size = n(),
            sd_date = sd(mean_box/sqrt(n())))
            


data %>% 
  group_by(year, area) %>% 
  summarise(n_obs = n())

data %>% 
  group_by(box) %>% 
  filter(row_number() == 1) %>% 
  group_by(year, area) %>% 
  summarise(n_obs = n())

##
##
##### initial visualisation & data summary #####
##
##
ggplot(data = data, aes(x = inc_start_aprildays, y = activity_onset.1)) +
  geom_point() +
  theme_bw() +
  facet_grid(~site) +
  stat_smooth(method = "lm", formula = y~poly(x,2))

ggplot(data = data, aes(x = day_before_hatch, y = activity_onset.1)) +
  geom_point() +
  theme_bw() +
  facet_grid(~area) +
  stat_smooth(method = "lm", formula = y~poly(x,2))


ggplot(data = data, aes(y = activity_onset.1, x = date_aprildays, group = box)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  facet_grid(~area) +
  stat_smooth(method = "lm", formula = y~poly(x,1), se = F)
##
##
##### Model for Onset of activity - absolute time #####
##
##
summary(data$activity_onset.1)
# model fit
model_onset <- lmer(activity_onset.1 ~ 
                      area:poly(inc_start_aprildays,2)[,2] + 
                      area:poly(inc_start_aprildays,2)[,1] +  
                      
                      area:poly(day_before_hatch,2)[,2] + 
                      area:poly(day_before_hatch,2)[,1] +  
                      
                      poly(inc_start_aprildays,2)[,2] +
                      poly(inc_start_aprildays,2)[,1] +
                      poly(day_before_hatch,2)[,2] +
                      poly(day_before_hatch,2)[,1] +
                      area + 
                      meantemp +
                      clutch_size + 
                      (1|site) +
                      (1|box) +
                      (1|year), 
                    REML = F,
                    na.action = "na.fail",
                    data= data) #full model
summary(model_onset)

# model diagnostics
plot(model_onset)
hist(residuals(model_onset))
shapiro.test(residuals(model_onset))


lmerTest::rand(model_onset)
#repeatability
as.numeric(summary(model_onset)$varcor[1]) / 
  (as.numeric(summary(model_onset)$varcor[1]) + 
     summary(model_onset)$sigma^2)

##
##
##### Model selection LRT #####
##
##
drop1(model_onset, test = "Chisq")
m1 <- update(model_onset, . ~ . - clutch_size)

#1
drop1(m1, test = "Chisq")
m2 <- update(m1, . ~ . - poly(inc_start_aprildays, 2)[, 2]:area)

#2
drop1(m2, test = "Chisq")
m3 <- update(m2, . ~ . - poly(inc_start_aprildays, 2)[, 2])

#3
drop1(m3, test = "Chisq")
m4 <- update(m3, . ~ . - poly(day_before_hatch, 2)[, 2]:area)

#4
drop1(m4, test = "Chisq")
m5 <- update(m4, . ~ . - meantemp)

#5
drop1(m5, test = "Chisq")


#mam
m_mam <- m5
drop1(m_mam, test = "Chisq")
summary(m_mam)

CI_onset <- confint(m_mam, 
                    level = 0.95, 
                    method = "boot", 
                    nsim = 500, 
                    boot.type = "norm")


# p for non significant terms
m_clutch_size <- update(m_mam, . ~ . + clutch_size)
anova(m_clutch_size, m_mam, test = "Chisq")

m_int <- update(m_mam, . ~ . + area:poly(day_before_hatch, 2)[, 2])
anova(m_int, m_mam, test = "Chisq")

m_aprildays2 <- update(m_mam, . ~ . + poly(inc_start_aprildays, 2)[, 2])
m_aprildays3 <- update(m_mam, . ~ . + poly(inc_start_aprildays, 2)[, 2] + poly(inc_start_aprildays, 2)[, 2]:area)
anova(m_aprildays2, m_mam, test = "Chisq")
anova(m_aprildays3, m_aprildays2, test = "Chisq")

m_meantemp <- update(m_mam, . ~ . + meantemp)
anova(m_meantemp, m_mam, test = "Chisq")



##
## Top model
top_model_onset <- lmer(activity_onset.1 ~ 
                          area:inc_start_aprildays + 
                          area:day_before_hatch +
                          
                          I(day_before_hatch^2) +
                          day_before_hatch +
                          inc_start_aprildays +
                          area + 
                          (1|site) +
                          (1|box) +
                          (1|year),
                        REML = T,
                        na.action = "na.fail",
                        data=data)
summary(top_model_onset)
drop1(top_model_onset, test = "Chisq")

lmerTest::rand(top_model_onset)
#repeatability
as.numeric(summary(top_model_onset)$varcor[1]) / 
  ((as.numeric(summary(top_model_onset)$varcor[1]) + 
     (as.numeric(summary(top_model_onset)$varcor[2])) +
        (as.numeric(summary(top_model_onset)$varcor[3])) +
     summary(top_model_onset)$sigma^2))


rep_absolute_onset <- rpt(activity_onset.1 ~ 
             area:inc_start_aprildays + 
             area:day_before_hatch +
             
             I(day_before_hatch^2) +
             day_before_hatch +
             inc_start_aprildays +
             area + 
             (1|site) +
             (1|box) +
             (1|year), 
           grname = "box", 
           data = data, 
           datatype = "Gaussian", 
           nboot = 1000, 
           npermut = 0)

##
## predictions plot
df_pred <- expand.grid(day_before_hatch = seq(min(data$day_before_hatch), 
                                              max(data$day_before_hatch), 1),
                       area = c("City", "Forest"),
                       inc_start_aprildays = seq(min(data$inc_start_aprildays), 
                                             max(data$inc_start_aprildays), 1))
df_pred$prediction <- predict(top_model_onset, df_pred, re.form = NA)

# plot only data in range
data %>% 
  group_by(area) %>% 
  summarise(min_chr = min(inc_start_aprildays),
            max_chr = max(inc_start_aprildays))

remove_city <- which((df_pred$inc_start_aprildays < 20 | 
                        df_pred$inc_start_aprildays > 44) & 
                       df_pred$area == "City")
remove_forest <- which((df_pred$inc_start_aprildays < 28 | 
                          df_pred$inc_start_aprildays > 53) & 
                         df_pred$area == "Forest")
df_pred <- df_pred[-c(remove_city, remove_forest),]

# need to adapt code to new R version
# SE for mean predictions
mm <- model.matrix(~ area:inc_start_aprildays + 
                     area:day_before_hatch +
                     
                     I(day_before_hatch^2) +
                     day_before_hatch +
                     inc_start_aprildays +
                     area,
                   data = df_pred)

pvar1 <- diag(mm %*% tcrossprod(vcov(top_model_onset),mm))
cmult <- 1 ## 1 SE
df_pred <- data.frame(
  df_pred
  , plo = df_pred$prediction-cmult*sqrt(pvar1)
  , phi = df_pred$prediction+cmult*sqrt(pvar1)
)


library(chron) 
min_time <- min(times(strftime(dmy_hm(data$first_offbout),"%H:%M:%S", tz = "UTC")))
max_time <- max(times(strftime(dmy_hm(data$first_offbout),"%H:%M:%S", tz = "UTC")))

hm <- merge(0:23, seq(0, 30, by = 30))

interval_15min <- data.frame('INTERVAL' = chron(time = paste(hm$x, ':', hm$y, ':', 0)))
interval_15min <- data.frame('INTERVAL' = interval_15min[order(interval_15min$INTERVAL), ])
labels_time <- substr(interval_15min[!interval_15min$INTERVAL <= min_time & !interval_15min$INTERVAL >= max_time,],
                      start = 0, 
                      stop = 5)


## plot hatching
onset_hatching <- ggplot(data = data, 
                         aes(x = day_before_hatch, 
                             y = activity_onset.1,
                             fill = area,
                             color = area)) +
  geom_point(position = position_jitterdodge(dodge.width = 0.5,
                                             jitter.width = 0.25),
             alpha = 0.4,
             size = 1.25,
             shape = 21, 
             color = "black") +
  theme_bw() +
  #stat_smooth(method = "lm", formula = y ~ poly(x,2)) +
  theme(legend.position = "top",
        legend.text = element_text("Arial", size = 10),
        panel.grid = element_blank(),
        axis.title = element_text("Arial", size = 10),
        axis.text = element_text("Arial", size = 10)) +
  geom_errorbar(data = df_pred %>% 
                  group_by(area, day_before_hatch) %>% 
                  summarise(mean_plo = mean(plo),
                            mean_phi = mean(phi),
                            mean_pred = mean(prediction)), 
                aes(ymin = mean_plo, 
                    ymax = mean_phi, 
                    y = mean_pred),
                width = 0,
                size = 1.5,
                color = "black",
                position = position_dodge(width = 0.5)) +
  geom_point(data = df_pred %>% 
               group_by(area, day_before_hatch) %>% 
               summarise(mean_plo = mean(plo),
                         mean_phi = mean(phi),
                         mean_pred = mean(prediction)), 
             aes(y = mean_pred), 
             size = 3, 
             shape = 21,
             color = "black",
             position = position_dodge(width = 0.5)) +
  labs(x = "Days before hatching", 
       y = "Absolute onset of activity") +
  scale_x_continuous(breaks = -15:-1, labels = 15:1) +
  scale_y_continuous(breaks = seq(240,420, 30),
                     labels = labels_time) +
  scale_fill_manual(name = "", labels = c("City", "Forest"), 
                    values = c("#af8dc3", "#7fbf7b")) +
  scale_color_manual(name = "", labels = c("City", "Forest"), 
                     values = c("#af8dc3", "#7fbf7b"))


ggsave(filename = "./plots/1a_Onset of activity/ABS_onset_pred_hatching.jpeg", 
       plot = onset_hatching, 
       device = "jpeg", 
       units = "mm",
       width = 125, 
       height = 125)  


# date from April 1

# checking sunrise times
x <- (data$activity_onset.1/60) - (data$diff_sunrise/60) 
y <- data$sunrise_inc_start_dec

plot(x, y)
abline(b = 1, a = 0)



onset_plot_abs_date <- ggplot(data = data, 
                              aes(x = inc_start_aprildays, 
                                  y = activity_onset.1,
                                  fill = area,
                                  color = area)) +
  geom_line(data = data,
            aes(y = (sunrise_inc_start_dec)*60, color = "black"),
            size = 1,
            alpha = 0.75,
            linetype = 2) +
  geom_point(alpha = 0.25,
             size = 1.25,
             shape = 21, 
             color = "black") +
  theme_bw() +
  facet_wrap(~area) +
  #stat_smooth(method = "lm", formula = y ~ poly(x,2)) +
  theme(legend.position = "none",
        legend.text = element_text("Arial", size = 10),
        strip.background = element_blank(),
        strip.text = element_text("Arial", size = 12),
        panel.grid = element_blank(),
        axis.title = element_text("Arial", size = 10),
        axis.text = element_text("Arial", size = 10)) +
  geom_ribbon(data = df_pred %>% 
                group_by(area, inc_start_aprildays) %>% 
                summarise(mean_plo = mean(plo),
                          mean_phi = mean(phi),
                          mean_pred = mean(prediction)), 
              aes(ymin = mean_plo, 
                  ymax = mean_phi, 
                  y = mean_pred),
              color = NA,
              alpha = 0.5) +
  geom_line(data = df_pred %>% 
              group_by(area, inc_start_aprildays) %>% 
              summarise(mean_plo = mean(plo),
                        mean_phi = mean(phi),
                        mean_pred = mean(prediction)), 
            aes(y = mean_pred), 
            size = 1.5) +
  labs(x = "Incubation start date (days after April 1)", 
       y = "Time of activity onset") +
  scale_y_continuous(breaks = seq(240,420, 30),
                   #  trans = "reverse",
                     labels = labels_time) + 
  scale_fill_manual(name = "", labels = c("City", "Forest"), 
                    values = c("#af8dc3", "#7fbf7b")) +
  scale_color_manual(name = "", labels = c(NA, "City", "Forest"), 
                     values = c("black", "#af8dc3", "#7fbf7b"))


ggsave(filename = "./plots/1a_Onset of activity/ABS_onset_plot_abs_date_pred.jpeg", 
       plot = onset_plot_abs_date, 
       device = "jpeg", 
       units = "mm",
       width = 135, 
       height = 100)  




