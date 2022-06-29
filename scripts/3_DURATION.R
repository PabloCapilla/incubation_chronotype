##
##### libraries #####
##
pacman::p_load(lubridate, lme4, MuMIn, dplyr, ggplot2, ggridges, extrafont)
loadfonts()


##
##### functions #####
##
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}



##
##
##### Data sets #####
##
##
data <- read.csv("./data/Full_incubation_reproduction_data_ROBYN_&_CIARA.csv")
data$box <- substr(as.character(data$box),1,nchar(as.character(data$box))-5) # correct format for box ID
data$last_onbout <- ifelse(data$last_onbout == "", NA, data$last_onbout)
data$activity_end <- hm(format(dmy_hm(data$last_onbout), format = "%H:%M"))
data <- left_join(x = data, 
                  y = data %>% 
                    group_by(year, box) %>% 
                    summarise(mean_days_april = mean(date_aprildays, na.rm = T)), 
                  by = c("year", "box")) %>% 
  filter(day_before_hatch < 0) # very little data on the day of hatching

head(data)

## the observations below were removed from Robyn's and Ciara's analyses. I assume there are 
## reasons for that, so I remove them too
box_remove <- c("27KG", "245SAL","32KG","12CASH","425CASH","425CASH","244SAL", "128PEN", "29KG","151FS")
recording_date_remove <- c("2018-05-14","2018-05-23","2018-05-21",
                           "2018-05-16","2018-05-16","2018-05-23",
                           "2018-05-23","2016-05-18","2017-05-05",
                           "2016-05-13")
remove_df<- data.frame(box = box_remove, recording_date = recording_date_remove)
remove_df$box <- as.character(remove_df$box)
remove_df$recording_date <- as.character(remove_df$recording)


##
## IMPORTANT, REPORT THESE REMOVALS OF DATA AND WHY
# Ciara's code for Robyn's data set
data$recording_date <- as.character(data$recording_date)
data_trimmed <- anti_join(data, remove_df, by=c("box","recording_date"))
data_trimmed<- data_trimmed[!is.na(data_trimmed$diff_sunset),]
data_trimmed<- data_trimmed %>% filter(diff_sunset > -300)



data_trimmed$duration_act <- as.numeric(as.difftime(data_trimmed$activity_end - hms(data_trimmed$activity_onset)),
                                        units = "hours")

# checking new and old calculations of duration
plot(data_trimmed$duration_act, data_trimmed$act_window)



# sample size
data_trimmed %>% 
  group_by(year, area) %>% 
  summarise(n_obs = n())

data_trimmed %>% 
  group_by(box, year, area) %>% 
  filter(row_number() == 1) %>% 
  group_by(year, area) %>% 
  summarise(n_obs = n())


##
##
##### Preliminary plots ####
##
##
ggplot(data = data_trimmed, aes(x = inc_start_aprildays, y = duration_act, color = area)) +
  facet_wrap(~area) +
  geom_point() +
  stat_smooth(method = "lm")

ggplot(data = data_trimmed, aes(x = day_before_hatch, y = duration_act, color = area)) +
  facet_wrap(~area) +
  geom_point() +
  stat_smooth(method = "lm")



##
##
##### Models #####
##
##

# Robyn's model
day_length_model <- lmer(duration_act ~ 
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
                         na.action = "na.fail",
                         REML = F,                 
                         data=data_trimmed) 

summary(day_length_model)

# model diagnostics
plot(day_length_model)
hist(residuals(day_length_model))


##
##
##### Model selection LRT #####
##
##
drop1(day_length_model, test = "Chisq")
m1 <- update(day_length_model, . ~ . - area:poly(inc_start_aprildays, 2)[, 2])

#1
drop1(m1, test = "Chisq")
m2 <- update(m1, . ~ . - poly(inc_start_aprildays, 2)[, 2])

#2
drop1(m2, test = "Chisq")
m3 <- update(m2, . ~ . - poly(day_before_hatch, 2)[, 2]:area)

#3
drop1(m3, test = "Chisq")
m4 <- update(m3, . ~ . - meantemp)

#4
drop1(m4, test = "Chisq")
m5 <- update(m4, . ~ . - clutch_size)

#5
drop1(m5, test = "Chisq")
m6 <- update(m5, . ~ . - poly(day_before_hatch, 2)[, 2])

#6
drop1(m6, test = "Chisq")

#mam
m_mam <- m6
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

m_int <- update(m_mam, . ~ . + poly(day_before_hatch, 2)[, 2])
anova(m_int, m_mam, test = "Chisq")

m_aprildays2 <- update(m_mam, . ~ . + poly(inc_start_aprildays, 2)[, 2])
anova(m_aprildays2, m_mam, test = "Chisq")

m_meantemp <- update(m_mam, . ~ . + meantemp)
anova(m_meantemp, m_mam, test = "Chisq")

m_hatch_int1 <- update(m_mam, . ~ . + poly(day_before_hatch, 2)[, 2] + poly(day_before_hatch, 2)[, 2]:area)
m_hatch_int2 <- update(m_mam, . ~ . + poly(day_before_hatch, 2)[, 2])
anova(m_hatch_int1, m_hatch_int2, test = "Chisq")

m_days_int1 <- update(m_mam, . ~ . + poly(inc_start_aprildays, 2)[, 2] + poly(inc_start_aprildays, 2)[, 2]:area)
m_days_int2 <- update(m_mam, . ~ . + poly(inc_start_aprildays, 2)[, 2])
anova(m_days_int1, m_days_int2, test = "Chisq")




##
##
##### Plot model predictions #####
##
##

# final model
day_length_final <- lmer(duration_act ~ 
                             area:inc_start_aprildays +  
                             area:day_before_hatch +  
                             inc_start_aprildays +
                             day_before_hatch +
                             area + 
                           (1|site) +
                           (1|box) +
                           (1|year),
                       REML = T,                 
                       data=data_trimmed)
r.squaredGLMM(day_length_final)
summary(day_length_final)

##
## repeatability
rep_duration <- rpt(duration_act ~ 
                             area:inc_start_aprildays +  
                             area:day_before_hatch +  
                             inc_start_aprildays +
                             day_before_hatch +
                             area + 
                             (1|site) +
                             (1|box) +
                             (1|year), 
                           grname = "box", 
                           data = data_trimmed, 
                           datatype = "Gaussian", 
                           nboot = 1000, 
                           npermut = 5)


CI_onset <- confint(day_length_final, 
                    level = 0.95, 
                    method = "boot", 
                    nsim = 500, 
                    boot.type = "norm")



df_pred <- expand.grid(day_before_hatch = seq(min(data_trimmed$day_before_hatch), 
                                              max(data_trimmed$day_before_hatch), 1),
                       area = c("City", "Forest"),
                       inc_start_aprildays = seq(min(data_trimmed$inc_start_aprildays), 
                                                 max(data_trimmed$inc_start_aprildays), 1))
df_pred$prediction <- predict(day_length_final, df_pred, re.form = NA)

# need to adapt code to new R version
# SE for mean predictions
mm <- model.matrix(~ area:inc_start_aprildays +  
                     area:day_before_hatch +  
                     inc_start_aprildays +
                     day_before_hatch +
                     area,
                   data = df_pred)

pvar1 <- diag(mm %*% tcrossprod(vcov(day_length_final),mm))
cmult <- 1 ## 1 SE
df_pred <- data.frame(
  df_pred
  , plo = df_pred$prediction-cmult*sqrt(pvar1)
  , phi = df_pred$prediction+cmult*sqrt(pvar1)
)

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



duration_hatch_plot <- ggplot(data = data_trimmed, 
                              aes(x = day_before_hatch, 
                                  y = duration_act,
                                  fill = area,
                                  color = area)) +
  geom_point(position = position_jitterdodge(dodge.width = 0.5,
                                             jitter.width = 0.25),
             alpha = 0.4,
             size = 1.25,
             shape = 21, 
             color = "black") +
  theme_bw() +
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
             size = 3.5, 
             shape = 21,
             color = "black",
             position = position_dodge(width = 0.5)) +
  labs(x = "Days before hatching", 
       y = "Duration of activity (hours)") +
  scale_x_continuous(breaks = -15:-1, labels = 15:1) +
  scale_fill_manual(name = "", labels = c("City", "Forest"), 
                    values = c("#af8dc3", "#7fbf7b")) +
  scale_color_manual(name = "", labels = c("City", "Forest"), 
                     values = c("#af8dc3", "#7fbf7b"))


ggsave(filename = "./plots/duration_hatch.jpeg",
       plot = duration_hatch_plot, 
       device = "jpeg", 
       units = "mm",
       width = 125, 
       height = 125)  


##
## plot for days from April 1
duration_abs_date <- ggplot(data = data_trimmed, 
                              aes(x = inc_start_aprildays, 
                                  y = duration_act,
                                  fill = area,
                                  color = area)) +
  geom_point(alpha = 0.25,
             size = 1.25,
             shape = 21, 
             color = "black") +
  theme_bw() +
  facet_wrap(~area) +
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
       y = "Duration of activity (hours)") +
  scale_fill_manual(name = "", labels = c("City", "Forest"), 
                    values = c("#af8dc3", "#7fbf7b")) +
  scale_color_manual(name = "", labels = c("City", "Forest"), 
                     values = c("#af8dc3", "#7fbf7b")) 


ggsave(filename = "./plots/Duration/duration_april1.jpeg", 
       plot = duration_abs_date, 
       device = "jpeg", 
       units = "mm",
       width = 135, 
       height = 100)  






