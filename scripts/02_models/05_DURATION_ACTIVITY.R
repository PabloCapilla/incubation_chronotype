###
###
#' 
#' Script for:
#' Reproductive fitness is associated with female chronotype in a songbird
#' Womack, et al. 
#' Preprint: 10.1101/2022.07.01.498449v1
#' 
#' Latest update: 2022-07-27
#' 
###
###

# Clear memory to make sure there are not files loaded that could cause problems
rm(list=ls())

##
##
##### Script aim: #####
#' Analysis of duration of activity
#' 
##
##

##
##
##### libraries #####
##
##
pacman::p_load(openxlsx, 
               lubridate, dplyr, tidyr, rptR,
               lme4, performance,
               ggplot2, extrafont)
loadfonts()

##
##
##### data #####
##
##
data00 <- readRDS("./data/data_incubation.RDS")
head(data00)


##
## filtering data for end of activity
##

## days were females were disturbed before end of activity. Therefore, observations removed from dataset
box_remove <- c("27KG", "245SAL","32KG","12CASH","425CASH","425CASH","244SAL", "128PEN", "29KG","151FS")
recording_date_remove <- c("2018-05-14","2018-05-23","2018-05-21",
                           "2018-05-16","2018-05-16","2018-05-23",
                           "2018-05-23","2016-05-18","2017-05-05",
                           "2016-05-13")
remove_df<- data.frame(box = box_remove, recording_date = recording_date_remove)
remove_df$box <- as.character(remove_df$box)
remove_df$recording_date <- as.character(remove_df$recording)

data00$recording_date <- as.character(data00$recording_date)
data <- anti_join(data00, remove_df, by=c("box","recording_date"))
data<- data[!is.na(data$activity_end_relative),]
data<- data %>% filter(activity_end_relative > -300)

## duration of activity in min
data$duration_act <- as.numeric(data$activity_end_abs - data$activity_onset_abs)

##
##
##### Models for duration of activity #####
##
##

# model fit
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
                         data=data) 
summary(day_length_model)

# model diagnostics
check_model(day_length_model, panel = F) # not too bad overall

# significance of random effects 
lmerTest::rand(day_length_model)

# maternal repeatability of relative onset of activity
as.numeric(summary(day_length_model)$varcor[1]) / 
  (as.numeric(summary(day_length_model)$varcor[1]) + 
     as.numeric(summary(day_length_model)$varcor[2]) +
     as.numeric(summary(day_length_model)$varcor[3]) +
     summary(day_length_model)$sigma^2)

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

# Final model
m_mam <- m6
summary(m_mam)


##
##
##### Likelihood-ratio test results #####
##
##

# for effects in final model
drop1(m_mam, test = "Chisq")

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
##### Final model results #####
##
##

# final model
day_length_final <- lmer(duration_act/60 ~ 
                           area:inc_start_aprildays +  
                           area:day_before_hatch +  
                           inc_start_aprildays +
                           day_before_hatch +
                           area + 
                           (1|site) +
                           (1|box) +
                           (1|year),
                         REML = T,                 
                         data=data)
summary(day_length_final)

# 95%CIs
CI_onset <- confint(day_length_final, 
                    level = 0.95, 
                    method = "boot", 
                    nsim = 500, 
                    boot.type = "norm")

# maternal repeatability of absolute end of activity
as.numeric(summary(day_length_final)$varcor[1]) / 
  ((as.numeric(summary(day_length_final)$varcor[1]) + 
      (as.numeric(summary(day_length_final)$varcor[2])) +
      (as.numeric(summary(day_length_final)$varcor[3])) +
      summary(day_length_final)$sigma^2))

# maternal repeatability with 95%CI of duration of active day
rep_duration <- rpt(duration_act/60 ~ 
                             area:inc_start_aprildays +  
                             area:day_before_hatch +  
                             inc_start_aprildays +
                             day_before_hatch +
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
##
##### Plot model predictions #####
##
##

# new dataframe to predict
df_pred <- expand.grid(day_before_hatch = seq(min(data$day_before_hatch), 
                                              max(data$day_before_hatch), 1),
                       area = c("City", "Forest"),
                       inc_start_aprildays = seq(min(data$inc_start_aprildays), 
                                                 max(data$inc_start_aprildays), 1))
df_pred$prediction <- predict(day_length_final, df_pred, re.form = NA)

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


## 
## plot days to hatching - duration
duration_hatch_plot <- ggplot(data = data, 
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

##
## plots not included in manuscript
#ggsave(filename = "./plots/duration_hatch.jpeg",
#       plot = duration_hatch_plot, 
#       device = "jpeg", 
#       units = "mm",
#       width = 125, 
#       height = 125)  


##
## plot for days from April 1 - Figure S2
duration_abs_date <- ggplot(data = data, 
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


ggsave(filename = "./plots/Figure S2.jpeg", 
       plot = duration_abs_date, 
       device = "jpeg", 
       units = "mm",
       width = 135, 
       height = 100)  






