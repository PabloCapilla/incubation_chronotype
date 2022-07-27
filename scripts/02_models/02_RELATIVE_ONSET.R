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
#' Analysis of relative onset of activity
#' 
##
##

##
##
##### libraries #####
##
##
pacman::p_load(openxlsx, 
               lubridate, dplyr, tidyr,
               lme4, performance, rptR,
               ggplot2, extrafont)
loadfonts()

##
##
##### data #####
##
##
data <- readRDS("./data/data_incubation.RDS")
head(data)



##
##
##### Sample sizes and data summary #####
##
##
head(data)
nrow(data) # total observations
length(unique(data$box)) # number of nest-boxes included

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

# number of nest-boxes per site (not clutches)
data %>% 
  group_by(site, box) %>% 
  filter(row_number() == 1) %>% 
  group_by(site) %>% 
  summarise(n_boxes = n())

# mean start of incubation per habitat
data %>% 
  group_by(year, area, box) %>% 
  summarise(mean_box = mean(inc_start_aprildays)) %>% 
  group_by(area) %>% 
  summarise(mean_date = mean(mean_box),
            n_size = n(),
            sd_date = sd(mean_box/sqrt(n())))

##
##### models for relative onset of activity #####
##
##
model_relative_onset <- lmer(activity_onset_relative ~ 
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
summary(model_relative_onset)

# model diagnostics
check_model(model_relative_onset, panel = F) # not too bad overall


# maternal repeatability of absolute onset of activity
as.numeric(summary(model_relative_onset)$varcor[1]) / 
  (as.numeric(summary(model_relative_onset)$varcor[1]) +
     as.numeric(summary(model_relative_onset)$varcor[2]) + 
     as.numeric(summary(model_relative_onset)$varcor[3]) +
     summary(model_relative_onset)$sigma^2)



##
##
##### Model selection LRT #####
##
##
drop1(model_relative_onset, test = "Chisq")
m1 <- update(model_relative_onset, . ~ . - clutch_size)

#1
drop1(m1, test = "Chisq")
m2 <- update(m1, . ~ . - poly(day_before_hatch, 2)[, 2]:area)

#2
drop1(m2, test = "Chisq")
m3 <- update(m2, . ~ . - meantemp)

#3
drop1(m3, test = "Chisq")
m4 <- update(m3, . ~ . - poly(inc_start_aprildays, 2)[, 2]:area)

#4
drop1(m4, test = "Chisq")
m5 <- update(m4, . ~ . - poly(inc_start_aprildays, 2)[, 2])

#mam
m_mam <- m5
summary(m_mam)

# 95%CIs
CI_onset <- confint(m_mam, 
                    level = 0.95, 
                    method = "boot", 
                    nsim = 500, 
                    boot.type = "norm")

##
##
##### Likelihood-ratio test results #####
##
##

# for effects in final model
drop1(m_mam, test = "Chisq")

# for effects not in final model
m_clutch_size <- update(m_mam, . ~ . + clutch_size)
anova(m_clutch_size, m_mam, test = "Chisq")

m_int <- update(m_mam, . ~ . + poly(day_before_hatch, 2)[, 2]:area)
anova(m_int, m_mam, test = "Chisq")

m_meantemp <- update(m_mam, . ~ . + meantemp)
anova(m_meantemp, m_mam, test = "Chisq")

m_aprildays2 <- update(m_mam, . ~ . + poly(inc_start_aprildays, 2)[, 2])
m_aprildays3 <- update(m_mam, . ~ . + poly(inc_start_aprildays, 2)[, 2] + poly(inc_start_aprildays, 2)[, 2]:area)
anova(m_aprildays2, m_mam, test = "Chisq")
anova(m_aprildays3, m_aprildays2, test = "Chisq")


##
##
##### Final model results #####
##
##
top_relative_onset <- lmer(activity_onset_relative ~ 
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
summary(top_relative_onset)

# maternal repeatability of relative onset of activity
as.numeric(summary(top_relative_onset)$varcor[1]) / 
  ((as.numeric(summary(top_relative_onset)$varcor[1]) + 
      (as.numeric(summary(top_relative_onset)$varcor[2])) +
      (as.numeric(summary(top_relative_onset)$varcor[3])) +
      summary(top_relative_onset)$sigma^2))

# maternal repeatability with 95%CI of relative onset of activity
rep_relative_onset <- rpt(activity_onset_relative ~ 
                            area:inc_start_aprildays + 
                            area:day_before_hatch +
                            
                            I(day_before_hatch^2) +
                            day_before_hatch +
                            inc_start_aprildays +
                            area + 
                            (1|site) +
                            (1|box) +
                            (1|year), 
                          grname = c("box", "Fixed"), 
                          data = data, 
                          datatype = "Gaussian", 
                          nboot = 1000, 
                          npermut = 0)
summary(rep_relative_onset)

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
df_pred$prediction <- predict(top_relative_onset, df_pred, re.form = NA)

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

# SE for mean predictions
mm <- model.matrix(~ area:inc_start_aprildays + 
                     area:day_before_hatch +
                     
                     I(day_before_hatch^2) +
                     day_before_hatch +
                     inc_start_aprildays +
                     area,
                   data = df_pred)
pvar1 <- diag(mm %*% tcrossprod(vcov(top_relative_onset),mm))
cmult <- 1 ## 1 SE
df_pred <- data.frame(
  df_pred
  , plo = df_pred$prediction-cmult*sqrt(pvar1)
  , phi = df_pred$prediction+cmult*sqrt(pvar1)
)




## plot hatching - relative to sunrise
relative_onset_hatching <- ggplot(data = data, 
                         aes(x = day_before_hatch, 
                             y = activity_onset_relative,
                             fill = area,
                             color = area)) +
  geom_hline(aes(yintercept=0), 
             colour="black", 
             linetype="dashed") +
  geom_point(position = position_jitterdodge(dodge.width = 0.5,
                                             jitter.width = 0.25),
             alpha = 0.40,
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
       y = "Onset of activity (minutes after sunrise)") +
  scale_x_continuous(breaks = -15:-1, labels = 15:1) +
  scale_fill_manual(name = "", labels = c("City", "Forest"), 
                    values = c("#af8dc3", "#7fbf7b")) +
  scale_color_manual(name = "", labels = c("City", "Forest"), 
                     values = c("#af8dc3", "#7fbf7b")) +
  geom_text(aes(-2, -15), 
            label = "Sunrise time", vjust = -1, 
            size = 4,
            color = "black") 

ggsave(filename = "./plots/Figure S4a.png", 
       plot = relative_onset_hatching, 
       device = "png", 
       units = "mm",
       width = 125, 
       height = 125)  


##
## plot for date from April 1
relative_onset_date <- ggplot(data = data, 
                              aes(x = inc_start_aprildays, 
                                  y = activity_onset_relative,
                                  fill = area,
                                  color = area)) +
  geom_hline(aes(yintercept=0), 
             colour="black", 
             linetype="dashed",
             size = 1,
             alpha = 0.75) +
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
       y = "Activity onset after sunrise (min)") +
  scale_fill_manual(name = "", labels = c("City", "Forest"), 
                    values = c("#af8dc3", "#7fbf7b")) +
  scale_color_manual(name = "", labels = c("City", "Forest"), 
                     values = c("#af8dc3", "#7fbf7b")) +
  geom_text(aes(48, -20), 
            label = "Sunrise time", vjust = -1, 
            size = 3.5,
            color = "black") 

ggsave(filename = "./plots/Figure 1cd.png", 
       plot = relative_onset_date, 
       device = "png", 
       units = "mm",
       width = 135, 
       height = 100)  
