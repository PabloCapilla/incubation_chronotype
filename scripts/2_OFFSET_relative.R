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


##
## sample size
nrow(data) # total data points

# data per site and year
data %>% 
  group_by(year, site) %>% 
  summarise(n_obs = n())

# number of boxes per site and year
data %>% 
  group_by(year, site, box) %>% 
  summarise(n_obs = n()) %>% 
  group_by(year, site) %>%
  summarise(n_box = n())



data_trimmed %>% 
  group_by(year, area) %>% 
  summarise(n_obs = n())

data_trimmed %>% 
  group_by(box, year, area) %>% 
  filter(row_number() == 1) %>% 
  group_by(year, area) %>% 
  summarise(n_obs = n())



##
##### Initial visualisation of data #####
##
ggplot(data = data_trimmed, aes(x = diff_sunset)) +
  geom_histogram() +
  theme_bw() +
  labs(x = "Activity onset - sunrise time", y = "Count")

ggplot(data = data_trimmed, aes(y = diff_sunset, x = day_before_hatch, color = area)) +
  geom_point() +
  stat_smooth(method = "lm") +
  theme_bw() +
  labs(x = "Days before incubation", 
       y = "End of activity (min after sunset)") +
  scale_fill_manual(name = "", labels = c("City", "Forest"), 
                    values = c("#3399CC", "#339900")) +
  scale_color_manual(name = "", labels = c("City", "Forest"), 
                     values = c("#3399CC", "#339900")) +
  geom_text(aes(-2, -15), 
            label = "Sunset time", vjust = -1, 
            size = 4,
            color = "black") 


ggplot(data = data_trimmed, 
       aes(y = diff_sunset, 
           x = mean_days_april, 
           color = area)) +
  geom_point() +
  stat_smooth(method = "lm") +
  theme_bw() 






ggplot(data = data_trimmed, aes(x = diff_sunset, y = factor(area))) +
  geom_density_ridges(scale = 0.8, 
                      alpha = 0.8, 
                      aes(fill = area)) +
  theme_bw() +
  scale_fill_manual(values = c("#8DA0CB", "#5AAE61")) +
  labs(x = "Last on-bout - sunset time", y = "Count")



ggplot(data = data_trimmed, aes(x = day_before_hatch, y = date_aprildays, group = box)) +
  geom_point() +
  geom_line() +
  facet_grid(year~area) +
  theme_bw()

ggplot(data = data_trimmed, aes(x = day_before_hatch, y = mean_days_april, group = box)) +
  geom_point() +
  geom_line() +
  facet_grid(year~area) +
  theme_bw()


##
##
##### Model for end of activity - relative to sunset #####
##
##

# model fit
model_end <- lmer(diff_sunset ~ 
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
                  data= data_trimmed) #full model
summary(model_end)

# model diagnostics
plot(model_end)
hist(residuals(model_end))

#repeatability
as.numeric(summary(model_end)$varcor[1]) / 
  (as.numeric(summary(model_end)$varcor[1]) + 
     summary(model_end)$sigma^2)

##
##
##### Model selection LRT #####
##
##
drop1(model_end, test = "Chisq")
m1 <- update(model_end, . ~ . - meantemp)

#1
drop1(m1, test = "Chisq")
m2 <- update(m1, . ~ . - area:poly(inc_start_aprildays, 2)[, 1])

#2
drop1(m2, test = "Chisq")
m3 <- update(m2, . ~ . - poly(inc_start_aprildays, 2)[, 2]:area)

#3
drop1(m3, test = "Chisq")
m4 <- update(m3, . ~ . - poly(inc_start_aprildays, 2)[, 2])

#4
drop1(m4, test = "Chisq")
m5 <- update(m4, . ~ . - poly(day_before_hatch, 2)[, 2]:area)

#5
drop1(m5, test = "Chisq")
m6 <- update(m5, . ~ .  - poly(day_before_hatch, 2)[, 2])

#6
drop1(m6, test = "Chisq")
m7 <- update(m6, . ~ . - clutch_size)

#8
drop1(m7, test = "Chisq")

#mam
m_mam <- m7
drop1(m_mam, test = "Chisq")
summary(m_mam)

CI_onset <- confint(m_mam, 
                    level = 0.95, 
                    method = "boot", 
                    nsim = 500, 
                    boot.type = "norm")

# p for non significant terms
m_int <- update(m_mam, . ~ . + area:poly(inc_start_aprildays, 2)[, 1])
anova(m_int, m_mam, test = "Chisq")

m_aprildays2 <- update(m_mam, . ~ . + poly(inc_start_aprildays, 2)[, 2])
anova(m_aprildays2, m_mam, test = "Chisq")

m_aprildays2_int1 <- update(m_mam, . ~ . + poly(inc_start_aprildays, 2)[, 2] + poly(inc_start_aprildays, 2)[, 2]:area)
m_aprildays2_int2 <- update(m_mam, . ~ . + poly(inc_start_aprildays, 2)[, 2])
anova(m_aprildays2_int1, m_aprildays2_int2, test = "Chisq")


m_hatch2 <- update(m_mam, . ~ . + poly(day_before_hatch, 2)[, 2])
anova(m_hatch2, m_mam, test = "Chisq")

m_hatch2_int1 <- update(m_mam, . ~ . + poly(day_before_hatch, 2)[, 2] + poly(day_before_hatch, 2)[, 2]:area)
m_hatch2_int2 <- update(m_mam, . ~ . + poly(day_before_hatch, 2)[, 2])
anova(m_hatch2_int1, m_hatch2_int2, test = "Chisq")


m_clutch_size <- update(m_mam, . ~ . + clutch_size)
anova(m_clutch_size, m_mam, test = "Chisq")

m_meantemp <- update(m_mam, . ~ . + meantemp)
anova(m_meantemp, m_mam, test = "Chisq")







##
## predictions plot

# top model for plot
model_end_plot <- lmer(diff_sunset ~ 
                         area:day_before_hatch +  
                         
                         day_before_hatch +
                         inc_start_aprildays +
                         area + 
                         (1|site) +
                         (1|box) +
                         (1|year), 
                       REML = T,
                       na.action = "na.fail",
                       data= data_trimmed) 
summary(model_end_plot)


##
## repeatability
rep_relative_offset <- rpt(diff_sunset ~ 
                             area:day_before_hatch +  
                             
                             day_before_hatch +
                             inc_start_aprildays +
                             area + 
                             (1|site) +
                             (1|box) +
                             (1|year), 
                           grname = "box", 
                           data = data_trimmed, 
                           datatype = "Gaussian", 
                           nboot = 1000, 
                           npermut = 0)


CI_onset <- confint(model_end_plot, 
                    level = 0.95, 
                    method = "boot", 
                    nsim = 500, 
                    boot.type = "norm")

df_pred <- expand.grid(day_before_hatch = seq(min(data$day_before_hatch), 
                                              max(data$day_before_hatch), 1),
                       inc_start_aprildays = seq(min(data$inc_start_aprildays), 
                                             max(data$inc_start_aprildays), 1),
                       area = c("City", "Forest"))

df_pred$prediction <- predict(model_end_plot, df_pred, re.form = NA)

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
mm <- model.matrix(~ area:day_before_hatch +  
                     
                     day_before_hatch +
                     inc_start_aprildays +
                     area,
                   data = df_pred)

## or newdat$distance <- mm %*% fixef(fm1)
pvar1 <- diag(mm %*% tcrossprod(vcov(model_end_plot),mm))
cmult <- 1 ## 1 SE
df_pred <- data.frame(
  df_pred
  , plo = df_pred$prediction-cmult*sqrt(pvar1)
  , phi = df_pred$prediction+cmult*sqrt(pvar1)
)

## plot end of activity averaging across April days
end_plot <- ggplot(data = data_trimmed, 
                   aes(x = day_before_hatch, 
                       y = diff_sunset,
                       fill = area,
                       color = area)) +
  geom_hline(aes(yintercept=0), 
             colour="black", 
             linetype="dashed") +
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
       y = "End of activity (minutes before sunset)") +
  scale_y_continuous(breaks = seq(-200, 0, 50), labels = seq(-200, 0, 50)*-1) +
  scale_x_continuous(breaks = -15:-1, labels = 15:1) +
  scale_fill_manual(name = "", labels = c("City", "Forest"), 
                    values = c("#af8dc3", "#7fbf7b")) +
  scale_color_manual(name = "", labels = c("City", "Forest"), 
                     values = c("#af8dc3", "#7fbf7b")) +
  geom_text(aes(-2, -30), 
            label = "Sunset time", vjust = -1, 
            size = 4,
            color = "black") 


ggsave(filename = "./plots/1b_End of activity/end_pred.jpeg", 
       plot = end_plot, 
       device = "jpeg", 
       units = "mm",
       width = 125, 
       height = 125)  



## plot end of activity averaging across days to hatching
end_plot_jan1 <- ggplot(data = data_trimmed, 
                        aes(x = inc_start_aprildays, 
                            y = diff_sunset,
                            fill = area,
                            color = area)) +
  theme_bw() +
  facet_wrap(~area) +
  theme(legend.position = "none",
        legend.text = element_text("Arial", size = 10),
        strip.background = element_blank(),
        strip.text = element_text("Arial", size = 12),
        panel.grid = element_blank(),
        axis.title = element_text("Arial", size = 10),
        axis.text = element_text("Arial", size = 10)) +
  geom_hline(aes(yintercept=0), 
             colour="black", 
             linetype="dashed",
             size = 1,
             alpha = 0.75) +
  geom_point(position = position_jitterdodge(dodge.width = 0.5,
                                             jitter.width = 0.25),
             alpha = 0.25,
             size = 1.25,
             shape = 21, 
             color = "black") +
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
       y = "End of activity (minutes before sunset)") +
  scale_y_continuous(breaks = seq(-200, 0, 50), labels = seq(-200, 0, 50)*-1) +
  scale_fill_manual(name = "", labels = c("City", "Forest"), 
                    values = c("#af8dc3", "#7fbf7b")) +
  scale_color_manual(name = "", labels = c("City", "Forest"), 
                     values = c("#af8dc3", "#7fbf7b")) +
  geom_text(aes(48, -33), 
            label = "Sunset time", vjust = -1, 
            size = 3.5,
            color = "black") 


ggsave(filename = "./plots/1b_End of activity/end_pred_april1.jpeg", 
       plot = end_plot_jan1, 
       device = "jpeg", 
       units = "mm",
       width = 135, 
       height = 100)  




















