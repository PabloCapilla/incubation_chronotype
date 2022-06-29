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



##
##
##### Preliminary visualisation of data #####
##
##
ggplot(data = data, aes(x = night_var)) +
  geom_histogram() +
  theme_bw() +
  labs(x = "Nighttime variation in temperature", y = "Count")


ggplot(data = data, aes(y = log(night_var), x = inc_start_aprildays, color = area)) +
  geom_point() +
  facet_wrap(year~area) +
  stat_smooth(method = "lm", formula = y ~ poly(x,2)) +
  theme_bw() +
  scale_fill_manual(name = "", labels = c("City", "Forest"), 
                    values = c("#3399CC", "#339900")) +
  scale_color_manual(name = "", labels = c("City", "Forest"), 
                     values = c("#3399CC", "#339900")) 

ggplot(data = data, aes(y = log(night_var), x = day_before_hatch, color = area)) +
  geom_point() +
  facet_wrap(year~area) +
  stat_smooth(method = "lm", formula = y ~ poly(x,2)) +
  theme_bw() +
  scale_fill_manual(name = "", labels = c("City", "Forest"), 
                    values = c("#3399CC", "#339900")) +
  scale_color_manual(name = "", labels = c("City", "Forest"), 
                     values = c("#3399CC", "#339900")) 

##
##
###### Models for night time variance (restlessness) #####
##
##
data_nightvar <- data %>% filter(!is.na(night_var))

# model fit
model_nightvar <- lmer(log(night_var) ~ 
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
                         as.factor(year) + 
                         (1|box), 
                  REML = T,
                  na.action = "na.fail",
                  data= data_nightvar) #full model

summary(model_nightvar)

# model diagnostics
plot(model_nightvar)
hist(residuals(model_nightvar))


##
##
##### Model selection LRT #####
##
##
drop1(model_nightvar, test = "Chisq")
m1 <- update(model_nightvar, . ~ . - clutch_size)

#1
drop1(m1, test = "Chisq")
m2 <- update(m1, . ~ . - poly(day_before_hatch, 2)[, 2]:area)

#2
drop1(m2, test = "Chisq")
m3 <- update(m2, . ~ . - poly(day_before_hatch, 2)[, 2])

#3
drop1(m3, test = "Chisq")
m4 <- update(m3, . ~ . - poly(day_before_hatch, 2)[, 1]:area)

#4
drop1(m4, test = "Chisq")


#mam
m_mam <- m4
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

m_hatch2 <- update(m_mam, . ~ . + poly(day_before_hatch, 2)[, 2])
anova(m_hatch2, m_mam, test = "Chisq")

m_int <- update(m_mam, . ~ . + poly(day_before_hatch, 2)[, 1]:area)
anova(m_int, m_mam, test = "Chisq")


##
## predictions plot

# top model for plot
model_var_plot <- lmer(log(night_var) ~ +
                         area:I(inc_start_aprildays^2) + 
                         area:inc_start_aprildays +  

                         I(inc_start_aprildays^2) +
                         inc_start_aprildays +
                         day_before_hatch +
                         area + 
                         meantemp +
                         as.factor(year) + 
                         (1|box),
                       REML = T,
                       na.action = "na.fail",
                       data= data_nightvar) 
summary(model_var_plot)
drop1(model_var_plot, test = "Chisq")

df_pred <- expand.grid(day_before_hatch = seq(min(data_nightvar$day_before_hatch), 
                                              max(data_nightvar$day_before_hatch), 1),
                       inc_start_aprildays = seq(min(data_nightvar$inc_start_aprildays), 
                                                max(data_nightvar$inc_start_aprildays), 1),
                       area = c("City", "Forest"),
                       year = c("2016", "2017", "2018"))
df_pred$meantemp <- ifelse(df_pred$area == "City", 
                           mean(data_nightvar$meantemp[data_nightvar$area == "City"]),
                           mean(data_nightvar$meantemp[data_nightvar$area == "Forest"]))

df_pred$prediction <- predict(model_var_plot, df_pred, re.form = NA)

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
mm <- model.matrix(~ area:I(inc_start_aprildays^2) + 
                     area:inc_start_aprildays +  
                     
                     I(inc_start_aprildays^2) +
                     inc_start_aprildays +
                     day_before_hatch +
                     area + 
                     meantemp +
                     as.factor(year),
                   data = df_pred)

## or newdat$distance <- mm %*% fixef(fm1)
pvar1 <- diag(mm %*% tcrossprod(vcov(model_var_plot),mm))
cmult <- 1 ## 1 SE
df_pred <- data.frame(
  df_pred
  , plo = df_pred$prediction-cmult*sqrt(pvar1)
  , phi = df_pred$prediction+cmult*sqrt(pvar1)
)

## 
##
##### Plots with model predictions #####
##
##

# days to hatching
varnight_plot <- ggplot(data = data_nightvar, 
                        aes(x = day_before_hatch, 
                            y = log(night_var),
                            fill = area,
                            color = area)) +
  geom_point(position = position_jitterdodge(dodge.width = 0.5,
                                             jitter.width = 0.25),
             alpha = 0.35,
             size = 1) +
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
       y = "Log Variance in Nighttime incubation temperature") +
  scale_x_continuous(breaks = -15:-1, labels = 15:1) +
  scale_fill_manual(name = "", labels = c("City", "Forest"), 
                    values = c("#3399CC", "#339900")) +
  scale_color_manual(name = "", labels = c("City", "Forest"), 
                     values = c("#3399CC", "#339900"))

ggsave(filename = "./plots/nightvar_hatch.jpeg", 
       plot = varnight_plot, 
       device = "jpeg", 
       units = "mm",
       width = 130, 
       height = 120)  


# days from April 1
varnight_april_date <- ggplot(data = data_nightvar, 
                              aes(x = inc_start_aprildays, 
                                  y = log(night_var),
                                  fill = area,
                                  color = area)) +
  geom_point(alpha = 0.35,
             size = 1) +
  theme_bw() +
  facet_wrap(~area) +
  #stat_smooth(method = "lm", formula = y ~ poly(x,2)) +
  theme(legend.position = "none",
        legend.text = element_text("Arial", size = 10),
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
  labs(x = "Days after April 1", 
       y = "Log Variance in Nighttime incubation temperature") +
  scale_fill_manual(name = "", labels = c("City", "Forest"), 
                    values = c("#3399CC", "#339900")) +
  scale_color_manual(name = "", labels = c("City", "Forest"), 
                     values = c("#3399CC", "#339900")) 


ggsave(filename = "./plots/nightvar_april.jpeg", 
       plot = varnight_april_date, 
       device = "jpeg", 
       units = "mm",
       width = 130, 
       height = 120)  





