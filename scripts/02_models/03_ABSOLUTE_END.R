###
###
#' 
#' Script for:
#' Reproductive fitness is associated with female chronotype in a songbird
#' Womack, et al. 
#' Preprint: 10.1101/2022.07.01.498449v1
#' 
#' Latest update: 2023-05-16
#' 
###
###

# Clear memory to make sure there are not files loaded that could cause problems
rm(list=ls())

##
##
##### Script aim: #####
#' Analysis of absolute end of activity
#' 
##
##

##
##
##### libraries #####
##
##
pacman::p_load(openxlsx, gtsummary, gt,
               lubridate, dplyr, tidyr,
               lme4, performance, rptR,
               ggplot2, extrafont)
loadfonts()
source("./scripts/FUNCTION_drop1_output.R")

#####

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

#####

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

#####

##
##
##### Models for absolute end of activity #####
##
##

# model fit
model_absolute_end <- lmer(activity_end_abs ~ 
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
                             (1|year) +
                             (1|site) +
                             (1|box),
                           REML = F,
                           na.action = "na.fail",
                           data= data) 
summary(model_absolute_end)

# model diagnostics
check_model(model_absolute_end, panel = T) # not too bad overall

#####

##
##
##### Are interactions significant? #####
##
##

## 1
drop1(model_absolute_end, test = "Chisq")

## interaction 'area:poly(inc_start_aprildays,2)[,2]'
anova(model_absolute_end, 
      update(model_absolute_end, 
             . ~ . - area:poly(inc_start_aprildays, 2)[, 2]),
      test = "LRT")

## interaction 'area:poly(inc_start_aprildays, 2)[, 1]'
anova(model_absolute_end, 
      update(model_absolute_end, 
             . ~ . - area:poly(inc_start_aprildays, 2)[, 1]),
      test = "LRT")

## interaction 'area:poly(day_before_hatch, 2)[, 2]'
anova(model_absolute_end, 
      update(model_absolute_end, 
             . ~ . - area:poly(day_before_hatch, 2)[, 2]),
      test = "LRT")



## removing both interactions 
full_model <- update(model_absolute_end, 
                     . ~ . 
                     - area:poly(inc_start_aprildays, 2)[, 2] 
                     - area:poly(inc_start_aprildays, 2)[, 1]
                     - area:poly(day_before_hatch, 2)[, 2])

drop1(full_model, test = "Chisq")
summary(full_model)


##
##
##### Final model results #####
##
##
summary(full_model)
lmerTest::rand(full_model)

# maternal repeatability of absolute onset of activity
as.numeric(summary(full_model)$varcor[1]) / 
  ((as.numeric(summary(full_model)$varcor[1]) + 
      (as.numeric(summary(full_model)$varcor[2])) +
      (as.numeric(summary(full_model)$varcor[3])) +
      summary(full_model)$sigma^2))

# maternal repeatability with 95%CI of absolute onset of activity
rep_absolute_end <- rpt(activity_end_abs ~ 
                          area:poly(day_before_hatch,2)[,1] +  
                          
                          poly(inc_start_aprildays,2)[,2] +
                          poly(inc_start_aprildays,2)[,1] +
                          poly(day_before_hatch,2)[,2] +
                          poly(day_before_hatch,2)[,1] +
                          
                          meantemp +
                          clutch_size + 
                          area + 
                          (1|year) +
                          (1|site) +
                          (1|box), 
                        grname = "box", 
                        data = data, 
                        datatype = "Gaussian", 
                        nboot = 1000, 
                        npermut = 0)
summary(rep_absolute_end)

#####

##
##
##### Table of results A2 #####
##
##

## base table
table_absolute_end00 <- full_model %>%
  tbl_regression(intercept = T,
                 label = list(
                   `(Intercept)` = "Intercept",
                   `poly(inc_start_aprildays, 2)[, 1]`  = "Incubation start date1", 
                   `poly(inc_start_aprildays, 2)[, 2]`  = "Incubation start date2",
                   `poly(day_before_hatch, 2)[, 1]` = "Days before hatching1",
                   `poly(day_before_hatch, 2)[, 2]` = "Days before hatching2",
                   `meantemp` = "Mean daily temperatures",
                   `clutch_size` = "Clutch size",
                   `area` = "Habitat",
                   `poly(day_before_hatch, 2)[, 1]:area` = "Days before hatching1 x Habitat"),
                 estimate_fun = ~ style_number(.x, digits = 2))

## add features
table_absolute_end <- table_absolute_end00 %>% 
  add_global_p(anova_fun = drop1_output) %>% 
  bold_p(t = 0.05) %>% 
  bold_labels() %>%
  italicize_levels() %>% 
  modify_table_body(fun = function(.){
    output <- left_join(x = .,
                        y = drop1_output(x=full_model) %>% 
                          dplyr::select(variable = term, Chisq=statistic, df),
                        by = "variable")
    output$df <- ifelse(output$row_type == "label",  output$df, NA)
    output$Chisq <- ifelse(output$row_type == "label",  output$Chisq, NA)
    return(output)
  }) %>% 
  modify_fmt_fun(c(Chisq) ~ function(x) style_number(x, digits = 2)) %>%
  modify_fmt_fun(c(std.error) ~ function(x) style_number(x, digits = 2)) %>%
  modify_fmt_fun(c(p.value) ~ function(x) style_number(x, digits = 3)) %>%
  modify_table_body(~.x %>% dplyr::relocate(p.value, .after = df)) %>% 
  modify_header(label ~ "**Fixed effect**") %>% 
  modify_header(std.error ~ "**SE**") %>%
  modify_header(estimate ~ "**Estimate**") %>%
  modify_header(df ~ "**df**") %>% 
  modify_header(Chisq ~ html("<b>&chi;<sup>2</sup></b>")) %>% 
  as_gt() %>% 
  opt_footnote_marks(marks = "LETTERS")

##
## save table
gtsave(table_absolute_end, "./tables/TABLE A2.html")

#####

##
##
##### Plot model predictions #####
##
##

# full model for predictions
full_model_predictions <- lmer(activity_end_abs ~ 
                                 area:day_before_hatch +  
                                 
                                 I(inc_start_aprildays^2) +
                                 inc_start_aprildays +
                                 I(day_before_hatch^2) +
                                 day_before_hatch +
                                 meantemp +
                                 clutch_size + 
                                 area + 
                                 (1|year) +
                                 (1|site) +
                                 (1|box), 
                               REML = F,
                               na.action = "na.fail",
                               data= data) 


# new dataframe to predict
df_pred <- expand.grid(day_before_hatch = seq(min(data$day_before_hatch), 
                                              max(data$day_before_hatch), 1),
                       inc_start_aprildays = seq(min(data$inc_start_aprildays), 
                                             max(data$inc_start_aprildays), 1),
                       meantemp = mean(data$meantemp),
                       
                       area = c("City", "Forest"))
df_pred$clutch_size <- NA
df_pred$clutch_size[df_pred$area == "Forest"] <- mean(data$clutch_size[data$area == "Forest"])
df_pred$clutch_size[df_pred$area == "City"] <- mean(data$clutch_size[data$area == "City"])
df_pred$prediction <- predict(full_model_predictions, df_pred, re.form = NA)

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
                     
                     I(inc_start_aprildays^2) +
                     inc_start_aprildays +
                     I(day_before_hatch^2) +
                     day_before_hatch +
                     meantemp +
                     clutch_size + 
                     area,
                   data = df_pred)

pvar1 <- diag(mm %*% tcrossprod(vcov(full_model_predictions),mm))
cmult <- 1 ## 1 SE
df_pred <- data.frame(
  df_pred
  , plo = df_pred$prediction-cmult*sqrt(pvar1)
  , phi = df_pred$prediction+cmult*sqrt(pvar1)
)

df_pred$area <- factor(x = ifelse(df_pred$area == "City", "Urban", "Forest"),
                       levels = c("Urban", "Forest"),
                       ordered = T)
data$area <- factor(x = ifelse(data$area == "City", "Urban", "Forest"),
                            levels = c("Urban", "Forest"),
                            ordered = T)

## labels for plot
library(chron) 
min_time <- min(times(strftime(dmy_hm(data$last_onbout),"%H:%M:%S", tz = "UTC")))
max_time <- max(times(strftime(dmy_hm(data$last_onbout),"%H:%M:%S", tz = "UTC")))

hm <- merge(0:23, seq(0, 30, by = 30))

interval_15min <- data.frame('INTERVAL' = chron(time = paste(hm$x, ':', hm$y, ':', 0)))
interval_15min <- data.frame('INTERVAL' = interval_15min[order(interval_15min$INTERVAL), ])
labels_time <- substr(interval_15min[!interval_15min$INTERVAL <= min_time & !interval_15min$INTERVAL >= max_time,],
                      start = 0, 
                      stop = 5)
labels_time <- sub(x = labels_time, pattern = ":", replacement = "")


## 
## plot days to hatching - absolute onset
absolute_end_hatching <- ggplot(data = data, 
                                aes(x = day_before_hatch, 
                                    y = activity_end_abs,
                                    fill = area,
                                    color = area)) +
  geom_point(position = position_jitterdodge(dodge.width = 0.5,
                                             jitter.width = 0.25),
             alpha = 0.4,
             size = 1.25,
             shape = 21, 
             color = "black") +
  theme_bw() +
  theme(legend.position = "none",
        legend.text = element_text("Arial", size = 15),
        panel.grid = element_blank(),
        axis.title = element_text("Arial", size = 15),
        axis.text = element_text("Arial", size = 13)) +
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
       y = "Clock time end of activity") +
  scale_x_continuous(breaks = -15:-1, labels = 15:1) +
  scale_y_continuous(breaks = seq(1050,1320, 30),
                     labels = labels_time) +
  scale_fill_manual(name = "", labels = c("Urban", "Forest"), 
                    values = c("#af8dc3", "#7fbf7b")) +
  scale_color_manual(name = "", labels = c("Urban", "Forest"), 
                     values = c("#af8dc3", "#7fbf7b"))

ggsave(filename = "./plots/Figure A4b.png", 
       plot = absolute_end_hatching, 
       device = "png", 
       units = "mm",
       width = 125, 
       height = 125)  

##
## plot for date from April 1
absolute_end_date <- ggplot(data = data, 
                            aes(x = inc_start_aprildays, 
                                y = activity_end_abs,
                                fill = area,
                                color = area)) +
  geom_line(data = data,
            aes(y = (1+sunset_inc_start_dec)*60, color = "black"),
            size = 1,
            alpha = 0.75,
            linetype = 2) +
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
                        mean_pred = mean(prediction)) %>% 
              mutate(area = factor(area, ordered = T, levels = c("Urban", "Forest"))), 
            aes(y = mean_pred, color = area), 
            size = 1.5) +
  labs(x = "Incubation start date (days after April 1)", 
       y = "Clock time end of activity") +
  scale_y_continuous(breaks = seq(1050,1320, 30),
                     labels = labels_time) +
  scale_fill_manual(name = "", labels = c("Urban", "Forest"), 
                    values = c("#af8dc3", "#7fbf7b")) +
  scale_color_manual(name = "", labels = c(NA, "Urban", "Forest"), 
                     values = c("black", "#7fbf7b", "#af8dc3"))


ggsave(filename = "./plots/Figure A2ab.png", 
       plot = absolute_end_date, 
       device = "png", 
       units = "mm",
       width = 135, 
       height = 100)  













