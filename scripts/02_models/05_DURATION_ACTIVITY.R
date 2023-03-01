###
###
#' 
#' Script for:
#' Reproductive fitness is associated with female chronotype in a songbird
#' Womack, et al. 
#' Preprint: 10.1101/2022.07.01.498449v1
#' 
#' Latest update: 2023-02-27
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
pacman::p_load(openxlsx, gtsummary, gt,
               lubridate, dplyr, tidyr, rptR,
               lme4, performance,
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

## duration of activity in min
data$duration_act <- as.numeric(data$activity_end_abs - data$activity_onset_abs)

#####

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
                           (1|year) +
                           (1|site) +
                           (1|box),
                         na.action = "na.fail",
                         REML = F,                 
                         data=data) 
summary(day_length_model)

# model diagnostics
check_model(day_length_model, panel = T) # not too bad overall

#####

##
##
##### Are interactions significant? #####
##
##

## 1
drop1(day_length_model, test = "Chisq")

## interaction 'area:poly(inc_start_aprildays,2)[,2]'
anova(day_length_model, update(day_length_model, 
                                   . ~ . - area:poly(inc_start_aprildays,2)[,2]),
      test = "LRT")

## interaction 'area:poly(inc_start_aprildays, 2)[, 1]'
anova(day_length_model, 
      update(day_length_model, 
             . ~ . - area:poly(inc_start_aprildays, 2)[, 1]),
      test = "LRT")

## interaction 'area:poly(day_before_hatch, 2)[, 2]'
anova(day_length_model, 
      update(day_length_model, 
             . ~ . - area:poly(day_before_hatch, 2)[, 2]),
      test = "LRT")

## removing both interactions 
full_model <- update(day_length_model, 
                     . ~ . 
                     - area:poly(inc_start_aprildays, 2)[, 2] 
                     - area:poly(inc_start_aprildays, 2)[, 1]
                     - area:poly(day_before_hatch, 2)[, 2])

drop1(full_model, test = "Chisq")
summary(full_model)

#####

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
rep_day_length <- rpt(duration_act ~ 
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
                      grname = "box", 
                      data = data, 
                      datatype = "Gaussian", 
                      nboot = 1000, 
                      npermut = 0)
summary(rep_day_length)

#####

##
##
##### Table of results S4 #####
##
##

## base table
table_day_length00 <- full_model %>%
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
table_day_length <- table_day_length00 %>% 
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
gtsave(table_day_length, "./tables/TABLE S4.html")

#####

##
##
##### Plot model predictions #####
##
##

# full model for predictions
full_model_predictions <- lmer(duration_act ~ 
                                 area:day_before_hatch +  
                                 
                                 I(inc_start_aprildays^2) +
                                 inc_start_aprildays +
                                 I(day_before_hatch^2) +
                                 day_before_hatch +
                                 meantemp +
                                 clutch_size + 
                                 area + 
                                 (1|year/site) +
                                 (1|box), 
                               REML = F,
                               na.action = "na.fail",
                               data= data) 



# new dataframe to predict
df_pred <- expand.grid(day_before_hatch = seq(min(data$day_before_hatch), 
                                              max(data$day_before_hatch), 1),
                       area = c("City", "Forest"),
                       meantemp = mean(data$meantemp),
                       inc_start_aprildays = seq(min(data$inc_start_aprildays), 
                                                 max(data$inc_start_aprildays), 1))
df_pred$clutch_size <- NA
df_pred$clutch_size[df_pred$area == "Forest"] <- mean(data$clutch_size[data$area == "Forest"])
df_pred$clutch_size[df_pred$area == "City"] <- mean(data$clutch_size[data$area == "City"])
df_pred$prediction <- predict(full_model_predictions, df_pred, re.form = NA)

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

#####


##
##
##### Plot days to hatching - duration #####
##
##
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

#####

##
##
##### Plot for days from April 1 - Figure S2 #####
##
##
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


ggsave(filename = "./plots/Figure S3.jpeg", 
       plot = duration_abs_date, 
       device = "jpeg", 
       units = "mm",
       width = 135, 
       height = 100)  

#####




