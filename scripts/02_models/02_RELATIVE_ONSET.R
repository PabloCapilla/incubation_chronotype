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
#' Analysis of relative onset of activity
#' 
##
##

##
##
##### libraries #####
##
##
pacman::p_load(openxlsx, gt, gtsummary,
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
data <- readRDS("./data/data_incubation.RDS")
head(data)

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

# breeding attempts per nestbox
ba_per_box <- data %>% 
  group_by(year, box) %>% 
  filter(row_number() == 1) %>% 
  group_by(box) %>% 
  summarise(count_ba = n())
table(ba_per_box$count_ba)


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
            sd_date = sd(mean_box),
            se_date = sd(mean_box/sqrt(n())))

##### 

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
                               (1|year)+
                               (1|site) +
                               (1|box), 
                             REML = F,
                             na.action = "na.fail",
                             data= data) #full model
summary(model_relative_onset)

# model diagnostics
check_model(model_relative_onset, panel = T) # not too bad overall

#####


##
##
##### Are interactions significant? #####
##
##

## 1
drop1(model_relative_onset, test = "Chisq")

## interaction 'area:poly(inc_start_aprildays,2)[,2]'
anova(model_relative_onset, update(model_relative_onset, 
                                   . ~ . - area:poly(inc_start_aprildays,2)[,2]),
      test = "LRT")

## interaction 'area:poly(day_before_hatch, 2)[, 2]'
anova(model_relative_onset, update(model_relative_onset, 
                                   . ~ . - area:poly(day_before_hatch, 2)[, 2]),
      test = "LRT")

## removing both interactions 
full_model <- update(model_relative_onset, 
                     . ~ . 
                     - area:poly(inc_start_aprildays,2)[,2] 
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

# maternal repeatability with 95%CI of relative onset of activity
rep_relative_onset <- rpt(activity_onset_relative ~ 
                            area:poly(inc_start_aprildays,2)[,1] +  
                            area:poly(day_before_hatch,2)[,1] +  
                            
                            poly(inc_start_aprildays,2)[,2] +
                            poly(inc_start_aprildays,2)[,1] +
                            poly(day_before_hatch,2)[,2] +
                            poly(day_before_hatch,2)[,1] +
                            
                            meantemp +
                            clutch_size + 
                            area + 
                            (1|year)+
                            (1|site) +
                            (1|box), 
                          grname = c("box"), 
                          data = data, 
                          datatype = "Gaussian", 
                          nboot = 1000, 
                          npermut = 0)
summary(rep_relative_onset)

#####

##
##
##### Table of results 2 #####
##
##

## base table
table_relative_onset00 <- full_model %>%
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
                   `poly(inc_start_aprildays, 2)[, 1]:areaForest` = "Incubation start date1 x Habitat",
                   `poly(day_before_hatch, 2)[, 1]:areaForest` = "Days before hatching1 x Habitat"),
                 estimate_fun = ~ style_number(.x, digits = 2))

## add features
table_relative_onset <- table_relative_onset00 %>% 
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
gtsave(table_relative_onset, "./tables/TABLE 2.html")

#####

##
##
##### Plot model predictions #####
##
##
full_model_predictions <- lmer(activity_onset_relative ~ 
                                 area:inc_start_aprildays +  
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
                               data= data) #full model
summary(full_model_predictions)

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
  theme(legend.position = "top",
        legend.text = element_text("Arial", size = 17),
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
             size = 3, 
             shape = 21,
             color = "black",
             position = position_dodge(width = 0.5)) +
  labs(x = "Days before hatching", 
       y = "Onset of activity (minutes after sunrise)") +
  scale_x_continuous(breaks = -15:-1, labels = 15:1) +
  scale_fill_manual(name = "", labels = c("Urban", "Forest"), 
                    values = c("#af8dc3", "#7fbf7b")) +
  scale_color_manual(name = "", labels = c("Urban", "Forest"), 
                     values = c("#af8dc3", "#7fbf7b")) +
  geom_text(aes(-2, -15), 
            label = "Sunrise time", vjust = -1, 
            size = 4,
            color = "black") 

ggsave(filename = "./plots/Figure A5a.png", 
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
