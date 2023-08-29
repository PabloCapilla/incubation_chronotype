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
#' Analysis of reproductive success, number of nestlings surviving to fledging
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
source("./scripts/FUNCTION_drop1_output.R")

#####

##
##
##### data #####
##
##
data00 <- readRDS("./data/data_incubation.RDS")
head(data00)

# hatching data in days after Jan 0
data00$hatching_date_jul <- yday(dmy(data00$hatching_date))

## calculate absolute female chronotype
chrono_box <- data00 %>%
  filter(!is.na(chronotype)) %>%
  group_by(year, box) %>% 
  summarise(chronotype_abs = mean(activity_onset_abs))

# add female chronotype to dataset
data <- data00 %>%
  filter(!is.na(chronotype)) %>%
  filter(!is.na(fledglings)) %>% 
  group_by(year, box) %>% 
  filter(row_number() == 1) %>% 
  mutate(failed_n = hatchlings - fledglings) %>% 
  left_join(., y = chrono_box, by = c("year", "box"))

#####

##
##
##### Models for number of fledglings #####
##
##
full_fledge_model <- glmer(fledglings ~
                             poly(hatching_date_jul,2)[,1] : area +
                             poly(hatching_date_jul,2)[,2] : area +
                             scale(chronotype):area+
                             
                             poly(hatching_date_jul,2)[,1] +
                             poly(hatching_date_jul,2)[,2] +
                             scale(chronotype)+
                             area+
                             clutch_size+
                             (1|year) +
                             (1|site) +
                             (1|box), 
                           family = "poisson",
                           control = glmerControl(optimizer = "bobyqa"),
                           na.action = "na.fail",
                           data=data)
summary(full_fledge_model)

# model diagnostics
residuals <- DHARMa::simulateResiduals(full_fledge_model, n = 1000, refit = F)
DHARMa::testDispersion(residuals)

#####

##
##
##### Are interactions significant? #####
##
##

## 1
drop1(full_fledge_model, test = "Chisq")

## interaction 'area:poly(inc_start_aprildays,2)[,2]'
anova(full_fledge_model, 
      update(full_fledge_model, 
             . ~ . - poly(hatching_date_jul,2)[,2] : area),
      test = "LRT")

## interaction 'area:poly(day_before_hatch, 2)[, 2]'
anova(full_fledge_model, 
      update(full_fledge_model, 
             . ~ . - poly(hatching_date_jul,2)[,1] : area),
      test = "LRT")

## interaction 'scale(chronotype):area'
anova(full_fledge_model, 
      update(full_fledge_model, 
             . ~ . - scale(chronotype):area),
      test = "LRT")


## removing both interactions 
full_model <- update(full_fledge_model, 
                     . ~ . 
                     - scale(chronotype):area 
                     - poly(hatching_date_jul,2)[,2] : area
                     - poly(hatching_date_jul,2)[,1] : area)

## comparison model without any interaction against initial
anova(full_fledge_model, 
      full_model,
      test = "LRT")

#####

##
##
##### Final model results #####
##
##
summary(full_model)
lmerTest::rand(full_model)

drop1(full_model, test = "Chisq")
summary(full_model)

##
## effect of absolute chronotype
drop1(update(full_model, .~.
             +scale(chronotype_abs)
             -scale(chronotype)), 
             test = "Chisq")

drop1(update(full_model, .~.
             +scale(chronotype_abs):area
             +scale(chronotype_abs)
             -scale(chronotype)), 
      test = "Chisq")

##
## final model with year as a fixed effect
full_fledge_model_year <- glmer(fledglings ~

                             poly(hatching_date_jul,2)[,1] +
                             poly(hatching_date_jul,2)[,2] +
                             scale(chronotype)+
                             area+
                             clutch_size+
                             as.factor(year) +
                             (1|site) +
                             (1|box), 
                           family = "poisson",
                           control = glmerControl(optimizer = "bobyqa"),
                           na.action = "na.fail",
                           data=data)
summary(full_fledge_model_year)
drop1(full_fledge_model_year, test = "Chisq")


#####

##
##
##### Table of results 3 #####
##
##

## base table
table_fledglings00 <- full_model %>%
  tbl_regression(intercept = T,
                 label = list(
                   `(Intercept)` = "Intercept",
                   `poly(hatching_date_jul, 2)[, 1]` = "Hatching date1",
                   `poly(hatching_date_jul, 2)[, 2]` = "Hatching date2",
                   `clutch_size` = "Clutch size",
                   `scale(chronotype)` = "Female chronotype",
                   `area` = "Habitat"),
                 estimate_fun = ~ style_number(.x, digits = 2))

## add features
table_fledglings <- table_fledglings00 %>% 
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
gtsave(table_fledglings, "./tables/TABLE 3.html")

#####

##
##
##### Plot model predictions #####
##
##

##
## model specification for predictions
full_model_predictions <- glmer(fledglings ~
                                  
                                  I(hatching_date_jul^2) +
                                  hatching_date_jul +
                                  chronotype+
                                  area+
                                  clutch_size+
                                  (1|year) +
                                  (1|site) +
                                  (1|box), 
                                family = "poisson",
                                control = glmerControl(optimizer = "bobyqa"),
                                na.action = "na.fail",
                                data=data)

# new data frame to predict
df_pred <- expand.grid(chronotype = seq(min(data$chronotype), 
                                        max(data$chronotype), 0.5),
                       area = c("City", "Forest"))
df_pred$clutch_size <- ifelse(df_pred$area == "City", 
                              mean(data$clutch_size[data$area == "City"]),
                              mean(data$clutch_size[data$area == "Forest"]))
df_pred$hatching_date_jul <- ifelse(df_pred$area == "City", 
                              mean(data$hatching_date_jul[data$area == "City"]),
                              mean(data$hatching_date_jul[data$area == "Forest"]))
df_pred$prediction_link <- predict(full_model_predictions,
                                   df_pred, 
                                   type = "link",
                                   re.form = NA)

# plot only data in range
data %>% 
  group_by(area) %>% 
  summarise(min_chr = min(chronotype),
            max_chr = max(chronotype))

remove_city <- which((df_pred$chronotype < -21 | 
                        df_pred$chronotype > 50) & 
                       df_pred$area == "City")
remove_forest <- which((df_pred$chronotype < -9.5 | 
                          df_pred$chronotype > 91.9) & 
                         df_pred$area == "Forest")
df_pred <- df_pred[-c(remove_city, remove_forest),]


# need to adapt code to new R version
# SE for mean predicitons
mm <- model.matrix(~ 
                     I(hatching_date_jul^2) +
                     hatching_date_jul +
                     chronotype+
                     area+
                     clutch_size,
                   data = df_pred)

## or newdat$distance <- mm %*% fixef(fm1)
pvar1 <- diag(mm %*% tcrossprod(vcov(full_model_predictions),mm))
cmult <- 1 ## 1 SE
df_pred <- data.frame(
  df_pred
  , plo_link = df_pred$prediction_link-cmult*sqrt(pvar1)
  , phi_link = df_pred$prediction_link+cmult*sqrt(pvar1)
)
df_pred$prediction <- exp(df_pred$prediction_link)
df_pred$plo <- exp(df_pred$plo_link)
df_pred$phi <- exp(df_pred$phi_link)


##
## Figure 2a
fledglings_plot <- ggplot(data = data, aes(x = chronotype, 
                                             y = fledglings,
                                             fill = area,
                                             color = area)) +
  geom_point(alpha = 0.5,
             size = 2,
             color = "black",
             aes(shape = area), 
             position = position_jitter(height = 0.05)) +
  #facet_grid(.~area) +
  theme_bw() +
  theme(legend.position = "top",
        #strip.background = element_blank(),
        #strip.text = element_text("Arial", size = 12),
        legend.text = element_text("Arial", size = 10),
        panel.grid = element_blank(),
        axis.title = element_text("Arial", size = 10),
        axis.text = element_text("Arial", size = 10)) +
  geom_ribbon(data = df_pred, aes(ymin = plo, 
                                  ymax = phi, 
                                  y = prediction),
              alpha = 0.25,
              color = NA,
              position = position_dodge(width = 0.5)) +
  geom_line(data = df_pred, 
            aes(y = prediction), 
            size = 1.5) +
  labs(x = expression(atop("Female chronotype", 
                   "(relative onset of activity [minutes after sunrise])")), 
    x = "Female chronotype",
    y = "Number of fledglings") +
  scale_shape_manual(name = "",
                     values = c(21,24),
                     labels = c("Urban", "Forest")) +  
  scale_fill_manual(name = "", labels = c("Urban", "Forest"), 
                    values = c("#af8dc3", "#7fbf7b")) +
  scale_color_manual(name = "", labels = c("Urban", "Forest"), 
                     values = c("#af8dc3", "#7fbf7b")) +
  scale_y_continuous(breaks = seq(0,10,2), labels = seq(0,10,2))


ggsave(filename = "./plots/Figure 2a.jpeg", 
       plot = fledglings_plot, 
       device = "jpeg", 
       units = "mm",
       width = 90, 
       height = 90)  

##
## simplified x lab
fledglings_plot_no_lab <- fledglings_plot +
  labs(x = "Female chronotype", 
       y = "Number of fledglings") 

ggsave(filename = "./plots/Figure 2a - simplified label.jpeg", 
       plot = fledglings_plot_no_lab, 
       device = "jpeg", 
       units = "mm",
       width = 90, 
       height = 90)  
#####

##
##
##### Analysis for total brood failure #####
##
##

##
## clean previous model objects
rm(list = c("m1","m2","m3","m4","m5","mam"))


data$total_failure <- ifelse(data$fledglings == 0, 1, 0) # binomial variable for full failure
full_failure <- glmer(total_failure ~

                        poly(hatching_date_jul,2)[,1] +
                        poly(hatching_date_jul,2)[,2] +
                        scale(chronotype)+
                        area+
                        clutch_size+
                        (1|site) +
                        (1|year) +
                        (1|box), 
                      family = "binomial",
                      na.action = "na.fail",
                      data=data)
summary(full_failure)
drop1(full_failure, test = "Chisq")

#####

##
##
##### Table of results A6 #####
##
##

## base table
table_fledglings_failure00 <- full_failure %>%
  tbl_regression(intercept = T,
                 label = list(
                   `(Intercept)` = "Intercept",
                   `poly(hatching_date_jul, 2)[, 1]` = "Hatching date1",
                   `poly(hatching_date_jul, 2)[, 2]` = "Hatching date2",
                   `clutch_size` = "Clutch size",
                   `scale(chronotype)` = "Female chronotype",
                   `area` = "Habitat"),
                 estimate_fun = ~ style_number(.x, digits = 2))

## add features
table_fledglings_failure <- table_fledglings_failure00 %>% 
  add_global_p(anova_fun = drop1_output) %>% 
  bold_p(t = 0.05) %>% 
  bold_labels() %>%
  italicize_levels() %>% 
  modify_table_body(fun = function(.){
    output <- left_join(x = .,
                        y = drop1_output(x=full_failure) %>% 
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
gtsave(table_fledglings_failure, "./tables/TABLE A6.html")

#####


##
##
##### Analysis for number of fledglings excluding total failures #####
##
##

# data excluding total failures
dataledglings <- data %>% 
  filter(fledglings != 0)

## model
full_fledglings_no0 <- lmer(fledglings ~

                          poly(hatching_date_jul,2)[,1] +
                          poly(hatching_date_jul,2)[,2] +
                          scale(chronotype)+
                          area+
                          clutch_size+
                          (1|year) +
                          (1|site) +
                          (1|box), 
                        REML = F,
                        na.action = "na.fail",
                        data=dataledglings)

##
## histogram
summary(full_fledglings_no0)
hist(residuals(full_fledglings_no0))
drop1(full_fledglings_no0, test = "Chisq")

##
## checking model convergence
check_model(full_fledglings_no0, panel = T)
normality <- check_normality(full_fledglings_no0)
hist <- plot(normality)
qqplot <- plot(normality, type = "qq")

## saving plots for reviewers
ggsave(filename = "./plots/fledging_model_normality.png", 
       plot = hist, 
       height = 100, 
       width = 100,
       device = "png", 
       units = "mm")

ggsave(filename = "./plots/fledging_model_normality_qq.png", 
       plot = qqplot, 
       height = 100, 
       width = 100,
       device = "png", 
       units = "mm")





##### 

##
##
##### Table of results A7 #####
##
##

## base table
table_fledglings_no000 <- full_fledglings_no0 %>%
  tbl_regression(intercept = T,
                 label = list(
                   `(Intercept)` = "Intercept",
                   `poly(hatching_date_jul, 2)[, 1]` = "Hatching date1",
                   `poly(hatching_date_jul, 2)[, 2]` = "Hatching date2",
                   `clutch_size` = "Clutch size",
                   `scale(chronotype)` = "Female chronotype",
                   `area` = "Habitat"),
                 estimate_fun = ~ style_number(.x, digits = 2))

## add features
table_fledglings_no0 <- table_fledglings_no000 %>% 
  add_global_p(anova_fun = drop1_output) %>% 
  bold_p(t = 0.05) %>% 
  bold_labels() %>%
  italicize_levels() %>% 
  modify_table_body(fun = function(.){
    output <- left_join(x = .,
                        y = drop1_output(x=full_fledglings_no0) %>% 
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
gtsave(table_fledglings_no0, "./tables/TABLE A7.html")

#####

##
##
##### Models for number of fledglings against absolute chronotype and Table A5 #####
##
##
full_fledge_abs_chr_model <- glmer(fledglings ~
                                     poly(hatching_date_jul,2)[,1] +
                                     poly(hatching_date_jul,2)[,2] +
                                     scale(chronotype_abs)+
                                     area+
                                     clutch_size+
                                     (1|year) +
                                     (1|site) +
                                     (1|box), 
                                   family = "poisson",
                                   control = glmerControl(optimizer = "bobyqa"),
                                   na.action = "na.fail",
                                   data=data)
summary(full_fledge_abs_chr_model)

##
## Table of results A5

## base table
table_fledglings_abs_chr00 <- full_fledge_abs_chr_model %>%
  tbl_regression(intercept = T,
                 label = list(
                   `(Intercept)` = "Intercept",
                   `poly(hatching_date_jul, 2)[, 1]` = "Hatching date1",
                   `poly(hatching_date_jul, 2)[, 2]` = "Hatching date2",
                   `clutch_size` = "Clutch size",
                   `scale(chronotype_abs)` = "Clock female chronotype",
                   `area` = "Habitat"),
                 estimate_fun = ~ style_number(.x, digits = 2))

## add features
table_fledglings_abs_chr <- table_fledglings_abs_chr00 %>% 
  add_global_p(anova_fun = drop1_output) %>% 
  bold_p(t = 0.05) %>% 
  bold_labels() %>%
  italicize_levels() %>% 
  modify_table_body(fun = function(.){
    output <- left_join(x = .,
                        y = drop1_output(x=full_fledge_abs_chr_model) %>% 
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
gtsave(table_fledglings_abs_chr, "./tables/TABLE A5.html")

#####

