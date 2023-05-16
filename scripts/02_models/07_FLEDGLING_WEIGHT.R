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
#' Analysis of fledgling weight
#' 
##
##

##
##
##### libraries #####
##
##
pacman::p_load(openxlsx, gt, gtsummary,
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
weight_box <- data00 %>%
  filter(!is.na(chronotype)) %>%
  group_by(year, box) %>% 
  summarise(chronotype_abs = mean(activity_onset_abs))

## calculate fledgling weight per nest
data_weight <- data00 %>%
  filter(!is.na(chronotype)) %>%
  filter(!is.na(avg_chickweight)) %>% 
  group_by(year, box) %>% 
  filter(row_number() == 1) %>% 
  left_join(., y = weight_box, by = c("year", "box"))

####

##
##
##### Model for nestling weight #####
##
##
full_weight <- lmer(avg_chickweight~ 
                      poly(hatching_date_jul,2)[,1] : area +
                      poly(hatching_date_jul,2)[,2] : area +
                      scale(chronotype):area+

                      poly(hatching_date_jul,2)[,1] +
                      poly(hatching_date_jul,2)[,2] +
                      scale(chronotype)+
                      area+
                      clutch_size+
                      (1|site) +
                      (1|year) +
                      (1|box), 
                    REML = F,
                    na.action = "na.fail",
                    data=data_weight)
summary(full_weight)

# model diagnostics
hist(residuals(full_weight))

##
##
##### Are interactions significant? #####
##
##

## 1
drop1(full_weight, test = "Chisq")

## interaction 'area:poly(inc_start_aprildays,2)[,2]'
anova(full_weight, 
      update(full_weight, 
             . ~ . - area:poly(hatching_date_jul, 2)[, 2]),
      test = "LRT")

## interaction 'scale(chronotype):area'
anova(full_weight, 
      update(full_weight, 
             . ~ . 
             - scale(chronotype):area),
      test = "LRT")


## removing both interactions 
full_model <- update(full_weight, 
                     . ~ . 
                     - scale(chronotype):area 
                     - poly(hatching_date_jul,2)[,2] : area)

## comparison model without any interaction against initial
anova(full_weight, 
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

#####

##
##
##### Table of results 4 #####
##
##

## base table
table_weight00 <- full_model %>%
  tbl_regression(intercept = T,
                 label = list(
                   `(Intercept)` = "Intercept",
                   `poly(hatching_date_jul, 2)[, 1]` = "Hatching date1",
                   `poly(hatching_date_jul, 2)[, 2]` = "Hatching date2",
                   `clutch_size` = "Clutch size",
                   `scale(chronotype)` = "Female chronotype",
                   `area` = "Habitat", 
                   `poly(hatching_date_jul, 2)[, 1]:areaForest` = "Hatching date1 x Habitat"),
                 estimate_fun = ~ style_number(.x, digits = 2))

## add features
table_weight <- table_weight00 %>% 
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
gtsave(table_weight, "./tables/TABLE 4.html")

#####

##
##
##### Plot model predictions #####
##
##

##
## model specification for predictions
full_model_predictions <- lmer(avg_chickweight~ 
                                 hatching_date_jul : area +

                                 hatching_date_jul +
                                 I(hatching_date_jul^2) +
                                 chronotype+
                                 area+
                                 clutch_size+
                                 (1|site) +
                                 (1|year) +
                                 (1|box), 
                               REML = F,
                               na.action = "na.fail",
                               data=data_weight)


# new dataframe to predict
df_pred <- expand.grid(hatching_date_jul = seq(min(data_weight$hatching_date_jul), 
                                               max(data_weight$hatching_date_jul), 1),
                       area = c("City", "Forest"))
df_pred$clutch_size <- ifelse(df_pred$area == "City", 
                              mean(data_weight$clutch_size[data_weight$area == "City"]),
                              mean(data_weight$clutch_size[data_weight$area == "Forest"]))
df_pred$chronotype <- ifelse(df_pred$area == "City", 
                              mean(data_weight$chronotype[data_weight$area == "City"]),
                              mean(data_weight$chronotype[data_weight$area == "Forest"]))
df_pred$prediction <- predict(full_model_predictions, df_pred, re.form = NA)

# plot only data in range
data_weight %>% 
  group_by(area) %>% 
  summarise(min_jul = min(hatching_date_jul),
            max_jul = max(hatching_date_jul))

remove_city <- which((df_pred$hatching_date_jul < 125 | 
                        df_pred$hatching_date_jul > 141) & 
                       df_pred$area == "City")
remove_forest <- which((df_pred$hatching_date_jul < 134 | 
                          df_pred$hatching_date_jul > 156) & 
                         df_pred$area == "Forest")
df_pred <- df_pred[-c(remove_city, remove_forest),]


# SE for mean predictions
mm <- model.matrix(~ hatching_date_jul : area +
                     
                     hatching_date_jul +
                     I(hatching_date_jul^2) +
                     chronotype+
                     area+
                     clutch_size,
                   data = df_pred)

## or newdat$distance <- mm %*% fixef(fm1)
pvar1 <- diag(mm %*% tcrossprod(vcov(full_model_predictions),mm))
cmult <- 1 ## 1 SE
df_pred <- data.frame(
  df_pred
  , plo = df_pred$prediction-cmult*sqrt(pvar1)
  , phi = df_pred$prediction+cmult*sqrt(pvar1)
)


# figure plot against hatching date
weight_plot <- ggplot(data = data_weight, aes(x = hatching_date_jul, 
                                               y = avg_chickweight,
                                               fill = area,
                                               color = area)) +
  geom_point(position = position_jitterdodge(dodge.width = 0.5,
                                             jitter.width = 0.25),
             alpha = 0.5,
             size = 2,
             aes(shape = area), 
             color = "black") +
  theme_bw() +
  theme(legend.position = "top",
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
  labs(x = "Hatching day (days after April 1)", 
       y = "Weight nestlings day 13 (g)") +
  scale_x_continuous(breaks = seq(125, 156, 5), labels = seq(35, 66, 5)) +
  scale_shape_manual(name = "",
                     values = c(21,24),
                     labels = c("Urban", "Forest")) + 
  scale_fill_manual(name = "", 
                    labels = c("Urban", "Forest"), 
                    values = c("#af8dc3", "#7fbf7b")) +
  scale_color_manual(name = "", 
                     labels = c("Urban", "Forest"), 
                     values = c("#af8dc3", "#7fbf7b")) 

ggsave(filename = "./plots/Figure 2b.jpeg", 
       plot = weight_plot, 
       device = "jpeg", 
       units = "mm",
       width = 90, 
       height = 90)  



















