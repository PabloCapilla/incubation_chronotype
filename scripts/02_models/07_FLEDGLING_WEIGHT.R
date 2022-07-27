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
#' Analysis of fledgling weight
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

##
##
##### Model for nestling weight #####
##
##
full_weight <- lmer(avg_chickweight~ 
                      poly(hatching_date_jul,2)[,1] : area +
                      poly(hatching_date_jul,2)[,2] : area +
                      chronotype:area+

                      poly(hatching_date_jul,2)[,1] +
                      poly(hatching_date_jul,2)[,2] +
                      chronotype+
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
##### Model selection LRT #####
##
##
drop1(full_weight, test = "Chisq")
m1 <- update(full_weight, . ~ . - chronotype:area)

#1
drop1(m1, test = "Chisq")
m2 <- update(m1, . ~ . - area:poly(hatching_date_jul, 2)[, 2])

#2
drop1(m2, test = "Chisq")
m3 <- update(m2, . ~ . - chronotype)

#3
drop1(m3, test = "Chisq")

# mam
m_mam <- m3
summary(m_mam)


##
##
##### Likelihood-ratio test results #####
##
##

# for effects in final model
drop1(m_mam, test = "Chisq")

# p for non significant terms
m_int <- update(m_mam, . ~ . + area:poly(hatching_date_jul, 2)[, 2])
anova(m_int, m_mam, test = "Chisq")

m_chr_int <- update(m_mam, . ~ . + chronotype + chronotype : area)
anova(m_chr_int, update(m_mam, .~.+chronotype), test = "Chisq")

m_chrono <- update(m_mam, . ~ . + chronotype)
anova(m_chrono, m_mam, test = "Chisq")



##
##
##### Final model results #####
##
##

# final model
top_model_weight <- lmer(avg_chickweight~
                           poly(hatching_date_jul,2,raw = T)[,1] : area +
                           poly(hatching_date_jul,2,raw = T)[,1] +
                           poly(hatching_date_jul,2,raw = T)[,2] +
                           area +
                           clutch_size +
                           (1|site) +
                           (1|year) +
                           (1|box),  
                         REML = T,
                         na.action = "na.fail",
                         data=data_weight)
summary(top_model_weight)
drop1(top_model_weight, test = "Chisq")

# 95%CIs
CI_weight <- confint(m_mam, 
                     level = 0.95, 
                     method = "boot", 
                     nsim = 500, 
                     boot.type = "norm")

##
##
##### Plot model predictions #####
##
##

# new dataframe to predict
df_pred <- expand.grid(hatching_date_jul = seq(min(data_weight$hatching_date_jul), 
                                               max(data_weight$hatching_date_jul), 1),
                       area = c("City", "Forest"))
df_pred$clutch_size <- ifelse(df_pred$area == "City", 
                              mean(data_weight$clutch_size[data_weight$area == "City"]),
                              mean(data_weight$clutch_size[data_weight$area == "Forest"]))
df_pred$prediction <- predict(top_model_weight, df_pred, re.form = NA)

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
mm <- model.matrix(~ poly(hatching_date_jul,2,raw = T)[,1] : area +
                     poly(hatching_date_jul,2,raw = T)[,1] +
                     poly(hatching_date_jul,2,raw = T)[,2] +
                     area +
                     clutch_size,
                   data = df_pred)

## or newdat$distance <- mm %*% fixef(fm1)
pvar1 <- diag(mm %*% tcrossprod(vcov(top_model_weight),mm))
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
                                  y = prediction,
                                  shape = area),
              alpha = 0.25,
              color = NA,
              position = position_dodge(width = 0.5)) +
  geom_line(data = df_pred, 
            aes(y = prediction, shape = area), 
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



















