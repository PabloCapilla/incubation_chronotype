###
###
#' 
#' Script for:
#' Reproductive fitness is associated with female chronotype in a songbird
#' Womack, et al. 
#' Preprint: 10.1101/2022.07.01.498449v1
#' 
#' Latest update: 2022-10-25
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
                             scale(clutch_size)+
                             (1|site) +
                             (1|year) +
                             (1|box), 
                           family = "poisson",
                           na.action = "na.fail",
                           data=data)
summary(full_fledge_model)
# model diagnostics
residuals <- DHARMa::simulateResiduals(full_fledge_model, n = 1000, refit = F)
DHARMa::testDispersion(residuals)

##
## full model without interactions (i.e., not doing model selection on single terms)
full_fledge_model_b <- glmer(fledglings ~

                             poly(hatching_date_jul,2)[,1] +
                             poly(hatching_date_jul,2)[,2] +
                             scale(chronotype)+
                             area+
                             scale(clutch_size)+
                             (1|site) +
                             (1|year) +
                             (1|box), 
                           family = "poisson",
                           na.action = "na.fail",
                           data=data)
summary(full_fledge_model_b)
drop1(full_fledge_model_b, test = "Chisq")

##
##
##### Model selection LRT #####
##
##
drop1(full_fledge_model, test = "Chisq")
m1 <- update(full_fledge_model, . ~ . - area:poly(hatching_date_jul, 2)[, 2])

#1
drop1(m1, test = "Chisq")
m2 <- update(m1, . ~ . - poly(hatching_date_jul, 2)[, 1]:area)

#2
drop1(m2, test = "Chisq")
m3 <- update(m2, . ~ . - scale(chronotype):area)

#3
drop1(m3, test = "Chisq")
m4 <- update(m3, . ~ . - poly(hatching_date_jul, 2)[, 2])

#4
drop1(m4, test = "Chisq")
m5 <- update(m4, . ~ . - poly(hatching_date_jul, 2)[, 1])

#5
drop1(m5, test = "Chisq")

# Final model
m_mam <- m5
summary(m_mam)

##
##
##### Likelihood-ratio test results #####
##
##

# for effects in final model
drop1(m_mam, test = "Chisq")

# p for non significant terms
m_int2 <- update(m_mam, . ~ . + 
                   poly(hatching_date_jul, 2)[, 1]+
                   poly(hatching_date_jul, 2)[, 2]+
                   area:poly(hatching_date_jul, 2)[, 2])
anova(m_int2, update(m_mam, .~.+
                       poly(hatching_date_jul, 2)[, 1]+
                       poly(hatching_date_jul, 2)[, 2]), test = "Chisq")

m_int1 <- update(m_mam, . ~ . + 
                   poly(hatching_date_jul, 2)[, 1]+
                   area:poly(hatching_date_jul, 2)[, 1])
anova(m_int1, update(m_mam, .~.+
                       poly(hatching_date_jul, 2)[, 1]), test = "Chisq")

m_chr_int <- update(m_mam, . ~ . + scale(chronotype) : area)
anova(m_chr_int, m_mam, test = "Chisq")

m_date2 <- update(m_mam, . ~ . + 
                    poly(hatching_date_jul, 2)[, 1]+
                    poly(hatching_date_jul, 2)[, 2])
anova(m_date2, update(m_mam, .~.+
                        poly(hatching_date_jul, 2)[, 1]), test = "Chisq")

m_date1 <- update(m_mam, . ~ . + 
                    poly(hatching_date_jul, 2)[, 1])
anova(m_date1, m_mam, test = "Chisq")

##
##
##### Final model results #####
##
##

# final model
top_fledge_model <- glmer(fledglings ~
                               chronotype+
                               area+
                               clutch_size +
                               (1|site) +
                               (1|box) +
                               (1|year),
                             na.action = "na.fail",
                             family = "poisson",
                             glmerControl(optimizer = "bobyqa"),
                             data=data)
summary(top_fledge_model)

# 95%CIs
CI_onset <- confint(top_fledge_model, 
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
df_pred <- expand.grid(chronotype = seq(min(data$chronotype), 
                                        max(data$chronotype), 0.5),
                       area = c("City", "Forest"))
df_pred$clutch_size <- ifelse(df_pred$area == "City", 
                              mean(data$clutch_size[data$area == "City"]),
                              mean(data$clutch_size[data$area == "Forest"]))
df_pred$prediction_link <- predict(top_fledge_model,
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
mm <- model.matrix(~ chronotype+
                     area+
                     clutch_size,
                   data = df_pred)

## or newdat$distance <- mm %*% fixef(fm1)
pvar1 <- diag(mm %*% tcrossprod(vcov(top_fledge_model),mm))
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
                     labels = c("City", "Forest")) +  
  scale_fill_manual(name = "", labels = c("City", "Forest"), 
                    values = c("#af8dc3", "#7fbf7b")) +
  scale_color_manual(name = "", labels = c("City", "Forest"), 
                     values = c("#af8dc3", "#7fbf7b")) +
  scale_y_continuous(breaks = seq(0,10,2), labels = seq(0,10,2))


ggsave(filename = "./plots/Figure 2a.jpeg", 
       plot = fledglings_plot, 
       device = "jpeg", 
       units = "mm",
       width = 90, 
       height = 90)  



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
                        poly(hatching_date_jul,2)[,1] : area +
                        poly(hatching_date_jul,2)[,2] : area +
                        scale(chronotype):area+
                        
                        poly(hatching_date_jul,2)[,1] +
                        poly(hatching_date_jul,2)[,2] +
                        scale(chronotype)+
                        area+
                        scale(clutch_size)+
                        (1|site) +
                        (1|year) +
                        (1|box), 
                      family = "binomial",
                      na.action = "na.fail",
                      data=data)
summary(full_failure)

##
##
##### Model selection LRT #####
##
##
drop1(full_failure, test = "Chisq")
ffm1 <- update(full_failure, . ~ . - scale(clutch_size))

#1
drop1(ffm1, test = "Chisq")
ffm2 <- update(ffm1, . ~ . - poly(hatching_date_jul, 2)[, 2]:area)

#2
drop1(ffm2, test = "Chisq")
ffm3 <- update(ffm2, . ~ . - poly(hatching_date_jul, 2)[, 2])

#3
drop1(ffm3, test = "Chisq")
ffm4 <- update(ffm3, . ~ . - scale(chronotype):area)

#4
drop1(ffm4, test = "Chisq")
ffm5 <- update(ffm4, . ~ . - poly(hatching_date_jul, 2)[, 1]:area)

#5
drop1(ffm5, test = "Chi")
ffm6 <- update(ffm5, . ~ . - poly(hatching_date_jul, 2)[, 1])

#5
drop1(ffm6, test = "Chisq")
ffm7 <- update(ffm6, . ~ . - scale(chronotype))

#6
drop1(ffm7, test = "Chisq")
ffm8 <- update(ffm7, . ~ . - area)


# mam
m_mam <- ffm8
drop1(m_mam, test = "Chisq")
summary(m_mam)


##
##
##### Likelihood-ratio test results #####
##
##
m_int2 <- update(m_mam, . ~ . + 
                   poly(hatching_date_jul, 2)[, 1]+
                   poly(hatching_date_jul, 2)[, 2]+
                   area:poly(hatching_date_jul, 2)[, 2]) # tested term
anova(m_int2, update(m_mam, .~.+
                       poly(hatching_date_jul, 2)[, 1]+
                       poly(hatching_date_jul, 2)[, 2]), test = "Chisq")

m_int1 <- update(m_mam, . ~ . + 
                   poly(hatching_date_jul, 2)[, 1]+
                   area:poly(hatching_date_jul, 2)[, 1])
anova(m_int1, update(m_mam, .~.+
                       poly(hatching_date_jul, 2)[, 1]), test = "Chisq")

m_chr_int <- update(m_mam, . ~ . + area + chronotype + chronotype : area)
anova(m_chr_int, update(m_mam, .~.+
                          area + chronotype), test = "Chisq")

m_chr <- update(m_mam, . ~ . + chronotype)
anova(m_chr, m_mam, test = "Chisq")

m_area <- update(m_mam, . ~ . + area)
anova(m_area, m_mam, test = "Chisq")

m_clutchsize <- update(m_mam, . ~ . + clutch_size)
anova(m_clutchsize, m_mam, test = "Chisq")

m_date2 <- update(m_mam, . ~ . + 
                    poly(hatching_date_jul, 2)[, 1]+
                    poly(hatching_date_jul, 2)[, 2])
anova(m_date2, update(m_mam, .~.+
                        poly(hatching_date_jul, 2)[, 1]), test = "Chisq")

m_date1 <- update(m_mam, . ~ . + 
                    poly(hatching_date_jul, 2)[, 1])
anova(m_date1, m_mam, test = "Chisq")


##
##
##### Analysis for number of fledglings excluding total failures #####
##
##

# data excluding total failures
dataledglings <- data %>% 
  filter(fledglings != 0)

## model
full_fledglings <- lmer(fledglings ~
                          poly(hatching_date_jul,2)[,1] : area +
                          poly(hatching_date_jul,2)[,2] : area +
                          scale(chronotype):area+
                          
                          poly(hatching_date_jul,2)[,1] +
                          poly(hatching_date_jul,2)[,2] +
                          scale(chronotype)+
                          area+
                          scale(clutch_size)+
                          (1|site) +
                          (1|year) +
                          (1|box), 
                        REML = F,
                        na.action = "na.fail",
                        data=dataledglings)
summary(full_fledglings)
hist(residuals(full_fledglings))
shapiro.test(residuals(full_fledglings))


##
## clean previous model objects
rm(list = c("ffm1","ffm2","ffm3","ffm4","ffm5", "ffm6", "ffm7", "ffm8", "mam"))

##
##
##### Model selection LRT #####
##
##
drop1(full_fledglings, test = "Chisq")
m1 <- update(full_fledglings, . ~ . - area:scale(chronotype))

#1
drop1(m1, test = "Chisq")
m2 <- update(m1, . ~ . - poly(hatching_date_jul, 2)[, 1]:area)

#2
drop1(m2, test = "Chisq")
m3 <- update(m2, . ~ . - poly(hatching_date_jul, 2)[, 2]:area)

#3
drop1(m3, test = "Chisq")
m4 <- update(m3, . ~ . - poly(hatching_date_jul, 2)[, 2])

#4
drop1(m4, test = "Chisq")
m5 <- update(m4, . ~ . - poly(hatching_date_jul, 2)[, 1])

#5
drop1(m5, test = "Chisq")
m6 <- update(m5, . ~ . - area)

# 6
drop1(m6, test = "Chisq")

# Final model
m_mam <- m6 #m6
drop1(m_mam, test = "Chisq")
summary(m_mam)


##
##
##### Likelihood-ratio test results #####
##
##

# p for non significant terms
m_int2 <- update(m_mam, . ~ . + 
                   area +
                   poly(hatching_date_jul, 2)[, 1]+
                   poly(hatching_date_jul, 2)[, 2]+
                   area:poly(hatching_date_jul, 2)[, 2])
anova(m_int2, update(m_mam, .~.+
                       area +
                       poly(hatching_date_jul, 2)[, 1]+
                       poly(hatching_date_jul, 2)[, 2]), test = "Chisq")

m_int1 <- update(m_mam, . ~ . + 
                   area +
                   poly(hatching_date_jul, 2)[, 1] +
                   area:poly(hatching_date_jul, 2)[, 1])
anova(m_int1, update(m_mam, . ~ . + 
                       area +
                       poly(hatching_date_jul, 2)[, 1]), test = "Chisq")

m_chr_int <- update(m_mam, . ~ . + area + scale(chronotype) : area)
anova(m_chr_int, update(m_mam, . ~ . + 
                          area), test = "Chisq")


m_date2 <- update(m_mam, . ~ . + 
                    poly(hatching_date_jul, 2)[, 1] +
                    poly(hatching_date_jul, 2)[, 2])
anova(m_date2, update(m_mam, . ~ . + 
                        poly(hatching_date_jul, 2)[, 1]), test = "Chisq")

m_date1 <- update(m_mam, . ~ . + 
                    poly(hatching_date_jul, 2)[, 1])
anova(m_date1, m_mam, test = "Chisq")

m_habitat <- update(m_mam, . ~ . + area)
anova(m_habitat, m_mam, test = "Chisq")



##
##
##### Final model results #####
##
##

# final model
top_model_fledglings <- lmer(fledglings ~
                               chronotype+
                               clutch_size +
                               (1|site) +
                               (1|year) +
                               (1|box), 
                             REML=F,
                             na.action = "na.fail",
                             data=dataledglings)
summary(top_model_fledglings)

# 95%CIs
CI_onset <- confint(top_model_fledglings, 
                    level = 0.95, 
                    method = "boot", 
                    nsim = 500, 
                    boot.type = "norm")

#