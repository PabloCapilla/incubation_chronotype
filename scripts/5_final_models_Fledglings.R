##
##### libraries #####
##
pacman::p_load(lubridate, lme4, MuMIn, dplyr, ggplot2, ggridges, extrafont, MCMCglmm)
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
data$hatching_date_jul <- yday(dmy(data$hatching_date))

head(data)

##
## data for model
chrono_box <- data %>%
  filter(!is.na(chronotype)) %>%
  group_by(year, box) %>% 
  summarise(chronotype_abs = mean(activity_onset.1))


data_f <- data %>%
  filter(!is.na(chronotype)) %>%
  filter(!is.na(fledglings)) %>% 
  group_by(year, box) %>% 
  filter(row_number() == 1) %>% 
  mutate(failed_n = hatchlings - fledglings) %>% 
  left_join(., y = chrono_box, by = c("year", "box"))


## sample sizes
data_f %>% 
  group_by(year, area) %>% 
  summarise(n_obs = n())

data_f %>% 
  group_by(box, year, area) %>% 
  filter(row_number() == 1) %>% 
  group_by(year, area) %>% 
  summarise(n_obs = n())


##
##
##### Raw data plot #####
##
##
ggplot(data = data_f, 
       aes(y = fledglings, 
           x = chronotype_abs,
           fill = area,
           color = area)) +
  geom_point() +
  stat_smooth(method = "glm", 
              method.args = c(family = "poisson"), 
              formula = y ~ x)


ggplot(data = data_f, 
       aes(y = hatching_date_jul, 
           x = chronotype,
           fill = area,
           color = area)) +
  geom_point() +
  stat_smooth(method = "glm", 
              method.args = c(family = "poisson"), 
              formula = y ~ x)


##
##
##### Model for nestling weight #####
##
##
full_f <- glmer(fledglings ~
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
                data=data_f)
summary(full_f)

# model diagnostics
plot(full_f)
hist(residuals(full_f))

residuals <- DHARMa::simulateResiduals(full_f, n = 1000, refit = F)
DHARMa::testUniformity(residuals)
DHARMa::testDispersion(residuals)

# aic model selection
model_f_aic <- dredge(full_f,
                           subset = dc("poly(hatching_date_jul, 2)[, 1]", 
                                       "poly(hatching_date_jul, 2)[, 2]"),
                           trace = 3, 
                           rank = "AIC")   # one can use 'BIC' as well if wanted
subset(model_f_aic, delta < 6)
subset(model_f_aic, delta < 6 & !nested(.)) 




##
##
##### Model selection LRT #####
##
##
drop1(full_f, test = "Chisq")
m1 <- update(full_f, . ~ . - area:poly(hatching_date_jul, 2)[, 2])

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



# mam
m_mam <- m5
drop1(m_mam, test = "Chisq")
summary(m_mam)

CI_onset <- confint(m_mam, 
                    level = 0.95, 
                    method = "boot", 
                    nsim = 500, 
                    boot.type = "norm")


##
## predictions plot
top_model_f_poisson <- glmer(fledglings ~
                           chronotype+
                           area+
                           clutch_size +
                             (1|site) +
                             (1|box) +
                             (1|year),
                         na.action = "na.fail",
                         family = "poisson",
                         glmerControl(optimizer = "bobyqa"),
                         data=data_f)
summary(top_model_f_poisson)

CI_onset <- confint(top_model_f_poisson, 
                    level = 0.95, 
                    method = "boot", 
                    nsim = 500, 
                    boot.type = "norm")








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
## predictions plot
df_pred <- expand.grid(chronotype = seq(min(data_f$chronotype), 
                                        max(data_f$chronotype), 0.5),
                       area = c("City", "Forest"))
df_pred$clutch_size <- ifelse(df_pred$area == "City", 
                              mean(data_f$clutch_size[data_f$area == "City"]),
                              mean(data_f$clutch_size[data_f$area == "Forest"]))
df_pred$prediction_link <- predict(top_model_f_poisson,
                                   df_pred, 
                                   type = "link",
                                   re.form = NA)

# plot only data in range
data_f %>% 
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
pvar1 <- diag(mm %*% tcrossprod(vcov(top_model_f_poisson),mm))
cmult <- 1 ## 1 SE
df_pred <- data.frame(
  df_pred
  , plo_link = df_pred$prediction_link-cmult*sqrt(pvar1)
  , phi_link = df_pred$prediction_link+cmult*sqrt(pvar1)
)
df_pred$prediction <- exp(df_pred$prediction_link)
df_pred$plo <- exp(df_pred$plo_link)
df_pred$phi <- exp(df_pred$phi_link)


fledglings_plot <- ggplot(data = data_f, aes(x = chronotype, 
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
  labs(#x = expression(atop("Female chronotype", 
            #               "(relative onset of activity [minutes after sunrise])")), 
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
             

ggsave(filename = "./plots/2_Reproduction/fledglings_pred_xaxis2.jpeg", 
       plot = fledglings_plot, 
       device = "jpeg", 
       units = "mm",
       width = 90, 
       height = 90)  




##
##
##### total brood failure analysis #####
##
##
data_f$total_failure <- ifelse(data_f$fledglings == 0, 1, 0)
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
                      data=data_f)
summary(full_failure)

# model diagnostics
plot(full_failure)
hist(residuals(full_failure))


# aic model selection
model_failure_aic <- dredge(full_failure,
                      subset = dc("poly(hatching_date_jul, 2)[, 1]", 
                                  "poly(hatching_date_jul, 2)[, 2]"),
                      trace = 3, 
                      rank = "AIC")   # one can use 'BIC' as well if wanted
subset(model_failure_aic, delta < 6)
subset(model_failure_aic, delta < 6 & !nested(.)) 




##
##
##### Model selection LRT #####
##
##
drop1(full_failure, test = "Chisq")
m1 <- update(full_failure, . ~ . - scale(clutch_size))

#1
drop1(m1, test = "Chisq")
m2 <- update(m1, . ~ . - poly(hatching_date_jul, 2)[, 2]:area)

#2
drop1(m2, test = "Chisq")
m3 <- update(m2, . ~ . - poly(hatching_date_jul, 2)[, 2])

#3
drop1(m3, test = "Chisq")
m4 <- update(m3, . ~ . - scale(chronotype):area)

#4
drop1(m4, test = "Chisq")
m5 <- update(m4, . ~ . - poly(hatching_date_jul, 2)[, 1]:area)

#5
drop1(m5, test = "Chi")
m6 <- update(m5, . ~ . - poly(hatching_date_jul, 2)[, 1])
#5
drop1(m6, test = "Chisq")
m7 <- update(m6, . ~ . - scale(chronotype))

#6
drop1(m7, test = "Chisq")
m8 <- update(m7, . ~ . - area)

summary




# mam
m_mam <- m9
drop1(m_mam, test = "Chisq")
summary(m_mam)


m_mam_glm <- glm(total_failure ~
                   as.factor(year),
                 family = "binomial",
                 na.action = "na.fail",
                 data=data_f)

CI_onset <- confint(m_mam_glm, 
                    level = 0.95, 
                    method = "boot", 
                    nsim = 500, 
                    boot.type = "norm")



# p for non significant terms
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
##### number of fledglings #####
##
##
data_fledglings <- data_f %>% 
  filter(fledglings != 0)
table(data_fledglings$fledglings)
hist(data_fledglings$fledglings)

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
                         data=data_fledglings)
summary(full_fledglings)
hist(residuals(full_fledglings))
shapiro.test(residuals(full_fledglings))



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

# if we drop hatching date
drop1(m6, test = "Chisq")


# mam
m_mam <- m6 #m6
drop1(m_mam, test = "Chisq")
summary(m_mam)



##
## predictions plot
top_model_fledglings <- lmer(fledglings ~
                               chronotype+
                               clutch_size +
                               (1|site) +
                               (1|year) +
                               (1|box), 
                             REML=F,
                             na.action = "na.fail",
                             data=data_fledglings)
summary(top_model_fledglings)
r.squaredGLMM(top_model_fledglings)
drop1(top_model_fledglings, test = "Chisq")

AIC(top_model_fledglings)

CI_onset <- confint(top_model_fledglings, 
                    level = 0.95, 
                    method = "boot", 
                    nsim = 500, 
                    boot.type = "norm")


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

m_chr_int <- update(m_mam, . ~ . + area + chronotype : area)
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


## plot
##
## predictions plot
df_pred <- expand.grid(chronotype = seq(min(data_fledglings$chronotype), 
                                        max(data_fledglings$chronotype), 0.5),
                       area = c("City", "Forest"))
df_pred$hatching_date_jul <- ifelse(df_pred$area == "City", 
                                    mean(data_fledglings$hatching_date_jul[data_fledglings$area == "City"]),
                                    mean(data_fledglings$hatching_date_jul[data_fledglings$area == "Forest"]))
df_pred$clutch_size <- ifelse(df_pred$area == "City", 
                                    mean(data_fledglings$clutch_size[data_fledglings$area == "City"]),
                                    mean(data_fledglings$clutch_size[data_fledglings$area == "Forest"]))


df_pred$prediction <- predict(top_model_fledglings,
                                   df_pred, 
                                   type = "response",
                                   re.form = NA)

# plot only data in range
data_fledglings %>% 
  group_by(area) %>% 
  summarise(min_chr = min(chronotype),
            max_chr = max(chronotype))

remove_city <- which((df_pred$chronotype < -21 | 
                        df_pred$chronotype > 50) & 
                       df_pred$area == "City")
remove_forest <- which((df_pred$chronotype < -9.5 | 
                          df_pred$chronotype > 71.6) & 
                         df_pred$area == "Forest")
df_pred <- df_pred[-c(remove_city, remove_forest),]


# need to adapt code to new R version
# SE for mean predicitons
mm <- model.matrix(~ #hatching_date_jul +
                     chronotype+
                     area+
                     clutch_size,
                   data = df_pred)

## or newdat$distance <- mm %*% fixef(fm1)
pvar1 <- diag(mm %*% tcrossprod(vcov(top_model_fledglings),mm))
cmult <- 1 ## 1 SE
df_pred <- data.frame(
  df_pred
  , plo = df_pred$prediction-cmult*sqrt(pvar1)
  , phi = df_pred$prediction+cmult*sqrt(pvar1)
)


fledglings_noZeros_plot <- ggplot(data = data_fledglings, 
                                  aes(x = chronotype, 
                                      y = fledglings,
                                      fill = area,
                                      color = area)) +
  geom_point(alpha = 0.5,
             size = 2.5) +
  facet_grid(.~area) +
  theme_bw() +
  theme(legend.position = "top",
        strip.background = element_blank(),
        strip.text = element_text("Arial", size = 12),
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
  labs(x = "Female's chronotype", 
       y = "Number of fledglings (excluding 0s)") +
  scale_fill_manual(name = "", labels = c("City", "Forest"), 
                    values = c("#3399CC", "#339900")) +
  scale_color_manual(name = "", labels = c("City", "Forest"), 
                     values = c("#3399CC", "#339900")) 


ggsave(filename = "./plots/fledglings_noZeros_pred.jpeg", 
       plot = fledglings_noZeros_plot, 
       device = "jpeg", 
       units = "mm",
       width = 150, 
       height = 120)  
