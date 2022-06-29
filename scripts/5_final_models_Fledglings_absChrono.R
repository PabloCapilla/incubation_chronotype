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
  group_by(year, box, area) %>% 
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
                  scale(chronotype_abs):area+

                  poly(hatching_date_jul,2)[,1] +
                  poly(hatching_date_jul,2)[,2] +
                  scale(chronotype_abs) +
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
m3 <- update(m2, . ~ . - scale(chronotype_abs):area)

#3
drop1(m3, test = "Chisq")
m4 <- update(m3, . ~ . - scale(chronotype_abs))

#4
drop1(m4, test = "Chisq")
m5 <- update(m4, . ~ . - poly(hatching_date_jul, 2)[, 2])

#5
drop1(m5, test = "Chisq")
m6 <- update(m5, . ~ . - poly(hatching_date_jul, 2)[, 1])


#6
drop1(m6, test = "Chisq")



# mam
m_mam <- m6
drop1(m_mam, test = "Chisq")
summary(m_mam)


##
## predictions plot
top_model_f_poisson <- glmer(fledglings ~
                               area+
                               clutch_size +
                               (1|site) +
                               (1|box) +
                               (1|year),
                             na.action = "na.fail",
                             family = "poisson",
                             glmerControl(optimizer = "bobyqa"),
                             data=data_f)



# p for non significant terms
#1
m_int2 <- update(m_mam, . ~ . + 
                   poly(hatching_date_jul, 2)[, 1]+
                   poly(hatching_date_jul, 2)[, 2]+
                   area:poly(hatching_date_jul, 2)[, 2])
anova(m_int2, update(m_mam, .~.+
                       poly(hatching_date_jul, 2)[, 1]+
                       poly(hatching_date_jul, 2)[, 2]), test = "Chisq")

#2
m_int1 <- update(m_mam, . ~ . + 
                   poly(hatching_date_jul, 2)[, 1]+
                   area:poly(hatching_date_jul, 2)[, 1])
anova(m_int1, update(m_mam, .~.+
                       poly(hatching_date_jul, 2)[, 1]), test = "Chisq")

#3
m_chr_int <- update(m_mam, . ~ . + scale(chronotype_abs) : area)
anova(m_chr_int, update(m_mam, .~.+
                          scale(chronotype_abs)), 
      test = "Chisq")

#4
m_date2 <- update(m_mam, . ~ . + 
                    poly(hatching_date_jul, 2)[, 1]+
                    poly(hatching_date_jul, 2)[, 2])
anova(m_date2, update(m_mam, .~.+
                        poly(hatching_date_jul, 2)[, 1]), test = "Chisq")

#5
m_date1 <- update(m_mam, . ~ . + 
                    poly(hatching_date_jul, 2)[, 1])
anova(m_date1, m_mam, test = "Chisq")

#6
m_chr1 <- update(m_mam, . ~ . + 
                    scale(chronotype_abs))
anova(m_chr1, m_mam, test = "Chisq")
