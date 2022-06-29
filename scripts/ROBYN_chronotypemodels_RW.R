##
##
##### Original 'chronotypemodels_RW.R' modified by PCL, JAN 2020 #####
##
##

##### libraries #####
library(ggplot2)
library(lme4)
library(lmtest)

library(dplyr)
library(lubridate)
#####


##
##
##### data for incubation #####
##
##
buttondata <- read.csv (file="./data/ROBYN_incubationdata.csv", 
                        na.strings=c("", "NA"))
head(buttondata)




##
## Is the code below really needed for the analysis??
##

#gsubset <- subset(buttondata, select=c(box, date_aprildays, year, area, 
#                                      clutch_size, day_before_hatch, 
#                                     diff_sunrise, diff_sunset, meantemp))

#remove missing rows, create index
#Aindex <- complete.cases(gsubset)
#table(Aindex)


#  subset[row,column]
Asubset <- buttondata #gsubset[Aindex,]
Asubset$year <- as.factor(Asubset$year) 
Asubset$diff_sunrise <- as.numeric(Asubset$diff_sunrise)
Asubset$day_before_hatch <- as.numeric(Asubset$day_before_hatch)
str(Asubset)



##
##
##### models onset of activity #####
##
##

## NOTE: if you want to use 'lrtest', you need to fit models using 'REML=F'


# global model - note the differnece in notation, see my e-mail to Ciara for a details on a warning that your models were producing
model1_onset <- lmer(diff_sunrise ~ 
                 area:poly(day_before_hatch,2, raw = TRUE)[,2] + # interaction area - quadratic 'day_before_hatch' 
                 area:poly(day_before_hatch,2,raw = TRUE)[,1] +  # interaction area - linear 'day_before_hatch' 
                 
                 poly(day_before_hatch,2,raw = TRUE)[,2] +       # quadratic 'day_before_hatch' 
                 poly(day_before_hatch,2,raw = TRUE)[,1] +       # linear 'day_before_hatch' 
                 meantemp + 
                 year + 
                 date_aprildays + 
                 clutch_size + 
                 area + 
                 (1|box), 
               REML = F,
               na.action = "na.fail",
               data=Asubset) #full model

summary(model1_onset)

## 'drop1' takes a model (e.g. 'model1_onset' in the case above) and makes every LRT comparison 
## removing each term one at a time. For example, if you look at the results that it produces, the statistic and 
## p-value associated with, say, 'meantemp' refer to a LRT test of 'model1_onset' against a similar nested model
## without the 'meantemp' variable. And it does so for all predictors! (so it'd save you a lot of time and code).

drop1(model1_onset, test = "Chisq")

## Look at the results and remove that predictor with highest p-value. In your case 'clutch_size', then create a model
## without 'clutch_size':


# model 2 - 'clutch size' removed
model2_onset <- lmer(diff_sunrise ~ 
                 area:poly(day_before_hatch,2, raw = TRUE)[,2] + 
                 area:poly(day_before_hatch,2, raw = TRUE)[,1] +  
                 
                 poly(day_before_hatch,2, raw = TRUE)[,2] +
                 poly(day_before_hatch,2, raw = TRUE)[,1] +
                 meantemp + 
                 year + 
                 date_aprildays + 
                 #clutch_size + 
                 area + 
                 (1|box),
               REML = F,
               data=Asubset) 

## Use 'drop1' again and repeat the process removing the predictor with highest p-value:
drop1(model2_onset, test = "Chisq")


# model 3 - 'area:poly(day_before_hatch,2)[,2]' removed
model3_onset <- lmer(diff_sunrise ~ 
                 #area:poly(day_before_hatch,2, raw = TRUE)[,2] + 
                 area:poly(day_before_hatch,2, raw = TRUE)[,1] +  
                 
                 poly(day_before_hatch,2, raw = TRUE)[,2] +
                 poly(day_before_hatch,2, raw = TRUE)[,1] +
                 meantemp + 
                 year + 
                 date_aprildays + 
                 #clutch_size + 
                 area + 
                 (1|box),
               REML = F,
               data=Asubset) 

## Use 'drop1' again and repeat the process removing the predictor with highest p-value:
drop1(model3_onset, test = "Chisq")


# model 4 - 'year' removed
model4_onset <- lmer(diff_sunrise ~ 
                 #area:poly(day_before_hatch,2, raw = TRUE)[,2] + 
                 area:poly(day_before_hatch,2, raw = TRUE)[,1] +  
                 
                 poly(day_before_hatch,2, raw = TRUE)[,2] +
                 poly(day_before_hatch,2, raw = TRUE)[,1] +
                 meantemp + 
                 #year + 
                 date_aprildays + 
                 #clutch_size + 
                 area + 
                 (1|box),
               REML = F,
               data=Asubset) 

## Use 'drop1' again and repeat the process removing the predictor with highest p-value:
drop1(model4_onset, test = "Chisq")


# model 5 - 'date_aprildays' removed
model5_onset <- lmer(diff_sunrise ~ 
                 #area:poly(day_before_hatch,2, raw = TRUE)[,2] + 
                 area:poly(day_before_hatch,2, raw = TRUE)[,1] +  
                 
                 poly(day_before_hatch,2, raw = TRUE)[,2] +
                 poly(day_before_hatch,2, raw = TRUE)[,1] +
                 meantemp + 
                 #year + 
                 #date_aprildays + 
                 #clutch_size + 
                 area + 
                 (1|box),
               REML = F,
               data=Asubset) 

## Use 'drop1' again and repeat the process removing the predictor with highest p-value:
drop1(model5_onset, test = "Chisq")


# model 6 - 'meantemp' removed
model6_onset <- lmer(diff_sunrise ~ 
                 #area:poly(day_before_hatch,2, raw = TRUE)[,2] + 
                 area:poly(day_before_hatch,2, raw = TRUE)[,1] +  
                 
                 poly(day_before_hatch,2, raw = TRUE)[,2] +
                 poly(day_before_hatch,2, raw = TRUE)[,1] +
                 #meantemp + 
                 #year + 
                 #date_aprildays + 
                 #clutch_size + 
                 area + 
                 (1|box),
               REML = F,
               data=Asubset) 

## Use 'drop1' again and repeat the process removing the predictor with highest p-value:
drop1(model6_onset, test = "Chisq")



# model 7 - 'poly(day_before_hatch, 2)[, 2]' removed
model7_onset <- lmer(diff_sunrise ~ 
                 #area:poly(day_before_hatch,2, raw = TRUE)[,2] + 
                 area:poly(day_before_hatch,2, raw = TRUE)[,1] +  
                 
                 #poly(day_before_hatch,2, raw = TRUE)[,2] +
                 poly(day_before_hatch,2, raw = TRUE)[,1] +
                 #meantemp + 
                 #year + 
                 #date_aprildays + 
                 #clutch_size + 
                 area + 
                 (1|box),
               REML = F,
               data=Asubset) 

## Use 'drop1' again and repeat the process removing the predictor with highest p-value:
drop1(model7_onset, test = "Chisq")


# model 8 - 'area:poly(day_before_hatch, 2)[, 1]' removed
model8_onset <- lmer(diff_sunrise ~ 
                 #area:poly(day_before_hatch,2, raw = TRUE)[,2] + 
                 #area:poly(day_before_hatch,2, raw = TRUE)[,1] +  
                 
                 #poly(day_before_hatch,2, raw = TRUE)[,2] +
                 poly(day_before_hatch,2, raw = TRUE)[,1] +
                 #meantemp + 
                 #year + 
                 #date_aprildays + 
                 #clutch_size + 
                 area + 
                 (1|box),
               REML = F,
               data=Asubset) 

## Use 'drop1' again and repeat the process removing the predictor with highest p-value:
drop1(model8_onset, test = "Chisq")

##
##
## 'model8_onset' is your final model because the removal of any of its predictors causes a significant change 
## in the power of the model (i.e. LR test for both predictors are siginificant as shown in the results of 'drop1')

##
## To calculate the significance of predictors not included in your final model, you added the predictor 
## of interest back to the final model and compare these two models via a LRT as before.

## You can do that for each preditor using 'anova':
anova(update(model8_onset, .~.+clutch_size, REML=F),   # final model + clutch size
      model8_onset,                            # final model
      test = "Chisq")
anova(update(model8_onset, .~.+date_aprildays, REML=F), # final model + date_aprildays
      model8_onset, 
      test = "Chisq")
anova(update(model8_onset, .~.+year, REML=F),           # final model + year
      model8_onset, 
      test = "Chisq")
anova(update(model8_onset, .~.+meantemp, REML=F),       # final model + meantemp
      model8_onset, 
      test = "Chisq")
anova(update(model8_onset, .~.+area:poly(day_before_hatch,2, raw = TRUE)[,1], REML=F), # final model + area:poly(day_before_hatch,2)[,1])
      model8_onset, 
      test = "Chisq")
anova(update(model8_onset, .~.+poly(day_before_hatch,2, raw = TRUE)[,2], REML=F),     # final model + poly(day_before_hatch,2)[,2])
      model8_onset, 
      test = "Chisq")



##
## To report model estimates, you can fit the model using REML=T
## Model 8 is the MAM
model8_onset_final <- update(model8_onset, .~.,REML=TRUE) # REML final estimates
summary(model8_onset_final)

##
##
##### Plotting model predictions - onset #####
##
##
summary(Asubset$day_before_hatch)

# new data table to generate predictions
data_predictions <- expand.grid(day_before_hatch = -17:0,
                                area = c("urban", "forest"))

# mean predictions
data_predictions$fit <- predict(model8_onset_final, newdata = data_predictions, re.form = NA)

# SE for mean predicitons
mm <- model.matrix(~ day_before_hatch + 
                     area,
                   data = data_predictions)

## or newdat$distance <- mm %*% fixef(fm1)
pvar1 <- diag(mm %*% tcrossprod(vcov(model8_onset_final),mm))
cmult <- 1 ## 1 SE
data_predictions <- data.frame(
  data_predictions
  , plo = data_predictions$fit-cmult*sqrt(pvar1)
  , phi = data_predictions$fit+cmult*sqrt(pvar1)
)

## plot
onsetplot <- ggplot(data = Asubset, 
                    aes(x=day_before_hatch, 
                        y=diff_sunrise, 
                        fill=area))+
  geom_point(position = position_dodge(width = -0.45), 
             alpha = 0.60, 
             color = "black", 
             shape = 21, 
             size = 0.75) +
  geom_errorbar(data = data_predictions,
                aes(x=day_before_hatch, 
                    y=fit, 
                    ymin = plo, ymax = phi),
                position = position_dodge(width = 1),
                width = 0,
                size = 1) +
  geom_point(data = data_predictions,
             aes(x=day_before_hatch, 
                 y=fit, 
                 fill=area),
             position = position_dodge(width = 1), 
             size = 2.5,
             color = "black", 
             shape = 21) +
  theme_bw() +
  labs(y = "Onset of activity (mins from sunrise)", x = "Days before hatching") +
  scale_fill_manual(name = "", labels = c("Forest", "City"), 
                    values = c("#339900","#3399CC")) 

## warning about position_dodge nothing to worry about
  

onset_plot_final <- onsetplot + 
  geom_hline(aes(yintercept=0), colour="black", linetype="dashed") +
  geom_text(aes( -2, -25, 
                 label = "Morning civil twilight", 
                 vjust = -1), 
            size = 4) + 
  theme(legend.position="top")


ggsave(filename = "./plots/onset_ROBYN.jpeg",
       plot = onset_plot_final, 
       device = "jpeg", 
       units = "mm", 
       width = 150, 
       height = 100)

##
##
##### models activity offset #####
##
##
# removing very large outliers
Asubset_offset <- Asubset %>% 
  filter(diff_sunset > -300)

## exactly the same approach as above can be used for the analysis of end of activity

# global model
model1_offset <- lmer(diff_sunset ~ 
                        area:poly(day_before_hatch,2, raw=T)[,2] + 
                        area:poly(day_before_hatch,2, raw=T)[,1] +  
                        
                        poly(day_before_hatch,2, raw=T)[,2] +
                        poly(day_before_hatch,2, raw=T)[,1] +
                        meantemp + 
                        year + 
                        date_aprildays + 
                        clutch_size + 
                        area + 
                        (1|box), 
                      REML = F,
                      na.action = "na.fail",
                      data=Asubset_offset) #full model
summary(model1_offset)
drop1(model1_offset, test = "Chisq")

# model 2 - 'meantemp' removed
model2_offset <- lmer(diff_sunset ~ 
                        area:poly(day_before_hatch,2, raw=T)[,2] + 
                        area:poly(day_before_hatch,2, raw=T)[,1] +  
                        
                        poly(day_before_hatch,2, raw=T)[,2] +
                        poly(day_before_hatch,2, raw=T)[,1] +
                        #meantemp + 
                        year + 
                        date_aprildays + 
                        clutch_size + 
                        area + 
                        (1|box), 
                      REML = F,
                      data=Asubset_offset) #full model
summary(model2_offset)
drop1(model2_offset, test = "Chisq")

# model 3 - 'year' removed
model3_offset <- lmer(diff_sunset ~ 
                        area:poly(day_before_hatch,2, raw=T)[,2] + 
                        area:poly(day_before_hatch,2, raw=T)[,1] +  
                        
                        poly(day_before_hatch,2, raw=T)[,2] +
                        poly(day_before_hatch,2, raw=T)[,1] +
                        #meantemp + 
                        #year + 
                        date_aprildays + 
                        clutch_size + 
                        area + 
                        (1|box), 
                      REML = F,
                      data=Asubset_offset) #full model
summary(model3_offset)
drop1(model3_offset, test = "Chisq")

# model 4 - 'area:poly(day_before_hatch,2)[,2]' removed
model4_offset <- lmer(diff_sunset ~ 
                        #area:poly(day_before_hatch,2, raw=T)[,2] + 
                        area:poly(day_before_hatch,2, raw=T)[,1] +  
                        
                        poly(day_before_hatch,2, raw=T)[,2] +
                        poly(day_before_hatch,2, raw=T)[,1] +
                        #meantemp + 
                        #year + 
                        date_aprildays + 
                        clutch_size + 
                        area + 
                        (1|box), 
                      REML = F,
                      data=Asubset_offset) #full model
summary(model4_offset)
drop1(model4_offset, test = "Chisq")

# model 5 - 'poly(day_before_hatch, 2)[, 2]' removed
model5_offset <- lmer(diff_sunset ~ 
                        #area:poly(day_before_hatch,2, raw=T)[,2] + 
                        area:poly(day_before_hatch,2, raw=T)[,1] +  
                        
                        #poly(day_before_hatch,2, raw=T)[,2] +
                        poly(day_before_hatch,2, raw=T)[,1] +
                        #meantemp + 
                        #year + 
                        date_aprildays + 
                        clutch_size + 
                        area + 
                        (1|box), 
                      REML = F,
                      data=Asubset_offset) #full model
summary(model5_offset)
drop1(model5_offset, test = "Chisq")

# model 6 - 'clutch_size' removed
model6_offset <- lmer(diff_sunset ~ 
                        #area:poly(day_before_hatch,2, raw=T)[,2] + 
                        area:poly(day_before_hatch,2, raw=T)[,1] +  
                        
                        #poly(day_before_hatch,2, raw=T)[,2] +
                        poly(day_before_hatch,2, raw=T)[,1] +
                        #meantemp + 
                        #year + 
                        date_aprildays + 
                        #clutch_size + 
                        area + 
                        (1|box), 
                      REML = F,
                      data=Asubset_offset) #full model
summary(model6_offset)
drop1(model6_offset, test = "Chisq")



## final model offset of activity
model6_offset_final <- lmer(diff_sunset ~ 
                              area:day_before_hatch +  
                              day_before_hatch +
                              date_aprildays + 
                              area + 
                              (1|box), 
                            REML = T,
                            data=Asubset_offset) #full model
summary(model6_offset_final)
drop1(model6_offset_final, test="Chisq")


## testing the significance level of predictors not included in MAM
anova(update(model6_offset, .~.+clutch_size), model6_offset, test = "Chisq")
anova(update(model6_offset, .~.+year), model6_offset, test = "Chisq")
anova(update(model6_offset, .~.+meantemp), model6_offset, test = "Chisq")
anova(update(model6_offset, .~.+poly(day_before_hatch,2)[,2]), model6_offset, test = "Chisq")


##
##
##### Plotting model predictions #####
##
##
summary(Asubset$day_before_hatch)

# new data table to generate predictions
data_predictions <- expand.grid(day_before_hatch = -17:0,
                                date_aprildays = mean(Asubset_offset$date_aprildays),
                                area = c("City", "Forest"))

# mean predictions
data_predictions$fit <- predict(model6_offset_final, newdata = data_predictions, re.form = NA)

# SE for mean predicitons
mm <- model.matrix(~ area:day_before_hatch +  
                     day_before_hatch +
                     date_aprildays + 
                     area,
                   data = data_predictions)

## or newdat$distance <- mm %*% fixef(fm1)
pvar1 <- diag(mm %*% tcrossprod(vcov(model6_offset_final),mm))
cmult <- 1 ## could use 1.96
data_predictions <- data.frame(
  data_predictions
  , plo = data_predictions$fit-cmult*sqrt(pvar1)
  , phi = data_predictions$fit+cmult*sqrt(pvar1)
)





offsetplot <- ggplot(data = Asubset_offset, 
                    aes(x=day_before_hatch, 
                        y=diff_sunset, 
                        fill=area))+
  geom_point(position = position_dodge(width = -0.45), 
             alpha = 0.60, 
             color = "black", 
             shape = 21,
             size = 0.75) +
  geom_errorbar(data = data_predictions,
                aes(x=day_before_hatch, 
                    y=fit, 
                    ymin = plo, ymax = phi),
                position = position_dodge(width = 1),
                width = 0,
                size = 1) +
  geom_point(data = data_predictions,
             aes(x=day_before_hatch, 
                 y=fit, 
                 fill=area),
             position = position_dodge(width = 1), 
             size = 2.5,
             color = "black", 
             shape = 21) +
  theme_bw() +
  labs(y = "Offset of activity (mins from sunrise)", x = "Days before hatching") +
  scale_fill_manual(name = "", labels = c("Forest", "City"), 
                    values = c("#339900","#3399CC")) 

offsetplot_final <- offsetplot + 
  geom_hline(aes(yintercept=0), 
             colour="black", 
             linetype="dashed") +
geom_text(aes( -13, 0, 
               label = "Evening civil twilight", vjust = -1), 
          size = 4) + 
  theme(legend.position="top")



ggsave(filename = "./plots/offset_ROBYN.jpeg",
       plot = offsetplot_final, 
       device = "jpeg", 
       units = "mm", 
       width = 150, 
       height = 100)


##
##
##### models percentage of time in #####
##
##

hist(Asubset$percentage_in)
# global model
model1_act <- lmer(percentage_in ~ 
                       area:poly(day_before_hatch,2, raw=T)[,2] + 
                       area:poly(day_before_hatch,2, raw=T)[,1] +  
                       
                       poly(day_before_hatch,2, raw=T)[,2] +
                       poly(day_before_hatch,2, raw=T)[,1] +
                       meantemp + 
                       year + 
                       date_aprildays + 
                       clutch_size + 
                       area + 
                       (1|box), 
                     REML = F,
                     na.action = "na.fail",
                     data=Asubset_offset) #full model

summary(model1_act)
drop1(model1_act, test = "Chisq")

# model 2 - "area:poly(day_before_hatch, 2)[, 1]" removed
model2_act <- lmer(percentage_in ~ 
                     area:poly(day_before_hatch,2, raw=T)[,2] + 
                     #area:poly(day_before_hatch,2, raw=T)[,1] +  
                     
                     poly(day_before_hatch,2, raw=T)[,2] +
                     poly(day_before_hatch,2, raw=T)[,1] +
                     meantemp + 
                     year + 
                     date_aprildays + 
                     clutch_size + 
                     area + 
                     (1|box), 
                   REML = F,
                   na.action = "na.fail",
                   data=Asubset_offset) 
summary(model2_act)
drop1(model2_act, test = "Chisq")

# model 3 - "area:poly(day_before_hatch, 2)[, 2]" removed
model3_act <- lmer(percentage_in ~ 
                     #area:poly(day_before_hatch,2)[,2] + 
                     #area:poly(day_before_hatch,2)[,1] +  
                     
                     poly(day_before_hatch,2, raw=T)[,2] +
                     poly(day_before_hatch,2, raw=T)[,1] +
                     meantemp + 
                     year + 
                     date_aprildays + 
                     clutch_size + 
                     area + 
                     (1|box), 
                   REML = F,
                   na.action = "na.fail",
                   data=Asubset_offset) 
summary(model3_act)
drop1(model3_act, test = "Chisq")

# model 4 - "meantemp" removed
model4_act <- lmer(percentage_in ~ 
                     #area:poly(day_before_hatch,2, raw=T)[,2] + 
                     #area:poly(day_before_hatch,2, raw=T)[,1] +  
                     
                     poly(day_before_hatch,2, raw=T)[,2] +
                     poly(day_before_hatch,2, raw=T)[,1] +
                     #meantemp + 
                     year + 
                     date_aprildays + 
                     clutch_size + 
                     area + 
                     (1|box), 
                   REML = F,
                   na.action = "na.fail",
                   data=Asubset_offset) 
summary(model4_act)
drop1(model4_act, test = "Chisq") # removing 'area' causes a model convergency problem


## final model percentage of activity
model4_act_final <- lmer(percentage_in ~ 
                           poly(day_before_hatch,2, raw=T)[,2] +
                           poly(day_before_hatch,2, raw=T)[,1] +
                           year + 
                           date_aprildays + 
                           clutch_size + 
                           area + 
                           (1|box), 
                         REML = F,
                         na.action = "na.fail",
                         data=Asubset_offset) 
summary(model4_act_final)
drop1(model4_act_final, test="Chisq")


## testing the significance level of predictors not included in MAM
anova(update(model4_act, .~.+meantemp), model4_act, test = "Chisq")
anova(update(model4_act, .~.+area:poly(day_before_hatch,2, raw=T)[,1]), model4_act, test = "Chisq")
anova(update(model4_act, .~.+area:poly(day_before_hatch,2, raw=T)[,2]), model4_act, test = "Chisq")

anova(update(model4_act, .~.+area:year), model4_act, test = "Chisq")


##
##
##### Plotting model predictions percentage of activity #####
##
##

# new data table to generate predictions
data_predictions <- expand.grid(day_before_hatch = -17:0,
                                year = c("2016", "2017", "2018"),
                                clutch_size = mean(Asubset_offset$clutch_size),
                                date_aprildays = mean(Asubset_offset$date_aprildays),
                                area = c("urban", "forest"))

# mean predictions
data_predictions$fit <- predict(model4_act_final, newdata = data_predictions, re.form = NA)

# SE for mean predicitons
mm <- model.matrix(~ I(day_before_hatch^2) +
                     day_before_hatch + 
                     year + 
                     date_aprildays + 
                     clutch_size + 
                     area,
                   data = data_predictions)

## or newdat$distance <- mm %*% fixef(fm1)
pvar1 <- diag(mm %*% tcrossprod(vcov(model4_act_final),mm))
cmult <- 1 ## could use 1.96
data_predictions <- data.frame(
  data_predictions
  , plo = data_predictions$fit-cmult*sqrt(pvar1)
  , phi = data_predictions$fit+cmult*sqrt(pvar1)
)



percinplot_final <- ggplot(data = Asubset_offset, 
                     aes(x=day_before_hatch, 
                         y=percentage_in, 
                         fill=area))+
  geom_point(position = position_dodge(width = -0.45), 
             alpha = 0.60, 
             color = "black", 
             shape = 21,
             size = 0.75) +
  geom_errorbar(data = data_predictions,
                aes(x=day_before_hatch, 
                    y=fit, 
                    ymin = plo, ymax = phi),
                position = position_dodge(width = 1),
                width = 0,
                size = 1) +
  geom_point(data = data_predictions,
             aes(x=day_before_hatch, 
                 y=fit, 
                 fill=area),
             position = position_dodge(width = 1), 
             size = 2.5,
             color = "black", 
             shape = 21) +
  facet_grid(.~year) +
  theme_bw() +
  theme(legend.position="top") +
  labs(y = "Percentage time in nest box (%)", x = "Days before hatching") +
  scale_fill_manual(name = "", labels = c("Forest", "City"), 
                    values = c("#339900","#3399CC")) 




ggsave(filename = "./plots/perc_act_ROBYN.jpeg",
       plot = percinplot_final, 
       device = "jpeg", 
       units = "mm", 
       width = 150, 
       height = 100)




##
##
##### models duration bird's day #####
##
##
 
Asubset_duration <- Asubset %>% 
  filter(birds_day > 650)

hist(Asubset_duration$birds_day)

# global model
model1_duration <- lmer(birds_day ~ 
                     area:poly(day_before_hatch,2, raw=T)[,2] + 
                     area:poly(day_before_hatch,2, raw=T)[,1] +  
                     
                     poly(day_before_hatch,2, raw=T)[,2] +
                     poly(day_before_hatch,2, raw=T)[,1] +
                     meantemp + 
                     year + 
                     date_aprildays + 
                     clutch_size + 
                     area + 
                     (1|box), 
                   REML = F,
                   na.action = "na.fail",
                   data=Asubset_duration) #full model

summary(model1_duration)
drop1(model1_duration, test = "Chisq")

# model 2 - 'area:poly(day_before_hatch,2, raw=T)[,1]' removed
model2_duration <- lmer(birds_day ~ 
                          area:poly(day_before_hatch,2, raw=T)[,2] + 
                          #area:poly(day_before_hatch,2, raw=T)[,1] +  
                          
                          poly(day_before_hatch,2, raw=T)[,2] +
                          poly(day_before_hatch,2, raw=T)[,1] +
                          meantemp + 
                          year + 
                          date_aprildays + 
                          clutch_size + 
                          area + 
                          (1|box), 
                        REML = F,
                        na.action = "na.fail",
                        data=Asubset_duration) 

summary(model2_duration)
drop1(model2_duration, test = "Chisq")

# model 3 - 'date_aprildays' removed
model3_duration <- lmer(birds_day ~ 
                          area:poly(day_before_hatch,2, raw=T)[,2] + 
                          #area:poly(day_before_hatch,2, raw=T)[,1] +  
                          
                          poly(day_before_hatch,2, raw=T)[,2] +
                          poly(day_before_hatch,2, raw=T)[,1] +
                          meantemp + 
                          year + 
                          #date_aprildays + 
                          clutch_size + 
                          area + 
                          (1|box), 
                        REML = F,
                        na.action = "na.fail",
                        data=Asubset_duration) 

summary(model3_duration)
drop1(model3_duration, test = "Chisq")

# model 4 - 'year' removed
model4_duration <- lmer(birds_day ~ 
                          area:poly(day_before_hatch,2, raw=T)[,2] + 
                          #area:poly(day_before_hatch,2, raw=T)[,1] +  
                          
                          poly(day_before_hatch,2, raw=T)[,2] +
                          poly(day_before_hatch,2, raw=T)[,1] +
                          meantemp + 
                          #year + 
                          #date_aprildays + 
                          clutch_size + 
                          area + 
                          (1|box), 
                        REML = F,
                        na.action = "na.fail",
                        data=Asubset_duration) 

summary(model4_duration)
drop1(model4_duration, test = "Chisq")

# model 5 - 'meantemp' removed
model5_duration <- lmer(birds_day ~ 
                          area:poly(day_before_hatch,2, raw=T)[,2] + 
                          #area:poly(day_before_hatch,2, raw=T)[,1] +  
                          
                          poly(day_before_hatch,2, raw=T)[,2] +
                          poly(day_before_hatch,2, raw=T)[,1] +
                          #meantemp + 
                          #year + 
                          #date_aprildays + 
                          clutch_size + 
                          area + 
                          (1|box), 
                        REML = F,
                        na.action = "na.fail",
                        data=Asubset_duration) 

summary(model5_duration)
drop1(model5_duration, test = "Chisq")

# model 6 - 'clutch_size' removed
model6_duration <- lmer(birds_day ~ 
                          area:poly(day_before_hatch,2, raw=T)[,2] + 
                          #area:poly(day_before_hatch,2, raw=T)[,1] +  
                          
                          poly(day_before_hatch,2, raw=T)[,2] +
                          poly(day_before_hatch,2, raw=T)[,1] +
                          #meantemp + 
                          #year + 
                          #date_aprildays + 
                          #clutch_size + 
                          area + 
                          (1|box), 
                        REML = F,
                        na.action = "na.fail",
                        data=Asubset_duration) 

summary(model6_duration)
drop1(model6_duration, test = "Chisq")


## final mode
model6_duration_final <- update(model6_duration, REML=T)

## testing the significance level of predictors not included in MAM
anova(update(model6_duration, .~.+clutch_size), model6_duration, test = "Chisq")
anova(update(model6_duration, .~.+date_aprildays), model6_duration, test = "Chisq")
anova(update(model6_duration, .~.+year), model6_duration, test = "Chisq")
anova(update(model6_duration, .~.+meantemp), model6_duration, test = "Chisq")
anova(update(model6_duration, .~.+area:poly(day_before_hatch,2, raw=T)[,1]), model6_duration, test = "Chisq")



##
##
##### Plotting model predictions percentage of activity #####
##
##

# new data table to generate predictions
data_predictions <- expand.grid(day_before_hatch = -17:0,
                                area = c("urban", "forest"))

# mean predictions
data_predictions$fit <- predict(model6_duration_final, newdata = data_predictions, re.form = NA)

# SE for mean predicitons
mm <- model.matrix(~ area:poly(day_before_hatch,2, raw=T)[,2] + 
                     poly(day_before_hatch,2, raw=T)[,2] +
                     poly(day_before_hatch,2, raw=T)[,1] +
                     area,
                   data = data_predictions)

## or newdat$distance <- mm %*% fixef(fm1)
pvar1 <- diag(mm %*% tcrossprod(vcov(model6_duration_final),mm))
cmult <- 1 ## could use 1.96
data_predictions <- data.frame(
  data_predictions
  , plo = data_predictions$fit-cmult*sqrt(pvar1)
  , phi = data_predictions$fit+cmult*sqrt(pvar1)
)


birdsdayplot_final <- ggplot(data = Asubset_duration, 
                       aes(x=day_before_hatch, 
                           y=birds_day, 
                           fill=area))+
  geom_point(position = position_dodge(width = -0.45), 
             alpha = 0.60, 
             color = "black", 
             shape = 21,
             size = 0.75) +
  geom_errorbar(data = data_predictions,
                aes(x=day_before_hatch, 
                    y=fit, 
                    ymin = plo, ymax = phi),
                position = position_dodge(width = 1),
                width = 0,
                size = 1) +
  geom_point(data = data_predictions,
             aes(x=day_before_hatch, 
                 y=fit, 
                 fill=area),
             position = position_dodge(width = 1), 
             size = 2.5,
             color = "black", 
             shape = 21) +
  theme_bw() +
  theme(legend.position="top") +
  theme(legend.position="top") +
  labs(y = "Duration of active day (in minutes)", x = "Days before hatching") +
  scale_fill_manual(name = "", labels = c("Forest", "City"), 
                    values = c("#339900","#3399CC")) 



ggsave(filename = "./plots/day_duration_ROBYN.jpeg",
       plot = birdsdayplot_final, 
       device = "jpeg", 
       units = "mm", 
       width = 150, 
       height = 100)





