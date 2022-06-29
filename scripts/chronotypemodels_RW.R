library(ggplot2)
library(lme4)
library(lmtest)

#### activity onset ###
buttondata <- read.csv (file="C:\\Users\\Robyn\\Dropbox\\iButton Chronotype Chapter\\incubationdata.csv", na.strings=c("", "NA"))

gsubset <- subset(buttondata, select=c(box, date_aprildays, year, area, clutch_size, day_before_hatch, diff_sunrise, meantemp))
#remove missing rows, create index
Aindex <- complete.cases(gsubset)
#  subset[row,column]
Asubset <- gsubset[Aindex,]
Asubset$year <- as.factor(Asubset$year) 
Asubset$diff_sunrise <- as.numeric(Asubset$diff_sunrise)
Asubset$day_before_hatch <- as.numeric(Asubset$day_before_hatch)
str(Asubset)

model1 <- lmer(diff_sunrise~(1|box)+
                 area + 
                 area:poly(day_before_hatch,2) + 
                 area:poly(day_before_hatch,1) +  
                 day_before_hatch+ 
                 poly(day_before_hatch,2) + 
                 clutch_size+
                 meantemp+
                 date_aprildays+
                 year, 
               data=Asubset) #full model
summary(model1)
drop1(model1, test = "Chisq")

model2 <- lmer(diff_sunrise~(1|box)+
                 area + 
                 area:poly(day_before_hatch,1) +  
                 day_before_hatch+ 
                 poly(day_before_hatch,2) + 
                 clutch_size+
                 meantemp+
                 date_aprildays+
                 year, 
               data=Asubset) #full model

lrtest(model1,model2) #keep area:poly(day_before_hatch,2)


model3 <- lmer(diff_sunrise~(1|box)+area + area:day_before_hatch^2 +  day_before_hatch+ poly(day_before_hatch,2) + clutch_size+meantemp+date_aprildays+year, data=Asubset) #full model
lrtest(model1,model3) # keep linear interaction
model4 <- lmer(diff_sunrise~(1|box)+area + area:poly(day_before_hatch,2) + area:poly(day_before_hatch,1) +  day_before_hatch+ poly(day_before_hatch,2) + meantemp+date_aprildays+year, data=Asubset) #full model
lrtest(model1,model4) #drop clutch size
model5 <- lmer(diff_sunrise~(1|box)+area + area:poly(day_before_hatch,2) + area:poly(day_before_hatch,1) +  day_before_hatch+ poly(day_before_hatch,2) + date_aprildays+year, data=Asubset) #full model
lrtest(model4,model5)  # drop mean temp
model6 <- lmer(diff_sunrise~(1|box)+area + area:poly(day_before_hatch,2) + area:poly(day_before_hatch,1) +  day_before_hatch+ poly(day_before_hatch,2) +year, data=Asubset) #full model
lrtest(model5,model6) 
model7 <- lmer(diff_sunrise~(1|box)+area + area:poly(day_before_hatch,2) + area:poly(day_before_hatch,1) +  day_before_hatch+ poly(day_before_hatch,2), data=Asubset) #full model
lrtest(model6,model7) #keep year
summary(model6)
write.csv( summary(model6)$coefficients, "C:\\Users\\Robyn\\Dropbox\\iButton Chronotype Chapter\\activityonset_take2.csv" )


#predictions#
predict <- predict(model6, newdata = Asubset)
summary(predict)
write.csv(predict, "C:\\Users\\Robyn\\Dropbox\\iButton Chronotype Chapter\\onsetpredictions.csv" )
write.csv(Asubset, "C:\\Users\\Robyn\\Dropbox\\iButton Chronotype Chapter\\onsetsubset.csv" )

#plotting onset
buttondata <- read.csv (file="C:\\Users\\Robyn\\Dropbox\\iButton Chronotype Chapter\\onsetsubset.csv", na.strings=c("", "NA"))

gsubset <- subset(buttondata, select=c(box, year, area, day_before_hatch, diff_sunrise, predictmeans))
#remove missing rows, create index
Aindex <- complete.cases(gsubset)
#  subset[row,column]
Asubset <- gsubset[Aindex,]
Asubset$year <- as.factor(Asubset$year) 
Asubset$diff_sunrise <- as.numeric(Asubset$diff_sunrise)
Asubset$predictmeans <- as.numeric(Asubset$predictmean)
str(Asubset)

theme_set(theme_bw(base_size = 15))
onsetplot <- ggplot(data = Asubset, 
                    aes(x=day_before_hatch, 
                        y=diff_sunrise, 
                        fill=as.factor(area)))+
  geom_point(position = position_dodge(width = 0.45), 
             alpha = 0.60, 
             color = "black", 
             shape = 21) +
  theme_bw() +
  labs(y = "Onset of activity (mins from sunrise)", x = "Days before hatching") +
  scale_fill_manual(name = "", labels = c("City", "Forest"), 
                    values = c("#3399CC", "#339900")) +
  geom_point(aes(y=predictmeans, 
                 colour=as.factor(area)),
             position = position_dodge(width = 0.85), 
             size = 3.7,
             color = "black", 
             shape = 21)

onsetplot + geom_hline(aes(yintercept=0), colour="black", linetype="dashed") +
  geom_text(aes( -2, -25, label = "Morning civil twilight", vjust = -1), size = 4) + theme(legend.position="top")



#### activity offset ###




buttondata <- read.csv (file="C:\\Users\\Robyn\\Dropbox\\iButton Chronotype Chapter\\incubationdata.csv", na.strings=c("", "NA"))

gsubset <- subset(buttondata, select=c(box, date_aprildays, year, area, clutch_size, day_before_hatch, diff_sunset, meantemp))
#remove missing rows, create index
Aindex <- complete.cases(gsubset)
#  subset[row,column]
Asubset <- gsubset[Aindex,]
Asubset$year <- as.factor(Asubset$year) 
Asubset$diff_sunset <- as.numeric(Asubset$diff_sunset)
Asubset$day_before_hatch <- as.numeric(Asubset$day_before_hatch)
str(Asubset)

model1 <- lmer(diff_sunset~(1|box)+area + area:poly(day_before_hatch,2) + area:poly(day_before_hatch,1) +  day_before_hatch+ poly(day_before_hatch,2) + clutch_size+meantemp+date_aprildays+year, data=Asubset) #full model
summary(model1)
drop1(model1, test = "Chisq")
model2 <- lmer(diff_sunset~(1|box)+area + area:poly(day_before_hatch,1) +  day_before_hatch+ poly(day_before_hatch,2) + clutch_size+meantemp+date_aprildays+year, data=Asubset) #full model
lrtest(model1,model2) #keep interaction
model3 <- lmer(diff_sunset~(1|box)+area + area:day_before_hatch^2 +  day_before_hatch+ poly(day_before_hatch,2) + clutch_size+meantemp+date_aprildays+year, data=Asubset) #full model
lrtest(model1,model3) # keep interaction
model4 <- lmer(diff_sunset~(1|box)+area + area:poly(day_before_hatch,2) + area:poly(day_before_hatch,1) +  day_before_hatch+ poly(day_before_hatch,2)+meantemp+date_aprildays+year, data=Asubset) #full model
lrtest(model1,model4) #keep clutch size
model5 <- lmer(diff_sunset~(1|box)+area + area:poly(day_before_hatch,2) + area:poly(day_before_hatch,1) +  day_before_hatch+ poly(day_before_hatch,2)+date_aprildays+year, data=Asubset) #full model
lrtest(model1,model5) #drop mean temp
model6 <- lmer(diff_sunset~(1|box)+area + area:poly(day_before_hatch,2) + area:poly(day_before_hatch,1) +  day_before_hatch+ poly(day_before_hatch,2)+year, data=Asubset) #full model
lrtest(model5,model6)
model7 <- lmer(diff_sunset~(1|box)+area + area:poly(day_before_hatch,2) + area:poly(day_before_hatch,1) +  day_before_hatch+ poly(day_before_hatch,2)+date_aprildays, data=Asubset) #full model
lrtest(model5,model7) #keep year
summary(model5)
write.csv( summary(model5)$coefficients, "C:\\Users\\Robyn\\Dropbox\\iButton Chronotype Chapter\\activityoffset_take2.csv" )

#predictions#
predict <- predict(model5, newdata = Asubset)
summary(predict)
write.csv(predict, "C:\\Users\\Robyn\\Dropbox\\iButton Chronotype Chapter\\offsetpredictions.csv" )
write.csv(Asubset, "C:\\Users\\Robyn\\Dropbox\\iButton Chronotype Chapter\\offsetsubset.csv" )


#plotting offset
buttondata <- read.csv (file="C:\\Users\\Robyn\\Dropbox\\iButton Chronotype Chapter\\offsetsubset.csv", na.strings=c("", "NA"))

gsubset <- subset(buttondata, select=c(box, year, area, day_before_hatch, diff_sunset, predictmeans))
#remove missing rows, create index
Aindex <- complete.cases(gsubset)
#  subset[row,column]
Asubset <- gsubset[Aindex,]
Asubset$year <- as.factor(Asubset$year) 
Asubset$diff_sunset <- as.numeric(Asubset$diff_sunset)
Asubset$predictmeans <- as.numeric(Asubset$predictmeans)
str(Asubset)

theme_set(theme_bw(base_size = 15))
offsetplot <- ggplot(data = Asubset, 
                    aes(x=day_before_hatch, 
                        y=diff_sunset, 
                        fill=as.factor(area)))+
  geom_point(position = position_dodge(width = 0.45), 
             alpha = 0.60, 
             color = "black", 
             shape = 21) +
  theme_bw() +
  labs(y = "Offset of activity (mins from sunset)", x = "Days before hatching") +
  scale_fill_manual(name = "", labels = c("City", "Forest"), 
                    values = c("#3399CC", "#339900")) +
  geom_point(aes(y=predictmeans, 
                 colour=as.factor(area)),
             position = position_dodge(width = 0.85), 
             size = 3.7,
             color = "black", 
             shape = 21)
  

offsetplot + geom_hline(aes(yintercept=0), colour="black", linetype="dashed") +
geom_text(aes( -13, 0, label = "Evening civil twilight", vjust = -1), size = 4) + theme(legend.position="top")



#percentage time in box #

buttondata <- read.csv (file="C:\\Users\\Robyn\\Dropbox\\iButton Chronotype Chapter\\incubationdata.csv", na.strings=c("", "NA"))

gsubset <- subset(buttondata, select=c(box, date_aprildays, year, area, clutch_size, day_before_hatch, percentage_in, meantemp))
#remove missing rows, create index
Aindex <- complete.cases(gsubset)
#  subset[row,column]
Asubset <- gsubset[Aindex,]
Asubset$year <- as.factor(Asubset$year) 
Asubset$percentage_in <- as.integer(Asubset$percentage_in)
Asubset$day_before_hatch <- as.numeric(Asubset$day_before_hatch)
str(Asubset)
plot(Asubset$percentage_in~Asubset$day_before_hatch)

model1 <- lmer(percentage_in~(1|box)+area + area:poly(day_before_hatch,2) + area:poly(day_before_hatch,1) +  day_before_hatch+ poly(day_before_hatch,2) + clutch_size+meantemp+date_aprildays+year, data=Asubset) #full model
summary(model1)
drop1(model1, test = "Chisq")
model2 <- lmer(percentage_in~(1|box)+area +  area:poly(day_before_hatch,1) +  day_before_hatch+ poly(day_before_hatch,2) + clutch_size+meantemp+date_aprildays+year, data=Asubset) #full model
lrtest(model1,model2)
model3 <- lmer(percentage_in~(1|box)+area + area:day_before_hatch^2 +  day_before_hatch+ poly(day_before_hatch,2) + clutch_size+meantemp+date_aprildays+year, data=Asubset) #full model
lrtest(model1,model3) #keep linear
summary(model2)
model4 <- lmer(percentage_in~(1|box)+area + area:poly(day_before_hatch,2) + area:poly(day_before_hatch,1) +  day_before_hatch+ poly(day_before_hatch,2)+meantemp+date_aprildays+year, data=Asubset) #full model
lrtest(model1,model4) #drop clutch size
model5 <- lmer(percentage_in~(1|box)+area + area:poly(day_before_hatch,2) + area:poly(day_before_hatch,1) +  day_before_hatch+ poly(day_before_hatch,2)+date_aprildays+year, data=Asubset) #full model
lrtest(model4,model5) #drop mean temp
model6 <- lmer(percentage_in~(1|box)+area + area:poly(day_before_hatch,2) + area:poly(day_before_hatch,1) +  day_before_hatch+ poly(day_before_hatch,2)+year, data=Asubset) #full model
lrtest(model5,model6) #drop date
model7 <- lmer(percentage_in~(1|box)+area + area:poly(day_before_hatch,2) + area:poly(day_before_hatch,1) +  day_before_hatch+ poly(day_before_hatch,2), data=Asubset) #full model
lrtest(model6,model7) #keep year
summary(model6)
write.csv( summary(model6)$coefficients, "C:\\Users\\Robyn\\Dropbox\\iButton Chronotype Chapter\\percentagein_take2.csv" )


predict <- predict(model6, newdata = Asubset)
summary(predict)
write.csv(predict, "C:\\Users\\Robyn\\Dropbox\\iButton Chronotype Chapter\\percentageinpredictions.csv" )
write.csv(Asubset, "C:\\Users\\Robyn\\Dropbox\\iButton Chronotype Chapter\\percentageinsubset.csv" )


#plotting percentage in#
buttondata <- read.csv (file="C:\\Users\\Robyn\\Dropbox\\iButton Chronotype Chapter\\percentageinsubset.csv", na.strings=c("", "NA"))

gsubset <- subset(buttondata, select=c(box, year, area, day_before_hatch, percentage_in, predictmeans))
#remove missing rows, create index
Aindex <- complete.cases(gsubset)
#  subset[row,column]
Asubset <- gsubset[Aindex,]
Asubset$year <- as.factor(Asubset$year) 
Asubset$percentage_in <- as.numeric(Asubset$percentage_in)
Asubset$predictmeans <- as.numeric(Asubset$predictmeans)
str(Asubset)

theme_set(theme_bw(base_size = 15))
percinplot <- ggplot(data = Asubset, 
                     aes(x=day_before_hatch, 
                         y=percentage_in, 
                         fill=as.factor(area)))+
  geom_point(position = position_dodge(width = 0.45), 
             alpha = 0.60, 
             color = "black", 
             shape = 21) +
  theme_bw() +
  labs(y = "Percentage time in nest box (%)", x = "Days before hatching") +
  scale_fill_manual(name = "", labels = c("City", "Forest"), 
                    values = c("#3399CC", "#339900")) +
  geom_point(aes(y=predictmeans, 
                 colour=as.factor(area)),
             position = position_dodge(width = 0.85), 
             size = 3.7,
             color = "black", 
             shape = 21)


percinplot + theme(legend.position="top")




## night time restlessness ##
buttondata <- read.csv (file="C:\\Users\\Robyn\\Dropbox\\iButton Chronotype Chapter\\incubationdata.csv", na.strings=c("", "NA"))

gsubset <- subset(buttondata, select=c(box, date_aprildays, year, area, clutch_size, day_before_hatch, night_var, light_atnight, meantemp))
#remove missing rows, create index
Aindex <- complete.cases(gsubset)
#  subset[row,column]
Asubset <- gsubset[Aindex,]
Asubset$year <- as.factor(Asubset$year) 
Asubset$night_var <- as.numeric(Asubset$night_var) 

model1 <- lmer(night_var~(1|box)+area + area:poly(day_before_hatch,2) + area:poly(day_before_hatch,1) +  day_before_hatch+ poly(day_before_hatch,2) + clutch_size+meantemp+date_aprildays+year, data=Asubset) #full model
summary(model1)
drop1(model1, test = "Chisq")
model2 <- lmer(night_var~(1|box)+area +  area:poly(day_before_hatch,1) +  day_before_hatch+ poly(day_before_hatch,2) + clutch_size+meantemp+date_aprildays+year, data=Asubset) #full model
lrtest(model1,model2) #drop quadratic interaction
model3 <- lmer(night_var~(1|box)+area +  day_before_hatch+ poly(day_before_hatch,2) + clutch_size+meantemp+date_aprildays+year, data=Asubset) #full model
lrtest(model2,model3) #keep linear interaction
model3.5 <- lmer(night_var~(1|box)+area +  area:poly(day_before_hatch,1) +  day_before_hatch + clutch_size+meantemp+date_aprildays+year, data=Asubset) #full model
lrtest(model2,model3.5) #drop quadratic

model4 <- lmer(night_var~(1|box)+area +  area:poly(day_before_hatch,1) +  day_before_hatch+meantemp+date_aprildays+year, data=Asubset) #full model
lrtest(model2,model4)  #drop clutch size
model5 <- lmer(night_var~(1|box)+area +  area:poly(day_before_hatch,1) +  day_before_hatch+date_aprildays+year, data=Asubset) #full model
lrtest(model4,model5) #drop mean temp
model6 <- lmer(night_var~(1|box)+area +  area:poly(day_before_hatch,1) +  day_before_hatch+year, data=Asubset) #full model
lrtest(model5,model6) #drop date
model7 <- lmer(night_var~(1|box)+area +  area:poly(day_before_hatch,1) +  day_before_hatch, data=Asubset) #full model
lrtest(model6,model7) #keep year
summary(model6)
write.csv( summary(model6)$coefficients, "C:\\Users\\Robyn\\Dropbox\\iButton Chronotype Chapter\\nightvar_take2.csv" )

#plotting night var

predict <- predict(model6, newdata = Asubset)
summary(predict)
write.csv(predict, "C:\\Users\\Robyn\\Dropbox\\iButton Chronotype Chapter\\nightvarpredictions.csv" )
write.csv(Asubset, "C:\\Users\\Robyn\\Dropbox\\iButton Chronotype Chapter\\nightvarsubset.csv" )

buttondata <- read.csv (file="C:\\Users\\Robyn\\Dropbox\\iButton Chronotype Chapter\\nightvarsubset.csv", na.strings=c("", "NA"))

gsubset <- subset(buttondata, select=c(box, year, area, day_before_hatch, night_var, predictmeans))
#remove missing rows, create index
Aindex <- complete.cases(gsubset)
#  subset[row,column]
Asubset <- gsubset[Aindex,]
Asubset$year <- as.factor(Asubset$year) 
Asubset$night_var <- as.numeric(Asubset$night_var)
Asubset$predictmeans <- as.numeric(Asubset$predictmeans)
str(Asubset)

theme_set(theme_bw(base_size = 15))
nightvarplot <- ggplot(data = Asubset, 
                     aes(x=day_before_hatch, 
                         y=night_var, 
                         fill=as.factor(area)))+
  geom_point(position = position_dodge(width = 0.45), 
             alpha = 0.60, 
             color = "black", 
             shape = 21) +
  theme_bw() +
  labs(y = "Night time temperature variance (C)", x = "Days before hatching") +
  scale_fill_manual(name = "", labels = c("City", "Forest"), 
                    values = c("#3399CC", "#339900")) +
  geom_point(aes(y=predictmeans, 
                 colour=as.factor(area)),
             position = position_dodge(width = 0.85), 
             size = 3.7,
             color = "black", 
             shape = 21)


nightvarplot + theme(legend.position="top")


## birds daylength ##
buttondata <- read.csv (file="C:\\Users\\Robyn\\Dropbox\\iButton Chronotype Chapter\\incubationdata.csv", na.strings=c("", "NA"))

gsubset <- subset(buttondata, select=c(box, date_aprildays, year, area, clutch_size, day_before_hatch, birds_day, meantemp))
#remove missing rows, create index
Aindex <- complete.cases(gsubset)
#  subset[row,column]
Asubset <- gsubset[Aindex,]
Asubset$year <- as.factor(Asubset$year) 
Asubset$day_before_hatch <- as.numeric(Asubset$day_before_hatch) 
str(Asubset) 

model1 <- lmer(birds_day~(1|box)+area + area:poly(day_before_hatch,2) + area:poly(day_before_hatch,1) +  day_before_hatch+ poly(day_before_hatch,2) + clutch_size+meantemp+date_aprildays+year, data=Asubset) #full model
summary(model1)
drop1(model1, test = "Chisq")
model2 <- lmer(birds_day~(1|box)+area +  area:poly(day_before_hatch,1) +  day_before_hatch+ poly(day_before_hatch,2) + clutch_size+meantemp+date_aprildays+year, data=Asubset) #full model
lrtest(model1,model2) #keep quadratic
model3 <- lmer(birds_day~(1|box)+area + area:day_before_hatch^2 + day_before_hatch+ poly(day_before_hatch,2) + clutch_size+meantemp+date_aprildays+year, data=Asubset) #full model
lrtest(model1,model3) #keep linear
model4 <- lmer(birds_day~(1|box)+area + area:poly(day_before_hatch,2) + area:poly(day_before_hatch,1) +  day_before_hatch+ poly(day_before_hatch,2) +meantemp+date_aprildays+year, data=Asubset) #full model
lrtest(model1,model4) #drop clutch size
model5 <- lmer(birds_day~(1|box)+area + area:poly(day_before_hatch,2) + area:poly(day_before_hatch,1) +  day_before_hatch+ poly(day_before_hatch,2) +date_aprildays+year, data=Asubset) #full model
lrtest(model4,model5)  #drop mean temp
model6 <- lmer(birds_day~(1|box)+area + area:poly(day_before_hatch,2) + area:poly(day_before_hatch,1) +  day_before_hatch+ poly(day_before_hatch,2) +year, data=Asubset) #full model
lrtest(model5,model6) #drop date
model7 <- lmer(birds_day~(1|box)+area + area:poly(day_before_hatch,2) + area:poly(day_before_hatch,1) +  day_before_hatch+ poly(day_before_hatch,2), data=Asubset) #full model
lrtest(model6,model7) #keep year
write.csv( summary(model6)$coefficients, "C:\\Users\\Robyn\\Dropbox\\iButton Chronotype Chapter\\daylength_take2.csv" )

plot(Asubset$birds_day~Asubset$area)
hist(Asubset$birds_day)

summary(model6)
#predictions#
predict <- predict(model6, newdata = Asubset)
summary(predict)
write.csv(predict, "C:\\Users\\Robyn\\Dropbox\\iButton Chronotype Chapter\\birdsdaypredictions.csv" )
write.csv(Asubset, "C:\\Users\\Robyn\\Dropbox\\iButton Chronotype Chapter\\birdsdaysubset.csv" )


#plotting day length
buttondata <- read.csv (file="C:\\Users\\Robyn\\Dropbox\\iButton Chronotype Chapter\\birdsdaysubset.csv", na.strings=c("", "NA"))

gsubset <- subset(buttondata, select=c(box, year, area, day_before_hatch, birds_day, predictmeans))
#remove missing rows, create index
Aindex <- complete.cases(gsubset)
#  subset[row,column]
Asubset <- gsubset[Aindex,]
Asubset$year <- as.factor(Asubset$year) 
Asubset$birds_day <- as.numeric(Asubset$birds_day)
Asubset$predictmeans <- as.numeric(Asubset$predictmeans)
str(Asubset)

theme_set(theme_bw(base_size = 15))
birdsdayplot <- ggplot(data = Asubset, 
                     aes(x=day_before_hatch, 
                         y=birds_day, 
                         fill=as.factor(area)))+
  geom_point(position = position_dodge(width = 0.45), 
             alpha = 0.60, 
             color = "black", 
             shape = 21) +
  theme_bw() +
  labs(y = "Duration of active day (in minutes)", x = "Days before hatching") +
  scale_fill_manual(name = "", labels = c("City", "Forest"), 
                    values = c("#3399CC", "#339900")) +
  geom_point(aes(y=predictmeans, 
                 colour=as.factor(area)),
             position = position_dodge(width = 0.85), 
             size = 3.7,
             color = "black", 
             shape = 21)


birdsdayplot + theme(legend.position="top")



###fitness models chronotype####
#chick weight#
buttondata <- read.csv (file="C:\\Users\\Robyn\\Dropbox\\iButton Chronotype Chapter\\chickfitness_ibuttons.csv", na.strings=c("", "NA"))
gsubset <- subset(buttondata, select=c(Chick, Box, Weight_d13, Area, Year, Brood_size, Hatch_date_aprildays, chronotype))
#remove missing rows, create index
Aindex <- complete.cases(gsubset)
#  subset[row,column]
Asubset <- gsubset[Aindex,]
str(Asubset)
Asubset$Year <- as.factor(Asubset$Year)

model1 <- lmer(Weight_d13~chronotype*Area+chronotype+Area+Year+Hatch_date_aprildays+Brood_size+(1|Box), data=Asubset)
summary(model1)
model2 <- lmer(Weight_d13~chronotype+Area+Year+Hatch_date_aprildays+Brood_size+(1|Box), data=Asubset)
lrtest(model1,model2) # keep interaction
model3 <- lmer(Weight_d13~chronotype*Area+chronotype+Area+Hatch_date_aprildays+Brood_size+(1|Box), data=Asubset)
lrtest(model1,model3) #drop year
model4 <- lmer(Weight_d13~chronotype*Area+chronotype+Area+Brood_size+(1|Box), data=Asubset)
lrtest(model3,model4) #keep hatch date
model5 <- lmer(Weight_d13~chronotype*Area+chronotype+Area+Hatch_date_aprildays+(1|Box), data=Asubset)
lrtest(model3,model5) # drop brood size
summary(model5)
write.csv( summary(model5)$coefficients, "C:\\Users\\Robyn\\Dropbox\\iButton Chronotype Chapter\\chickweight_take2.csv" )


#predictions#
predict <- predict(model5, newdata = Asubset)
summary(predict)
write.csv(predict, "C:\\Users\\Robyn\\Dropbox\\iButton Chronotype Chapter\\weightpredictions.csv" )
write.csv(Asubset, "C:\\Users\\Robyn\\Dropbox\\iButton Chronotype Chapter\\weightsubset.csv" )

#plotting weight #
buttondata <- read.csv (file="C:\\Users\\Robyn\\Dropbox\\iButton Chronotype Chapter\\weightsubset.csv", na.strings=c("", "NA"))

gsubset <- subset(buttondata, select=c(Box, Area, Weight_d13, chronotype, predictweight))
#remove missing rows, create index
Aindex <- complete.cases(gsubset)
#  subset[row,column]
Asubset <- gsubset[Aindex,]
Asubset$year <- as.factor(Asubset$year) 
Asubset$diff_sunrise <- as.numeric(Asubset$diff_sunrise)
Asubset$predictmeans <- as.numeric(Asubset$predictmean)
str(Asubset)

library(ggplot2)
theme_set(theme_bw(base_size = 15))
weightplot <- ggplot(data = Asubset, 
                    aes(x=chronotype, 
                        y=Weight_d13, 
                        fill=as.factor(Area)))+
  geom_point(position = position_dodge(width = 0.45), 
             alpha = 0.60, 
             color = "black", 
             shape = 21) +
  theme_bw() +
  labs(y = "Chick weight (g)", x = "Mean activity onset of mother relative to sunrise (chronotype)") +
  scale_fill_manual(name = "", labels = c("City", "Forest"), 
                    values = c("#3399CC", "#339900")) +
  geom_point(aes(y=predictweight, 
                 colour=as.factor(Area)),
             position = position_dodge(width = 0.85), 
             size = 3.7,
             color = "black", 
             shape = 21)

weightplot + theme(legend.position="top")



#old weight plot #
library(ggplot2)
theme_set(theme_bw(base_size = 15))
weightplot <- ggplot(data = Asubset, 
                       aes(x=chronotype, 
                           y=Weight_d13, 
                           fill=as.factor(Area)))+
  geom_point(position = position_dodge(width = 0.45), 
             alpha = 0.60, 
             color = "black", 
             shape = 21,
             size = 2) +
  theme_bw() +
  labs(y = "Chick weight (g)", x = "Mean activity onset of mother relative to sunrise (chronotype)") +
  scale_fill_manual(name = "", labels = c("City", "Forest"), 
                    values = c("#3399CC", "#339900"))


weightplot + theme(legend.position="top")



###hatch date
buttondata <- read.csv (file="C:\\Users\\Robyn\\Dropbox\\iButton Chronotype Chapter\\hatchsuccessdata.csv", na.strings=c("", "NA"))
gsubset <- subset(buttondata, select=c(Box, Area, clutch_size, Hatch_date, Year, chronotype, Hatched))
#remove missing rows, create index
Aindex <- complete.cases(gsubset)
#  subset[row,column]
Asubset <- gsubset[Aindex,]
Asubset$Year <- as.factor(Asubset$Year)
Asubset$Hatched<- as.factor(Asubset$Hatched)

model1 <- glmer(Hatched~Area+chronotype+(chronotype*Area)+clutch_size+Hatch_date+Year+(1|Box), data=Asubset, family=binomial)
summary(model1)
model1.5 <- glmer(Hatched~Area+chronotype+(chronotype*Area)+clutch_size+Year+(1|Box), data=Asubset, family=binomial)
lrtest(model1,model1.5) # hatch date not significant

model2 <- glmer(Hatched~Area+chronotype+clutch_size+Year+(1|Box), data=Asubset, family=binomial)
lrtest(model1,model2) #drop interaction
model3 <- glmer(Hatched~Area+chronotype+Year+(1|Box), data=Asubset, family=binomial)
lrtest(model2,model3) #drop clutch size
model4 <- glmer(Hatched~Area+chronotype+(1|Box), data=Asubset, family=binomial)
lrtest(model3,model4) #drop year
model5 <- glmer(Hatched~chronotype+(1|Box), data=Asubset, family=binomial)
lrtest(model4,model5) #drop area
model6 <- glmer(Hatched~Area+(1|Box), data=Asubset, family=binomial)
lrtest(model4,model6) #drop chronotype
summary(model4)
write.csv( summary(model4)$coefficients, "C:\\Users\\Robyn\\Dropbox\\iButton Chronotype Chapter\\hatchsuccess_take2.csv" )


#fledging success##
buttondata <- read.csv (file="C:\\Users\\Robyn\\Dropbox\\iButton Chronotype Chapter\\fledgesuccessdata.csv", na.strings=c("", "NA"))
gsubset <- subset(buttondata, select=c(Box, Area, Brood_size, Year, chronotype, Hatch_date, Fledged))
#remove missing rows, create index
Aindex <- complete.cases(gsubset)
#  subset[row,column]
Asubset <- gsubset[Aindex,]
str(Asubset)
Asubset$Year <- as.factor(Asubset$Year)
Asubset$Fledged<- as.factor(Asubset$Fledged)

model1 <- glmer(Fledged~(chronotype*Area)+Area+chronotype+Hatch_date+Brood_size+Year+(1|Box), data=Asubset, family=binomial)
summary(model1)
model2 <- glmer(Fledged~Area+chronotype+Hatch_date+Brood_size+Year+(1|Box), data=Asubset, family=binomial)
lrtest(model1,model2) #drop interaction
model3 <- glmer(Fledged~Area+chronotype+Hatch_date+Year+(1|Box), data=Asubset, family=binomial)
lrtest(model2,model3) #drop brood size
model4 <- glmer(Fledged~Area+chronotype+Hatch_date+(1|Box), data=Asubset, family=binomial)
lrtest(model3,model4) #drop year
model5 <- glmer(Fledged~Area+chronotype+(1|Box), data=Asubset, family=binomial)
lrtest(model4,model5) # drop hatch date
model6 <- glmer(Fledged~Area+(1|Box), data=Asubset, family=binomial)
lrtest(model5,model6) #drop chronotype
model7 <- glmer(Fledged~chronotype+(1|Box), data=Asubset, family=binomial)
lrtest(model5,model7) #drop area
summary(model5)
str(Asubset)
write.csv( summary(model5)$coefficients, "C:\\Users\\Robyn\\Dropbox\\iButton Chronotype Chapter\\fledgesuccess_take2.csv" )
