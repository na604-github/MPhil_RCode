### Analyses reported Chapter 4 #####
### Age-Related Behavioural Shifts in Flanged Male Orangutans:  
### Insights from Female Association Behaviour and Long Call Patterns

library(plyr)
library(glmmTMB)
library(multcomp)
library(DHARMa)
library(MuMIn)
library(ggplot2)
library(ggeffects)
library(effects)
library(cowplot)
library(car)
library(magick)
library(performance)
library(coxme)
library(tidyverse)

# import data set 
data_unclean <- read.csv("combined_data.csv",
                         stringsAsFactors = TRUE)

# Subset dataset to 
data_1 <- droplevels(subset(data, FollowType=="NN")) ## only NN-follows
data_2 <- subset(data_1, !(Morph == "ufm")) # Remove unflanged males
data_3 <- subset(data_2, !(Name == "flm(unid)")) # Remove unidentified flanged males
data_4 <- subset(data_3, !(Name == "Unknown Flanged Male")) # Remove unidentified flanged males
data_5 <- droplevels(subset(data_4, !is.na(FAI))) ## exclude days with unknown FAI
data_6 <-  droplevels(subset(data_5, !is.na(ObsTime))) ## exclude unknown Observation hours
data_7 <- data_6[data_6$Morph == "flm", ] ## exclude non-flanged males e.g. flanging
data <- droplevels(subset(data_7, !is.na(YearSinceFlanged))) ## exclude individuals not seen in unflanged morph


########## Long call models ###########

# 1. Number of long calls per follow day

#Null
lc_made_null <- glmmTMB(LCM ~ 1 + +(1|Name/Month/Year),
                        data=data, 
                        family=poisson(link = "log"))


summary(lc_made_null)

# Full
lc_made_model <- glmmTMB(LCM ~ YearSinceFlanged + FAI + femAssocH + total.rain + LCH +offset(ObsTime) +(1|Name/Month/Year),
                         data=data, 
                         family=poisson(link = "log"))


summary(lc_made_model)

# Check assumptions

#Linearity of Predictors
predicted_values <- predict(lc_made_model, type = "response")
plot(predicted_values, data$LCM, main = "Observed vs. Predicted", xlab = "Predicted", ylab = "Observed")
abline(0, 1)

#Overdispersion
dispersion_glmmTMB <- function(model){
  rdf <- df.residual(model)
  res_dev <- sum(resid(model, type = "pearson")^2)
  return(res_dev / rdf)
}

dispersion_value <- dispersion_glmmTMB(lc_made_model)
print(dispersion_value)


# Check model is significantly better than null

aic_null <- AIC(lc_made_null)
aic_full <- AIC(lc_made_model)
aic_difference <- aic_null - aic_full
aic_difference

anova(lc_made_null, lc_made_model, test="Chisq")

# 2. Number of long call pulses per long call
lc_nr_pulses_null <- glmmTMB(AvgPulseNr ~ 1 +(1|Name/Month/Year),
                             data=data, 
                             family=poisson(link = "log"))


summary(lc_nr_pulses_null)

lc_nr_pulses_full <- glmmTMB(AvgPulseNr ~ YearSinceFlanged + FAI + femAssocH + total.rain + LCH +offset(AvgDuration) +(1|Name/Year/Month),
                              data=data, 
                              family=poisson(link = "log"))


summary(lc_nr_pulses_full)

# Check assumptions

#Linearity of Predictors
predicted_values <- predict(lc_nr_pulses_full, type = "response")
plot(predicted_values, data$AvgPulseNr, main = "Observed vs. Predicted", xlab = "Predicted", ylab = "Observed")
abline(0, 1)

#Overdispersion
dispersion_glmmTMB <- function(model){
  rdf <- df.residual(model)
  res_dev <- sum(resid(model, type = "pearson")^2)
  return(res_dev / rdf)
}

dispersion_value <- dispersion_glmmTMB(lc_nr_pulses_full)
print(dispersion_value) #less than one, that's good


# Check model is significantly better than null

aic_null <- AIC(lc_nr_pulses_null)
aic_full <- AIC(lc_nr_pulses_full)
aic_difference <- aic_null - aic_full
aic_difference

anova(lc_nr_pulses_null, lc_nr_pulses_full, test="Chisq")

# 3. Long call duration

lc_duration_null_0 <- glmmTMB(AvgDuration ~ 1 + (1|Name/Month/Year),
                              data=data, family = gaussian())

lc_duration_full_0 <- glmmTMB(AvgDuration ~ YearSinceFlanged + FAI + femAssocH + LCH + total.rain + (1|Name/Month/Year),
                              data=data, family = gaussian()) # Singularity issues so reducing random effect structure

lc_duration_null <- glmmTMB(AvgDuration ~ 1 + (1|Name),
                            data=data, family = gaussian())

lc_duration_full <- glmmTMB(AvgDuration ~ YearSinceFlanged + FAI + femAssocH + LCH + total.rain + (1|Name),
                            data=data, family = gaussian())

summary(lc_duration_full)

# Check assumptions

#Linearity of Predictors
predicted_values <- predict(lc_duration_full, type = "response")
plot(predicted_values, data$AvgDuration, main = "Observed vs. Predicted", xlab = "Predicted", ylab = "Observed")
abline(0, 1)

#Overdispersion
dispersion_glmmTMB <- function(model){
  rdf <- df.residual(model)
  res_dev <- sum(resid(model, type = "pearson")^2)
  return(res_dev / rdf)
}

dispersion_value <- dispersion_glmmTMB(lc_duration_full)
print(dispersion_value)


# Check model is significantly better than null

aic_null <- AIC(lc_duration_null)
aic_full <- AIC(lc_duration_full)
aic_difference <- aic_null - aic_full
aic_difference

anova(lc_duration_null, lc_duration_full, test="Chisq")

########## Female association models ###########

# 4. Number of fertile females in association per follow day


# Null model
fertfemNr_null <- glmmTMB(NrFertFem ~ 1  + (1|Name/Year/Month),
                          data=data, 
                          family=poisson(link = "log"))
summary(fertfemNr_null)

# Full model
fertfemNr_full <- glmmTMB(NrFertFem ~ YearSinceFlanged + FAI + LCM + offset(ObsTime)+(1|Name/Year/Month),
                          data=data, 
                          family=poisson(link = "log"))

summary(fertfemNr_full)

# Check assumptions

#Linearity of Predictors
predicted_values <- predict(fertfemNr_full, type = "response")
plot(predicted_values, data$NrFertFem, main = "Observed vs. Predicted", xlab = "Predicted", ylab = "Observed")
abline(0, 1)

#Overdispersion
dispersion_glmmTMB <- function(model){
  rdf <- df.residual(model)
  res_dev <- sum(resid(model, type = "pearson")^2)
  return(res_dev / rdf)
}

dispersion_value <- dispersion_glmmTMB(fertfemNr_full)
print(dispersion_value)


# Check model is significantly better than null

aic_null <- AIC(fertfemNr_null)
aic_full <- AIC(fertfemNr_full)
aic_difference <- aic_null - aic_full
aic_difference

anova(fertfemNr_null, fertfemNr_full, test="Chisq")

# 5. Number of females in association (total)

# Null model
femNr_null <- glmmTMB(NrFemAssoc ~ 1  + (1|Name/Year/Month),
                      data=data, 
                      family=poisson(link = "log"))
summary(femNr_null)

# Full model
femNr_full <- glmmTMB(NrFemAssoc ~ YearSinceFlanged + FAI + LCM + offset(ObsTime)+(1|Name/Year/Month),
                      data=data, 
                      family=poisson(link = "log"))

summary(femNr_full)

# Check assumptions

#Linearity of Predictors
predicted_values <- predict(femNr_full, type = "response")
plot(predicted_values, data$femNr, main = "Observed vs. Predicted", xlab = "Predicted", ylab = "Observed")
abline(0, 1)

#Overdispersion
dispersion_glmmTMB <- function(model){
  rdf <- df.residual(model)
  res_dev <- sum(resid(model, type = "pearson")^2)
  return(res_dev / rdf)
}

dispersion_value <- dispersion_glmmTMB(femNr_full)
print(dispersion_value)


# Check model is significantly better than null

aic_null <- AIC(femNr_null)
aic_full <- AIC(femNr_full)
aic_difference <- aic_null - aic_full
aic_difference

anova(femNr_null, femNr_full, test="Chisq")


#6. Hours where females are in association

# Null model
femAssociation_null_0 <- glmmTMB(FemAssocH ~ 1  + (1|Name/Year/Month),
                                 data=data, 
                                 family=gaussian())
summary(femAssociation_null_0)


# Full model
femAssociation_full_0 <- glmmTMB(FemAssocH ~ YearSinceFlanged + FAI + LCM + offset(ObsTime)+(1|Name/Year/Month),
                                 data=data, 
                                 family=gaussian()) #Singularity issues, reducing number of random effects

# Null model
femAssociation_null <- glmmTMB(FemAssocH ~ 1  + (1|Name),
                               data=data, 
                               family=gaussian())
summary(femAssociation_null)


# Full model
femAssociation_full <- glmmTMB(FemAssocH ~ YearSinceFlanged + FAI + LCM + offset(ObsTime)+(1|Name),
                               data=data, 
                               family=gaussian())

summary(femAssociation_full)


# Check assumptions

#Linearity of Predictors
predicted_values <- predict(femAssociation_full, type = "response")
plot(predicted_values, data$FemAssocH, main = "Observed vs. Predicted", xlab = "Predicted", ylab = "Observed")
abline(0, 1)

#Overdispersion
dispersion_glmmTMB <- function(model){
  rdf <- df.residual(model)
  res_dev <- sum(resid(model, type = "pearson")^2)
  return(res_dev / rdf)
}

dispersion_value <- dispersion_glmmTMB(femAssociation_full)
print(dispersion_value)


# Check model is significantly better than null

aic_null <- AIC(femAssociation_null)
aic_full <- AIC(femAssociation_full)
aic_difference <- aic_null - aic_full
aic_difference

anova(femAssociation_null, femAssociation_full, test="Chisq")




########### Plots #################

#Fig 4.1.

obs_plot <- ggplot(data, aes(x = Date, y = Name)) +
  geom_point(aes(color = Name), shape = 4, size = 1.5, position = position_jitter(width = 0.1, height = 0)) +
  facet_wrap(~ Site, ncol = 2, scales = "free") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  theme_minimal() +
  labs(x = "Year", y = "Individual Name", 
       title = "Full day orangutan follows where the individual was previously observed in the\nunflanged morph") +
  theme(legend.position = "none")
#plot.title = element_text(hjust = 0.5),
#axis.text.x = element_text(size = 6),  # Adjust the x-axis font size
#axis.text.y = element_text(size = 5))  # Adjust the y-axis font size

obs_plot

ggsave("/Volumes/mawasHDD/MPhil/chapter_4/obs_plot.eps", plot = obs_plot, device = "eps", width = 8, height = 6)


#Fig 4.2.
femAssociation_predict <- ggpredict(femAssociation_full, terms = c("YearSinceFlanged")) ## model predictions of full model 


femAssociation_plot <- ggplot() +
  geom_line(data= femAssociation_predict,
            mapping = aes(x=x, y=conf.low), linetype="dashed", color="black") +
  geom_line(data= femAssociation_predict,
            mapping = aes(x=x, y=conf.high), linetype="dashed", color="black") +
  geom_line(data= femAssociation_predict,
            mapping = aes(x, predicted), size=1.2) +
  geom_jitter(data=data,
              mapping= aes(YearSinceFlanged, FemAssocH),
              color="#2596be",
              width=0.15, height=0.1) +
  theme_set(theme_cowplot()) +
  labs(x="Years since flanged",
       y="Hours of female association per\nfollow day",
       title = "Total hours in association with a female per full focal follow day") +
  theme(legend.box.background = element_rect(colour="black"),
        legend.box.margin = margin(2,4,2,4),
        title=element_text(size=12),
        legend.position = "none", # No legend
        plot.title.position = "plot",
        panel.grid.major = element_line(colour = "grey80")) +
  scale_colour_brewer(palette = "Set1") +
  scale_x_continuous(limits = c(0,9), breaks = seq(0,9,1)) # Scale axis 1 to 10, separation of one


femAssociation_plot

ggsave("/Volumes/mawasHDD/MPhil/chapter_4/femAssociation_plot.eps", plot = FemAssocPlot, device = "eps", width = 8, height = 4)


# Fig 4.3. 

lc_nr_pulses_full_predict <- ggpredict(lc_nr_pulses_full, terms = c("YearSinceFlanged")) ## model predictions of full model 

lc_pulses_plot <- ggplot()+
  geom_line(data= lc_nr_pulses_full_predict,
            mapping = aes(x=x, y=conf.low), linetype="dashed", color="black")+
  geom_line(data= lc_nr_pulses_full_predict,
            mapping = aes(x=x, y=conf.high), linetype="dashed", color="black")+
  geom_line(data= lc_nr_pulses_full_predict,
            mapping = aes(x, predicted), size=1.2)+
  geom_jitter(data,
              mapping= aes(YearSinceFlanged, AvgPulseNr),
              width=0.15, height=0.1, col = "#367EB8")+
  theme_set(theme_cowplot()) +
  theme(legend.box.background = element_rect(colour="black"),
        legend.box.margin = margin(2,4,2,4),
        title=element_text(size=12, hjust=0.5),
        legend.position = c(0.65,0.9),
        panel.grid.major = element_line(colour = "grey80"))+
  labs(x="Years since flanged",
       y="Average number of pulses per\nlong call per follow day",
       title = "Average number of pulses per\nlong call per follow day",
       colour = "Site") +
  scale_colour_brewer(palette = "Set1")+
  guides(col="none")+
  scale_x_continuous(limits = c(0,9), breaks = seq(0,9,1)) # Scale axis 1 to 10, separation of one

lc_pulses_plot

lc_duration_predict <- ggpredict(lc_duration_full, terms = c("YearSinceFlanged")) ## model predictions of full model 

lc_duration_plot <- ggplot()+
  geom_line(data= lc_duration_predict,
            mapping = aes(x=x, y=conf.low), linetype="dashed", color="black")+
  geom_line(data= lc_duration_predict,
            mapping = aes(x=x, y=conf.high), linetype="dashed", color="black")+
  geom_line(data= lc_duration_predict,
            mapping = aes(x, predicted), size=1.2)+
  geom_jitter(data,
              mapping= aes(YearSinceFlanged, AvgDuration),
              width=0.15, height=0.1, col = "#367EB8")+
  theme_set(theme_cowplot()) +
  theme(legend.box.background = element_rect(colour="black"),
        legend.box.margin = margin(2,4,2,4),
        title=element_text(size=12, hjust=0.5),
        legend.position = c(0.65,0.9),
        panel.grid.major = element_line(colour = "grey80"))+
  labs(x="Years since flanged",
       y="Average long call duration per\nfollow day (s)",
       title = "Average long call duration per\nfollow day",
       colour = "Site") +
  scale_colour_brewer(palette = "Set1")+
  guides(col="none")+
  scale_x_continuous(limits = c(0,9), breaks = seq(0,9,1)) + # Scale axis 1 to 10, separation of one 
  scale_y_continuous(limits = c(0,120), breaks = seq(0,120,10)) # Scale axis 1 to 10, separation of one

lc_duration_plot

combined_plot <- grid.arrange(lc_pulses_plot, lc_duration_plot, ncol = 2)

ggsave("/Volumes/mawasHDD/MPhil/chapter_4/combined_plot.eps", plot = combined_plot, device = "eps", width = 8, height = 4)
