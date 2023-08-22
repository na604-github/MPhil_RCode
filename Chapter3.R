### Chapter 3: Reductions in Orangutan Flange Size after Initial Development: Relationship to FAI, Antagonistic Interactions, and Age

# Load packages
setwd("/Volumes/mawasHDD/MPhil/behavioural_data")
library(ggplot2)
library(vegan)
library(Rlab)
library(DataExplorer)
library(lattice)
library(dplyr)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(reshape2)
library(car)
library(lme4)
library(lubridate)
library(nortest)
library(irr)
library(COMPoissonReg)
library(glm2)
library(lmerTest)
library(forcats)
library(lubridate)

# Load data

data <- read.csv("flangeandotherdata.csv", header = TRUE)
head(data)

flange_data <- subset(data, subset=="y")

# Clean and tidy data

### Turn year and month into factors

flange_data$Month <- as.factor(flange_data$Month)
flange_data$Name <- as.factor(flange_data$Name)
flange_data$Site <- as.factor(flange_data$Site)

flange_data_clean <- subset(flange_data, !is.na(Site) & !is.na(Name))
flange_data_clean <- droplevels(flange_data_clean)

# Create flange proportions & FAI

flange_data$FAI_average_3month <- as.numeric(flange_data$FAI_average_3month)

flange_data <- flange_data %>%
  group_by(Site) %>%
  mutate(zFAI_average_3month = scale(FAI_average_3month))

flange_data <- flange_data %>%
  group_by(Name) %>%
  mutate(Prop_eyeflange = relative_eyeflange_to_eyeflange / max(relative_eyeflange_to_eyeflange, na.rm = TRUE))

flange_data <- flange_data %>%
  group_by(Name) %>%
  mutate(Prop_flangearea = relative_flange_size / max(relative_flange_size, na.rm = TRUE))

flange_data_subset <- flange_data %>%
  select(Name, Site, FAI, FAI_lag_1month, FAI_lag_2month, FAI_lag_3month, FAI_average_3month, relative_flange_size, relative_eyeflange_to_eyeflange, YearSinceFirstSeenTab, Prop_eyeflange, Prop_flangearea, )

flange_data$FemAssocH <- as.numeric(flange_data$FemAssocH)

flange_data <- flange_data %>%
  group_by(Site) %>%
  mutate(zAssocH = scale(FemAssocH))


# Further subsets
flange_suaq <- subset(flange_data, field_site=="suaq")
flange_tuanan <- subset(flange_data, field_site=="tuanan")
flange_data_headon <- subset(flange_data=="y")


###### Analysis #######

# Model 1
  model_1_null <- lmer(Prop_eyeflange ~ 1 + (1 | Name), data = flange_data_subset)

  model_1 <- lmer(Prop_eyeflange ~ YearSinceFirstSeenTab+zFAI_average_3month*Site + (1 | Name), data = flange_data_subset)
  summary(model_1)

  # Check assumptions

  # Plot residuals vs. fitted values
  fitted_vals <- fitted(model_1)
  resid_vals <- residuals(model_1)
  plot(fitted_vals, resid_vals, main="Residuals vs. Fitted", xlab="Fitted values", ylab="Residuals")
  abline(h=0, col="red")

  # Normality of residuals
  qqnorm(resid_vals)
  qqline(resid_vals)

  # Multicolinearity
  vif(lm(Prop_eyeflange ~ YearSinceFirstSeenTab + zFAI_average_3month*Site, data = flange_data_subset))

  # Check model is significantly better than null

  aic_null <- AIC(model_1_null)
  aic_full <- AIC(model_1)
  aic_difference <- aic_null - aic_full
  aic_difference
  
  anova(model_1_null, model_1, test="Chisq")
  

# Model 2
  model_2_null <- lmer(Prop_flangearea ~ 1 + (1 | Name), data = flange_data_subset)
  
  model_2 <- lmer(Prop_flangearea ~ YearSinceFirstSeenTab+zFAI_average_3month*Site + (1 | Name), data = flange_data_subset)
  summary(model_2)
  
  # Check assumptions
  
  # Plot residuals vs. fitted values
  fitted_vals <- fitted(model_2)
  resid_vals <- residuals(model_2)
  plot(fitted_vals, resid_vals, main="Residuals vs. Fitted", xlab="Fitted values", ylab="Residuals")
  abline(h=0, col="red")
  
  # Normality of residuals
  qqnorm(resid_vals)
  qqline(resid_vals)
  
  # Multicolinearity
  vif(lm(Prop_flangearea ~ YearSinceFirstSeenTab + zFAI_average_3month*Site, data = flange_data_subset))
  
  # Check model is significantly better than null
  
  aic_null <- AIC(model_2_null)
  aic_full <- AIC(model_2)
  aic_difference <- aic_null - aic_full
  aic_difference
  
  anova(model_2_null, model_2, test="Chisq")
  
# Model 3
  model_3_null <- lmer(Prop_eyeflange ~ 1 + (1 | Name), data = flange_data_subset) # Same null as before
  
  model_3_full <- lmer(Prop_eyeflange ~ YearSinceFirstSeenTab+zFAI_average_3month+Site + (1 | Name), data = flange_data_subset)
  summary(model_3_full)
  
  # Check assumptions
  
  # Plot residuals vs. fitted values
  fitted_vals <- fitted(model_3)
  resid_vals <- residuals(model_3)
  plot(fitted_vals, resid_vals, main="Residuals vs. Fitted", xlab="Fitted values", ylab="Residuals")
  abline(h=0, col="red")
  
  # Normality of residuals
  qqnorm(resid_vals)
  qqline(resid_vals)
  
  # Multicolinearity
  vif(lm(Prop_eyeflange ~ YearSinceFirstSeenTab + zFAI_average_3month+Site, data = flange_data_subset))
  
  # Check model is significantly better than null
  
  aic_null <- AIC(model_3_null)
  aic_full <- AIC(model_3_full)
  aic_difference <- aic_null - aic_full
  aic_difference
  
  anova(model_3_null, model_3_full, test="Chisq")
  
# Model 4
  
  model_4_null <- lmer(Prop_flangearea ~ 1 + (1 | Name), data = flange_data_subset) # Same as model 2
  
  model_4_full <- lmer(Prop_flangearea ~ YearSinceFirstSeenTab+zFAI_average_3month+Site + (1 | Name), data = flange_data_subset)
  summary(model_4_full)
  
  # Check assumptions
  
  # Plot residuals vs. fitted values
  fitted_vals <- fitted(model_4_full)
  resid_vals <- residuals(model_4_full)
  plot(fitted_vals, resid_vals, main="Residuals vs. Fitted", xlab="Fitted values", ylab="Residuals")
  abline(h=0, col="red")
  
  # Normality of residuals
  qqnorm(resid_vals)
  qqline(resid_vals)
  
  # Multicolinearity
  vif(lm(Prop_flangearea ~ YearSinceFirstSeenTab + zFAI_average_3month+Site, data = flange_data_subset))
  
  # Check model is significantly better than null
  
  aic_null <- AIC(model_4_null)
  aic_full <- AIC(model_4_full)
  aic_difference <- aic_null - aic_full
  aic_difference
  
  anova(model_4_null, model_4, test="Chisq")
  
#### Tuanan only models #####
  
  # Model 5
  model_5_null <- lmer(Prop_eyeflange ~ 1 + (1 | Name), data = flange_tuanan)
  
  model_5_full <- lmer(Prop_eyeflange ~ YearSinceFirstSeenTab+antag_cumulative*zFAI_average_3month + (1 | Name), data = flange_tuanan)
  summary(model_5_full)
  
  # Check assumptions
  
  # Plot residuals vs. fitted values
  fitted_vals <- fitted(model_5_full)
  resid_vals <- residuals(model_5_full)
  plot(fitted_vals, resid_vals, main="Residuals vs. Fitted", xlab="Fitted values", ylab="Residuals")
  abline(h=0, col="red")
  
  # Normality of residuals
  qqnorm(resid_vals)
  qqline(resid_vals)
  
  # Multicolinearity
  vif(lm(Prop_eyeflange ~ YearSinceFirstSeenTab+antag_cumulative*zFAI_average_3month, data = flange_tuanan))
  
  # Check model is significantly better than null
  
  aic_null <- AIC(model_5_null)
  aic_full <- AIC(model_5_full)
  aic_difference <- aic_null - aic_full
  aic_difference
  
  anova(model_5_null, model_5_full, test="Chisq")
  
  
  # Model 6
  model_6_null <- lmer(Prop_flangearea ~ 1 + (1 | Name), data = flange_tuanan)
  
  model_6_full <- lmer(Prop_flangearea ~ YearSinceFirstSeenTab+antag_cumulative*zFAI_average_3month + (1 | Name), data = flange_tuanan)
  summary(model_6_full)
  
  # Check assumptions
  
  # Plot residuals vs. fitted values
  fitted_vals <- fitted(model_6_full)
  resid_vals <- residuals(model_6_full)
  plot(fitted_vals, resid_vals, main="Residuals vs. Fitted", xlab="Fitted values", ylab="Residuals")
  abline(h=0, col="red")
  
  # Normality of residuals
  qqnorm(resid_vals)
  qqline(resid_vals)
  
  # Multicolinearity
  vif(lm(Prop_flangearea ~ YearSinceFirstSeenTab+antag_cumulative*zFAI_average_3month + (1 | Name), data = flange_tuanan))
  
  # Check model is significantly better than null
  
  aic_null <- AIC(model_6_null)
  aic_full <- AIC(model_6_full)
  aic_difference <- aic_null - aic_full
  aic_difference
  
  anova(model_6_null, model_6_full, test="Chisq")
  
  # Model 7
  model_7_null <- lmer(Prop_eyeflange ~ 1 + (1 | Name), data = flange_tuanan) # Same null as before
  
  model_7_full <- lmer(Prop_eyeflange ~ YearSinceFirstSeenTab+antag_cumulative+zFAI_average_3month + (1 | Name), data = flange_tuanan)
  summary(model_7_full)
  
  # Check assumptions
  
  # Plot residuals vs. fitted values
  fitted_vals <- fitted(model_7_full)
  resid_vals <- residuals(model_7_full)
  plot(fitted_vals, resid_vals, main="Residuals vs. Fitted", xlab="Fitted values", ylab="Residuals")
  abline(h=0, col="red")
  
  # Normality of residuals
  qqnorm(resid_vals)
  qqline(resid_vals)
  
  # Multicolinearity
  vif(lm(Prop_eyeflange ~ YearSinceFirstSeenTab+antag_cumulative+zFAI_average_3month + (1 | Name), data = flange_tuanan))
  
  # Check model is significantly better than null
  
  aic_null <- AIC(model_7_null)
  aic_full <- AIC(model_7_full)
  aic_difference <- aic_null - aic_full
  aic_difference
  
  anova(model_7_null, model_7_full, test="Chisq")
  
  # Model 8
  
  model_8_null <- lmer(Prop_flangearea ~ 1 + (1 | Name), data = flange_tuanan) # Same as model 2
  
  model_8_full <- lmer(Prop_flangearea ~ YearSinceFirstSeenTab+antag_cumulative+zFAI_average_3month + (1 | Name), data = flange_tuanan)
  summary(model_8_full)
  
  # Check assumptions
  
  # Plot residuals vs. fitted values
  fitted_vals <- fitted(model_8_full)
  resid_vals <- residuals(model_8_full)
  plot(fitted_vals, resid_vals, main="Residuals vs. Fitted", xlab="Fitted values", ylab="Residuals")
  abline(h=0, col="red")
  
  # Normality of residuals
  qqnorm(resid_vals)
  qqline(resid_vals)
  
  # Multicolinearity
  vif(lm(Prop_flangearea ~ YearSinceFirstSeenTab+antag_cumulative+zFAI_average_3month + (1 | Name)))
  
  # Check model is significantly better than null
  
  aic_null <- AIC(model_8_null)
  aic_full <- AIC(model_8_full)
  aic_difference <- aic_null - aic_full
  aic_difference
  
  anova(model_8_null, model_8_full, test="Chisq")
  
# Model 9

library(glmmTMB) # Loading to account for zero-inflation in behavioural data

  model_9_null_0 <- glmmTMB(LCM ~ 1 +offset(ObsTime) +(1|Name/Month/Year),
                           data=flange_data, 
                           family=poisson(link = "log"))
  model_9_full_0 <- glmmTMB(LCM ~ Prop_eyeflange + zFAI + Site + LCH + zAssocH + total.rain + offset(ObsTime) +(1|Name/Month/Year),
                           data=flange_data, 
                           family=poisson(link = "log"))
  summary(model_9_full_0)
  
  # Check assumptions
  
  overdisp_ratio <- as.numeric(sum(resid(model_9_full_0, type = "pearson")^2) / model_9_full_0$df.residual)
  

  # Overdispersion, going to use negative binomial instead
  
  model_9_null <- glmmTMB(LCM ~ 1 + offset(ObsTime) + (1|Name/Month/Year),
                             data=flange_data, 
                             family=nbinom2(link = "log"))
  model_9_full <- glmmTMB(LCM ~ Prop_eyeflange + zFAI + Site + LCH + zAssocH + total.rain + offset(ObsTime) +(1|Name/Month/Year),
                             data=flange_data, 
                             family=nbinom2(link = "log"))
  summary(model_9_full)
  
  # Test for overdispersion
  
  plot(resid(model_9_full, type = "pearson") ~ fitted(model_9_full))
  abline(h = 0, col = "red")
  
  library(lmtest)
  lrtest(model_9_full, model_9_full_0) # Significantly better fit
  
  # Check model is significantly better than null
  
  aic_null <- AIC(model_9_null)
  aic_full <- AIC(model_9_full)
  aic_difference <- aic_null - aic_full
  aic_difference
  
  anova(model_9_null, model_9_full, test="Chisq")
  
# Model 10
  
  model_10_null_0 <- glmmTMB(AvgPulseNr ~ 1 +offset(LCM) +(1|Name/Month/Year),
                            data=flange_data, 
                            family=poisson(link = "log"))
  model_10_full_0 <- glmmTMB(AvgPulseNr ~ Prop_eyeflange + zFAI + Site + LCH + zAssocH + total.rain + offset(LCM) +(1|Name/Month/Year),
                            data=flange_data, 
                            family=poisson(link = "log"))
  summary(model_10_full_0) # Singularity issues, reducing random effect structure

  model_10_null <- glmmTMB(AvgPulseNr ~ 1 + offset(LCM) + (1|Name/Month/Year),
                           data=flange_data, 
                           family=poisson(link = "log"))
  model_10_full <- glmmTMB(AvgPulseNr ~ Prop_eyeflange + zFAI + Site + LCH + zAssocH + total.rain + offset(LCM) +(1|Name/Month/Year),
                           data=flange_data, 
                           family=poisson(link = "log"))
  summary(model_10_full)
  
  # Check assumptions
  
  overdisp_ratio <- as.numeric(sum(resid(model_10_full_0, type = "pearson")^2) / model_10_full$df.residual)
  
  # Overdispersion, going to use negative binomial instead
  
  model_10_null <- glmmTMB(AvgPulseNr ~ 1 + offset(LCM) + (1|Name/Month/Year),
                          data=flange_data, 
                          family=nbinom2(link = "log"))
  model_10_full <- glmmTMB(AvgPulseNr ~ Prop_eyeflange + zFAI + Site + LCH + zAssocH + total.rain + offset(LCM) +(1|Name/Month/Year),
                          data=flange_data, 
                          family=nbinom2(link = "log"))
  summary(model_10_full)
  
  
  # Test for overdispersion
  
  plot(resid(model_10_full, type = "pearson") ~ fitted(model_9_full))
  abline(h = 0, col = "red")
  
  lrtest(model_10_full, model_10_full_0) # Significantly better fit
  
  
  # Check model is significantly better than null
  
  aic_null <- AIC(model_10_null)
  aic_full <- AIC(model_10_full)
  aic_difference <- aic_null - aic_full
  aic_difference
  
  anova(model_10_null, model_10_full, test="Chisq")
  
# Model 11
  
  model_11_null <- glmmTMB(TotalCop ~ 1 + offset(FemAssocH) +(1|Name/Month/Year),
                             data=flange_data, 
                             family=poisson(link = "log"))
  model_11_full <- glmmTMB(TotalCop ~ Prop_eyeflange + zFAI + Site + LCM + offset(FemAssocH) +(1|Name/Month/Year),
                             data=flange_data, 
                             family=poisson(link = "log"))
  summary(model_11_full)
  
  # Check assumptions
  
  overdisp_ratio <- as.numeric(sum(resid(model_11_full, type = "pearson")^2) / model_11_full$df.residual)
  plot(resid(model_11_full, type = "pearson") ~ fitted(model_11_full))
  abline(h = 0, col = "red")
  
  # Check model is significantly better than null
  
  aic_null <- AIC(model_11_null)
  aic_full <- AIC(model_11_full)
  aic_difference <- aic_null - aic_full
  aic_difference
  
  anova(model_11_null, model_11_full, test="Chisq")
  
# Model 12
  
  model_12_null <- lmer(AvgDuration ~ 1 + (1|Name/Month/Year),
                                   data=flange_data)
  model_12_full <- lmer(AvgDuration ~ Prop_eyeflange + zFAI + Site + zAssocH + LCH + total.rain + (1|Name/Month/Year),
                                   data=flange_data)
  summary(model_12_full)
  
  # Check assumptions
  
  # Plot residuals vs. fitted values
  fitted_vals <- fitted(model_12_full)
  resid_vals <- residuals(model_12_full)
  plot(fitted_vals, resid_vals, main="Residuals vs. Fitted", xlab="Fitted values", ylab="Residuals")
  abline(h=0, col="red")
  
  # Normality of residuals
  qqnorm(resid_vals)
  qqline(resid_vals)
  
  # Multicolinearity
  vif(lm(AvgDuration ~ Prop_eyeflange + zFAI + Site + zAssocH + LCH + `total rain` + (1|Name/Month/Year), data=flange_data))
  
  # Check model is significantly better than null
  
  aic_null <- AIC(model_12_null)
  aic_full <- AIC(model_12_full)
  aic_difference <- aic_null - aic_full
  aic_difference
  
  anova(model_12_null, model_12_full, test="Chisq")
  
  
# Model 13
  
  model_13_null <- lmer(NrFemAssoc ~ 1 + (1|Name/Month/Year),
                        data=flange_data)
  model_13_full <- lmer(NrFemAssoc ~ Prop_eyeflange + zFAI + Site + LCM + (1|Name/Month/Year),
                        data=flange_data)
  summary(model_13_full)
  
  # Check assumptions
  
  # Plot residuals vs. fitted values
  fitted_vals <- fitted(model_13_full)
  resid_vals <- residuals(model_13_full)
  plot(fitted_vals, resid_vals, main="Residuals vs. Fitted", xlab="Fitted values", ylab="Residuals")
  abline(h=0, col="red")
  
  # Normality of residuals
  qqnorm(resid_vals)
  qqline(resid_vals)
  
  # Multicolinearity
  vif(lm(NrFemAssoc ~ Prop_eyeflange + zFAI + Site + LCM + (1|Name/Month/Year), data=flange_data))
  
  # Check model is significantly better than null
  
  aic_null <- AIC(model_13_null)
  aic_full <- AIC(model_13_full)
  aic_difference <- aic_null - aic_full
  aic_difference
  
  anova(model_13_null, model_13_full, test="Chisq")

######## Plots ###########

# Fig 3.2.

flange_data <- flange_data %>% arrange(Site, Name) # Order by site and name for better plot
flange_data$Date <- as.Date(flange_data$Date, format = "%d/%m/%Y") # Convert the Date column to Date class
flange_data$Name <- as.character(flange_data$Name) # Convert Name to a character vector
flange_data <- flange_data[order(flange_data$Name), ] # Sort the data frame by Name in alphabetical order
flange_data$Name <- factor(flange_data$Name, levels = unique(flange_data$Name)) # Convert Name back to an ordered factor variable


sightings_plot <- ggplot(flange_data, aes(x = Date, y = Name)) +
  geom_point(aes(color = Name), shape = 4, size = 2) +
  facet_grid(rows = vars(Site), scales = "free", space = "free") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  theme_minimal() +
  labs(x = "Year", y = "Individual Name", 
       title = "Orangutan follows where flange measurements are available") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))

ggsave("/Volumes/G-DRIVE/MPhil/behavioural_data/sightings_plot.eps", plot = sightings_plot, device = "eps", width = 8, height = 5)

# Fig. 3.3

fai_suaq <- read.csv("fai_suaq.csv", header = TRUE)
fai_tuanan <- read.csv("fai_tuanan.csv", header = TRUE)

fai_suaq_plot <- ggplot(data = fai_suaq, aes(x = FAIMonth, y = FAIFruit, color = factor(FAIYear))) +
  geom_point() +
  geom_smooth(se = FALSE) +
  #geom_smooth(aes(color = NULL), se = FALSE) + # Removed color from aes() in geom_smooth()
  theme_minimal() +
  labs(x = "Month", y = "Fruit Availability Index (FAI)", title = "Suaq FAI by Month and Year with Trend Lines") +
  scale_color_discrete(name = "Year") +
  scale_x_continuous(breaks = 1:12, labels = month.abb) + # Added this line to change month labels
  theme(legend.position = "bottom") +
  scale_y_continuous(limits = c(0, 20), breaks = c(0, 5, 10, 15, 20))

fai_tuanan_plot <- ggplot(data = fai_tuanan, aes(x = month, y = FAI, color = factor(year))) +
  geom_point() +
  geom_smooth(se = FALSE) +
  #geom_smooth(aes(color = NULL), se = FALSE) + # Removed color from aes() in geom_smooth()
  theme_minimal() +
  labs(x = "Month", y = "Fruit Availability Index (FAI)", title = "Tuanan FAI by Month and Year with Trend Lines") +
  scale_color_discrete(name = "Year") +
  scale_x_continuous(breaks = 1:12, labels = month.abb) + # Added this line to change month labels
  theme(legend.position = "bottom") +
  scale_y_continuous(limits = c(0, 20), breaks = c(0, 5, 10, 15, 20))

combined_plot <- grid.arrange(fai_suaq_plot, fai_tuanan_plot, ncol = 2)

ggsave("/Volumes/mawasHDD/MPhil/chapter_3/fai_plot.eps", plot = combined_plot, device = "eps", width = 8, height = 4)

# Fig. 3.4

model_3_predict <- ggpredict(model_3, terms = c("YearSinceFirstSeenTab")) ## model predictions of full model 

model_3_predict_plot <- ggplot()+
  geom_line(data= model_3_predict,
            mapping = aes(x=x, y=conf.low), linetype="dashed", color="black")+
  geom_line(data= model_3_predict,
            mapping = aes(x=x, y=conf.high), linetype="dashed", color="black")+
  geom_line(data= model_3_predict,
            mapping = aes(x=x, y=predicted), size=1.2)+
  geom_jitter(data = subset_merged_data,
              mapping= aes(YearSinceFirstSeenTab, Prop_eyeflange, colour = Site, shape = Site),
              width=0.15, height=0.001) +  
  theme_set(theme_cowplot())+
  theme(legend.box.background = element_rect(colour="black"),
        legend.box.margin = margin(2,4,2,4),
        title=element_text(size=12),
        legend.position = c(0, 0),    # position updated
        legend.justification = c(0, 0),  # aligns the legend box with the plot corners
        panel.grid.major = element_line(colour = "grey80"))+
  labs(x="Year Since First Observed",
       y="Proportion of Maximum Observed Flange Width",
       title = "Change in proportion of flange width by relative age",
       colour = "Site",
       shape = "Site") +
  scale_colour_brewer(palette = "Set1") +
  scale_y_continuous(limits = c(NA, 1))

model_3_predict_plot

ggsave("/Volumes/mawasHDD/MPhil/behavioural_data/flange_width.eps", plot = model_2_predict_plot, device = "eps", width = 6, height = 4)
