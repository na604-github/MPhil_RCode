# Chapter 2 

## Load relevant packages
# Loading commonly used packages and setting working directory
#setwd("/Volumes/LaCie/MPhil/analysis")
setwd("/Volumes/G-Drive/MPhil/analysis")
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
 

## Load data
 
flange <- read.csv("flange_measurements.csv", header = TRUE) # Import dataset into R
flange_hq <- subset(flange, head_on=="y")

replicates <- subset(flange_hq, replicate=="y")
originals <- subset(flange_hq, replicate=="n")

summary(flange_hq)
 

## Normality plots

ggplot(flange_hq, aes(x = H1_signed)) +
  geom_freqpoly(binwidth = 0.1)
ggplot(flange_hq, aes(x = H2_signed)) +
  geom_freqpoly(binwidth = 0.1)
ggplot(flange_hq, aes(x = H3_signed)) +
  geom_freqpoly(binwidth = 0.1)
ggplot(flange_hq, aes(x = V1_signed)) +
  geom_freqpoly(binwidth = 0.1)
ggplot(flange_hq, aes(x = V2_signed)) +
  geom_freqpoly(binwidth = 0.1)
ggplot(flange_hq, aes(x = V3_signed)) +
  geom_freqpoly(binwidth = 0.1)

ggplot(flange_hq, aes(x = innereye)) +
  

## Repeatability tests

#A high correlation indicates repetability of measurements.

flange_repeat <- read.csv("flange_measurements_repeatability.csv", header = TRUE) # Import dataset into R
summary(flange_repeat)
flange_a <- subset(flange_repeat, photo=="a")
flange_b <- subset(flange_repeat, photo=="b")

cor.test(flange_a$H1,flange_b$H1)
cor.test(flange_a$H2,flange_b$H2)
cor.test(flange_a$H3,flange_b$H3)
cor.test(flange_a$V1,flange_b$V1)
cor.test(flange_a$V2,flange_b$V2)
cor.test(flange_a$V3,flange_b$V3)


cor.test(originals$H1,replicates$H1)
cor.test(originals$H2,replicates$H2)
cor.test(originals$H3,replicates$H3)
cor.test(originals$V1,replicates$V1)
cor.test(originals$V2,replicates$V2)
cor.test(originals$V3,replicates$V3)

H1_ICC <- data.frame(A=c(8.666666667,2.666666667,17.33333333,9,10,11,8,13.5,3.333333333,6,11.83333333,2.666666667,6,2,6.333333333,8,1.5,1.166666667,7,8,10,6,7),
                     B=c(7.333333333,1,15.5,9.333333333,7.5, 13.83333333,5.5,11.33333333,5.333333333,8.333333333,16,0,5.333333333,2.5,5.833333333,8.666666667,2.666666667,1.166666667,7,8,11.33333333,6,7.5))

H2_ICC <- data.frame(A=c(14.66666667,2.166666667,15.33333333,7.5,9,14,8,11,10.66666667,3,12.83333333,0.666666667,2,1.5,6.333333333,4,1.5,1.666666667,3,4,14,2,8),
                     B=c(11.33333333,1,17,9.333333333,9,18.33333333,8.5,9.333333333,7.333333333,3.333333333,12,1.5,2.666666667,1,6.333333333,6.666666667,2.666666667,1.666666667,4,4,13.33333333,2,7))

H3_ICC <- data.frame(A=c(23.33333333,4.833333333,32.66666667,16.5,19,25,16,24.5,7.333333333,9,24.66666667,3.333333333,4,3.5,12.66666667,12,3,2.833333333,10,12,24,8,15),
                     B=c(18.66666667,2,32.5,18.66666667,16.5,32.16666667,14,20.66666667,12.66666667,11.66666667,28,1.5,2.666666667,3.5,12.16666667,15.33333333,5.333333333,2.833333333,3,12,24.66666667,8,14.5))

V1_ICC <- data.frame(A=c(20,3,8,6,10,6,0,12,14,0,6,4,20,7,8,4,3,0,0,0,4,8,2),
                     B=c(16,3,12,4,6,7,9,4,4,6,4,9,4,9,2,8,2,2,0,0,0,4,6))

V2_ICC <- data.frame(A=c(4,2,12,3,2,0,3,6,8,0,3,4,4,0,4,0,1,0,0,0,0,0,0),
                     B=c(0,0,6,4,3,3,0,16,8,0,4,6,4,0,8,0,0,0,0,0,0,0,0))

V3_ICC <- data.frame(A=c(4,3,8,3,1,6,0,1,6,2,4,0,8,1,7,0,0,0,0,0,0,0,0),
                     B=c(4,3,3,8,3,8,3,8,4,3,4,6,0,3,15,0,0,0,0,0,0,0,0))


#Imaging error correlates
 
cor.test(flange_a$imaging_error,flange_a$focal_length)
cor.test(flange_a$imaging_error,flange_a$innereye_to_innereye_length)
cor.test(flange_a$imaging_error,flange_a$innereye_to_innereye_length)
cor.test(flange_b$raw_FA,flange_a$raw_FA)

cor.test(flange_a$innereye_to_innereye_length,flange_a$composite_FA)
flange_a$focal_length <- as.numeric(flange_a$focal_length)



## Intrclass correlation coefficient tests for FA

flange_repeat <- read.csv("flange_measurements_repeatability.csv", header = TRUE) # Import dataset into R
ICC_H1 <- read.csv("ICC_H1.csv", header = TRUE) # Import dataset into R

icc(H1_ICC, model = "two", type = "consistency", unit = "single")
icc(H2_ICC, model = "two", type = "consistency", unit = "single")
icc(H3_ICC, model = "two", type = "consistency", unit = "single")
icc(V1_ICC, model = "two", type = "consistency", unit = "single")
icc(V2_ICC, model = "two", type = "consistency", unit = "single")
icc(V3_ICC, model = "two", type = "consistency", unit = "single")



## T-tests to test for directional assymmetry

#Deviance from normal distribution with a mean of 0 indicates directional assymetry
H1_ttest <- t.test(flange_repeat$H1_signed, mu = 0)
H1_ttest

H2_ttest <- t.test(flange_repeat$H2_signed, mu = 0)
H2_ttest

H3_ttest <- t.test(flange_repeat$H3_signed, mu = 0)
H3_ttest

V1_ttest <- t.test(flange_repeat$V1_signed, mu = 0)
V1_ttest

V2_ttest <- t.test(flange_repeat$V2_signed, mu = 0)
V2_ttest

V3_ttest <- t.test(flange_repeat$V3_signed, mu = 0)
V3_ttest


## Shaprio tests for normality


#Deviances from normality indicates not fluctuating assyemtry
shapiro.test(flange_repeat$H1_signed)
shapiro.test(flange_repeat$H2_signed)
shapiro.test(flange_repeat$H3_signed)
shapiro.test(flange_repeat$V1_signed)
shapiro.test(flange_repeat$V2_signed)
shapiro.test(flange_repeat$V3_signed)


## ANOVA for fluctuating assymetry
 
#Provided the interaction between name and side is significant, with an F value of >1 we have a measurable effect

two.way_H1 <-aov(H1_left_relative_FA~replicate*side*name, data = flange)
summary(two.way_H1)
two.way_H2 <-aov(H2 ~replicate*H2_left_relative_FA*name, data = flange)
summary(two.way_H2_hq)
two.way_H3 <-aov(H3 ~replicate*H3_left_relative_FA*name, data = flange)
summary(two.way_H3)

#VFA - all data
two.way_V1 <-aov(V1~replicate*H1_left_relative_FA*name, data = originals)
summary(two.way_V1)
two.way_V2 <-aov(V2 ~replicate*H2_left_relative_FA*name, data = flange)
summary(two.way_H2_hq)
two.way_V3 <-aov(V3 ~replicate*H3_left_relative_FA*name, data = flange)
summary(two.way_H3)

#HFA - hq data
two.way_length <-aov(H1~replicate*H1_left_relative_FA*name, data = flange)
summary(two.way_length)
two.way_length <-aov(H2~replicate*H2_left_relative_FA*name, data = flange)
summary(two.way_length)
two.way_length <-aov(H3~replicate*H3_left_relative_FA*name, data = flange_hq)
summary(two.way_length)

#HFA - repeatability data
two.way_H1_repeat <-aov(H1~replicate*H1_left_relative_FA*name, data = flange)
summary(two.way_H1_repeat)

two.way_H2_repeat <-aov(H2~filename*H2_left_relative_FA*photo, data = flange_repeat)
summary(two.way_H2_repeat)

summary(flange_repeat)


#2 way field site 
two.way_fieldsite <-aov(raw_FA~field_site*name, data = flange_hq)
summary(two.way_fieldsite)

#2 way flange measurements
two.way_flange_nose_repeat <-aov(noseflange_to_noseflange~replicate*left_relative_nose_flange*name, data = flange)
two.way_flange_repeat <-aov(eyeflange_to_eyeflange~replicate*left_relative_eyeflange_inner*name, data = flange)

summary(two.way_flange_nose_repeat)


df1 <- 26 # degrees of freedom for I*S (name * H1_left_relative_FA)
df2 <- 224 # degrees of freedom for Residuals
F_value <- 4.79

p_value <- 1 - pf(F_value, df1, df2)
p_value


H1_by_site <- aov(H1~field_site, data = originals)
summary(H1_by_site)

H1_mixed <- lmer(H1 ~ field_site + (1|name), data=flange)
anova(H1_mixed)



##Flange FA

flange_measures <- subset(flange_repeat, flange_subset=="y")

head(flange_measures)

flange_measures$unsignedFA__topright_to_bottomleft_x <- as.numeric(flange_measures$unsignedFA__topright_to_bottomleft_x)
flange_measures$unsignedFA_topleft_to_bottomright <- as.numeric(flange_measures$unsignedFA_topleft_to_bottomright)
flange_measures$unsignedFA_eyeflange_to_eyeflange <- as.numeric(flange_measures$unsignedFA_eyeflange_to_eyeflange)
flange_measures$unsignedFA_noseflange_to_nose_flange <- as.numeric(flange_measures$unsignedFA_noseflange_to_nose_flange)

flange_measures$signedFA__topright_to_bottomleft_x <- as.numeric(flange_measures$signedFA__topright_to_bottomleft_x)
flange_measures$signedFA_topleft_to_bottomright <- as.numeric(flange_measures$signedFA_topleft_to_bottomright)
flange_measures$signedFA_eyeflange_to_eyeflange <- as.numeric(flange_measures$signedFA_eyeflange_to_eyeflange)
flange_measures$signedFA_noseflange_to_nose_flange <- as.numeric(flange_measures$signedFA_noseflange_to_nose_flange)

flange_measures_a <- subset(flange_measures, photo=="a")
flange_measures_b <- subset(flange_measures, photo=="b")

cor.test(flange_measures$composite_FA, flange_measures$imaging_error)
cor.test(flange_measures$imaging_error, flange_measures$resolution)


cor.test(flange_measures_a$unsignedFA_topleft_to_bottomright, flange_measures_b$unsignedFA_topleft_to_bottomright)
cor.test(flange_measures_a$unsignedFA__topright_to_bottomleft_x, flange_measures_b$unsignedFA__topright_to_bottomleft_x)
cor.test(flange_measures_a$unsignedFA_eyeflange_to_eyeflange, flange_measures_b$unsignedFA_eyeflange_to_eyeflange)
cor.test(flange_measures_a$unsignedFA_noseflange_to_nose_flange, flange_measures_b$unsignedFA_noseflange_to_nose_flange)


length(flange_measures_a$unsignedFA__topright_to_bottomleft_x)
length(flange_measures_b$unsignedFA__topright_to_bottomleft_x)
summary(flange_measures)



## T-test for flange measurements

H1_ttest <- t.test(flange_measures$signedFA_noseflange_to_noseflange, mu = 0)
H1_ttest


##ICC for flange measurements


D4_ICC <- data.frame(A=c(522,106,186,92,53.83333333,222,265.5,244.3333333,81,168,94.83333333,127.5,47.33333333,196),
                     B=c(471.3333333,89.5,277,194,144,273.3333333,193.3333333,334.6666667,170.3333333,163,156,325.8333333,48,116.6666667))

D5_ICC <- data.frame(A=c(226,25.5,72,53,44.83333333,111,99.5,196.3333333,43,160,75.83333333,33,5.333333333,82),
                     B=c(155.3333333,85,109,84,76.5,139.3333333,67.33333333,94.66666667,102.3333333,54,94,113.8333333,13.5,66.66666667))

D6_ICC <- data.frame(A=c(172,33.5,120,46,28.16666667,135,90.5,184.6666667,5,144,41.66666667,42.5,7.666666667,56),
                     B=c(132.6666667,57.5,104.5,62,13.5,48.16666667,72.66666667,97.33333333,73.66666667,49,1,42.16666667,16,58.33333333))

D7_ICC <- data.frame(A=c(27.33333333,13.66666667,151.3333333,6,17,5,19,22.33333333,18,12,1,11.33333333,0.166666667,6),
                     B=c(113.3333333,0.5,218,22.66666667,19.5,40.33333333,26.66666667,4.333333333,7.666666667,17.33333333,0.5,98.16666667,8.666666667,11))

icc(D4_ICC, model = "two", type = "consistency", unit = "single")
icc(D5_ICC, model = "two", type = "consistency", unit = "single")
icc(D6_ICC, model = "two", type = "consistency", unit = "single")
icc(D7_ICC, model = "two", type = "consistency", unit = "single")


##Shaprio tests

shapiro.test(flange_measures$signedFA_topleft_to_bottomright)
shapiro.test(flange_measures$signedFA__topright_to_bottomleft_x)
shapiro.test(flange_measures$signedFA_eyeflange_to_eyeflange)
shapiro.test(flange_measures$signedFA_noseflange_to_nose_flange)


## T tests
D4_ttest <- t.test(flange_measures$signedFA_topleft_to_bottomright, mu = 0)
D4_ttest

D5_ttest <- t.test(flange_measures$signedFA__topright_to_bottomleft_x, mu = 0)
D5_ttest

D6_ttest <- t.test(flange_measures$signedFA_eyeflange_to_eyeflange, mu = 0)
D6_ttest

D7_ttest <- t.test(flange_measures$signedFA_noseflange_to_nose_flange, mu = 0)
D7_ttest


## FA correlations with other FA metrics

flange_measures$unsignedFA_topleft_to_bottomright <- as.numeric(flange_measures$unsignedFA_topleft_to_bottomright)
flange_measures$unsignedFA__topright_to_bottomleft_x <- as.numeric(flange_measures$unsignedFA__topright_to_bottomleft_x)
flange_measures$unsignedFA_eyeflange_to_eyeflange <- as.numeric(flange_measures$unsignedFA_eyeflange_to_eyeflange)
flange_measures$unsignedFA_noseflange_to_nose_flange <- as.numeric(flange_measures$unsignedFA_noseflange_to_nose_flange)

cor.test(flange_measures$unsignedFA_eyeflange_to_eyeflange, flange_measures$unsignedFA_noseflange_to_nose_flange)



cleaned_flange_data <- subset(flange_data, !is.na(relative_eyeflange_to_eyeflange))
repeated_measures_ANOVA <- lme(relative_eyeflange_to_eyeflange ~ Year, random = ~1 | Name, data = cleaned_flange_data)
summary(repeated_measures_ANOVA)
 

# Reshape the data to wide format
# Run the Friedman test
friedman_result <- friedman.test(relative_eyeflange_to_eyeflange ~ Year | Name, data = cleaned_flange_data)

# Check the results
print(friedman_result)

