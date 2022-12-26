## ---------------------------
##
## Script name: Mussel shell temperature for all infaunal invertebrate species experiments
##
## Purpose of script: Does mussel shell temperature differ between infested and non-infested 
##                    mussel beds, and does this effect varies with time and between dates?
##
## Author: Dr.Guy F. Sutton and Alexia Dievart
##
## Date Created: 2022-11-09
## Dates Updated: 2022-11-14; 2022-11-21; 2022-11-26
##
## Copyright (c) Guy Sutton 2022
## Email: g.sutton@ru.ac.za
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------

###############################################################################################
# Section: Session setup ----------------------------------------------------------------------
###############################################################################################

# Set working directory
setwd("D:/Little Bull/Etudes/phD - Rhodes/1.2_HOT LIMPETS/STATS/alexia_mussels")

# Load packages
library(pacman)
pacman::p_load(tidyverse,
               DHARMa,
               lme4,
               mgcv,
               gratia,
               ggplot2,
               ggtext,
               tidyr,
               dplyr,
               MuMIn,
               glue,
               glmmTMB,
               ggeffects)

# Set default ggplot theme
theme_set(theme_classic() +
            theme(panel.border = element_rect(colour = "black", fill = NA),
                  axis.text = element_text(colour = "black"),
                  axis.text.x = element_text(angle = 0, hjust = 0),
                  axis.title.x = element_text(margin = unit(c(2, 0, 0, 0), "mm")),
                  axis.title.y = element_text(margin = unit(c(0, 4, 0, 0), "mm")),
                  legend.position = "none"))



###############################################################################################
# Section: Load data --------------------------------------------------------------------------
###############################################################################################

# Load full dataset
data <- readr::read_csv2("./infrared_data.csv")
head(data)
dplyr::glimpse(data)

# Select only rows where species == "Perna perna"
mussels <- data[data$species == "Perna perna",]

# Create a unique identifier ID for each replicate, because of repeated measurements 
mussels <- mussels %>%
  dplyr::group_by(date, replicate, specimen) %>%
  tidyr::unite(col = "id", date, replicate, specimen, sep = "_", remove = FALSE) %>%
  dplyr::mutate(id = as.factor(id)) %>% 
  dplyr::ungroup() %>%
  dplyr::mutate(infestation = dplyr::if_else(infestation == "infested", 
                                             "Infested",
                                             "Non-infested"))
dplyr::glimpse(mussels)

# Convert some variables into other types 
mussels$time <- as.numeric(mussels$time)
mussels$shell_temp <- as.numeric(mussels$shell_temp)
mussels$specimen <- as.factor(mussels$specimen)
mussels$infestation <- as.factor(mussels$infestation)
mussels$replicate <- as.factor(mussels$replicate)
mussels$date <- as.factor(mussels$date)
dplyr::glimpse(mussels)



###############################################################################################
# Section: Visualize data ---------------------------------------------------------------------
###############################################################################################

# Summarise mussel shell temperatures by date, time and infestation status
plot_mussels <- mussels %>%
  dplyr::group_by(date, time, infestation) %>%
  dplyr::summarise(
    shell.mean = mean(shell_temp),
    shell.sd   = sd(shell_temp),
    n = n(),
    shell.se   = shell.sd / sqrt(n)
  )
head(plot_mussels)

# Plot for every single experimental date (n = 6) 
my_colors <- c("darkgrey", "brown3")
my_colors1 <- c("grey86", "rosybrown2")
labels <- c(inf = "Infested", noninf = "Non-infested")
breaks <- c("inf", "noninf")

plot_mussels1 <- plot_mussels %>%
  ggplot(data = .,
         aes(
           x = time,
           y = shell.mean,
           group = infestation, 
           color = infestation
         )) +
  geom_ribbon(
    aes(ymin = shell.mean - shell.se, 
        ymax = shell.mean + shell.se,
        fill = infestation),
    alpha = 0.5,
    linetype = 0) +
  geom_errorbar(
      aes(ymin = shell.mean - shell.se, 
          ymax = shell.mean + shell.se),
     width = 5,
     position=position_dodge(.05)
  ) +
  geom_line(lwd = 1, linetype = "solid") +
  geom_point(size = 2) +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors1) +
  scale_x_continuous(
    breaks = seq(0, 90, by = 30),
    limits = c(-5, 95)
  ) +
  scale_y_continuous(
    breaks = seq(20, 50, by = 5),
    limits = c(19, 50)
  ) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = c(0.92, 0.1), 
        legend.background = element_rect(size = 0.5, linetype = "solid", colour = "black"),
        legend.title = element_text(face = "bold", size = 22), legend.title.align = 0.5,
        legend.text = element_text(size=20)) +
  theme(axis.text = element_text(size=20),
        axis.text.x = element_text(hjust = 0.5),
        axis.title = element_text(size=22),
        strip.text.x = element_text(size = 20)) +
  # Changing the legend title
  labs(colour = "Infestation",
       fill = "Infestation",
       x = "\nTime (min)",
       y = "Shell temperature (°C)\n") +
  facet_wrap(~ date, nrow = 1)
plot_mussels1

# Create a data frame with labels for each facet (i.e., experimental date)
plot_labels <- data.frame( 
  date = c("2021-04-13", "2021-04-14", "2021-04-15", "2021-04-18", "2022-02-09", "2022-02-10"),
  label = c("UV index = 7<br>Wind speed = 7.88-8.69 km/h",
            "UV index = 6.5<br>Wind speed = 13.68 km/h",
            "UV index = 6<br>Wind speed = 16.90-19.79 km/h",
            "UV index = 5.5<br>Wind speed = 19.47-21.24 km/h",
            "UV index = 13<br>Wind speed = 19.47-21.24 km/h",
            "UV index = 12<br>Wind speed = 14.00-15-45 km/h"),
  infestation = "Infested")
View(plot_labels)

# Plot for each experimental date with the corresponding environmental variables
plot_mussels2 <- plot_mussels1 +
  geom_richtext(x = 45, y = 50, aes(label = label), data = plot_labels, colour = "black",
            size = 5, fill = NA, label.colour = NA)
plot_mussels2



###############################################################################################
# Section 1: ROBOMUSSELS without DATE EFFECT --------------------------------------------------
###############################################################################################

###############################################################################################
# Section 1.1: Model fitting without DATE EFFECT ==============================================
###############################################################################################

###############################################################################################
# Section 1.1.1: Model #1 - Gaussian GLM ######################################################
###############################################################################################

# Here, we model shell temperature as a function of infestation status (infested/non-infested),
# and time since the start of the experiment:
# - Does not account for repeated measurements of the same mussel at different time intervals 
mod1 <- glmmTMB::glmmTMB(
  # Response variable 
  shell_temp ~
    # Fixed effects 
    infestation + time, 
  data = mussels, 
  family = gaussian(link = "identity")
)

# Check model fit 
DHARMa::simulateResiduals(fittedModel = mod1, plot = T)
# Deviation not significant (KS test, p = 0.56) but quantile deviations detected.

# Not a good fit. Unsurprisingly, as we have pseudoreplication.



###############################################################################################
# Section 1.1.2: Model #2 - Gaussian LMM ######################################################
###############################################################################################

# Here, we model shell temperature as a function of infestation status (infested/non-infested),
# and time since the start of the experiment, specifying a random intercept term for each 
# individual mussel to account for repeated measurements over time. 
mod2 <- glmmTMB::glmmTMB(
  # Response variable 
  shell_temp ~
    # Fixed effects 
    infestation + time +
    # Random effects
    (1| id), 
  data = mussels, 
  family = gaussian(link = "identity")
)

# Check model fit 
DHARMa::simulateResiduals(fittedModel = mod2, plot = T)
# Deviation significant (KS test, p < 0.05) and quantile deviations detected.

# A bit worse.



###############################################################################################
# Section 1.1.3: Model #3 - Gaussian LMM ######################################################
###############################################################################################

# Here, we model shell temperature as a function of infestation status (infested/non-infested),
# and time since the start of the experiment, specifying a random slope term for each 
# individual mussel to account for repeated measurements over time. 
# - Differs from model #2 in that we allow the infestation effect to vary over time. 
mod3 <- glmmTMB::glmmTMB(
  # Response variable 
  shell_temp ~
    # Fixed effects 
    infestation * time +
    # Random effects
    (1 | id), 
  data = mussels, 
  family = gaussian(link = "identity")
)

# Check model fit 
DHARMa::simulateResiduals(fittedModel = mod3, plot = T)
# Deviation significant (KS test, p < 0.01) and quantile deviations significant.

# Check residuals vs each fixed effect to diagnose issue with model fit 
simulationOutput <- DHARMa::simulateResiduals(fittedModel = mod3, plot = F) 
plotResiduals(simulationOutput, mussels$infestation)  # within-group deviations significant
plotResiduals(simulationOutput, mussels$time) # within-group deviations significant

# The problem appears to be with the variance between infestation levels and over time. 
# We might need a model that specifically deals with temporal autocorrelation. 



###############################################################################################
# Section 1.1.4: Model #4 - Gaussian LMM - dispersion term for 'time' #########################
###############################################################################################

# Here, we model shell temperature as a function of infestation status (infested/non-infested),
# and time since the start of the experiment, specifying a random slope term for each 
# individual mussel to account for repeated measurements over time. 
# - Differs from model #2 in that we allow the infestation effect to vary over time. 
mod4 <- glmmTMB::glmmTMB(
  # Response variable 
  shell_temp ~
    # Fixed effects 
    infestation * time +
    # Random effects
    (1 | id), 
  dispformula = ~ 1 + time, 
  data = mussels, 
  family = gaussian(link = "identity")
)

# Check model fit 
DHARMa::simulateResiduals(fittedModel = mod4, plot = T)
# Deviation not significant (KS test, p = 0.12) but quantile deviations detected.

# Check residuals vs each fixed effect to diagnose issue with model fit 
simulationOutput <- DHARMa::simulateResiduals(fittedModel = mod4, plot = F)
plotResiduals(simulationOutput, mussels$infestation) # within-group deviation close to significant
plotResiduals(simulationOutput, mussels$time) # within-group deviation significant

# The problem appears to be with the variance over time.



###############################################################################################
# Section 1.1.5: Model #5 - Generalised additive model ########################################
###############################################################################################

# Fit GAM model #5.1
# - Linear effects of time and infestation 
m1_gam <- gam(shell_temp ~ 
                # Fixed effects (linear terms)
                infestation + time + 
                # Random effect for id (smoothed - penalised)
                s(id, bs = 're'),
              data = mussels, 
              method = 'REML')
summary(m1_gam)

# Fit GAM model #5.2
# - Linear effects of infestation and non-linear effect for time 
m2_gam <- gam(shell_temp ~ 
                # Fixed effects (linear terms)
                infestation + 
                # Non-linear terms - cubic regression spline with 5 knots
                # - Allows the pieces between time intervals to hinge 
                s(time, bs = "cr", k = 4) + 
                # Random effect for id (smoothed - penalised)
                s(id, bs = 're'),
              data = mussels, 
              method = 'REML')
summary(m2_gam)

# Fit GAM model #5.3
# - Interaction term: Allows the non-linear effect of time to vary between the different
#   levels of infestation, with a simple random intercept smoother  
m3_gam <- gam(shell_temp ~ 
                # Fixed effects (linear terms)
                infestation + 
                # Non-linear terms - cubic regression spline with 5 knots
                # - Allows the pieces between time intervals to hinge 
                s(time, bs = "cr", by = infestation, k = 4) + 
                # Random effect for id (smoothed - penalised)
                s(id, bs = 're'),
              data = mussels, 
              method = 'REML')
summary(m3_gam)

# Fit GAM model #5.4
# - Interaction term: Allows the non-linear effect of time to vary between the different
#   levels of infestation
#   - Fit a nested random intercept term for each mussel (id) nested within the mussel_bed
# m4_gam <- gam(shell_temp ~ 
                # Fixed effects (linear terms)
#               infestation + 
                # Non-linear terms - cubic regression spline with 5 knots
                # - Allows the pieces between time intervals to hinge 
#               s(time, bs = "cr", by = infestation, k = 4) + 
                # Random effect for id (smoothed - penalised)
#               s(id, bs = 're', by = mussel_bed),
#               data = mussels, 
#               method = 'REML')
# summary(m4_gam) 
# Error in gam(shell_temp ~ infestation + s(time, bs = "cr", by = infestation,  : 
# Model has more coefficients than data

## GAM model #5.4 captures the structure of the experimental design and tests the hypothesis
## that we are interested in. Because all shell temperatures of mussels were pooled together, 
## there is more coefficients (id) than observations. The GAM model is not capable to run the
## analysis in this context.


# Check model fit for m3_gam
# - We want to know whether the model captured the data structure properly
gratia::appraise(m3_gam, method = "simulate")
# This model captured the data structure properly
## This model does not account for the variance between mussels and between mussel beds.

# Compare the models based on AICc
AICc(m1_gam)
AICc(m2_gam)
AICc(m3_gam) # Best model (with k = 4)

# Information to report in the Supplementary Materials
summary(m3_gam) # for estimate and std error
anova(m3_gam, test = 'F') # for F-value and p-value
# Warning message:
# In anova.gam(m3_gam, test = "F") : test argument ignored


###############################################################################################
# Section 1.2: Model inference without DATE EFFECT =============================================
###############################################################################################

#################################################
# Test for treatment effect of infestation status 
#################################################

# Perform Wald-like test (Wood 2013a,b)
anova(m3_gam, test = "Chisq") 
# Warning message:
# In anova.gam(m3_gam, test = "Chisq") : test argument ignored

# There is evidence for an effect of infestation status on shell temperature 
# (F = 496.6, p < 0.001).



#################################################
# Test for treatment effect of infestation status to vary over time 
#################################################

# Fit null model without the time-varying parameter for infestation effect 
null_gam <- gam(shell_temp ~ 
                  # Fixed effects (linear terms)
                  infestation + 
                  # Non-linear terms - cubic regression spline with 5 knots
                  # - Allows the pieces between time intervals to hinge 
                  s(time, bs = "cr", k = 4) + 
                  # Random effect for id (smoothed - penalised)
                  s(id, bs = 're'),
                data = mussels, 
                method = 'REML')
summary(null_gam)

# Perform Wald's test (Wood 2013a,b)
anova(null_gam, m3_gam, test = "Chisq") 
# There is evidence for the effect of infestation to vary with time 
# (Chisq = 226.28, p < 0.001).



###############################################################################################
# Section 1.3: Plot model predictions without DATE EFFECT ======================================
###############################################################################################

# Create a set of data to predict over (using m3_gam)
new_data <- tidyr::expand(mussels, nesting(id, infestation),
                           time = unique(time))
head(new_data)

# Extract predictions from the model 
best_mod_pred <- bind_cols(new_data,
                           as.data.frame(predict(m3_gam, newdata = new_data,
                                                 se.fit = TRUE))) 
head(best_mod_pred)

# Extract marginal effects (average effect)
preds <- ggeffects::ggemmeans(
  m3_gam,
  terms = c("infestation", "time [0:90 by = 1"),
  type = "fe",
  interval = "confidence"
) %>%
  as.data.frame() %>%
  dplyr::mutate(group = readr::parse_integer(as.character(group))) %>%
  dplyr::rename(
    infestation = x,
    time = group
  )
head(preds)

# Make plot (all experimental dates pooled)
# - Each dashed line represents an individual mussel
# - The thick line represents the marginal (average/mean) effect
# - This conveys the variation in shell temperatures over time 
ggplot(data = best_mod_pred) +
  geom_line(aes(
    x = time,
    y = fit,
    group = id,
    colour = infestation
  ),
  linetype = "dashed"
  ) +
  geom_line(data = preds,
            aes(
              x = time,
              y = predicted,
              colour = infestation
            ),
            size = 2
   ) + 
  facet_wrap(~ infestation) +
  geom_point(data = mussels, aes(x = time, y = shell_temp, colour = infestation),
             size = 0.75) +
  labs(
    x = "Exposure time (minutes)",
    y = "Shell temperature (degrees C)"
  ) + 
  theme_bw()

# Make plot (all experimental dates pooled)
# - Each dot line represents an individual mussel
# - The line represents the marginal (average/mean) effect
# - The ribbon represents the confidence interval predicted by the model,
# Here: narrow interval because the model does not account for the nested effect (random mussels w/in mussel bed)
# - This conveys the variation in shell temperatures over time  
ggplot(data = preds,
       aes(
         x = time,
         y = predicted,
         colour = infestation,
         fill = infestation
       )) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = 0.2) +
  # This overlays each individual mussels raw data as points
  # Change to geom_line to plot each mussels data as a curve
  geom_point(
    data = best_mod_pred,
    aes(
      x = time,
      y = fit,
      group = id,
      colour = infestation
    )) +
  facet_wrap( ~ infestation, nrow = 1) +
  theme(legend.position = "right"
  ) +
  labs(
    x = "Exposure time (minutes)",
    y = "Shell temperature (degrees C)",
    fill = "Infestation status"
  ) +
  guides(colour = "none")



###############################################################################################
# Section 2: ROBOMUSSELS with DATE EFFECT ------------------------------------------------------
###############################################################################################

###############################################################################################
# Section 2.1: Model fitting with DATE EFFECT =================================================
###############################################################################################

###############################################################################################
# Section 2.1.1: Model #5 - Generalised additive model ########################################
###############################################################################################

# Fit GAM #5.3 without the date effect
# - Interaction term: Allows the non-linear effect of time to vary between the different
#   levels of infestation, with a simple random intercept smoother  
m3_gam <- gam(shell_temp ~ 
                # Fixed effects (linear terms)
                infestation + 
                # Non-linear terms - cubic regression spline with 5 knots
                # - Allows the pieces between time intervals to hinge 
                s(time, bs = "cr", by = infestation, k = 4) + 
                # Random effect for id (smoothed - penalised)
                s(id, bs = 're'),
              data = mussels, 
              method = 'REML')
summary(m3_gam)
## This model does not account for nested effect (random mussels within mussel bed)

# Fit GAM #5.5 with the date effect 
m5_gam <- gam(shell_temp ~ 
                # Fixed effects (linear terms)
                infestation + date + 
                # Non-linear terms - cubic regression spline with 5 knots
                # - Allows the pieces between time intervals to hinge 
                s(time, bs = "cr", by = infestation, k = 4) + 
                # Random effect for id (smoothed - penalised)
                s(id, bs = 're'),
              data = mussels, 
              method = 'REML')
summary(m5_gam)

#  Check model fits
gratia::appraise(m3_gam, method = "simulate")
gratia::appraise(m5_gam, method = "simulate")
## Models look okay. 

# Compare the models based on AICc
AICc(m3_gam)
AICc(m5_gam)  # Best model (with k = 4)
# This is a first quick check if including 'date' in the model improved fit.
# - There is an increase in model fit, indicated by the lower AICc when including 'date' in the model.



###############################################################################################
# Section 2.2: Model fitting with DATE EFFECT =================================================
###############################################################################################

# Formal hypothesis test to see if the 'date' effect is statistically significant 
# - We do this using a Wald's test (Wood 2013) which compares the two models, 
#   without (m3_gam) and with (m5_gam) the 'date' as a fixed effect
anova(m3_gam, m5_gam, test = "F")

## There is an effect of the 'date' on the temperature of mussels (F = 1.93, p = 0.02 *) 



###############################################################################################
# Section 2.3: Model fitting with DATE EFFECT =================================================
###############################################################################################

# Create a set of data to predict over (using m5_gam)
new_data1 <- tidyr::expand(mussels, nesting(id, date, infestation),
                          time = unique(time))
head(new_data1)

# Extract predictions from the model 
best_mod_pred1 <- bind_cols(new_data1,
                           as.data.frame(predict(m5_gam, newdata = new_data1,
                                                 se.fit = TRUE))) 
View(best_mod_pred1)

# For a plot that reads better, check plot_mussels2
plot_mussels2

# Make plot (with facetting for each experimental date)
ggplot(best_mod_pred1,
       aes(
         x = time,
         y = fit,
         group = id,
         colour = infestation
       ))+
  geom_line() +
  facet_wrap(~ infestation + date) +
  geom_point(
    data = mussels, 
    aes(
      x = time, 
      y = shell_temp, 
      colour = infestation, 
      shape = date
    ),
    size = 0.75) +
  labs(
    x = "Exposure time (minutes)",
    y = "Shell temperature (degrees C)"
  )

# Make plot (without the facetting)
ggplot(best_mod_pred1,
       aes(
         x = time,
         y = fit,
         group = id,
         colour = date
       ))+
  geom_line() +
  facet_wrap(~ infestation) +
  labs(
    x = "Exposure time (minutes)",
    y = "Shell temperature (degrees C)"
  ) +
  theme(legend.position = "right")
