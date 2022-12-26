## ---------------------------
##
## Script name: Infaunal shell temperature for Oxystele antoni and Helcion pectunculus experiments
##
## Purpose of script: Do snail and limpet (dark) shell temperature differ between infested and non-infested 
##                    mussel beds, and does this effect varies with time and between days?
##
## Author: Dr.Guy F. Sutton and Alexia Dievart
##
## Date Created: 2022-11-14
## Dates Updated: 2022-11-15
##
## Copyright (c) Guy Sutton 2022
## Email: g.sutton@ru.ac.za
##
## ---------------------------
##
## Notes:
##   Section 1: Only on snail and limpet shell temperature
##   Section 2: Addition of mussel shell temperature for the corresponding dates,
##              with species as fixed parametric factor in models
##
## ---------------------------

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

# Customize colours for plots
my_colors <- c("darkgrey", "darkorange")

###############################################################################################
# Section: Load data --------------------------------------------------------------------------
###############################################################################################

# Load full dataset
data <- readr::read_csv2("./infrared_data.csv")
head(data)
View(data)
dplyr::glimpse(data)

# Select only rows with relevant dates : 2021-04-13 and 2021-04-15
infauna_mussels <- data[data$date == c("2021-04-13", "2021-04-15"),]
View(infauna_mussels)

# Create a unique identifier ID for each replicate, because of repeated measurements 
infauna_mussels <- infauna_mussels %>%
  dplyr::group_by(date, replicate, specimen, species) %>%
  tidyr::unite(col = "id", date, replicate, specimen, species, sep = "_", remove = FALSE) %>%
  dplyr::mutate(id = as.factor(id)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(infestation = dplyr::if_else(infestation == "infested", 
                                             "Infested",
                                             "Non-infested"))
dplyr::glimpse(infauna_mussels)
View(infauna_mussels)

# Convert some variables into other types 
infauna_mussels$time <- as.numeric(infauna_mussels$time)
infauna_mussels$shell_temp <- as.numeric(infauna_mussels$shell_temp)
factor <- c("specimen", "infestation", "replicate", "date", "species")
infauna_mussels[,factor] <- lapply(infauna_mussels[,factor], factor)
dplyr::glimpse(infauna_mussels)

# Select only rows with snail and limpet shell temperatures
infauna <- infauna_mussels[infauna_mussels$species != "Perna perna",]
infauna$species <- as.factor(infauna$species)
View(infauna)

# Select corresponding mussel shell temperature
mussels <- infauna_mussels[infauna_mussels$species == "Perna perna",]
mussels$date <- as.factor(mussels$date)
View(mussels)

###############################################################################################
# Section: Visualize data ---------------------------------------------------------------------
###############################################################################################

# Summarise limpet, snail AND mussel shell temperature by date, time and infestation class
plot_infauna_mussels <- infauna_mussels %>%
  dplyr::group_by(date, time, infestation, species) %>%
  dplyr::summarise(
    shell.mean = mean(shell_temp),
    shell.sd   = sd(shell_temp),
    n = n(),
    shell.se   = shell.sd / sqrt(n)
  )
head(plot_infauna_mussels)

# Plot of limpet AND mussel shell temperatures for both days
lines <- c("Perna perna" = "dotted", "Oxystele antoni" = "solid", "Helcion pectunculus" = "longdash")

plot_infauna2 <- plot_infauna_mussels %>%
  ggplot(data = .,
         aes(
           x = time,
           y = shell.mean, 
           color = infestation
         )) +
  geom_errorbar(
    aes(ymin = shell.mean - shell.se, 
        ymax = shell.mean + shell.se),
    width = 1,
    size = 1) +
  geom_line(aes(linetype = species, color = infestation), lwd=1) +
  scale_color_manual(values = my_colors) +
  scale_linetype_manual(values = lines, labels = expression(italic("Perna perna"), 
                                                            italic("Oxystele antoni"), 
                                                            italic("Helcion pectunculus"))) +
  scale_x_continuous(
    breaks = seq(0, 90, by = 15),
    limits = c(-1, 91)
  ) +
  scale_y_continuous(
    breaks = seq(20, 50, by = 5),
    limits = c(19, 50)
  ) +
  ylab("Shell temperature (°C)\n") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "right", legend.key.width = unit(1, "cm"),
        legend.title = element_text(face = "bold"), legend.text.align = 0) +
  # Changing the legend title
  labs(colour = "Infestation", linetype = "Species",
       x = "Time (min)") +
  facet_wrap(~ date, nrow = 1) +
  ggtitle("\n")
plot_infauna2

###############################################################################################
# Section 1: LIMPET + SNAIL -------------------------------------------------------------------
###############################################################################################

###############################################################################################
# Section 1.1: Model fitting without DATE EFFECT ==============================================
###############################################################################################

###############################################################################################
# Section 1.1.1: Model #1 - Gaussian GLM ######################################################
###############################################################################################

# Here, we model shell temperature as a function of infestation status (infested/non-infested),
# species (HELPEC/OXYANT), and time since the start of the experiment:
# - Does not account for repeated measurements of the same invertebrates at different time intervals 
mod1 <- glmmTMB::glmmTMB(
  # Response variable 
  shell_temp ~
    # Fixed effects 
    infestation + species + time, 
  data = infauna, 
  family = gaussian(link = "identity")
)

# Check model fit 
DHARMa::simulateResiduals(fittedModel = mod1, plot = T)
# Deviation not significant (KS test, p = 0.1) but quantile deviations detected.

# Not the best fit. Unsurprisingly, as we have pseudoreplication.

###############################################################################################
# Section 1.1.2: Model #2 - Gaussian LMM ######################################################
###############################################################################################

# Here, we model shell temperature as a function of infestation status (infested/non-infested),
# species (HELPEC/OXYANT) and time since the start of the experiment, 
# specifying a random intercept term for each individual invertebrate to account 
# for repeated measurements over time. 
mod2 <- glmmTMB::glmmTMB(
  # Response variable 
  shell_temp ~
    # Fixed effects 
    infestation + species + time +
    # Random effects
    (1| id), 
  data = infauna, 
  family = gaussian(link = "identity")
)

# Check model fit 
DHARMa::simulateResiduals(fittedModel = mod2, plot = T)
# Deviation not significant (KS test, p = 0.07) but quantile deviations detected.

# A bit better, but not much.

###############################################################################################
# Section 1.1.3: Model #3 - Gaussian LMM ######################################################
###############################################################################################

# Here, we model shell temperature as a function of infestation status (infested/non-infested),
# species (HELPEC/OXYANT) and time since the start of the experiment, specifying a random slope term 
# for each individual invertebrate to account for repeated measurements over time. 
# - Differs from model #2 in that we allow the infestation effect to vary over time. 
mod3 <- glmmTMB::glmmTMB(
  # Response variable 
  shell_temp ~
    # Fixed effects 
    infestation * time + species +
    # Random effects
    (1 | id), 
  data = infauna, 
  family = gaussian(link = "identity")
)

# Check model fit 
DHARMa::simulateResiduals(fittedModel = mod3, plot = T)
# Deviation not significant (KS test, p = 0.1) but quantile deviations significant.

# Check residuals vs each fixed effect to diagnose issue with model fit 
simulationOutput <- DHARMa::simulateResiduals(fittedModel = mod3, plot = F) 
plotResiduals(simulationOutput, infauna$infestation)  # within-group deviations not significant
plotResiduals(simulationOutput, infauna$time) # within-group deviations significant

# The problem appears to be with the variance over time. 
# We might need a model that specifically deals with temporal autocorrelation.

###############################################################################################
# Section 1.1.4: Model #4 - Gaussian LMM - dispersion term for 'time' #########################
###############################################################################################

# Here, we model shell temperature as a function of infestation status (infested/non-infested),
# species (HELPEC/OXYANT) and time since the start of the experiment, specifying a random slope 
# term for each individual invertebrate to account for repeated measurements over time. 
# - Differs from model #2 in that we allow the infestation effect to vary over time. 
mod4 <- glmmTMB::glmmTMB(
  # Response variable 
  shell_temp ~
    # Fixed effects 
    infestation * time + species +
    # Random effects
    (1 | id), 
  dispformula = ~ 1 + time, 
  data = infauna, 
  family = gaussian(link = "identity")
)

# Check model fit 
DHARMa::simulateResiduals(fittedModel = mod4, plot = T)
# Deviation significant (KS test, p < 0.05) and quantile deviations detected.

# Check residuals vs each fixed effect to diagnose issue with model fit 
simulationOutput <- DHARMa::simulateResiduals(fittedModel = mod4, plot = F)
plotResiduals(simulationOutput, infauna$infestation) # within-group deviation significant
plotResiduals(simulationOutput, infauna$time) # within-group deviation significant

# The problem appears to be with the variance over time.

###############################################################################################
# Section 1.1.5: Model #5 - Generalised additive model  #######################################
###############################################################################################

# Fit GAM #5.1
# - Linear effects of time and infestation 
m1_gam <- gam(shell_temp ~ 
                # Fixed effects (linear terms)
                infestation + species + time + 
                # Random effect for id (smoothed - penalised)
                s(id, bs = 're'),
              data = infauna, 
              method = 'REML')
summary(m1_gam)

# Fit GAM #5.2
# - Linear effects of infestation and non-linear effect for time 
m2_gam <- gam(shell_temp ~ 
                # Fixed effects (linear terms)
                infestation + species +
                # Fixed effect (non-linear terms) - cubic regression spline with 4 knots
                # - Allows the pieces between time intervals to hinge 
                s(time, bs = "cr", k = 4) + 
                # Random effect for id (smoothed - penalised)
                s(id, bs = 're'),
              data = infauna, 
              method = 'REML')
summary(m2_gam)

# Fit GAM #5.3
# - Interaction term: Allows the non-linear effect of time to vary between the different
#   levels of infestation, with a simple random intercept smoother  
m3_gam <- gam(shell_temp ~ 
                # Fixed effects (linear terms)
                infestation + species +
                # Fixed effect (non-linear terms) - cubic regression spline with 5 knots
                # - Allows the pieces between time intervals to hinge 
                s(time, bs = "cr", by = infestation, k = 4) + 
                # Random effect for id (smoothed - penalised)
                s(id, bs = 're'),
              data = infauna, 
              method = 'REML')
summary(m3_gam)

# Check model fit for m3_gam (with k = 4)
# - We want to know whether the model captured the data structure properly
gratia::appraise(m3_gam, method = "simulate")
# Not the best fit, but not too bad.

# Test whether allowing for time specific effects to vary between infestation levels improves model fit 
anova(m2_gam, m3_gam, test = "Chisq")
# There is NO statistical difference between m2_gam and m3_gam (Chisq = 4.13, p = 0.9)

# Compare the models based on AICc
AICc(m1_gam)
AICc(m2_gam) # Best model without varying time effect between infestation levels.
AICc(m3_gam) 



###############################################################################################
# Section 1.1: Model inference without DATE EFFECT ============================================
###############################################################################################

# I chose to continue with m3_gam (even if the AICc is higher than m2_gam)

#################################################
# Test for treatment effect of infestation status 
#################################################

# Perform Wald-like test (Wood 2013a,b)
anova(m3_gam, test = "F") 
# There is no evidence for an effect of infestation or invertebrate species on shell temperature
# (infestation: F = 2.183, p = 0.142; species: F = 0.059, p = 0.808). 

#################################################
# Test for treatment effect of infestation status to vary over time 
#################################################

# Fit null model without the time-varying parameter for infestation effect 
null_gam <- gam(shell_temp ~ 
                  # Fixed effects (linear terms)
                  infestation + species +
                  # Non-linear terms - cubic regression spline with 5 knots
                  # - Allows the pieces between time intervals to hinge 
                  s(time, bs = "cr", k = 4) + 
                  # Random effect for id (smoothed - penalised)
                  s(id, bs = 're'),
                data = infauna, 
                method = 'REML')
summary(null_gam)

# Perform Wald's test (Wood 2013a,b)
anova(null_gam, m3_gam, test = "F")
# There is NO evidence for the effect of infestation to vary with time 
# (F = 0.2365, p = 0.8725).



###############################################################################################
# Section 1.1: Plot model predictions without DATE EFFECT =====================================
###############################################################################################

# Create a set of data to predict over (using m3_gam)
new_data <- tidyr::expand(infauna, nesting(id, infestation, species),
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
  terms = c("infestation", "time [0:90 by = 1", "species"),
  type = "fe",
  interval = "confidence"
) %>%
  as.data.frame() %>%
  dplyr::mutate(group = readr::parse_integer(as.character(group))) %>%
  dplyr::rename(
    infestation = x,
    time = group,
    species = facet
  )
head(preds)

# Make plot (all experimental dates pooled)
# - Each dashed line represents an individual invertebrate at the given time point
# - The thick line represents the marginal (average/mean) effect
# - This conveys the variation in shell temperatures over time 
ggplot(data = best_mod_pred) +
  geom_line(aes(
    x = time,
    y = fit,
    group = id,
    colour = infestation,
    linetype = species
  )) +
  geom_line(data = preds,
            aes(
              x = time,
              y = predicted,
              colour = infestation
            ),
            size = 2
  ) + 
  scale_color_manual(values = my_colors) +
  facet_wrap(~ infestation + species) +
  #geom_point(data = mussels, aes(x = time, y = shell.temp, colour = infestation),
  #           size = 0.75) +
  labs(
    x = "Exposure time (minutes)",
    y = "Shell temperature (degrees C)"
  )

# Make plot (all experimental dates pooled)
# - Each dot line represents an individual invertebrate
# - The line represents the marginal (average/mean) effect
# - The ribbon represents the confidence interval predicted by the model,
# Here: narrow interval because the model does not account for the nested effect (random invertebrates w/in mussel bed)
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
              alpha = 0.2, fill = c("darkgrey", "darkgrey", "darkgrey", "darkgrey","darkgrey",
                                    "darkgrey", "darkgrey", "darkgrey", "darkorange", 
                                    "darkorange", "darkorange", "darkorange", "darkorange",
                                    "darkorange", "darkorange", "darkorange")) +
  scale_color_manual(values = my_colors) +
  # This overlays each individual invertebrate raw data as points
  # Change to geom_line to plot each invertebrate data as a curve
  geom_point(
    data = best_mod_pred,
    aes(
      x = time,
      y = fit,
      group = id,
      colour = infestation
    )) +
  facet_wrap( ~ infestation + species, nrow = 2) +
  theme(legend.position = "right"
  ) +
  labs(
    x = "Exposure time (minutes)",
    y = "Shell temperature (degrees C)",
    fill = "Infestation status"
  ) 



###############################################################################################
# Section 1.2: Model fitting with DATE EFFECT =================================================
###############################################################################################

# Fit model without the date effect
m3_gam <- gam(shell_temp ~ 
                # Fixed effects (linear terms)
                infestation  + species +
                # Non-linear terms - cubic regression spline with 4 knots
                # - Allows the pieces between time intervals to hinge 
                s(time, bs = "cr", by = infestation, k = 4) + 
                # Random effect for id (smoothed - penalised)
                s(id, bs = 're'),
              data = infauna, 
              method = 'REML')
summary(m3_gam)

# Fit model with the date effect 
m5_gam <- gam(shell_temp ~ 
                # Fixed effects (linear terms)
                infestation + species + date + 
                # Non-linear terms - cubic regression spline with 5 knots
                # - Allows the pieces between time intervals to hinge 
                s(time, bs = "cr", by = infestation, k = 4) + 
                # Random effect for id (smoothed - penalised)
                s(id, bs = 're'),
              data = infauna, 
              method = 'REML')
summary(m5_gam)

#  Check model fits
gratia::appraise(m3_gam, method = "simulate")
gratia::appraise(m5_gam, method = "simulate")
# Not the best fit, but not too bad.

# Compare the models based on AICc
AICc(m3_gam)
AICc(m5_gam) # Best model with date effect (k = 4)

# This is a first quick check if including 'date' in the model improved fit
# - There is an increase in model fit, indicated by the lower AICc
#   value for the model with the date effect, when including 'date'.

###############################################################################################
# Section 1.2: Model inference with DATE EFFECT ================================================
###############################################################################################

# Formal hypothesis test to see if the 'date' effect is statistically significant 
# - We do this using a Wald's test (Wood 2013) which compares the two models, 
#   without (m3_gam) and with (m5_gam), the 'date' fixed effect
anova(m3_gam, m5_gam, test = "F")
# Similar to the AICc analysis above, there is a significant effect of date on
# invertebrate shell temperature (F = 2.64, p < 0.05 *)



###############################################################################################
# Section 1.2: Plot model predictions with DATE EFFECT =========================================
###############################################################################################

# Create a set of data to predict over (using m4_gam)
new_data1 <- tidyr::expand(infauna, nesting(id, infestation, species, date),
                           time = unique(time))
View(new_data1)


# Extract predictions from the model 
best_mod_pred1 <- bind_cols(new_data1,
                            as.data.frame(predict(m4_gam, newdata = new_data1,
                                                  se.fit = TRUE))) 
head(best_mod_pred1)

# Make plot
ggplot(best_mod_pred1,
       aes(
         x = time,
         y = fit,
         group = id,
         colour = infestation
         
       ))+
  geom_line() +
  facet_wrap(~ infestation + date + species, nrow = 2) +
  scale_color_manual(values = my_colors) +
  geom_point(
    data = infauna, 
    aes(
      x = time, 
      y = shell_temp, 
      colour = infestation
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
  facet_wrap(~ infestation + species) +
  labs(
    x = "Exposure time (minutes)",
    y = "Shell temperature (degrees C)"
  ) +
  theme(legend.position = "right")



###############################################################################################
# Section 2: LIMPET and SNAIL vs ROBOMUSSELS --------------------------------------------------
###############################################################################################

###############################################################################################
# Section 2.1: Model fitting without DATE EFFECT ==============================================
###############################################################################################

###############################################################################################
# Section 2.1.1: Model #5 - Generalised additive model  #######################################
###############################################################################################

# Fit GAM #5.1
# - Linear effects of time and infestation 
m1_gam <- gam(shell_temp ~ 
                # Fixed effects (linear terms)
                infestation + time + species +
                # Random effect for id (smoothed - penalised)
                s(id, bs = 're'),
              data = infauna_mussels, 
              method = 'REML')
summary(m1_gam)

# Fit GAM #5.2
# - Linear effects of infestation and non-linear effect for time 
m2_gam <- gam(shell_temp ~ 
                # Fixed effects (linear terms)
                infestation + species +
                # Fixed effect (non-linear terms) - cubic regression spline with 4 knots
                # - Allows the pieces between time intervals to hinge 
                s(time, bs = "cr", k = 4) + 
                # Random effect for id (smoothed - penalised)
                s(id, bs = 're'),
              data = infauna_mussels, 
              method = 'REML')
summary(m2_gam)

# Fit GAM #5.3
# - Interaction term: Allows the non-linear effect of time to vary between the different
#   levels of infestation, with a simple random intercept smoother  
m3_gam <- gam(shell_temp ~ 
                # Fixed effects (linear terms)
                infestation + species +
                # Fixed effect (non-linear terms) - cubic regression spline with 4 knots
                # - Allows the pieces between time intervals to hinge 
                s(time, bs = "cr", by = infestation, k = 4) + 
                # Random effect for id (smoothed - penalised)
                s(id, bs = 're'),
              data = infauna_mussels, 
              method = 'REML')
summary(m3_gam)

# Check model fit for m3_gam (with k = 4)
# - We want to know whether the model captured the data structure properly
gratia::appraise(m3_gam, method = "simulate")
# Model fit is not bad.

# Test whether allowing for time specific effects to vary between infestation levels improves model fit 
anova(m2_gam, m3_gam, test = "Chisq")
# There is NO statistical difference between m2 and m3_gam (Chisq = 2.24, p = 0.9)

# Compare the models based on AICc
AICc(m1_gam)
AICc(m2_gam) # Best model without the variation of time with infestation levels.
AICc(m3_gam)



##############################################################################################
# Section 2.1: Model inference without DATE EFFECT =============================================
###############################################################################################

# I chose to continue with m3_gam

#################################################
# Test for treatment effect of infestation status 
#################################################

# Perform Wald-like test (Wood 2013a,b)
anova(m3_gam, test = "F") 
# There is NO evidence for an effect of infestation status and invertebrate species on shell temperature 
# (infestation: F = 2.182, p = 0.141; species: F = 0.604, p =0.547). 

#################################################
# Test for treatment effect of infestation status to vary over time 
#################################################

# Fit null model without the time-varying parameter for infestation effect 
null_gam <- gam(shell_temp ~ 
                  # Fixed effects (linear terms)
                  infestation + species + 
                  # Non-linear terms - cubic regression spline with 4 knots
                  # - Allows the pieces between time intervals to hinge 
                  s(time, bs = "cr", k = 4) + 
                  # Random effect for id (smoothed - penalised)
                  s(id, bs = 're'),
                data = infauna_mussels, 
                method = 'REML')
summary(null_gam)

# Perform Wald's test (Wood 2013a,b)
anova(null_gam, m3_gam, test = "Chisq")
# There is NO evidence for the effect of infestation to vary with time 
# (F = 0.1189, p = 0.9498). 



###############################################################################################
# Section 2.1: Plot model predictions without DATE EFFECT ======================================
###############################################################################################

# Create a set of data to predict over (using m3_gam)
new_data <- tidyr::expand(infauna_mussels, nesting(id, infestation, species),
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
  terms = c("infestation", "time [0:90 by = 1", "species"),
  type = "fe",
  interval = "confidence"
) %>%
  as.data.frame() %>%
  dplyr::mutate(group = readr::parse_integer(as.character(group))) %>%
  dplyr::rename(
    infestation = x,
    time = group,
    species = facet
  )
head(preds)

# Make plot (all experimental dates pooled)
# - Each dashed line represents an individual invertebrate or mussel at the given time point
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
  scale_color_manual(values = my_colors) +
  facet_wrap(~ infestation + species) +
  #geom_point(data = mussels, aes(x = time, y = shell.temp, colour = infestation),
  #           size = 0.75) +
  labs(
    x = "Exposure time (minutes)",
    y = "Shell temperature (degrees C)"
  )

# Nice plot showing that mussel, snail and limpet shell temperatures was never different on both experimental dates.
ggplot(data = preds,
       aes(
         x = time,
         y = predicted,
         colour = infestation,
         fill = infestation
       )) +
  geom_line(aes(group = species, linetype = species, colour = infestation)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, group = species, colour = infestation),
              alpha = 0.2, fill = c("darkgrey", "darkgrey", "darkgrey", "darkgrey", "darkgrey", "darkgrey",
                                    "darkgrey", "darkgrey", "darkgrey", "darkgrey", "darkgrey", "darkgrey",
                                    "darkorange", "darkorange", "darkorange", "darkorange",
                                    "darkorange", "darkorange", "darkorange", "darkorange",
                                    "darkorange", "darkorange", "darkorange", "darkorange")) +
  scale_color_manual(values = my_colors) +
  # This overlays each individual invertebrate raw data as points
  # Change to geom_line to plot each invertebrate data as a curve
   geom_point(
    data = best_mod_pred,
    aes(
      x = time,
      y = fit,
      group = id,
      colour = infestation,
      shape = species
    )) +
  facet_wrap( ~ infestation, nrow = 1) +
  theme(legend.position = "right"
  ) +
  labs(
    x = "Exposure time (minutes)",
    y = "Shell temperature (degrees C)",
    fill = "Infestation status"
  ) 



###############################################################################################
# Section 2.2: Model fitting with DATE EFFECT ==================================================
###############################################################################################

# Fit model without the date effect
m3_gam <- gam(shell_temp ~ 
                # Fixed effects (linear terms)
                infestation  + species +
                # Non-linear terms - cubic regression spline with 4 knots
                # - Allows the pieces between time intervals to hinge 
                s(time, bs = "cr", by = infestation, k = 4) + 
                # Random effect for id (smoothed - penalised)
                s(id, bs = 're'),
              data = infauna_mussels, 
              method = 'REML')
summary(m3_gam)

# Fit model with the date effect 
m5_gam <- gam(shell_temp ~ 
                # Fixed effects (linear terms)
                infestation + date + species +
                # Non-linear terms - cubic regression spline with 4 knots
                # - Allows the pieces between time intervals to hinge 
                s(time, bs = "cr", by = infestation, k = 4) + 
                # Random effect for id (smoothed - penalised)
                s(id, bs = 're'),
              data = infauna_mussels, 
              method = 'REML')
summary(m5_gam)

#  Check model fits
gratia::appraise(m3_gam, method = "simulate")
gratia::appraise(m5_gam, method = "simulate")
# Models do not fit too bad.

# Compare the models based on AICc
AICc(m3_gam)
AICc(m5_gam) # Best model (with k = 4)

# This is a first quick check if including 'date' in the model improved fit
# - There is a slight increase in model fit, indicated by the lower AICc
#   value for the model without the date effect, when including 'date'.



###############################################################################################
# Section 2.2: Model inference with DATE EFFECT ================================================
###############################################################################################

# Formal hypothesis test to see if the 'date' effect is statistically significant 
# - We do this using a Wald's test (Wood 2013) which compares the two models, 
#   without (m3_gam) and with (m5_gam), the 'date' fixed effect
anova(m3_gam, m5_gam, test = "F")
# Similar to the AICc analysis above, the effect of date is conclusive. 
# There is statistical evidence for an effect of 'date' on shell temperatures
# (F = 1.71, p < 0.05 *)



###############################################################################################
# Section 2.2: Plot model predictions with DATE EFFECT =========================================
###############################################################################################

# Create a set of data to predict over (using m5_gam)
new_data1 <- tidyr::expand(infauna_mussels, nesting(id, infestation, date, species),
                           time = unique(time))
View(new_data1)


# Extract predictions from the model 
best_mod_pred1 <- bind_cols(new_data1,
                            as.data.frame(predict(m5_gam, newdata = new_data1,
                                                  se.fit = TRUE))) 
head(best_mod_pred1)

# Make plot
ggplot(best_mod_pred1,
       aes(
         x = time,
         y = fit,
         group = id,
         colour = infestation
       ))+
  geom_line(aes(group = id, linetype = species, colour = infestation)) +
  facet_wrap(~ infestation + date) +
  scale_color_manual(values = my_colors) +
  geom_point(
    data = infauna_mussels, 
    aes(
      x = time, 
      y = shell_temp, 
      colour = infestation,
      shape = species
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
  geom_line(aes(group = id, linetype = species)) +
  facet_wrap(~ infestation) +
  labs(
    x = "Exposure time (minutes)",
    y = "Shell temperature (degrees C)"
  ) +
  theme(legend.position = "right")