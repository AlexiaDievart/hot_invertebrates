## ---------------------------
##
## Script name: Desiccation stress for infaunal species on mussel beds with different infestation status
##
## Purpose of script: Do infaunal species lose less water and less fast on infested mussel beds
##                    rather than on non-infested mussel beds?
##
## Author: Dr.Guy F. Sutton and Alexia Dievart
##
## Date Created: 2022-11-17
## Dates Updated: 2022-11-27
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
data <- read.csv("./desic_data1.csv", dec = ",", header = T, sep = ";")
head(data)
dplyr::glimpse(data)

# Create a unique identifier ID for each replicate, because of repeated measurements 
data <- data %>%
  dplyr::group_by(date, species, replicate, specimen) %>%
  tidyr::unite(col = "id", date, species, replicate, specimen, sep = "_", remove = FALSE) %>%
  dplyr::mutate(id = as.factor(id)) %>% 
  dplyr::ungroup() %>%
  dplyr::mutate(infestation = dplyr::if_else(infestation == "infested", 
                                             "Infested",
                                             "Non-infested"))
dplyr::glimpse(data)

# Convert some variables into other types 
data$date <- as.factor(data$date)
data$time <- as.numeric(data$time)
data$species <- as.factor(data$species)
data$infestation <- as.factor(data$infestation)
data$replicate <- as.factor(data$replicate)
data$specimen <- as.factor(data$specimen)
data$wet_wgt <- as.numeric(data$wet_wgt)
dplyr::glimpse(data)



###############################################################################################
# Section: Visualise data ---------------------------------------------------------------------
###############################################################################################

# Summarise wet weight by date, species, time and infestation class (only 1 day)
plot_data <- data %>%
  dplyr::group_by(date, species, time, infestation) %>%
  dplyr::summarise(
    wgt.mean = mean(wet_wgt),
    wgt.sd   = sd(wet_wgt),
    n = n(),
    wgt.se   = wgt.sd / sqrt(n)
  )
head(plot_data)
View(plot_data)

# Make plot 
my_colors <- c("darkgrey", "darkorange")

plot_data1 <- plot_data %>%
  ggplot(data = .,
         aes(
           x = time,
           y = wgt.mean,
           group = infestation, 
           color = infestation
         )) +
  geom_errorbar(
    aes(ymin = wgt.mean - wgt.se, 
        ymax = wgt.mean + wgt.se),
    width = 1,
    size = 1) +
  geom_line(lwd = 1, linetype = "solid") +
  scale_color_manual(values = my_colors) +
  scale_x_continuous(
    breaks = seq(0, 90, by = 15),
    limits = c(-1, 91)
  ) +
  scale_y_continuous(
    breaks = seq(0, 5500, by = 500),
    limits = c(0, 5500)
  ) +
  ylab("Wet weight (mg)\n") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "right") +
  # Changing the legend title
  labs(colour = "Infestation",
       x = "Time (min)") +
  facet_wrap(~ date + species)
plot_data1

# There is no clear difference in wet weight between organisms on infested and non-infested mussel beds,
# which is true for all species investigated.



###############################################################################################
# Section: Model fitting without DAY EFFECT ---------------------------------------------------
###############################################################################################

######################################################
# Model #5: Generalised additive model 
######################################################

# Fit GAM #5.2
# - Linear effects of infestation and species and non-linear effect for time 
m2_gam <- gam(log(wet_wgt) ~ 
                # Fixed effects (linear terms)
                infestation + species + 
                # Non-linear terms - cubic regression spline with 5 knots
                # - Allows the pieces between time intervals to hinge 
                s(time, bs = "cr", k = 5) + 
                # Random effect for id (smoothed - penalised)
                s(id, bs = 're'),
              data = data, 
              method = 'REML')
summary(m2_gam)

# Fit Gam #5.3
# - Linear effects of infestation and species and non-linear effect for time 
m3_gam <- gam(log(wet_wgt) ~ 
                # Fixed effects (linear terms)
                infestation * species + 
                # Non-linear terms - cubic regression spline with 5 knots
                # - Allows the pieces between time intervals to hinge 
                s(time, bs = "cr", k = 5) + 
                # Random effect for id (smoothed - penalised)
                s(id, bs = 're'),
              data = data, 
              method = 'REML')
summary(m3_gam)

# Fit GAM #5.4
m4_gam <- gam(log(wet_wgt) ~ 
                # Fixed effects (linear terms)
                infestation * species + 
                # Non-linear terms - cubic regression spline with 5 knots
                # - Allows the pieces between time intervals to hinge 
                s(time, bs = "cr", by = infestation, k = 5) + 
                # Random effect for id (smoothed - penalised)
                s(id, bs = 're'),
              data = data, 
              method = 'REML')
summary(m4_gam)

# Check it log-transformation is important for the wet weight
m5_gam <- gam(wet_wgt ~ 
                # Fixed effects (linear terms)
                infestation * species + 
                # Non-linear terms - cubic regression spline with 5 knots
                # - Allows the pieces between time intervals to hinge 
                s(time, bs = "cr", by = infestation, k = 5) + 
                # Random effect for id (smoothed - penalised)
                s(id, bs = 're'),
              data = data, 
              method = 'REML')

# Check model fit for fitted GAM's 
# - We want to know whether the model captured the data structure properly
gratia::appraise(m2_gam, method = "simulate")
gratia::appraise(m3_gam, method = "simulate")
gratia::appraise(m4_gam, method = "simulate")
gratia::appraise(m5_gam, method = "simulate")

# They all seem fine to me once we log transform the `wet_weight` variable 

# Compare the models based on AICc
AICc(m2_gam)
AICc(m3_gam) # Best model with k = 5
AICc(m4_gam)



###############################################################################################
# Section: Model inference --------------------------------------------------------------------
###############################################################################################

#################################################
# Test for treatment effect of infestation status and species 
#################################################

# Perform F-tests to compare parametric terms 
anova(m4_gam, test = "F")
# There is NO evidence for an effect of infestation status on wet weight (F = 1.015, p = 0.314),
# but there is evidence for an effect of species on wet weight (F = 239.74, p < 0.001). 
# There is no interaction between species and infestation.

#################################################
# Test for treatment effect of species to vary over time 
#################################################

mnull <- gam(log(wet_wgt) ~ 
               # Fixed effects (linear terms)
               infestation + species + 
               # Non-linear terms - cubic regression spline with 5 knots
               # - Allows the pieces between time intervals to hinge 
               s(time, bs = "cr", by = infestation, k = 5) + 
               # Random effect for id (smoothed - penalised)
               s(id, bs = 're'),
             data = data, 
             method = 'REML')

# Perform Wald's test (Wood 2013a,b)
anova(m3_gam, m4_gam, test = "Chisq") 
# There is NO evidence for the wet weigth differences between species to vary with time 
# (Chisq = 0.37, p > 0.05). 



###############################################################################################
# Section: Plot model predictions -------------------------------------------------------------
###############################################################################################

# Create a set of data to predict over (using m3_gam)
new_data <- tidyr::expand(data, nesting(id, infestation, species),
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
  terms = c("species", "infestation", "time [0:90 by = 1"),
  type = "fe",
  interval = "confidence",
  back.transform = TRUE) %>%
  as.data.frame() %>%
  dplyr::rename(
    infestation = group,
    species = x,
    time = facet
  ) %>%
  dplyr::mutate(time = as.double(as.character(time)))
head(preds)

# Make plot
# - Each dashed line represents an individual organism
# - The thick line represents the marginal (average/mean) effect
# - This conveys the variation in wet weight over time 
preds %>%
  ggplot(data = ., aes(
    x = time,
    y = predicted,
    colour = infestation,
    fill = infestation
  )) +
  geom_line() +
  # Remove geom_ribbon to get rid of the messy confidence bands 
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  facet_wrap(~ species, nrow = 1) +
  labs(
    x = "Exposure time (minutes)",
    y = "Wet weight (mg)"
  ) +
  theme_bw()