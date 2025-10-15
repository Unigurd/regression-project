library(dplyr)
library(ggplot2)
library(survival)
library(tidyr)

setwd("C:/Users/SIMON/Desktop/GitHub/regression-project")
surv_data <- read.csv("training_data.csv") %>%
  filter(is.na(CVD) == F) # make sure there are no NAs in CVD (there are in fact none)
# I do this analysis based on the values for individuals in the baseline examination!  
# Below is basically Sigurd's code for the model_input.csv data but it keeps all variables
surv_data_final <- surv_data %>%
  tibble() |>
  mutate(LDLC=NULL, HDLC=NULL) |>
  filter(
    PERIOD == 1, # most importantly!
    HEARTRTE < 200,
    TOTCHOL < 600,
    SYSBP < 250,
  ) |>
  drop_na()

# check distribution of TIMECVD
time_dens <- plot(density(surv_data$TIMECVD))
time_dens <-  plot(density(surv_data_final$TIMECVD))
# plot the CDF
time_cdf <- ecdf(surv_data$TIMECVD)
plot(time_cdf) # I don't understand why this doesn't reach 1?
plot(ecdf(surv_data_final$TIMECVD))

# should we maybe consider TIMECVD only for those who actually got CVD?
# create a subset of people with CVD
CVD_positive <- subset(surv_data, CVD == 1) 
CVD_positive_final <- subset(surv_data_final, CVD == 1)
# distribution and CDF
time_dens_pos <- density(CVD_positive$TIMECVD)
time_cdf_pos <- ecdf(CVD_positive$TIMECVD)
plot(time_dens_pos)
plot(time_cdf_pos)

time_dens_pos_final <- density(CVD_positive_final$TIMECVD)
time_cdf_pos_final <- ecdf(CVD_positive_final$TIMECVD)
plot(time_dens_pos_final)
plot(time_cdf_pos_final)

# CDF is pretty much a linear line: "flow" of CVD-diagnosis is constant over time 
# with small bump at a couple of days after first examination

# survival functions: survfit
# our time variable is TIME and the event is CVD
# since we have right-rensored data I added type = right. I am not sure what this changes
# Is it right that we don't use the TIMECVD variable here? 
# I copy the structure here from RwR chap. 11 (fig. 11.2)
surv_arg <- Surv(time = surv_data$TIME, event = surv_data$CVD, type = "right")
surv_arg_final <- Surv(time = surv_data_final$TIMECVD, event = surv_data_final$CVD, type = "right")

# Kaplan-Meier version (default)
survfit_KM <- survfit(surv_arg ~ 1) # why ~ 1 ?
survfit_KM_final <- survfit(surv_arg_final ~ 1)
# Nelson-Aalen version
survfit_NA <- survfit(surv_arg ~ 1,
                      stype = 2)

survfit_NA_final <- survfit(surv_arg_final ~ 1,
                            stype = 2)
# Nelson-Aalen with Fleming-Harrington correction
survfit_NA_FH <- survfit(surv_arg ~ 1,
                      stype = 2,
                      ctype = 2)
survfit_NA_FH_final <- survfit(surv_arg_final ~ 1,
                         stype = 2,
                         ctype = 2)

surv_plot <- summary(survfit_KM, data.frame = TRUE) %>%
  ggplot(aes(time, surv)) +
  geom_step(data = summary(survfit_NA, data.frame = TRUE), color = "blue") +
  geom_step(data = summary(survfit_NA_FH, data.frame = TRUE), color = "red")  +
  geom_step() +
  ylim(0,1) +
  xlab("Time (days)") + 
  ylab("Probability of not having CVD")
surv_plot

# The result is almost identical for all specifications (KM, NA, NA_FH)

# This graph seems odd. There are clear steps for when people drop out
# Why does the probability decrease so drastically after ~4,000 days?
# I'm stupid: I used the whole training sample where I should only have used those that appear in the first examination
# Do it again with this limited sample AND with TIMECVD as the time variable as it should be:

surv_plot_final <- summary(survfit_KM_final, data.frame = TRUE) %>%
  ggplot(aes(time, surv)) +
  geom_step(data = summary(survfit_NA_final, data.frame = TRUE), color = "blue") +
  geom_step(data = summary(survfit_NA_FH_final, data.frame = TRUE), color = "red")  +
  geom_step() +
  ylim(0,1) +
  xlab("Time (days)") + 
  ylab("Probability of not having CVD")
surv_plot_final


# This gives a blank grid. Maybe because I use TIME and not TIMECVD?
# I use TIMECVD now

surv_arg_CVD <- Surv(time = surv_data_final$TIMECVD, event = surv_data_final$CVD, type = "right")
survfit_CVD <- survfit(surv_arg_CVD ~ 1)
surv_plot_CVD <- summary(survfit_CVD, data.frame = TRUE) %>%
  ggplot(aes(time, surv)) +
  geom_step()
surv_plot_CVD

# This looks like something we would want to see! The survival function is decreasing at a constant rate
# The survival prob. at the end of the sollow-up periods is about 70 %



