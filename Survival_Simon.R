# look to towards the end of this file for analysis that actually makes sense (the survreg part)
library(dplyr)
library(ggplot2)
library(survival)
library(tidyr)
library(broom)
library(ggbeeswarm)
library(ggsurvfit)
library(flexsurv)
library(survminer)

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
  mutate(Prev_sick = ifelse(PREVCHD == 1 | PREVAP == 1 | PREVMI == 1 | PREVSTRK == 1 | PREVHYP == 1 , 1, 0)) |>
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
  geom_step() +
  ylim(0,1) +
  xlab("Time (days)") + 
  ylab("Probability of not having CVD")
surv_plot_CVD

# This looks like something we would want to see! The survival function is decreasing at a constant rate
# The survival prob. at the end of the sollow-up periods is about 70 %

# What happens if we treat TIMECVD as uncensored?
# Will use this version as baseline
surv_arg_uncensored <- Surv(time = surv_data_final$TIMECVD, event = surv_data_final$CVD)
survfit_uncensored <- survfit(surv_arg_uncensored ~ 1)
surv_plot_uncensored <- summary(survfit_uncensored, data.frame = TRUE) %>%
  ggplot(aes(time+0.5, surv)) + # use time +0.5 here as in RwR
  geom_step() +
  ylim(0,1) +
  xlab("Time (days)") + 
  ylab("Probability of not having CVD")
surv_plot_uncensored
# no noticeable change

# Maybe there are interesting results when you condition on having had a coronary heart disease?
survfit_PREVCHD <- survfit(surv_arg_CVD ~ surv_data_final$PREVCHD)
surv_plot_CHD <- summary(survfit_PREVCHD, data.frame = TRUE) %>%
  ggplot(aes(time, surv, colour = )) +
  geom_step() +
  ylim(0,1) +
  xlab("Time (days)") + 
  ylab("Probability of not having CVD given prev. CHD")
surv_plot_CHD

# this plot is weird. 

# What about sex?
survfit_sex <- survfit(surv_arg_CVD ~ surv_data_final$SEX)
surv_plot_sex <- summary(survfit_sex, data.frame = TRUE) %>%
  ggplot(aes(time, surv)) +
  geom_step() +
  ylim(0,1) +
  xlab("Time (days)") + 
  ylab("Probability of not having CVD given sex")
surv_plot_sex

# what does this tell me?

# do a regression explaining impact of previous diseases
surv_reg_prev <- survreg(
  Surv(TIMECVD + 0.5 , CVD == 1) ~ Prev_sick,
  data = surv_data_final
)
surv_reg_prev

# below I do the Cox regression model BUT we actually CANNOT use it because the two key assumptions are 
# non-informative censoring
# proportional hazards: this is not satisfied as we showed in theory part
# we can adjust for this by adding tt(age) apparently, or rstpm2: stpm2

surv_reg_cox <- coxph(
  Surv(TIMECVD + 0.5, CVD == 1) ~ Prev_sick,
  data = surv_data_final
)
surv_reg_cox

baseline_cox <- survfit(surv_reg_cox) 
summary(baseline_cox, data.frame = TRUE) |>
  ggplot(aes(time + 0.5, surv)) + 
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "gray90") +
  geom_step() +
  ylim(0, 1) +
  xlab("time (months)") + 
  ylab("survival probability") 

# create predictor data:
pred_frame <- tibble(
  data = 1:2,
  Prev_sick = factor(1:2, labels = levels(as.factor(surv_data_final$Prev_sick))),
  strata = factor(1:2, labels = paste("Prev_sick=", Prev_sick, sep = ""))
)

  survfit(surv_reg_cox, newdata = pred_frame) |>
    summary(data.frame = TRUE) |>
    left_join(pred_frame) |>
    ggplot(aes(time, surv, color = strata)) + 
    geom_line(data = fit, linetype = 2) +
    geom_step() +
    ylim(0, 1) +
    xlab("time (months)") + 
    ylab("survival probability") +
    theme(legend.position = "top")
  
  
  augment(surv_reg_cox, data = surv_data_final) |> 
    ggplot(aes(Prev_sick, .resid)) + 
    ggbeeswarm::geom_beeswarm(alpha = 0.4) + 
    geom_point(
      stat = "summary", 
      fun = mean, 
      color = "blue", 
      size = 2.5
    ) +
    geom_errorbar(
      stat = "summary", 
      fun.data = mean_se, 
      fun.args = list(mult = 1.96),
      color = "blue", 
      width = 0.15, 
      linewidth = 0.7
    )

  #---------
  # new approach using ggsurvfit (from https://www.emilyzabor.com/survival-analysis-in-r.html)
  s1 <- survfit(Surv(TIMECVD, CVD) ~ 1, data = surv_data_final)
str(s1) 

survfit2(Surv(TIMECVD, CVD) ~ 1, data = surv_data_final) %>%
  ggsurvfit() +
  labs(
    x = "Days",
    y = "Overall survival probability"
  ) + 
  add_confidence_interval() +
  add_risktable()

summary(survfit(Surv(TIMECVD, CVD) ~ 1, data = surv_data_final), times = 3652.5)

#-----------------------------
# Proper approach using survreg
# start off by using same regressors as in my extended model from logmodel_firsteps
# add +0.5 to survival time to account for log(0) problem
reg_log <- survreg(formula = Surv(TIMECVD +0.1 , CVD) ~ .,
        data = surv_data_final, 
        dist = "loglogistiC") # can try this also with dist = one of the following:
# “extreme”, “logistic”, “gaussian”, “weibull”, “exponential”, “rayleigh”, “loggaussian”, “lognormal”, “loglogistic”, “t”

# residual plots
p1 <- augment(reg_log, data = surv_data_final, type.predict = "expected") |>
  survfit(Surv(.fitted, CVD == 1) ~ 1, data = _) |>
  summary(data.frame = TRUE) |>
  ggplot(aes(time, cumhaz)) + 
  geom_abline(slope = 1, color = "blue") +
  geom_point() +
  xlab("Cox-Snell residuals") + 
  ylab("Cumulative hazard") +


p2 <- augment(reg_log) |> 
  ggplot(aes(.fitted, .resid)) +
  geom_point(alpha = 0.4) + 
  geom_smooth() + 
  xlab("Fitted values") +
  ylab("Martingale residuals")

gridExtra::grid.arrange(p1, p2, ncol = 2)

# AGE vs. resids
augment(reg_log, data = surv_data_final) %>%
  ggplot(aes(AGE, .resid)) +
  geom_point() +
  geom_smooth() +
  geom_errorbar(
    stat = "summary", 
    fun.data = mean_se, 
    fun.args = list(mult = 1.96),
    color = "blue", 
    width = 0.15, 
    linewidth = 0.7
  )

# Cox-Snell residuals look like log?
# From chap. 11.5.1:
# n practice we plot the estimated cumulative hazards against the Cox-Snell residuals,
# which should fall on a straight line with slope 1

# Some model suggestions:
form <- Surv(TIMECVD + 0.1, CVD) ~ AGE + SEX + TOTCHOL + SYSBP + PREVCHD + CIGPDAY + DIABETES + BMI

fits <- list(
  weibull = flexsurvreg(form, data = surv_data_final, dist = "weibull"),
  lognormal = flexsurvreg(form, data = surv_data_final, dist = "lnorm"),
  loglogistic = flexsurvreg(form, data = surv_data_final, dist = "llogis"),
  gengamma = flexsurvreg(form, data = surv_data_final, dist = "gengamma")
)

