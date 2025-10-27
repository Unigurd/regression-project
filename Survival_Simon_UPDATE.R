## This file is the practical survival regression part of the assignment
## I start it from scratch again because my initial try was pathetic
library(dplyr)
library(ggplot2)
library(survival)
library(tidyr)
library(broom)
library(survminer)
library(riskCommunicator)
library(purrr)

data("framingham")
framingham_clean <- as_tibble(framingham) |>
  filter(PERIOD == 1) |>
  mutate(
    X         = NULL,
    ## GLUCOSE   = log(GLUCOSE),
    ## BMI       = log(BMI)
    ## HEARTRE   = log(HEARTRTE)
    SEX = SEX - 1,  # 0=Male, 1=Female
    SEX = factor(SEX, labels = c("Male", "Female")),
    educ      = factor(educ, levels = c(1, 2, 3, 4)),
    TIMECVD_YEARS = TIMECVD / 365.25
  )

#--------------------------------
# Cox proportional hazards model
#--------------------------------

# We want to compare the model to the binary regression model, so we use the same covariates
cox_form1  <- Surv(TIMECVD, CVD) ~ SEX + AGE + SYSBP + TOTCHOL +
  PREVMI + DIABETES + PREVHYP + log(HEARTRTE) + CIGPDAY +
  log(BMI) + educ + log(GLUCOSE) + PREVCHD + PREVAP + BPMEDS

cox <- coxph(cox_form1, data = framingham_clean)

# coefficients and significance:
cox_tidy <- broom::tidy(cox, exponentiate = TRUE, conf.int = TRUE)
knitr::kable(cox_tidy %>% select(term, estimate, conf.low, conf.high, p.value),
             digits=3, caption="Cox PH model: hazard ratios (HR), 95% CI, p-value")

# test proportional hazards assumption of Cox model
cox.zph(cox)

# different estimators of survival functions
# Kaplan Meier
framingham_KM <- survfit(
  Surv(TIMECVD, CVD) ~ 1,
  data = framingham_clean
)

# Nelson Aalen
framingham_NA <- survfit(
  Surv(TIMECVD, CVD) ~ 1,
  data = framingham_clean,
  stype = 2
)

# Nelson Aalen with Fleming-Harrington correction for ties
framingham_NA_FH <- survfit(
  Surv(TIMECVD, CVD) ~ 1,
  data = framingham_clean,
  stype = 2,
  ctype = 2
)

# plot all of these together
survival_plot <- summary(framingham_KM, data.frame = TRUE) |> 
  ggplot(aes(time, surv)) + 
  geom_step(data = summary(framingham_NA, data.frame = TRUE), color = "blue") +
  geom_step(data = summary(framingham_NA_FH, data.frame = TRUE), color = "red") +
  geom_step() +
  ylim(0, 1) +
  xlab("time (months)") + 
  ylab("probability of not having CVD")

# zoom in to detect how they differ
survival_plot_zoom <- survival_plot + coord_cartesian(xlim = c(2500, 2520), ylim = c(0.905, 0.915))
# difference is barely detectable and of no practical relevance

# We can see how the survival plot depends on sex
framingham_KM_sex <- survfit(
  Surv(TIMECVD, CVD) ~ SEX,
  data = framingham_clean
)

ggsurvplot(
  fit = framingham_KM_sex,
  data = framingham_clean,
  xlab = "Time (months)",
  ylab = "Probability of not having CVD",
  legend.title = "Sex",
  legend = "top"
)
# men are more likely to get CVD at any point in time and the hazard ratio seems to become bigger with time

# ------------------------------
# Model Diagnosis
#-------------------------------

# Cox-Snell residuals
Cox_Snell_resid <- augment(cox, data = framingham_clean, type.predict = "expected") |>
    survfit(Surv(.fitted, CVD) ~ 1, data = _) |>
    summary(data.frame = TRUE) |>
    ggplot(aes(time, cumhaz)) + 
    geom_abline(slope = 1, color = "blue") +
    geom_point() +
    xlab("Cox-Snell residuals") + 
    ylab("Cumulative hazard")
# This looks pretty!

# Martingale residuals
Martingale_resid <- augment(cox) |> 
  ggplot(aes(.fitted, .resid)) +
  geom_point() + 
  geom_smooth() + 
  xlab("Fitted values") +
  ylab("Martingale residuals")
# This one has approx. mean 0

surv_residual_plot <- gridExtra::grid.arrange(Cox_Snell_resid, Martingale_resid, ncol = 2)

# Deviance residuals
Deviance_resid <- augment(cox, type.residuals = "deviance") |>
  ggplot(aes(.fitted, .resid)) +
  geom_point(alpha = 0.5) +
  geom_smooth() +
  xlab("Fitted values") +
  ylab("Deviance residuals") +
  theme_minimal()

# Check distribution of deviance residuals
resid_data <- augment(cox, type.residuals = "deviance")

ggplot(resid_data, aes(x = .resid)) +
  geom_histogram(bins = 30, color = "white") +
  theme_minimal() +
  xlab("Deviance residuals")
# they are skewed (even though not terribly)

#------------------------------------
# Plot Martingale residuals against continuous variables to detect possible misspecifications
resid_data <- augment(cox)

# Select only continuous covariates to check nonlinearity
continuous_vars <- c("AGE", "SYSBP", "TOTCHOL", "log(HEARTRTE)", "CIGPDAY", "log(BMI)", "log(GLUCOSE)")

# Create a long-format data frame for faceting
resid_long <- resid_data |>
  select(all_of(continuous_vars), .resid) |>
  pivot_longer(cols = all_of(continuous_vars), names_to = "variable", values_to = "value")

# Plot Martingale residuals vs each covariate
ggplot(resid_long, aes(x = value, y = .resid)) +
  geom_point(alpha = 0.3) +
  geom_smooth(color = "blue") +
  facet_wrap(~ variable, scales = "free_x") +
  theme_minimal() +
  labs(
    x = "Covariate value",
    y = "Martingale residuals",
    title = "Martingale residuals by covariate",
    subtitle = "Curvature indicates potential nonlinearity"
  )

### can do the same with deviance residuals
resid_data <- augment(cox, type.residuals = "deviance")

# Create a long-format data frame for faceting
resid_long <- resid_data |>
  select(all_of(continuous_vars), .resid) |>
  pivot_longer(cols = all_of(continuous_vars), names_to = "variable", values_to = "value")

# Plot Martingale residuals vs each covariate
ggplot(resid_long, aes(x = value, y = .resid)) +
  geom_point(alpha = 0.3) +
  geom_smooth(color = "blue") +
  facet_wrap(~ variable, scales = "free_x") +
  theme_minimal() +
  labs(
    x = "Covariate value",
    y = "Deviance residuals",
    title = "Deviance residuals by covariate"
  )
# only log(Glucose) slightly wavy but still ok





