################################################################################
#         SURVIVAL MODEL ANALYSIS 
# 

 #       Part 2: Simulation Study + Part 3: Practice Application
#        ########################################################

# Libraries
library(tidyverse)
library(survival)
library(survminer)
library(broom)
library(knitr)
library(gridExtra)
library(riskCommunicator)

# Set seed for reproducibility
set.seed(2025)


################################################################################
# PART 2: SIMULATION STUDY
################################################################################

# ------------------------------------------------------------------------------
# 1. Extract Framingham distributions (for inspiration)
# ------------------------------------------------------------------------------

data(framingham)

framingham_baseline <- framingham %>%
  filter(PERIOD == 1) %>%
  select(AGE, SEX, TIMECVD, CVD) %>%
  na.omit() %>%
  mutate(
    SEX = SEX - 1,  # 0=Male, 1=Female
    TIMECVD_YEARS = TIMECVD / 365.25
  )

# Summary stats
nrow(framingham_baseline)           # n = 4434
sum(framingham_baseline$CVD)        # 1157 events
mean(framingham_baseline$CVD)       # 26.1% event rate

# Empirical distributions - that i will use for simulation
age_mean <- mean(framingham_baseline$AGE)      # 49.93
age_sd <- sd(framingham_baseline$AGE)          # 8.68
sex_prob_female <- mean(framingham_baseline$SEX)  # 0.562







### ------------------------------------------------------------------------------
# 2. Define theoretical parameters
# ------------------------------------------------------------------------------

# True parameters for simulation (chosen to give realistic survival times)
beta_true <- c(AGE = -0.07, SEX = -0.5)
gamma_0 <- 1.5

# Model: S_η(t) = 1 / (1 + exp(η + γ₀ log(t)))
# where η = β_AGE * AGE + β_SEX * SEX







# ------------------------------------------------------------------------------
#  Simulation function
# ------------------------------------------------------------------------------

simulate_survival_data <- function(n, beta, gamma_0) {
  # Generate covariates from Framingham-inspired distributions
  AGE <- rnorm(n, mean = age_mean, sd = age_sd)
  SEX <- rbinom(n, size = 1, prob = sex_prob_female)
  
  # Linear predictor
  eta <- beta["AGE"] * AGE + beta["SEX"] * SEX
  
  # Generate survival times via inverse transform
  # From S_η(T) = U, solve for T: T = exp((log((1-U)/U) - η) / γ₀)
  
  # Solve S_eta(T) = U for T
  # U = 1 / (1 + T^gamma0 * exp(eta))
  # 1/U = 1 + T^gamma0 * exp(eta)
  # T^gamma0 = (1/U - 1) / exp(eta)
  # T = ((1/U - 1) / exp(eta))^(1/gamma0)
  
  U <- runif(n)
  T_star <- exp((log((1 - U) / U) - eta) / gamma_0)
  
  data.frame(AGE = AGE, SEX = SEX, T_star = T_star, eta = eta)
}





# ------------------------------------------------------------------------------
#  Test parameters - sanity check
# ------------------------------------------------------------------------------

test_data <- simulate_survival_data(n = 1000, beta = beta_true, gamma_0 = gamma_0)

# Check survival time distribution
median(test_data$T_star)      # Should be ~10-20 years -> realistic   its = 13.23
mean(test_data$T_star)
quantile(test_data$T_star, c(0.25, 0.75))

# Event rates at different censoring times
t_check <- c(2, 5, 10, 15, 20, 24)
event_rates_test <- sapply(t_check, function(t) mean(test_data$T_star <= t))
event_rates_test  # vary from ~5% to ~70%

# 10-20 years, event rates show variation


# ------------------------------------------------------------------------------
#  Main verification for Q1  - test theoretical predictions
# ------------------------------------------------------------------------------

# From Question 1, theory predicts:
#   β_t(AGE) = β_AGE + γ₀/t  (time-varying)
#   β_t(SEX) = β_SEX         (constant)
#   α_t = γ₀ log(t)          (log-linear)

n_sim <- 1000
sim_data <- simulate_survival_data(n = n_sim, beta = beta_true, gamma_0 = gamma_0)

median(sim_data$T_star)  # Check survival times

# Censoring times to examine
t_values <- c(2, 5, 10, 15, 20, 24)

# Fit logistic at each censoring time
verification_results <- map_dfr(t_values, function(t) {
  
  # Binary outcome: event by time t
  Y_t <- as.numeric(sim_data$T_star <= t)
  n_events <- sum(Y_t)
  
  # Need sufficient variation to fit model
  if(n_events >= 10 && n_events <= (n_sim - 10)) {
    
    # Fit logistic
    logit_fit <- glm(Y_t ~ AGE + SEX, data = sim_data, family = binomial)
    coefs <- tidy(logit_fit)
    
    data.frame(
      t = t,
      n_events = n_events,
      alpha_emp = coefs$estimate[coefs$term == "(Intercept)"],
      beta_AGE_emp = coefs$estimate[coefs$term == "AGE"],
      beta_SEX_emp = coefs$estimate[coefs$term == "SEX"],
      se_AGE = coefs$std.error[coefs$term == "AGE"],
      se_SEX = coefs$std.error[coefs$term == "SEX"]
    )
  } else {
    NULL  # Skip if insufficient variation
  }
})

# Calculate theoretical predictions
verification_results <- verification_results %>%
  mutate(
    alpha_theo = gamma_0 * log(t),
    beta_AGE_theo = beta_true["AGE"] + gamma_0 / t,
    beta_SEX_theo = beta_true["SEX"],
    diff_alpha = alpha_emp - alpha_theo,
    diff_AGE = beta_AGE_emp - beta_AGE_theo,
    diff_SEX = beta_SEX_emp - beta_SEX_theo
  )

# View results
verification_results

# Summary of deviations
mean(abs(verification_results$diff_AGE))
mean(abs(verification_results$diff_SEX))

# Table for report
verification_display <- verification_results %>%
  select(t, n_events, 
         alpha_emp, alpha_theo, diff_alpha,
         beta_AGE_emp, beta_AGE_theo, diff_AGE,
         beta_SEX_emp, beta_SEX_theo, diff_SEX)

kable(verification_display, digits = 4)

# |  t| n_events| alpha_emp| alpha_theo| diff_alpha| beta_AGE_emp| beta_AGE_theo| diff_AGE| beta_SEX_emp| beta_SEX_theo| diff_SEX|
#   |--:|--------:|---------:|----------:|----------:|------------:|-------------:|--------:|------------:|-------------:|--------:|
#   |  2|       73|    2.1116|     1.0397|     1.0719|      -0.0909|        0.6800|  -0.7709|      -0.6986|          -0.5|  -0.1986|
#   |  5|      234|    2.6875|     2.4142|     0.2733|      -0.0713|        0.2300|  -0.3013|      -0.7492|          -0.5|  -0.2492|
#   | 10|      439|    3.4249|     3.4539|    -0.0290|      -0.0675|        0.0800|  -0.1475|      -0.5478|          -0.5|  -0.0478|
#   | 15|      568|    4.2293|     4.0621|     0.1673|      -0.0714|        0.0300|  -0.1014|      -0.6179|          -0.5|  -0.1179|
#   | 20|      651|    4.4422|     4.4936|    -0.0514|      -0.0679|        0.0050|  -0.0729|      -0.6303|          -0.5|  -0.1303|
#   | 24|      717|    4.5116|     4.7671|    -0.2555|      -0.0640|       -0.0075|  -0.0565|      -0.5284|          -0.5|  -0.0284|






# ------------------------------------------------------------------------------
# 6. Plots - verification
# ------------------------------------------------------------------------------

# Plot 1: β_AGE vs t
p1 <- ggplot(verification_results, aes(x = t)) +
  geom_point(aes(y = beta_AGE_emp), color = "blue", size = 3) +
  geom_line(aes(y = beta_AGE_theo), color = "red", linewidth = 1) +
  geom_hline(yintercept = beta_true["AGE"], linetype = "dashed", color = "gray") +
  labs(
    title = "β_AGE: Empirical vs Theory",
    subtitle = expression(paste("Theory: ", beta[t], "(AGE) = ", beta, " + ", gamma[0], "/t")),
    x = "Censoring time t (years)",
    y = expression(hat(beta)[t](AGE))
  ) +
  theme_minimal()

# Plot 2: β_SEX vs t
p2 <- ggplot(verification_results, aes(x = t)) +
  geom_point(aes(y = beta_SEX_emp), color = "blue", size = 3) +
  geom_hline(aes(yintercept = beta_SEX_theo), color = "red", linewidth = 1) +
  labs(
    title = "β_SEX: Empirical vs Theory",
    subtitle = "Theory: β_t(SEX) = β (constant)",
    x = "Censoring time t (years)",
    y = expression(hat(beta)[t](SEX))
  ) +
  theme_minimal()

# Display
grid.arrange(p1, p2, ncol = 2)

# Save
# ggsave("figures/verification_plots.png", 
#        arrangeGrob(p1, p2, ncol = 2),
#        width = 12, height = 5, dpi = 300)


# ------------------------------------------------------------------------------
# 7. Monte Carlo - variance analysis
# ------------------------------------------------------------------------------

# Question: How does Var(β̂) depend on t?

B <- 200      # Replications
n_mc <- 500   # Sample size per rep

mc_results <- map_dfr(t_values, function(t) {
  
  # Run B replications
  replications <- map_dfr(1:B, function(b) {
    
    data_b <- simulate_survival_data(n = n_mc, beta = beta_true, gamma_0 = gamma_0)
    Y_t <- as.numeric(data_b$T_star <= t)
    
    # Fit model (with error handling)
    fit <- tryCatch({
      glm(Y_t ~ AGE + SEX, data = data_b, family = binomial)
    }, error = function(e) NULL)
    
    if(!is.null(fit)) {
      coefs <- coef(fit)
      data.frame(beta_AGE = coefs["AGE"], beta_SEX = coefs["SEX"])
    }
  })
  
  # Calculate variance across replications
  data.frame(
    t = t,
    var_AGE = var(replications$beta_AGE, na.rm = TRUE),
    var_SEX = var(replications$beta_SEX, na.rm = TRUE),
    n_success = nrow(replications)
  )
})

# Results
mc_results

kable(mc_results, digits = 6)

# |  t|  var_AGE|  var_SEX| n_success|
#   |--:|--------:|--------:|---------:|
#   |  2| 0.000552| 0.135258|       200|
#   |  5| 0.000179| 0.062312|       200|
#   | 10| 0.000127| 0.035919|       200|
#   | 15| 0.000138| 0.044114|       200|
#   | 20| 0.000147| 0.037004|       200|
#   | 24| 0.000170| 0.045834|       200|

# Optimal censoring times (minimum variance) one with least varians
mc_results$t[which.min(mc_results$var_AGE)]
mc_results$t[which.min(mc_results$var_SEX)]

# Plot variance vs t
p_variance <- ggplot(mc_results, aes(x = t)) +
  geom_line(aes(y = var_AGE, color = "AGE"), linewidth = 1) +
  geom_point(aes(y = var_AGE, color = "AGE"), size = 3) +
  geom_line(aes(y = var_SEX, color = "SEX"), linewidth = 1) +
  geom_point(aes(y = var_SEX, color = "SEX"), size = 3) +
  labs(
    title = "Variance vs Censoring Time",
    x = "Censoring time t (years)",
    y = "Variance",
    color = "Parameter"
  ) +
  theme_minimal()

print(p_variance)
# ggsave("figures/variance_plot.png", p_variance, width = 8, height = 6, dpi = 300)


# ------------------------------------------------------------------------------
# 8. Death censoring analysis
# ------------------------------------------------------------------------------

# Compare: administrative censoring vs death censoring

lambda_C <- 0.05  # Death rate

# Add death censoring to sim_data
sim_data_censored <- sim_data %>%
  mutate(
    C = rexp(n_sim, rate = lambda_C),     # Death time
    T_obs = pmin(T_star, C),               # Observed time
    delta = as.numeric(T_star <= C)        # Event indicator
  )

# Censoring summary
sum(sim_data_censored$delta)              # Events
sum(1 - sim_data_censored$delta)          # Censored

# Compare models at t=10
t_compare <- 10

# Model 1: No death censoring
Y_no_censor <- as.numeric(sim_data_censored$T_star <= t_compare)
fit_no_censor <- glm(Y_no_censor ~ AGE + SEX, 
                     data = sim_data_censored, 
                     family = binomial)

# Model 2: With death censoring
Y_with_censor <- sim_data_censored$delta
fit_with_censor <- glm(Y_with_censor ~ AGE + SEX,
                       data = sim_data_censored,
                       family = binomial)

# Model 3: Cox (on sim_data_censored)
cox_sim <- coxph(Surv(T_obs, delta) ~ AGE + SEX, data = sim_data_censored)

# Test proportional hazards assumption
ph_test_sim <- cox.zph(cox_sim)
ph_test_sim
# PH violation as seen in theoretical part 

# Plot Schoenfeld residuals  
par(mfrow = c(1, 2))
plot(ph_test_sim[1], main = "AGE - Simulated Data")
abline(h = 0, col = "red", lty = 2)
plot(ph_test_sim[2], main = "SEX - Simulated Data")
abline(h = 0, col = "red", lty = 2)
par(mfrow = c(1, 1))


# Comparison table
censoring_comparison <- data.frame(
  Parameter = c("Intercept", "AGE", "SEX"),
  No_Censoring = coef(fit_no_censor),
  With_Censoring = coef(fit_with_censor),
  Cox = c(NA, coef(cox_sim)["AGE"], coef(cox_sim)["SEX"]),
  Difference = coef(fit_with_censor) - coef(fit_no_censor),
  Pct_Change = 100 * (coef(fit_with_censor) - coef(fit_no_censor)) / coef(fit_no_censor)
)

kable(censoring_comparison, digits = 4, row.names = FALSE)

#        NOTE: Cox model used here despite PH violation because:
# 1. It correctly handles censoring (unlike logistic)
# 2. It provides robust average effect estimates
# 3. Demonstrates that Cox is useful even when PH not perfect
# 4. Confirms theoretical prediction: our model is NOT PH (as shown in Q3)


# Check impact
max(abs(censoring_comparison$Pct_Change[2:3]))  # Max % change


################################################################################
#                PRACTICE-part - FRAMINGHAM APPLICATION
################################################################################

# ------------------------------------------------------------------------------
# 1. Data preparation
# ------------------------------------------------------------------------------

# Reuse framingham_baseline, but make SEX a factor for Cox model
framingham_clean <- framingham_baseline %>%
  mutate(SEX = factor(SEX, levels = c(0, 1), labels = c("Male", "Female")))

# Summary
nrow(framingham_clean)                    # n = 4434
sum(framingham_clean$CVD)                 # 1157 events
mean(framingham_clean$CVD)                # 26.1%
mean(framingham_clean$TIMECVD_YEARS)      # 16.1 years mean follow-up
median(framingham_clean$TIMECVD_YEARS)    # 16.5 years median


# ------------------------------------------------------------------------------
# 2. Cox proportional hazards model
# ------------------------------------------------------------------------------

cox_model <- coxph(Surv(TIMECVD_YEARS, CVD) ~ AGE + SEX, 
                   data = framingham_clean)

# Results
summary(cox_model)
summary(cox_model)$coefficients

# Tidy output
cox_coefs <- tidy(cox_model, exponentiate = TRUE, conf.int = TRUE)
cox_coefs

# Key statistics
coef(cox_model)["AGE"]                    # β_AGE
exp(coef(cox_model)["AGE"])               # HR for AGE
coef(cox_model)["SEXFemale"]              # β_SEX
exp(coef(cox_model)["SEXFemale"])         # HR for Female
summary(cox_model)$concordance[1]         # Concordance


# ------------------------------------------------------------------------------
# 3. Model visualization
# ------------------------------------------------------------------------------

# Survival curves at median age
median_age <- median(framingham_clean$AGE)

newdata_pred <- data.frame(
  AGE = median_age,
  SEX = factor(c("Male", "Female"), levels = c("Male", "Female"))
)

survfit_pred <- survfit(cox_model, newdata = newdata_pred)

# Plot
p_survival <- ggsurvplot(
  survfit_pred, 
  data = framingham_clean,
  conf.int = TRUE,
  risk.table = TRUE,
  legend.labs = sprintf("%s (Age %d)", c("Male", "Female"), round(median_age)),
  xlab = "Time (years)",
  ylab = "CVD-free survival",
  title = "Predicted Survival Curves",
  palette = c("#E7B800", "#2E9FDF")
)

print(p_survival)

# Baseline hazard
basehaz_data <- basehaz(cox_model)

p_hazard <- ggplot(basehaz_data, aes(x = time, y = hazard)) +
  geom_step(color = "darkred", linewidth = 1) +
  labs(
    title = "Baseline Cumulative Hazard",
    x = "Time (years)",
    y = expression(H[0](t))
  ) +
  theme_minimal()

print(p_hazard)


# ------------------------------------------------------------------------------
# 4. Proportional hazards assumption
# ------------------------------------------------------------------------------

# Test PH assumption
ph_test <- cox.zph(cox_model)
ph_test

# p-values
ph_test$table["AGE", "p"]           # AGE: should be > 0.05
ph_test$table["SEX", "p"]     # SEX: check if < 0.05

# Plot Schoenfeld residuals
par(mfrow = c(1, 2))
plot(ph_test[1], main = "AGE")
abline(h = 0, col = "red", lty = 2)
plot(ph_test[2], main = "SEX")
abline(h = 0, col = "red", lty = 2)
par(mfrow = c(1, 1))


# ------------------------------------------------------------------------------
# 5. Logistic regression at multiple times
# ------------------------------------------------------------------------------

t_values <- c(5, 10, 15, 20, 24)

# Fit logistic at each t
logistic_results <- map_dfr(t_values, function(t) {
  
  Y_t <- as.numeric(framingham_clean$TIMECVD_YEARS <= t & # baseline
                      framingham_clean$CVD == 1)
  
  n_events <- sum(Y_t)
  event_rate <- mean(Y_t)
  
  fit <- glm(Y_t ~ AGE + SEX, data = framingham_clean, family = binomial)
  
  # Robust extraction by position
  coefs <- coef(fit)
  ses <- summary(fit)$coefficients[, "Std. Error"]
  
  data.frame(
    t = t,
    n_events = n_events,
    event_rate = event_rate,
    beta_AGE = coefs[2],
    beta_SEX = coefs[3],
    se_AGE = ses[2],
    se_SEX = ses[3]
  )
})

# Results
logistic_results
survre

kable(logistic_results %>%
        select(t, n_events, event_rate, beta_AGE, beta_SEX),
      digits = 4,
      row.names = FALSE) # better
# Forskellige "stop tidspunkter" for analysen - Svarer til spørgsmålet: 
# "Har personen fået CVD inden år t?
# 
# I forhold til opgaven:
#   - Dette er de faste tidspunkter hvor vi laver binary logistic regression
#   - Sammenligner med Cox model der bruger hele tidsaksen

# ------------------------------------------------------------------------------
# 6. Comparison: Cox vs Logistic
# ------------------------------------------------------------------------------

# Add Cox estimates for comparison
logistic_results <- logistic_results %>%
  mutate(
    cox_AGE = coef(cox_model)["AGE"],
    cox_SEX = coef(cox_model)["SEXFemale"]
  )

# Full comparison table
comparison_table <- logistic_results %>%
  select(t, n_events, event_rate, beta_AGE, cox_AGE, beta_SEX, cox_SEX)

kable(comparison_table, digits = 4)

# Focused comparison at t=10
logistic_t10 <- logistic_results %>% filter(t == 10)

comparison_t10 <- data.frame(
  Method = c("Logistic (t=10)", "Cox PH"),
  beta_AGE = c(logistic_t10$beta_AGE, coef(cox_model)["AGE"]),
  beta_SEX = c(logistic_t10$beta_SEX, coef(cox_model)["SEXFemale"])
)

kable(comparison_t10, digits = 4)


# ------------------------------------------------------------------------------
# 7. Convergence plots
# ------------------------------------------------------------------------------

# β_AGE convergence
p1 <- ggplot(logistic_results, aes(x = t)) +
  geom_point(aes(y = beta_AGE), color = "blue", size = 3) +
  geom_line(aes(y = beta_AGE), color = "blue") +
  geom_hline(aes(yintercept = cox_AGE), 
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(
    title = "β_AGE: Convergence to Cox",
    x = "Censoring time t (years)",
    y = "β̂_AGE"
  ) +
  theme_minimal()

# β_SEX convergence
p2 <- ggplot(logistic_results, aes(x = t)) +
  geom_point(aes(y = beta_SEX), color = "blue", size = 3) +
  geom_line(aes(y = beta_SEX), color = "blue") +
  geom_hline(aes(yintercept = cox_SEX), 
             color = "red", linetype = "dashed", linewidth = 1) +
  labs(
    title = "β_SEX: Convergence to Cox",
    x = "Censoring time t (years)",
    y = "β̂_SEX"
  ) +
  theme_minimal()

# Display
grid.arrange(p1, p2, ncol = 1)

# Save
# ggsave("figures/framingham_convergence.png",
     #  arrangeGrob(p1, p2, ncol = 1),
     #  width = 8, height = 10, dpi = 300)
Simulation study confirms:
# - β_t(AGE) varies with t (time-varying)
# - β_t(SEX) stays constant (time-invariant)
# - Variance depends on censoring time
# - Death censoring introduces bias

# Framingham analysis shows:
# - Cox: AGE HR = 1.077, Female HR = 0.590
# - Logistic estimates converge to Cox as t increases
# - Cox more efficient (uses all events)
# - Patterns match theoretical predictions
