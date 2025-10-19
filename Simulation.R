### Define the Simulation Functions (bug-free)
# Define parameters
set.seed(2025)
n <- 2000
beta_true <- c(age = 0.03, smoker = 0.40)  # true regression coefficients
gamma0_val <- 0.5                           # baseline time parameter
times_vec <- c(1000, 2000, 4000, 8766)      # censoring times
R <- 500                                    # Monte Carlo repetitions
# Function: generate event times T* (no recursion issues)
generate_Tstar <- function(eta_values, gamma_val) {
  # Random uniform U(0,1)
  u <- runif(length(eta_values))
  # Inversion of logistic survival model
  tstar <- ((1 / u - 1) * exp(-eta_values))^(1 / gamma_val)
  return(tstar)
}
### Run a Single Simulation
single_sim <- function(n, beta, gamma_val, times) {
  X_age <- rnorm(n, mean = 50, sd = 10)
  X_smoker <- rbinom(n, 1, 0.43)
  eta_vals <- beta["age"] * X_age + beta["smoker"] * X_smoker
  Tstar <- generate_Tstar(eta_vals, gamma_val)
    res_list <- list()
  for (tc in times) {
    Yt <- as.integer(Tstar <= tc)
    df_sim <- data.frame(Yt, X_age, X_smoker)
    fit <- glm(Yt ~ X_age + X_smoker, family = binomial, data = df_sim)
    res_list[[as.character(tc)]] <- coef(fit)
  }
  return(res_list)
}
demo <- single_sim(n = n, beta = beta_true, gamma_val = gamma0_val, times = times_vec)
demo

### Monte Carlo Simulation Across Different Censoring Times
mc_repeat_fixed_t <- function(R, n, beta, gamma_val, times) {
  k <- length(times)
  out_age <- matrix(NA, nrow = R, ncol = k)
  out_smoker <- matrix(NA, nrow = R, ncol = k)
    for (r in seq_len(R)) {
    X_age <- rnorm(n, mean = 50, sd = 10)
    X_smoker <- rbinom(n, 1, 0.43)
    eta_vals <- beta["age"] * X_age + beta["smoker"] * X_smoker
    Tstar <- generate_Tstar(eta_vals, gamma_val)
        for (j in seq_along(times)) {
      tc <- times[j]
      Yt <- as.integer(Tstar <= tc)
      df_sim <- data.frame(Yt, X_age, X_smoker)
      fit <- glm(Yt ~ X_age + X_smoker, family = binomial, data = df_sim)
      out_age[r, j] <- coef(fit)["X_age"]
      out_smoker[r, j] <- coef(fit)["X_smoker"]
    }
  }
  return(list(age = out_age, smoker = out_smoker))
}
# Run Monte Carlo
set.seed(2025)
mc_out <- mc_repeat_fixed_t(R, n, beta_true, gamma0_val, times_vec)

### Summarize Results and Plot
summary_table <- data.frame(
  time = times_vec,
  age_mean = apply(mc_out$age, 2, mean),
  age_sd = apply(mc_out$age, 2, sd),
  smoker_mean = apply(mc_out$smoker, 2, mean),
  smoker_sd = apply(mc_out$smoker, 2, sd)
)
print(summary_table)
# Plot: Standard Error vs Time
plot_df <- tibble(
  time = rep(summary_table$time, 2),
  sd = c(summary_table$age_sd, summary_table$smoker_sd),
  term = rep(c("age", "smoker"), each = nrow(summary_table))
)
p <- ggplot(plot_df, aes(x = time, y = sd, color = term)) +
  geom_line(size = 1) + geom_point(size = 2) +
  labs(
    title = "Monte Carlo SD of Coefficients vs Censoring Time",
    x = "Censoring time t",
    y = "Monte Carlo SD (empirical SE)"
  ) +
  theme_minimal()
print(p)

###Add Death Censoring Scenarios
## Independent death censoring scenario
simulate_death_censoring <- function(R, n, beta, gamma_val, lambda_death = 1/12000) {
  out_age <- numeric(R)
  out_smoker <- numeric(R)
    for (r in seq_len(R)) {
    X_age <- rnorm(n, mean = 50, sd = 10)
    X_smoker <- rbinom(n, 1, 0.43)
    eta_vals <- beta["age"] * X_age + beta["smoker"] * X_smoker
    Tstar <- generate_Tstar(eta_vals, gamma_val)
       # Death censoring time ~ Exp(lambda_death)
    Cdeath <- rexp(n, rate = lambda_death)
        # Observed outcome: event occurs before death
    Yobs <- as.integer(Tstar <= Cdeath)
        fit <- glm(Yobs ~ X_age + X_smoker, family = binomial)
    out_age[r] <- coef(fit)["X_age"]
    out_smoker[r] <- coef(fit)["X_smoker"]
  }
   summary_death <- data.frame(
    scenario = "Independent death censoring",
    age_mean = mean(out_age),
    age_sd = sd(out_age),
    smoker_mean = mean(out_smoker),
    smoker_sd = sd(out_smoker)
  )
  return(summary_death)
}
# Run it
set.seed(2025)
summary_death <- simulate_death_censoring(R = R, n = n, beta = beta_true, gamma_val = gamma0_val)
print(summary_death)

## Informative death censoring scenario
simulate_informative_death <- function(R, n, beta, gamma_val, base_rate = 1/20000, k = 0.02) {
  out_age <- numeric(R)
  out_smoker <- numeric(R)
    for (r in seq_len(R)) {
    X_age <- rnorm(n, mean = 50, sd = 10)
    X_smoker <- rbinom(n, 1, 0.43)
    eta_vals <- beta["age"] * X_age + beta["smoker"] * X_smoker
    Tstar <- generate_Tstar(eta_vals, gamma_val)
        # Death rate depends on age → informative censoring
    death_rate <- base_rate * exp(k * (X_age - mean(X_age)))
    Cdeath <- rexp(n, rate = death_rate)
        Yobs <- as.integer(Tstar <= Cdeath)
        fit <- glm(Yobs ~ X_age + X_smoker, family = binomial)
    out_age[r] <- coef(fit)["X_age"]
    out_smoker[r] <- coef(fit)["X_smoker"]
  }
   summary_info <- data.frame(
    scenario = "Informative death censoring",
    age_mean = mean(out_age),
    age_sd = sd(out_age),
    smoker_mean = mean(out_smoker),
    smoker_sd = sd(out_smoker)
  )
   return(summary_info)
}
# Run it
set.seed(2025)
summary_info <- simulate_informative_death(R = R, n = n, beta = beta_true, gamma_val = gamma0_val)
print(summary_info)

### Combine fixed-time summary with death censoring results
summary_all <- dplyr::bind_rows(
  summary_table %>% mutate(scenario = paste0("Fixed t=", time)) %>% 
    select(scenario, age_mean, age_sd, smoker_mean, smoker_sd),
  summary_death,
  summary_info
)
print(summary_all)
# Save for report
write.csv(summary_all, "simulation_results_all.csv", row.names = FALSE)




