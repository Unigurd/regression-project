# Simon's modelling approaches
library(dplyr)
library(broom)
library(ggplot2)
library(patchwork)

# read in the data
setwd("C:/Users/SIMON/Desktop/GitHub/regression-project")
data <- read.csv("model_input.csv")

# here is again the nice correlation plot between variables included in this data we all use for modelling
# Compute Spearman correlation
cp <- cor(
  data.matrix(data), 
  use = "complete.obs", 
  method = "spearman"
)

# Visualize with hierarchical clustering
corrplot::corrplot(
  cp, 
  diag = FALSE, 
  order = "hclust", 
  addrect = 12,        # draws rectangles around clusters
  tl.srt = 45, 
  tl.col = "black", 
  tl.cex = 0.8
)

# -------
# Model assessment
# -------
# create a model assessment function
model_assess <- function(model){
  summary(model)
  print(exp(coef(model))) # Extracting the resulting odds-ratios
  
  tidy(model, exponentiate = TRUE, conf.int = TRUE) %>% # visualizing it
    ggplot(aes(x = term, y = estimate, ymin = conf.low, ymax = conf.high)) +
    geom_pointrange() +
    coord_flip() +
    labs(y = "Odds Ratio", x = "Predictor")
}

# Model 1
# do a simple log-model containing the most promising variables (highest correlation with CVD)
# omitting some variables where multicollinearity might be a problem (PREVMI, PREVAP, PREVHYP)
simple_form <- CVD ~ AGE + SEX + TOTCHOL + SYSBP + PREVCHD
log_model_simple <- glm(simple_form, family = binomial(link = "logit"), data = data)
model_assess(log_model_simple)

# Model 2
# adding some more variabls: CIGPDAY, DIBETES, BMI
extended_form <- CVD ~ AGE + SEX + TOTCHOL + SYSBP + PREVCHD + CIGPDAY + DIABETES + BMI
log_model_extended <- glm(extended_form, family = binomial(link = "logit"), data = data)
model_assess(log_model_extended)

extended_HRTR

### Comparing model 1 and model 2
# since the models are nested, we can use anova (likelihood ratio test) (AIC comparison leads to the same result)
anova(log_model_simple, log_model_extended, test = "Chisq")
# this shows that the extension makes the model fit significantly better
# Insights from the coefficients:
# PREVCHD, DIABETES, and SEX are strongest predictors of CVD
# The other variables are also significant, but seem to be less meaningful

# Model 3
# use all the variables available
super_extended_form <- CVD ~ .
log_model_super_extended <- glm(super_extended_form, family = binomial(link = "logit"), data = data)
model_assess(log_model_super_extended) 
# graph is useless because of the big CI of PREVMI: There are probably too few observations, so estimator is very uncertain
# similar for CVD in this case --> do not include these two variables 

# recommend to stick with extended model!

#------
# Assessment of model fit
#-----
augment(log_model_extended, type.residuals = "pearson")
# Estimated dispersion parameter
summary(log_model_extended)$dispersion
summary(log_model_simple)$dispersion
summary(log_model_super_extended)$dispersion
# ----- residual plots -------

# deviance residual plot (SIMPLE)
simple_diag <- augment(log_model_simple, type.predict = "response") # important to specify type.predict = repsonse to deliver probability outcomes
simple_diag$.dev.res <- residuals(log_model_simple, type = "deviance")
ggplot(simple_diag, aes(.fitted, .dev.res)) +
  geom_point() +
  geom_smooth()

# alternative plot with same result
plot(fitted(log_model_simple),
     residuals(log_model_simple, type = "deviance"),
     xlab = "Fitted values (Predicted probability of CVD)",
     ylab = "Deviance residuals",
     main = "Deviance Residuals vs Fitted Values")
abline(h = 0, lty = 2, col = "red")

# deviance residual plot (EXTENDED)
extended_diag <- augment(log_model_extended, type.predict = "response") # important to specify type.predict = repsonse to deliver probability outcomes
extended_diag$.dev.res <- residuals(log_model_extended, type = "deviance")
ggplot(extended_diag, aes(.fitted, .dev.res)) +
  geom_point() +
  geom_smooth()

# Age vs residuals: skewed
extended_diag %>% ggplot(aes(AGE, .resid))+
  geom_point()+
  geom_smooth()+
  ylab("Pearson residuals")

# variance vs fitted values
extended_diag %>%
  mutate(.fitted_group = cut(.fitted, quantile(.fitted, probs = seq(5, 95, 1) / 100)))

extended_diag %>% group_by(.fitted) |> # this one does not work for some reason
  summarize(.fitted_local = mean(.fitted), .var_local = var(.resid)) |>
  na.omit() |> 
  ggplot(aes(.fitted_local, .var_local)) + 
  geom_point() + 
  geom_hline(yintercept = 1, linetype = 2) +
  geom_smooth() + 
  xlab("fitted values") +
  ylab("variance")

# confidence intervals for all estimators
tidy(log_model_extended)[-1, "estimate"] |> 
  cbind(confint.lm(log_model_extended)[-1, ]) |>
  knitr::kable(digits = 2)

group_p <- function(y, p, n_probs = 100) {
  tibble(
    y = y,
    p = p,
    p_group = cut(p, quantile(p, probs = seq(0, n_probs, 1) / n_probs), include.lowest = TRUE)
  ) |>
    group_by(p_group) |>
    summarize(
      p_local = mean(p), 
      obs_local = mean(y), 
      n = n(), 
      se_local = sqrt(p_local * (1 - p_local) / n()),
    )
}

# ---- SUper extended ----

# deviance residual plot (SUPER EXTENDED)
super_extended_diag <- augment(log_model_super_extended, type.predict = "response") # important to specify type.predict = repsonse to deliver probability outcomes
super_extended_diag$.dev.res <- residuals(log_model_super_extended, type = "deviance")
ggplot(super_extended_diag, aes(.fitted, .dev.res)) +
  geom_point() +
  geom_smooth()

# deviance (black) and pearson (red) residual plot 
super_extended_diag$.pearson <- residuals(log_model_super_extended, type = "pearson")
 ggplot(super_extended_diag, aes(.fitted, .dev.res)) +
 geom_point() +
 geom_point(aes(y = .pearson), color = "red") +
 geom_smooth(color = "black", se = FALSE) +
 geom_smooth(aes(y = .pearson), color = "red", se = FALSE)
 
 # ------ simulating new responses from fitted model -----
 # simple
 B <- 8
  p_simple <- list()
  for (b in 1:B) {
  simple_diag$y <- simulate(log_model_simple)[, 1]
  sim_glm_tmp <- glm(simple_form, data = simple_diag, family = binomial(link = "logit"))
  sim_tmp_diag <- augment(sim_glm_tmp, type.residuals = "pearson", type.predict = "response")
  p_simple[[b]] <- ggplot(sim_tmp_diag, aes(.fitted, .resid)) +
      geom_point() +
      geom_smooth()
  }

  # extended
  B <- 8
  p_extended <- list()
  for (b in 1:B) {
    simple_diag$y <- simulate(log_model_extended)[, 1]
    sim_glm_tmp <- glm(extended_form, data = extended_diag, family = binomial(link = "logit"))
    sim_tmp_diag <- augment(sim_glm_tmp, type.residuals = "pearson", type.predict = "response")
    p_extended[[b]] <- ggplot(sim_tmp_diag, aes(.fitted, .resid)) +
      geom_point() +
      geom_smooth()
  }
  
  # super extended
  B <- 8
  p_super_extended <- list()
  for (b in 1:B) {
    simple_diag$y <- simulate(log_model_super_extended)[, 1]
    sim_glm_tmp <- glm(super_extended_form, data = super_extended_diag, family = binomial(link = "logit"))
    sim_tmp_diag <- augment(sim_glm_tmp, type.residuals = "pearson", type.predict = "response")
    p_super_extended[[b]] <- ggplot(sim_tmp_diag, aes(.fitted, .resid)) +
      geom_point() +
      geom_smooth()
  }
  
  
  # Combine all plots in a list 
  combined_plot_simple <- wrap_plots(p_simple, nrow = 2, ncol = 4)
  combined_plot_extended <- wrap_plots(p_extended, nrow = 2, ncol = 4)
  combined_plot_super_extended <- wrap_plots(p_super_extended, nrow = 2, ncol = 4)
  combined_plot_simple
  combined_plot_extended
  combined_plot_super_extended
  
  # Interpretation:
  # there is something very wrong with the super extended model: Immense residuals
  # 
  #####
  # Pearson training error
  residuals(log_model_extended, type = "pearson")^2 %>% mean()
  
  pearson <- function(y, p) (y - p)^2 / (p * (1 - p))
   cv <- function(form, data, K = 10, M = 20) {
     n <- nrow(data)
     response <- all.vars(form)[1]
     loss <- vector("list", M)
     mu_hat <- numeric(n)
     for(b in 1:M) {
       group <- sample(rep(1:K, length.out = n))
       for(i in 1:K) {
         model <- glm(form, data = data[group != i, ], family = "binomial")
         mu_hat[group == i] <-
           predict(model, newdata = data[group == i, ], type = "response")
         }
       loss[[b]] <- pearson(data[, "CVD"], mu_hat)
       }
     unlist(loss) |> mean()
   }
   cv(log_model_simple, data = data)
   cv(log_model_extended, data = data)
   cv(log_model_super_extended, data = data)
   # result: again confirmation that super_extended model is very bad
   # cv is a bit higher in extended model
   
   ### Estimating AUC
   log_model_extended_aug <- augment(log_model_extended, data = data)
   eta1_ext <- filter(log_model_extended_aug, CVD == 1)$.fitted
   eta0_ext <- filter(log_model_extended_aug, CVD == 0)$.fitted
   AUC_ext <- wilcox.test(eta1_ext, eta0_ext)$statistic / (length(eta1_ext) * length(eta0_ext))
   AUC_ext
   
   log_model_simple_aug <- augment(log_model_simple, data = data)
   eta1_simple <- filter(log_model_simple_aug, CVD == 1)$.fitted
   eta0_simple <- filter(log_model_simple_aug, CVD == 0)$.fitted
   AUC_simple <- wilcox.test(eta1_simple, eta0_simple)$statistic / (length(eta1_simple) * length(eta0_simple))
   AUC_simple
   
   log_model_super_ext_aug <-  augment(log_model_super_extended, data = data)
   eta1_super_ext <- filter(log_model_super_ext_aug, CVD == 1)$.fitted
   eta0_super_ext <- filter(log_model_super_ext_aug, CVD == 0)$.fitted
   AUC_super_ext <- wilcox.test(eta1_super_ext, eta0_super_ext)$statistic / (length(eta1_super_ext) * length(eta0_super_ext))
   AUC_super_ext
   # how do I fit AUC using true etas???
   # AUC is higher for extended version (and even higher for super_extended version??????????????)
  
  
  # Cross-validation only makes sense in the context of splines (??) which we don't have: we don't have anything to tune on
  # below is a skeleton for how to do it (taken from book / slides) in case someone arrives with a specification including splines
   cv <- function(form, data, K = 10, M = 20) {
    n <- nrow(data)
    response <- all.vars(form)[1]
    loss <- vector("list", M)
    mu_hat <- numeric(n)
    for(b in 1:M) {
      # Generating the random subdivision into groups
      group <- sample(rep(1:K, length.out = n))
      for(i in 1:K) {
        model <- lm(form, data = data[group != i, ])
        mu_hat[group == i] <- predict(model, newdata = data[group == i, ])
      }
      loss[[b]] <- (data[, response] - mu_hat)^2
    }
    unlist(loss) |> mean()
  }
  
  set.seed(42)
  cv_extended <- tibble(
    cv = Vectorize(\(d) cv(extended_form, data)),
    mse = Vectorize(\(d) lm(extended_form, data) |> 
                      (\(x) residuals(x)^2)() |> mean())
  )
  
  # calculating training error
  train_err <- function(form, data) { 
    model <- lm(form, data = data)
    residuals(model)^2 |> mean()
  }
  
  val_err <- function(form, data, newdata) { 
    response <- all.vars(form)[1]
    model <- lm(form, data = data)
    (newdata[, response] - predict(model, newdata = newdata))^2 |> unlist() |> mean()  
  }
  
  ##----------------
  
  group_p <- function(y, p, n_probs = 100) {
    tibble(
      y = y,
      p = p,
      p_group = cut(p, quantile(p, probs = seq(0, n_probs, 1) / n_probs), include.lowest = TRUE)
    ) |>
      group_by(p_group) |>
      summarize(
        p_local = mean(p), 
        obs_local = mean(y), 
        n = n(), 
        se_local = sqrt(p_local * (1 - p_local) / n()),
      )
  }
  
  cal_plot <- ggplot(mapping = aes(p_local, obs_local)) + 
    geom_point() + 
    geom_line(aes(y = p_local + 1.96 * se_local), color = "red") + 
    geom_line(aes(y = p_local - 1.96 * se_local), color = "red") +
    geom_abline(slope = 1, intercept = 0) + 
    geom_smooth() + 
    xlab("predicted probability") +
    ylab("observed probability")
  
  cal_plot %+% group_p(data$CVD, fitted(log_model_extended)) + 
    coord_cartesian(xlim = c(0, 0.35), ylim = c(0, 0.35))  
  
  cal_plot %+% group_p(data$CVD, fitted(log_model_simple)) + 
    coord_cartesian(xlim = c(0, 0.35), ylim = c(0, 0.35)) 
  
  cal_plot %+% group_p(data$CVD, fitted(log_model_super_extended)) + 
    coord_cartesian(xlim = c(0, 0.35), ylim = c(0, 0.35)) 
  