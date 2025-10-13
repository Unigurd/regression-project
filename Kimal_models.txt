library(tidyverse); library(broom); library(car); library(pROC)
library(ResourceSelection); library(rsample); library(glmnet)
library(recipes); library(yardstick); library(vip); library(splines)
library(rms); library(boot); library(ggplot2); library(dplyr)

set.seed(5055)
df <- read_csv("C:/Users/MARCMIRATECHGLOBAL/Desktop/CPH/model_input.csv")
###1) Clean & type conversions (exact names)
  df <- df %>%
  mutate(
    SEX = factor(SEX, levels = c(1, 2), labels = c("Male", "Female")),
    educ = factor(educ, levels = c(1, 2, 3, 4)),
    DIABETES = factor(DIABETES, levels = c(0, 1), labels = c("No", "Yes")),
    BPMEDS = factor(BPMEDS, levels = c(0, 1), labels = c("No", "Yes")),
    PREVCHD = factor(PREVCHD, levels = c(0, 1), labels = c("No", "Yes")),
    PREVAP = factor(PREVAP, levels = c(0, 1), labels = c("No", "Yes")),
    PREVMI = factor(PREVMI, levels = c(0, 1), labels = c("No", "Yes")),
    PREVHYP = factor(PREVHYP, levels = c(0, 1), labels = c("No", "Yes"))
      )

# Marginal distribution for continuous predictors
# Define the continuous variables
vars <- c("AGE","TOTCHOL","SYSBP","BMI","HEARTRTE","GLUCOSE")
# Reshape data to long format for easy plotting
df_long <- df %>%
  select(all_of(vars)) %>%
  pivot_longer(cols = everything(), names_to = "Variable", values_to = "Value")
# Plot histograms with density overlay
ggplot(df_long, aes(x = Value)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "skyblue", color = "white", alpha = 0.7) +
  geom_density(color = "red", linewidth = 1) +
  facet_wrap(~ Variable, scales = "free", ncol = 3) +
  theme_minimal(base_size = 13) +
  labs(title = "Marginal Distributions of Continuous Variables",
       x = "Value", y = "Density")
# Transformations for GLUCOSE and CIGPDAY
# Apply log  transformations for GLUCOSE
df <- df %>%
  mutate(
    log_GLUCOSE = log1p(GLUCOSE),  # handles GLUCOSE skewness
    SMOKE = factor(ifelse(CIGPDAY > 0, "Yes", "No"), levels = c("No","Yes"))
  )
# Check
str(df)
summary(df)
#Correlation of only the continuous variables
cont_vars <- c("AGE","TOTCHOL","SYSBP","BMI","HEARTRTE"," log_GLUCOSE ")
cor_matrix <- cor(df[, cont_vars], use = "complete.obs")
# Plot a correlation heatmap
library(corrplot)
corrplot(cor_matrix, method = "color", type = "upper",
         tl.col = "black", tl.srt = 45, addCoef.col = "black",
         number.cex = 0.7)

###2) Define candidate model(s)
# MODEL 1: All variables included
model1 <- glm(CVD ~ ., data = df, family = "binomial")
summary(model1)
# Collinearity test 
vif(model1)
#After vif for collinearity, PREVCHD and PREVAP are taken out of the model 
# MODEL 2: All variables except PREVCHD and PREVAP
model1 <- glm(CVD ~ AGE + SEX + SYSBP + TOTCHOL + SMOKE + BMI+ DIABETES
+ BPMEDS+ HEARTRTE + log_GLUCOSE + educ + PREVMI + PREVHYP, data = df, family = "binomial")
summary(model2)
# Exponentiated coefficients (Odds Ratios) and 95% CIs
or_tab <- broom::tidy(model1) %>%
  mutate(OR = exp(estimate),
         OR_low = exp(estimate - 1.96*std.error),
         OR_high = exp(estimate + 1.96*std.error)) %>%
  select(term, estimate, std.error, OR, OR_low, OR_high, p.value)
or_tab
#Model3: Interaction variables considered
model2 <- glm(CVD ~ AGE * SEX + SYSBP * BPMEDS + BMI * DIABETES + TOTCHOL + CIGPDAY + HEARTRTE + log_GLUCOSE + educ + PREVMI + PREVHYP, data = df, family = "binomial")
summary(model3)
# Model 4: Model with Spline
dd <- datadist(df)
options(datadist = "dd")
model4  <- glm(CVD ~ rcs(AGE, 4) + SEX + rcs(SYSBP, 4) + BMI + rcs(TOTCHOL, 4) + log_GLUCOSE + DIABETES + BPMEDS + PREVHYP + rcs(HEARTRTE, 4) + educ, data = df, family = "binomial")
anova(model4)

###3) Goodness-of-fit, residuals & calibration
# Hosmer-Lemeshow
hoslem <- ResourceSelection::hoslem.test(df$CVD, fitted(model1), g = 10)
hoslem
# Standardized residuals vs fitted
df <- df %>% mutate(pred = fitted(model1),
                            std_resid = rstandard(model1))
ggplot(df, aes(x = pred, y = std_resid)) +
  geom_point(alpha = 0.6) + geom_hline(yintercept = 0, linetype = "dashed") +
  geom_smooth(se = FALSE, color = "red") + theme_minimal()
# Calibration by deciles
cal_df <- df %>% mutate(dec = ntile(pred, 10)) %>%
  group_by(dec) %>% summarise(mean_pred = mean(pred), obs = mean(CVD), n = n())
ggplot(cal_df, aes(x = mean_pred, y = obs)) + geom_point() + geom_line() +
  geom_abline(intercept = 0, slope = 1, color = "red") + theme_minimal()
###4) Discrimination (AUC), Brier, and threshold metrics
# AUC
roc_obj <- roc(df$CVD, df$pred)
auc_val <- auc(roc_obj); auc_val
ci_auc <- ci.auc(roc_obj)
ci_auc
# Brier score
brier <- mean((df$CVD - df$pred)^2); brier
# Confusion matrix at threshold 0.5
df <- df %>% mutate(pred_class = as.integer(pred > 0.5))
table(Observed = df$CVD, Predicted = df$pred_class)
# Sens/Spec at optimal threshold (Youden)
coords(roc_obj, "best", ret = c("threshold","sensitivity","specificity"))

###5) Model selection & penalized regression (LASSO)
# Prepare design matrix
x <- model.matrix(~ AGE + SEX + SYSBP + TOTCHOL + SMOKE + BMI+ DIABETES
+ BPMEDS+ HEARTRTE + log_GLUCOSE + educ + PREVMI + PREVHYP, data = df)[,-1]
y <- df$CVD
set.seed(2025)
cv_lasso <- cv.glmnet(x, y, family = "binomial", alpha = 1, nfolds = 10)
plot(cv_lasso)
coef(cv_lasso, s = "lambda.min")    # coefficients at lambda.min
coef(cv_lasso, s = "lambda.1se")    # simpler model at lambda.1se
# Evaluate lasso chosen model on same data (or use CV predictions)
pred_lasso <- predict(cv_lasso, newx = x, s = "lambda.min", type = "response")
auc(pROC::roc(y, as.numeric(pred_lasso)))

