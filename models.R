#.  Make some models and seehow well it fits the data.

# Loading data (framingham - cleaned dataset)

model_data <- framingham

# Convert categorical variables to factors WITH MEANINGFUL LABELS
model_data <- model_data %>%
  mutate(
    # SEX: 1 for Male,   2 for = Female
    SEX = factor(SEX, levels = c(1, 2), labels = c("Male", "Female")),
    
    # educ: 1-4  (No change here - just for indsight)
    educ = factor(educ, levels = c(1, 2, 3, 4)),
    
    # Binære variable: 0 = No, 1 = Yes
    DIABETES = factor(DIABETES, levels = c(0, 1), labels = c("No", "Yes")),
    BPMEDS = factor(BPMEDS, levels = c(0, 1), labels = c("No", "Yes")),
    PREVCHD = factor(PREVCHD, levels = c(0, 1), labels = c("No", "Yes")),
    PREVAP = factor(PREVAP, levels = c(0, 1), labels = c("No", "Yes")),
    PREVMI = factor(PREVMI, levels = c(0, 1), labels = c("No", "Yes")),
    PREVHYP = factor(PREVHYP, levels = c(0, 1), labels = c("No", "Yes"))
  )

                  ############# family = binomial(link = logit) #############
# ==================================================================================================
# MODEL 1

model_Bin_logit_1 <- glm(CVD ~ AGE + SEX + SYSBP + TOTCHOL + CIGPDAY + BMI,
              data = model_data, family = "binomial")
summary(model_Bin_logit_1)

min(fitted(model_Bin_logit_2))
max(fitted(model_Bin_logit_2))

anova(model_Bin_logit_2)
anova(model_Bin_logit_2, model_G, model_Bin_logit_2, model_negBin_logit_2)
# MODEL 2:
# All variables in cleanes dataset
model_Bin_logit_2 <- glm(CVD ~ ., data = model_data, family = "binomial")
summary(model_Bin_logit_2)

pearson_resid_b_l_model1 <- resid(model_Bin_logit_1, type = "pearson")
deviance_resid_b_l_model2  <- resid(model_Bin_logit_2, type = "deviance")

pearson_resid_b_l_model1
deviance_resid_b_l_model2

std_resid_b_l_modelp1 <- rstandard(model_Bin_logit_1)

# 4x4 diagnostic plot matrix
par(mfrow = c(2, 2))

# Plot 1: Residuals vs fitted values
plot(fitted(model_Bin_logit_1), pearson_resid_b_l_model1, 
     main = "Residuals vs Fitted",
     xlab = "Fitted Values", ylab = "Pearson Residuals")
abline(h = 0, col = "red")

# Plot 2: Q-Q plot of residuals
qqnorm(std_resid_b_l_modelp1, main = "Q-Q Plot of Standardized Residuals")
qqline(std_resid_b_l_modelp1, col = "red")

# Plot 3: Scale-Location plot
plot(fitted(model_Bin_logit_1), sqrt(abs(std_resid_b_l_modelp1)),
     main = "Scale-Location Plot",
     xlab = "Fitted Values", ylab = "√|Standardized Residuals|")

# Plot 4: Residuals vs leverage
plot(hatvalues(model_Bin_logit_1), std_resid_b_l_modelp1,
     main = "Residuals vs Leverage",
     xlab = "Leverage", ylab = "Standardized Residuals")
abline(h = 0, col = "red")

par(mfrow = c(1, 1))


# ==================================================================================================
                ############# family = poisson(link = log) #############.  - Bare for sjovt. -> det giver ikke god mening. 
# ==================================================================================================
# Model 1

model_Poi_log_1 <- glm(CVD ~ AGE + SEX + SYSBP + TOTCHOL + CIGPDAY + BMI,
                   data = model_data, family = poisson("log"))
summary(model_Poi_log_1)

# model 2
model_Poi_log_2 <- glm(CVD ~ ., data = model_data, family = poisson(link = "log"))
summary(model_Poi_log_2)

# ==================================================================================================
                    #############  NegativBinomial.  #############
# ==================================================================================================
# model 1

model_negBin_logit_1 <- glm.nb(CVD ~ AGE + SEX + SYSBP + TOTCHOL + CIGPDAY + BMI,
                            data = model_data)
summary(model_negBin_logit_1)

# Model 2
model_negBin_logit_2 <- glm.nb(CVD ~ ., data = model_data)
summary(model_negBin_logit_2)

# ==================================================================================================
                  ############# family = binomial(link = logit) #############
# ==================================================================================================

# Model 1

library(rpart)
tree_model <- rpart(
  CVD ~ SEX + TOTCHOL + AGE + SYSBP + CIGPDAY + BMI + 
    DIABETES + HEARTRTE + PREVMI + PREVHYP,
  data = model_data,
  method = "class"  # For classification
)

# Visualisér træet
install.packages("rpart.plot")
library(rpart.plot)
rpart.plot(tree_model)




alternative_models <- list(
  binomial = model_Bin_logit_2,
  poisson = model_Poi_log_2,
  negbinomial = model_negBin_logit_2
)
summary(alternative_models)

performance_comparison <- data.frame(
  Model = names(alternative_models),
  AIC = sapply(alternative_models[1:3], AIC),  # Kun for GLM modeller
  Observations = nrow(model_data)
)
performance_comparison









# =============================================================================================================================




# Variable med p > 0.05: BPMEDS, GLUCOSE, educ
improved_model <- glm(
  CVD ~ SEX + TOTCHOL + AGE + SYSBP + CIGPDAY + BMI + 
    DIABETES + HEARTRTE + PREVCHD + PREVAP + PREVMI + PREVHYP,
  data = model_data,
  family = "binomial"
)

# Sammenligning
AIC(model2)
AIC(improved_model)
AIC(model2) - AIC(improved_model)

# Check om forbedringen er signifikant
lr_test <- anova(improved_model, model2, test = "LRT")

# Performance sammenligning
improved_auc <- calculate_auc(improved_model, model_data)
auc_model2
improved_auc

# VIF for den forbedrede model

vif_improved <- vif(improved_model)
print(vif_improved)
vif(model_Bin_logit_2)
print(vif(model_Bin_logit_2))
summary(improved_model)
*# ================================================================








    ### Playing with some models.... ultimately hybrid as last model.


####################################### Cleaned dataset period 1 (from sigurd)

library(tibble)
library(broom)
library(readr)
library(skimr)
library(tidyverse)
library(dplyr)
library(tidyr)
library(hexbin)
library(knitr)

library(conflicted)
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")


data(framingham, package = "riskCommunicator")

framingham_data <- framingham

set.seed(8088)
framingham <- read.csv("training_data.csv") |>
  tibble() |>
  mutate(LDLC=NULL, HDLC=NULL) |>
  filter(
    PERIOD == 1,
    HEARTRTE < 200,
    TOTCHOL < 600,
    SYSBP < 250,
  ) |>
  select(
    CVD, SEX, TOTCHOL, AGE, SYSBP, CIGPDAY, BMI, DIABETES, BPMEDS, HEARTRTE, GLUCOSE, educ, PREVCHD, PREVAP, PREVMI, PREVHYP,
  ) |>
  drop_na()
write.csv(framingham, "model_input.csv")

## Excluded:
## - X
## - RANDID
## - TIME
## - PERIOD
## - DIABP      Correlated with SYSBP
## - CURSMOKE   Correlated with CIGPDAY
## - HDLC       Only measured in third checkup
## - LDLC       Only measured in third checkup
## - PREVSTRK   Too few stroke cases
## - DEATH
## - ANGINA
## - HOSPMI
## - MI_FCHD
## - ANYCHD
## - STROKE
## - HYPERTEN
## - TIMEAP
## - TIMEMI
## - TIMEMIFC
## - TIMECHD
## - TIMESTRK
## - TIMECVD
## - TIMEDTH
## - TIMEHYP

# ________________________________________________

library(tibble)
library(broom)
library(readr)
library(skimr)
library(dplyr)
library(tidyr)
library(hexbin)

framingham <- read.csv("training_data.csv") |>
  tibble() |>
  mutate(LDLC=NULL, HDLC=NULL) |>
  filter(
    PERIOD == 1,
    HEARTRTE < 200,
    TOTCHOL < 600,
    SYSBP < 250,
  ) |>
  select(
    CVD, SEX, TOTCHOL, AGE, SYSBP, CIGPDAY, BMI, DIABETES, BPMEDS, HEARTRTE, GLUCOSE, educ, PREVCHD, PREVAP, PREVMI, PREVHYP,
  ) |>
  drop_na()


# Loading data (framingham - cleaned dataset)

model_data <- framingham

# Convert categorical variables to factors WITH MEANINGFUL LABELS
model_data <- model_data %>%
  mutate(
    # SEX: 1 for Male,   2 for = Female
    SEX = factor(SEX, levels = c(1, 2), labels = c("Male", "Female")),
    
    # educ: 1-4  (No change here - just for indsight)
    educ = factor(educ, levels = c(1, 2, 3, 4)),
    
    # Binære variable: 0 = No, 1 = Yes
    DIABETES = factor(DIABETES, levels = c(0, 1), labels = c("No", "Yes")),
    BPMEDS = factor(BPMEDS, levels = c(0, 1), labels = c("No", "Yes")),
    PREVCHD = factor(PREVCHD, levels = c(0, 1), labels = c("No", "Yes")),
    PREVAP = factor(PREVAP, levels = c(0, 1), labels = c("No", "Yes")),
    PREVMI = factor(PREVMI, levels = c(0, 1), labels = c("No", "Yes")),
    PREVHYP = factor(PREVHYP, levels = c(0, 1), labels = c("No", "Yes"))
  )

# Check
str(model_data)
summary(model_data)

# Models: 
#========#

# MODEL 1 ->  Core risk factors 
model1 <- glm(CVD ~ AGE + SEX + SYSBP + TOTCHOL + CIGPDAY + BMI,
              data = model_data, family = "binomial")

# MODEL 2: All variables in cleanes dataset
model2 <- glm(CVD ~ ., data = model_data, family = "binomial")

# MODEL 3: Strategic interactions
model3 <- glm(CVD ~ AGE * SEX + SYSBP * BPMEDS + BMI * DIABETES + 
                TOTCHOL + CIGPDAY + HEARTRTE + GLUCOSE + educ +
                PREVCHD + PREVAP + PREVMI + PREVHYP,
              data = model_data, family = "binomial")

# MODEL 4: Stepwise selection
full_model <- glm(CVD ~ ., data = model_data, family = "binomial")
model4 <- step(full_model, direction = "both", trace = 0) 



summary(model1)
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
#   (Intercept) -7.5675873  0.4360332 -17.356  < 2e-16 ***
#   AGE          0.0526834  0.0053384   9.869  < 2e-16 ***   Greater CVD risk with age
#   SEXFemale   -0.9473376  0.0909895 -10.412  < 2e-16 ***  Female -0.9473376 => exp(..) = OR ~ 0.39. (model 1)
#   SYSBP        0.0176290  0.0020176   8.737  < 2e-16 ***  Højere blodtryk = højere risiko
#   TOTCHOL      0.0039735  0.0009754   4.074 4.62e-05 ***  Higher bloodpreasure = higher risk 
#   CIGPDAY      0.0086287  0.0036143   2.387  0.01697 *    higher Cigs smoked = higher risk
#   BMI          0.0345457  0.0106345   3.248  0.00116 **   Higher BMI = højere risk



summary(model2)
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
#   (Intercept) -6.309615   0.587415 -10.741  < 2e-16 ***
#   SEXFemale   -0.850608   0.095244  -8.931  < 2e-16 ***
#   TOTCHOL      0.004480   0.001013   4.423 9.73e-06 ***
#   AGE          0.045069   0.005759   7.826 5.04e-15 ***
#   SYSBP        0.012559   0.002719   4.619 3.86e-06 ***
#   CIGPDAY      0.010519   0.003776   2.786  0.00533 ** 
#   BMI          0.031126   0.011102   2.804  0.00505 ** 
#   DIABETESYes  0.868612   0.290353   2.992  0.00278 **  Diabetes == higher risk 
#   BPMEDSYes    0.081543   0.215251   0.379  0.70482    
#   HEARTRTE    -0.010569   0.003753  -2.816  0.00486 ** 
#   GLUCOSE      0.003503   0.002209   1.586  0.11277    
#   educ2        0.102718   0.105834   0.971  0.33177    
#   educ3       -0.133954   0.132610  -1.010  0.31243    
#   educ4       -0.214585   0.147335  -1.456  0.14527    
#   PREVCHDYes   1.889604   1.174545   1.609  0.10766    
#   PREVAPYes   -1.707517   1.175243  -1.453  0.14625    
#   PREVMIYes    3.400092   0.811071   4.192 2.76e-05 *** Pre heart attack = Much higher risk
#   PREVHYPYes   0.365476   0.120896   3.023  0.00250 ** Pre hypertension = higher risk


summary(model3)
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
#   (Intercept)     -6.5655070  0.6714168  -9.779  < 2e-16 ***
#   AGE              0.0442246  0.0072127   6.131 8.71e-10 ***
#   SEXFemale       -0.9047751  0.5623119  -1.609  0.10761    
#   SYSBP            0.0149480  0.0029094   5.138 2.78e-07 ***
#   BPMEDSYes        3.1685924  1.3207558   2.399  0.01644 *  
#   BMI              0.0310096  0.0114298   2.713  0.00667 ** 
#   DIABETESYes      1.4181556  1.2319940   1.151  0.24969    
#   TOTCHOL          0.0044250  0.0010273   4.307 1.65e-05 ***
#   CIGPDAY          0.0106079  0.0037816   2.805  0.00503 ** 
#   HEARTRTE        -0.0103188  0.0037548  -2.748  0.00599 ** 
#   GLUCOSE          0.0035233  0.0021971   1.604  0.10880    
#   educ2            0.1015932  0.1059731   0.959  0.33773    
#   educ3           -0.1429201  0.1328034  -1.076  0.28185    
#   educ4           -0.2252495  0.1477707  -1.524  0.12743    
#   PREVCHDYes       1.8767361  1.1745111   1.598  0.11007    
#   PREVAPYes       -1.6767487  1.1750441  -1.427  0.15359    
#   PREVMIYes        3.3932949  0.8128147   4.175 2.98e-05 ***
#   PREVHYPYes       0.2945204  0.1245176   2.365  0.01802 *  
#   AGE:SEXFemale    0.0009491  0.0106747   0.089  0.92915    
#   SYSBP:BPMEDSYes -0.0186115  0.0078743  -2.364  0.01810 *  
#   BMI:DIABETESYes -0.0192936  0.0435137  -0.443  0.65748   

summary(model4)
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
#   (Intercept) -6.384992   0.565074 -11.299  < 2e-16 ***
#   SEXFemale   -0.832106   0.094145  -8.839  < 2e-16 ***
#   TOTCHOL      0.004378   0.001009   4.340 1.42e-05 ***
#   AGE          0.044555   0.005610   7.942 1.98e-15 ***
#   SYSBP        0.012977   0.002689   4.825 1.40e-06 ***
#   CIGPDAY      0.010789   0.003766   2.865  0.00417 ** 
#   BMI          0.031596   0.010984   2.877  0.00402 ** 
#   DIABETESYes  0.857811   0.289950   2.958  0.00309 ** 
#   HEARTRTE    -0.010122   0.003738  -2.708  0.00677 ** 
#   GLUCOSE      0.003541   0.002209   1.603  0.10901    
#   PREVCHDYes   1.890613   1.173362   1.611  0.10712    
#   PREVAPYes   -1.714938   1.174325  -1.460  0.14419    
#   PREVMIYes    3.360172   0.808247   4.157 3.22e-05 ***
#   PREVHYPYes   0.362866   0.120369   3.015  0.00257 ** 

# _________________________________


# All of these are HIGHLY SIGNIFICANT (p < 0.01)
# ' - AGE             Greater CVD risk with age
#   - SEXFemale       Female -0.9473376 => exp(..) = OR ~ 0.39. (model 1)
#   - SYSBP           Higher bloodpreasure = higher risk 
#   - TOTCHOL         Higher kolesterol = higher risk
#   - CIGPDAY         higher Cigs smoked = higher risk
#   - BMI             Higher BMI = højere risk
#   - PREVMIYes       Pre heart attack = Much higher risk
#   - PREVHYPYes      Pre hypertension = higher risk
#   - DIABETESYes     Diabetes == higher risk ' 


############ Model 1 (Core element/factor only few variables) ##########

#   -> shows and confirms the classical risk factor we'd expect to see
# ==========================================================================


# ############ model 2 (All variables in cleaned dataset)  ###########

#   -> New strong predictors: Diabetes, previous heart attack (PREVMI), 
#     and existing hypertension (PREVHYP) emerge as major risk factors

#   Interesting and suprisingly result: Higher heart rate (HEARTRTE) have (negative coefficient) appears protective, hm?

#   Not significant: Blood pressure meds alone, glucose levels, education, and some previous conditions don't add independent predictive value
# ==========================================================================


########### model 3: (interaction) 

#.  -> SYSBP:BPMEDSYes -0.0186115  0.0078743  -2.364  0.01810 * Significant interaction.
#.    taking blood pressure meds reduces risk from high blood pressure

# NON-significant interactions:
# AGE:SEXFemale    0.0009491     # No difference in age effect between genders
# BMI:DIABETESYes -0.0192936     # No interaction between BMI and diabetes

# ========================================================================

# Model 4: (Stepwise Selection, "both direction" )
#.    Algorithm chose to remove BPMEDS and educ but kept all other variables.
#    hence, it found most efficient combination of predictors without losing much predictive power.

??framingham


library(pROC)

# we have some regression models now, and want to see how well models has fitted the data.
# Area under Curve: 

calculate_auc <- function(model, data) {
  predictions <- predict(model, type = "response")
  roc_obj <- roc(data$CVD, predictions)
  return(auc(roc_obj))
}

auc_model1 <- calculate_auc(model1, model_data)
auc_model2 <- calculate_auc(model2, model_data)
auc_model3 <- calculate_auc(model3, model_data)
auc_model4 <- calculate_auc(model4, model_data)


auc_model1
auc_model2
auc_model3
auc_model4


# Saml AIC og BIC for alle modeller
model_comparison <- data.frame(
  Model = c("Core", "Comprehensive", "Interactions", "Stepwise"),
  AIC = c(AIC(model1), AIC(model2), AIC(model3), AIC(model4)),
  BIC = c(BIC(model1), BIC(model2), BIC(model3), BIC(model4))
)

model_comparison
# model_4 < model_2 < model_3 < model_1 
#  ==> Model 4 performs statistically best (Lower AIC than other models)

# AUC VALUES:
# Model 1 (Core):          0.7358
# Model 2 (Comprehensive): 0.7614
# Model 3 (Interactions):  0.7624  # HIGHEST AUC but with only small difference to model 2 and 4.
# Model 4 (Stepwise):      0.7596
# 
# # AIC VALUES:
# Model 1 (Core):          3519.4
# Model 2 (Comprehensive): 3373.3
# Model 3 (Interactions):  3373.6
# Model 4 (Stepwise):      3371.0  # LOWEST AIC

# ______________________________________________________

# What we have for now:


# OPTION A: Choose Model 3 (Interactions)
# Pros: Best predictive power, includes clinically important interaction
# Cons: Slightly more complex

# OPTION B: Choose Model 4 (Stepwise) 
# Pros: Most statistically efficient, simpler
# Cons: Slightly lower AUC, misses SYSBP*BPMEDS interaction

# OPTION C: Create HYBRID model where we use variables selected by model.4 as base!
# and we add significant interactions terms from model 3.
# namely: 
# SYSBP:BPMEDSYes -0.0186115  0.0078743  -2.364  0.01810 

model_G <- glm(CVD ~ SEX + TOTCHOL + AGE + SYSBP * BPMEDS + 
                     CIGPDAY + BMI + DIABETES + HEARTRTE + 
                     PREVMI + PREVHYP, 
                   data = model_data, family = "binomial")

auc_model_G <- calculate_auc(model_G, model_data)
aic_model_G <- AIC(model_G)

print(aic_model_G)
aic_model_G

# Comparison:

comparison <- data.frame(
  Model = c("Model 1", "Model 2", "Model 3", "Model 4", "Model_G"),
  AUC = c(auc_model1, auc_model2, auc_model3, auc_model4, auc_model_G),
  AIC = c(AIC(model1), AIC(model2), AIC(model3), AIC(model4), aic_model_G)
)
print(comparison)

print(comparison)
# Model       AUC      AIC
# 1 Model 1 0.7358162 3519.366
# 2 Model 2 0.7613803 3373.301
# 3 Model 3 0.7623736 3373.587
# 4 Model 4 0.7596172 3371.040
# 5 Model_G 0.7599056 3369.634. <- (Best_overall)

# This approach has successfully achieved:

# 1) Best statistical efficiency (lowest AIC)

# 2) Strong predictive power (AUC nearly identical to best model)
# 
# 3) Clinical relevance (includes important interaction we got from m3)
# 
# 4) Interpretability (removes confusing non-significant variables)

# ========================================

summary(model_G)
exp(confint(model_G))

#                      2.5 %        97.5 %
# (Intercept)     5.333582e-04  4.817659e-03  # Reference odds
# SEXFemale       3.582315e-01  5.183433e-01   Women: 64-66% LOWER odds vs men
# TOTCHOL         1.002322e+00  1.006295e+00   Each unit: 0.2-0.6% HIGHER odds
# AGE             1.034693e+00  1.057511e+00   Each year: 3.5 - 5.8% HIGHER odds 
# SYSBP           1.009771e+00  1.021297e+00   Each unit: 1.0 -2.1% HIGHER odds
# BPMEDSYes       1.538955e+00  2.763203e+02   Very wide CI - interpret with interaction
# CIGPDAY         1.003281e+00  1.018196e+00   Each cig: 0.3-1.8% HIGHER odds
# BMI             1.009575e+00  1.053989e+00   1-4.5% higher odds pr. unit
# DIABETESYes     2.018950e+00  5.044222e+00   Diabetes: 2.0-5.0x HIGHER odds
# HEARTRTE        9.833264e-01  9.978323e-01   Each unit: 1.2-1.7% LOWER odds
# PREVMIYes       2.120243e+01  4.276345e+02   Prev heart attack: 21-428x HIGHER odds
# PREVHYPYes      1.040886e+00  1.693874e+00   Prev hypertension: 1.0-1.7x HIGHER odds
# SYSBP:BPMEDSYes 9.673150e-01  9.977220e-01  -> blood pressur medics. reduces the risk associated with high SYSBP by 0.3-3.3%



# Final model (Model_G) odds ratios revealed several strong predictors of CVD
 # Most notably is that pre myocardial infarction conferred dramatically increased odds (OR = 21.2, 95% CI: 21.2-427.6)
 # While female sex was strongly decreased risk (OR = 0.44, 95% CI: 0.36-0.52)
 
# the significant SYSBP×BPMEDS interaction (OR = 0.98, 95% CI: 0.97-1.00) indicates blood pressure medication reduces heart disease risk.
 #  Konfidensintervallet (0.967-0.998) ligger UNDER 1.00
# → Interaktionen er statistisk signifikant (p < 0.05)


#Model 5 udvikledes gennem en two-staged tilgang: først identificerede stepwise selection de mest informative prædiktorer,
# hvorefter klinisk relevante interaktioner blev testet og inkluderet baseret på statistisk signifikans og faglig begrundelse


# _________________________________________________





# _________________________________________________



# Load required packages
library(car)    # For VIF and regression diagnostics
library(ggplot2) # For advanced plotting
library(pROC)   # For ROC curves

# Variance Inflation Factors
#
vif_results <- vif(model_G)
print("Variance Inflation Factors:")
print(vif_results)

# Interpretation: VIF > 5-10 indicates problematic multicollinearity

# 2. INFLUENTIAL OBSERVATIONS - Cook's Distance

cooks_d <- cooks.distance(model_G)

# Identify influential observations (Cook's D > 4/n)
n_obs <- nrow(model_data)
influential_cutoff <- 4/n_obs
influential_obs <- which(cooks_d > influential_cutoff)

 round(influential_cutoff, 4)
    length(influential_obs)

# RESIDUAL ANALYSIS

  
# Create comprehensive diagnostic plots
par(mfrow = c(2, 3))  # Set up 2x3 plot grid

# Plot 1: Residuals vs Fitted
plot(model_G, which = 1, main = "Residuals vs Fitted")

# Plot 2: Q-Q Plot for normality
plot(model_G, which = 2, main = "Q-Q Plot")

# Plot 3: Scale-Location Plot
plot(model_G, which = 3, main = "Scale-Location")

# Plot 4: Cook's Distance
plot(model_G, which = 4, main = "Cook's Distance")

# Plot 5: Residuals vs Leverage
plot(model_G, which = 5, main = "Residuals vs Leverage")

# Plot 6: Histogram of residuals
residuals <- resid(model_G, type = "pearson")
hist(residuals, main = "Distribution of Residuals", 
     xlab = "Pearson Residuals", col = "lightblue")

par(mfrow = c(1, 1))  # Reset plot layout

# 4. MODEL FIT STATISTICS

# goodness of fit
library(ResourceSelection)
hl_test <- hoslem.test(model_data$CVD, fitted(model_G), g = 10)

hl_test$statistic
hl_test$p.value

# PREDICTIVE PERFORMANCE

# AUC with confidence interval
roc_obj <- roc(model_data$CVD, predict(model_G, type = "response"))
auc_ci <- ci.auc(roc_obj)
auc(roc_obj)
auc_ci[1]
auc_ci[3]

# Classification table at different thresholds
predictions <- ifelse(predict(model_G, type = "response") > 0.5, 1, 0)
conf_matrix <- table(Predicted = predictions, Actual = model_data$CVD)
print(conf_matrix)
# ============================================
# =============================================
# ============================================

#_______INTERACTION EFFECT VISUALIZATION - CORRECTED VERSION

# Check the actual factor levels in your model

levels(model_data$BPMEDS)
levels(model_data$SEX)
levels(model_data$DIABETES)
levels(model_data$PREVMI)
levels(model_data$PREVHYP)

# Create plot data with CORRECT factor levels
interaction_plot_data <- expand.grid(
  SYSBP = seq(min(model_data$SYSBP), max(model_data$SYSBP), length = 50),
  BPMEDS = factor(c("No", "Yes"), levels = levels(model_data$BPMEDS)),  # Use actual levels
  AGE = median(model_data$AGE),
  SEX = factor("Male", levels = levels(model_data$SEX)),
  TOTCHOL = median(model_data$TOTCHOL),
  CIGPDAY = median(model_data$CIGPDAY),
  BMI = median(model_data$BMI),
  DIABETES = factor("No", levels = levels(model_data$DIABETES)),
  HEARTRTE = median(model_data$HEARTRTE),
  PREVMI = factor("No", levels = levels(model_data$PREVMI)),
  PREVHYP = factor("No", levels = levels(model_data$PREVHYP))
)

# Predict probabilities
interaction_plot_data$pred_prob <- predict(model_G, 
                                           newdata = interaction_plot_data, 
                                           type = "response")

# Create the interaction plot
library(ggplot2)
interaction_plot <- ggplot(interaction_plot_data, aes(x = SYSBP, y = pred_prob, color = BPMEDS)) +
  geom_line(size = 1.5) +
  labs(title = "SYSBP*BPMEDS Interaction Effect on CVD Risk",
       x = "Systolic Blood Pressure (mmHg)",
       y = "Predicted Probability of CVD",
       color = "Blood Pressure Medication") +
  theme_minimal() +
  scale_color_manual(values = c("No" = "red", "Yes" = "blue")) +
  theme(legend.position = "bottom")

print(interaction_plot)

# Alternative simple version if you still get errors:


# Simple approach using actual data
simple_plot_data <- model_data
simple_plot_data$pred_prob <- predict(model_G, type = "response")

ggplot(simple_plot_data, aes(x = SYSBP, y = pred_prob, color = BPMEDS)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "glm", method.args = list(family = "binomial"), se = FALSE) +
  labs(title = "SYSBP*BPMEDS Interaction (Actual Data)",
       x = "Systolic Blood Pressure",
       y = "Predicted CVD Probability",
       color = "BP Medication") +
  theme_minimal()





# Create simple grid with only the variables we need now
simple_grid <- expand.grid(
  SYSBP = seq(100, 200, length = 30),  # Reasonable SYSBP range
  BPMEDS = factor(c("No", "Yes")),
  SEX = factor("Male"),
  AGE = 55,  # Typical age
  TOTCHOL = 240,
  CIGPDAY = 0,
  BMI = 26,
  DIABETES = factor("No"),
  HEARTRTE = 75,
  PREVMI = factor("No"), 
  PREVHYP = factor("No")
)

# Ensure factor levels match the original data
levels(simple_grid$BPMEDS) <- levels(model_data$BPMEDS)
levels(simple_grid$SEX) <- levels(model_data$SEX)
levels(simple_grid$DIABETES) <- levels(model_data$DIABETES)
levels(simple_grid$PREVMI) <- levels(model_data$PREVMI)
levels(simple_grid$PREVHYP) <- levels(model_data$PREVHYP)

# Predict
simple_grid$pred_prob <- predict(model_G, newdata = simple_grid, type = "response")

# Plot
ggplot(simple_grid, aes(x = SYSBP, y = pred_prob, color = BPMEDS)) +
  geom_line(size = 1.5) +
  labs(title = "Effect of Blood Pressure Medication on CVD Risk",
       subtitle = "Higher SYSBP increases CVD risk, but medication reduces this effect",
       x = "Systolic Blood Pressure (mmHg)",
       y = "Predicted Probability of CVD",
       color = "Taking BP Medication") +
  theme_minimal() +
  scale_color_manual(values = c("No" = "red", "Yes" = "#3dd")) # tilfældigt tal!


# 7. FINAL SUMMARY

length(influential_obs)
hl_test$p.value
 round(auc(roc_obj), 3)




#################
##########
#############
##########












# ULTRA SIMPLE INTERACTION PLOT
# Create simple grid with only the variables we need
simple_grid <- expand.grid(
  SYSBP = seq(100, 200, length = 30),  # Reasonable SYSBP range
  BPMEDS = factor(c("No", "Yes")),
  SEX = factor("Male"),
  AGE = 55,  # Typical age for denne dsta
  TOTCHOL = 240,
  CIGPDAY = 0,
  BMI = 26,
  DIABETES = factor("No"),
  HEARTRTE = 75,
  PREVMI = factor("No"), 
  PREVHYP = factor("No")
)

# Ensure factor levels match the original data
levels(simple_grid$BPMEDS) <- levels(model_data$BPMEDS)
levels(simple_grid$SEX) <- levels(model_data$SEX)
levels(simple_grid$DIABETES) <- levels(model_data$DIABETES)
levels(simple_grid$PREVMI) <- levels(model_data$PREVMI)
levels(simple_grid$PREVHYP) <- levels(model_data$PREVHYP)

# Predict
simple_grid$pred_prob <- predict(model_G, newdata = simple_grid, type = "response")

# Plot
ggplot(simple_grid, aes(x = SYSBP, y = pred_prob, color = BPMEDS)) +
  geom_line(size = 1.5) +
  labs(title = "Effect of Blood Pressure Medication on CVD Risk",
       subtitle = "Higher SYSBP increases CVD risk, but medication reduces this effect",
       x = "Systolic Blood Pressure (mmHg)",
       y = "Predicted Probability of CVD",
       color = "Taking BP Medication") +
  theme_minimal() +
  scale_color_manual(values = c("No" = "purple", "Yes" = "white"))




















==========================================================






library(splines)








# Model 1: Original (lineær for kontinuerte variable)
model_original <- glm(
  CVD ~ AGE + SEX + TOTCHOL + SYSBP + PREV_CVD_EVENT_VMI_CHD + 
    CIGPDAY + DIABETES + BMI + HEARTRTE + PREVHYP,
  family = binomial(link = "logit"),
  data = training_data
)

# Model 2: Splines for kontinuerte variable (undtagen CIGPDAY)
model_splines <- glm(
  CVD ~ ns(AGE, df = 3) + SEX + ns(TOTCHOL, df = 3) + 
    ns(SYSBP, df = 3) + PREV_CVD_EVENT_VMI_CHD + 
    CIGPDAY + DIABETES + ns(BMI, df = 3) + 
    ns(HEARTRTE, df = 3) + PREVHYP,
  family = binomial(link = "logit"),
  data = training_data
)

####

# håndtering af CIGPDAY
log_model_splines_virknu <- glm(
  CVD ~ ns(AGE, df = 3) + SEX + ns(TOTCHOL, df = 3) + 
    ns(SYSBP, df = 3) + PREV_CVD_EVENT_VMI_CHD + 
    # Bedre håndtering af CIGPDAY:
    I(CIGPDAY == 0) + ns(CIGPDAY[CIGPDAY > 0], df = 2) +  # Separat behandling af 0 og >0
    DIABETES + ns(BMI, df = 3) + 
    ns(HEARTRTE, df = 3) + PREVHYP,
  family = binomial(link = "logit"),
  data = training_data
)




log_model_splines_categorical <- glm(
  CVD ~ ns(AGE, df = 3) + SEX + ns(TOTCHOL, df = 3) + 
    ns(SYSBP, df = 3) + PREV_CVD_EVENT_VMI_CHD + 
    smoking_cat + DIABETES + ns(BMI, df = 3) + 
    ns(HEARTRTE, df = 3) + PREVHYP,
  family = binomial(link = "logit"),
  data = training_data
)



library(splines)

# Forbered data: opret nye variabler for splines
training_data_prepared <- training_data %>%
  mutate(
    smoking_indicator = as.numeric(CIGPDAY > 0),  # 1 for rygere, 0 for ikke-rygere
    cigpday_positive = ifelse(CIGPDAY > 0, CIGPDAY, NA)  # Kun positive værdier
  )

#
log_model_splines_corrected <- glm(
  CVD ~ ns(AGE, df = 3) + SEX + ns(TOTCHOL, df = 3) + 
    ns(SYSBP, df = 3) + PREV_CVD_EVENT_VMI_CHD + 
    smoking_indicator + ns(cigpday_positive, df = 2) + 
    DIABETES + ns(BMI, df = 3) + 
    ns(HEARTRTE, df = 3) + PREVHYP,
  family = binomial(link = "logit"),
  data = training_data_prepared
)



log_model_splines_simple <- glm(
  CVD ~ ns(AGE, df = 3) + SEX + ns(TOTCHOL, df = 3) + 
    ns(SYSBP, df = 3) + PREV_CVD_EVENT_VMI_CHD + 
    # Simpel håndtering af CIGPDAY:
    I(CIGPDAY == 0) + I(log(CIGPDAY + 1)) + 
    DIABETES + ns(BMI, df = 3) + 
    ns(HEARTRTE, df = 3) + PREVHYP,
  family = binomial(link = "logit"),
  data = training_data
)





training_data <- training_data %>%
  mutate(
    smoking_cat = case_when(
      CIGPDAY == 0 ~ "Non_smoker",
      CIGPDAY <= 10 ~ "Light_1_10",
      CIGPDAY <= 20 ~ "Moderate_11_20", 
      TRUE ~ "Heavy_20_plus"
    ),
    smoking_cat = factor(smoking_cat)
  )



log_model_splines_final <- glm(
  CVD ~ ns(AGE, df = 3) + SEX + ns(TOTCHOL, df = 3) + 
    ns(SYSBP, df = 3) + PREV_CVD_EVENT_VMI_CHD + 
    smoking_cat + DIABETES + ns(BMI, df = 3) + 
    ns(HEARTRTE, df = 3) + PREVHYP,
  family = binomial(link = "logit"),
  data = training_data
)





# Tjek om modellen konvergerede
summary(log_model_splines_final)

}

# Sammenlign med original model
original_vs_splines <- anova(
  log_model_extend_2b_HRT_PPP_COMB_VMI_CHD, 
  log_model_splines_final, 
  test = "Chisq"
)

print(original_vs_splines)


# Forbered validation data
validation_data_prepared <- validation_data %>%
  mutate(
    smoking_cat = case_when(
      CIGPDAY == 0 ~ "Non_smoker",
      CIGPDAY <= 10 ~ "Light_1_10",
      CIGPDAY <= 20 ~ "Moderate_11_20",
      TRUE ~ "Heavy_20_plus"
    ),
    smoking_cat = factor(smoking_cat)
  )

# Forbered validation_data på samme måde som training_data
validation_data_prepared <- validation_data %>%
  mutate(
    PREV_CVD_EVENT_VMI_CHD = ifelse(PREVMI == 1 | PREVCHD == 1 | PREVAP == 1, 1, 0),
    smoking_cat = case_when(
      CIGPDAY == 0 ~ "Non_smoker",
      CIGPDAY <= 10 ~ "Light_1_10",
      CIGPDAY <= 20 ~ "Moderate_11_20",
      TRUE ~ "Heavy_20_plus"
    ),
    smoking_cat = factor(smoking_cat, levels = levels(training_data$smoking_cat))
  )

# Tjek  nødvendige variabler er på plads
names(validation_data_prepared)

# Forudsigelser og validering
splines_predictions <- predict(log_model_splines_final, 
                               newdata = validation_data_prepared, 
                               type = "response")

splines_auc <- roc(validation_data_prepared$CVD, splines_predictions)$auc

cat("Splines model validation AUC:", round(splines_auc, 3), "\n")

vars_required <- all.vars(formula(log_model_splines_final))
missing_var <- setdiff(vars_required, names(validation_data))
if(length(missing_var)>0){
  print(missing_var)
}




validation_data_prepared <- validation_data %>%
  mutate(
    PREV_CVD_EVENT_VMI_CHD = ifelse(PREVMI == 1 | PREVCHD == 1 | PREVAP == 1, 1, 0),
    smoking_cat = case_when(
      CIGPDAY == 0 ~ "Non_smoker",
      CIGPDAY <= 10 ~ "Light_1_10",
      CIGPDAY <= 20 ~ "Moderate_11_20",
      TRUE ~ "Heavy_20_plus"
    ),
    smoking_cat = factor(smoking_cat, levels = levels(training_data$smoking_cat))
  )

