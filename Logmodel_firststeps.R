# Simon's modelling approaches

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
  addrect = 6,        # draws rectangles around clusters
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
log_model_simple <- glm(CVD ~ AGE + SEX + TOTCHOL + SYSBP + PREVCHD, family = binomial(link = "logit"), data = data)
model_assess(log_model_simple)

# Model 2
# adding some more variabls: CIGPDAY, DIBETES, BMI
log_model_extend <- glm(CVD ~ AGE + SEX + TOTCHOL + SYSBP + PREVCHD + CIGPDAY + DIABETES + BMI, family = binomial(link = "logit"), data = data)
model_assess(log_model_extend)

### Comparing model 1 and model 2
# since the models are nested, we can use anova (likelihood ratio test) (AIC comparison leads to the same result)
anova(log_model_simple, log_model_extend, test = "Chisq")
# this shows that the extension makes the model fit significantly better
# Insights from the coefficients:
# PREVCHD, DIABETES, and SEX are strongest predictors of CVD
# The other variables are also significant, but seem to be less meaningful

# Model 3
# use all the variables available
log_model_super_extend <- glm(CVD ~ ., family = binomial(link = "logit"), data = data)
model_assess(log_model_super_extend) 
# graph is useless because of the big CI of PREVMI: There are probably too few observations, so estimator is very uncertain
# similar for CVD in this case --> do not include these two variables 

# recommend to stick with extended model!
