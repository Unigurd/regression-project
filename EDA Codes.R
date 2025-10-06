### Install & load packages
# Run install.packages(...)
install.packages(c("riskCommunicator", "tidyverse", "data.table",
                   "summarytools", "GGally", "naniar", "tableone",
                   "corrplot", "scales", "survival", "survminer", "mice"))
# Load libraries
library(riskCommunicator)
library(tidyverse)
library(data.table)
library(summarytools)
library(GGally)
library(naniar)
library(tableone)
library(corrplot)
library(scales)
library(survival)
library(survminer)
library(mice)

# #Set seed and create output folder (if not exists)
set.seed(2025)
dir.create("output", showWarnings = FALSE)

### Load framingham data
data("framingham", package = "riskCommunicator")  # loads object named 'framingham'
df <- as_tibble(framingham)   # convert to tibble for tidyverse
# Quick checks:
dim(df)       # number of rows and columns (should be ~4434 rows)
names(df)     # column names - inspect these to know variables available
glimpse(df)   # types and sample values
# Quick overview & variable inventory
dim(framingham)
# List all variable names
names(framingham)
# Glimpse at variable types and first few values
glimpse(framingham)
# Count missing values per variable
colSums(is.na(framingham))

# Basic cleaning — types & a few derived variables
df <- framingham %>%
  mutate(
    # Convert binary to factors with labels
    SEX       = factor(SEX, levels = c(0,1), labels = c("Female","Male")),
    CURSMOKE  = factor(CURSMOKE, levels = c(0,1), labels = c("No","Yes")),
    DIABETES  = factor(DIABETES, levels = c(0,1), labels = c("No","Yes")),
    BPMEDS    = factor(BPMEDS, levels = c(0,1), labels = c("No","Yes")),
    PREVCHD   = factor(PREVCHD, levels = c(0,1), labels = c("No","Yes")),
    PREVAP    = factor(PREVAP, levels = c(0,1), labels = c("No","Yes")),
    PREVMI    = factor(PREVMI, levels = c(0,1), labels = c("No","Yes")),
    PREVSTRK  = factor(PREVSTRK, levels = c(0,1), labels = c("No","Yes")),
    PREVHYP   = factor(PREVHYP, levels = c(0,1), labels = c("No","Yes")),
    DEATH     = factor(DEATH, levels = c(0,1), labels = c("Alive","Dead")),
    ANGINA    = factor(ANGINA, levels = c(0,1), labels = c("No","Yes")),
    HOSPMI    = factor(HOSPMI, levels = c(0,1), labels = c("No","Yes")),
    MI_FCHD   = factor(MI_FCHD, levels = c(0,1), labels = c("No","Yes")),
    ANYCHD    = factor(ANYCHD, levels = c(0,1), labels = c("No","Yes")),
    STROKE    = factor(STROKE, levels = c(0,1), labels = c("No","Yes")),
    CVD       = factor(CVD, levels = c(0,1), labels = c("No","Yes")),
    HYPERTEN  = factor(HYPERTEN, levels = c(0,1), labels = c("No","Yes"))
  )
# Quick check after conversion
glimpse(df)

### Missing data exploration
missing_summary <- tibble(
  variable = names(df),
  n_missing = colSums(is.na(df)),
  pct_missing = round(100 * colSums(is.na(df)) / nrow(df), 1)
) %>%
  arrange(desc(pct_missing))
print(missing_summary, n = 15)  # show top 15
# Visualize missingness pattern
library(naniar)
# Overall missingness heatmap
vis_miss(df)
# Mosaic view: where multiple vars are missing together
gg_miss_upset(df %>% select(TOTCHOL, CIGPDAY, BMI, GLUCOSE, HDLC, LDLC, educ))

##Univariate EDA (Numeric variables)
num_vars <- c("AGE","TOTCHOL","SYSBP","DIABP","BMI","HEARTRTE","GLUCOSE","TIMECVD")
summary(df[num_vars])
# Histograms
for (v in num_vars) {
  p <- ggplot(df, aes_string(x = v)) +
    geom_histogram(bins = 40, fill = "steelblue", color = "white") +
    theme_minimal() +
    ggtitle(paste("Distribution of", v))
  print(p)
}
# Select categorical variables of interest
cat_vars <- c("CURSMOKE","CVD","DIABETES","BPMEDS","PREVCHD","PREVAP",
              "PREVMI","PREVSTRK","PREVHYP","HYPERTEN","DEATH")
# Frequency tables and barplots
for (v in cat_vars) {
  tab <- df %>% count(!!sym(v)) %>% mutate(freq = n/sum(n))
  print(tab)
    p <- ggplot(tab, aes_string(x = v, y = "n")) +
    geom_col(fill = "steelblue") +
    ggtitle(paste("Counts of", v)) +
    theme_minimal()
  print(p)
}

### Explore outcome relationships (CVD vs predictors)
# 1. CVD rate by categorical predictors
# CVD by smoking
df %>% group_by(CURSMOKE, CVD) %>%
  summarise(n = n()) %>%
  mutate(prop = n/sum(n)) %>%
  print()
# CVD by diabetes
df %>% group_by(DIABETES, CVD) %>%
  summarise(n = n()) %>%
  mutate(prop = n/sum(n))
# CVD by hypertension
df %>% group_by(HYPERTEN, CVD) %>%
  summarise(n = n()) %>%
  mutate(prop = n/sum(n))
# 2. Continuous predictors vs CVD (boxplots)
# Age vs CVD
ggplot(df, aes(x = CVD, y = AGE, fill = CVD)) +
  geom_boxplot() +
  theme_minimal() +
  ggtitle("Age distribution by CVD outcome")
# Cholesterol vs CVD
ggplot(df, aes(x = CVD, y = TOTCHOL, fill = CVD)) +
  geom_boxplot() +
  theme_minimal() +
  ggtitle("Total cholesterol by CVD outcome")
# Systolic BP vs CVD
ggplot(df, aes(x = CVD, y = SYSBP, fill = CVD)) +
  geom_boxplot() +
  theme_minimal() +
  ggtitle("Systolic BP by CVD outcome")
# 3. Smoothed probability curves
# Logistic smooth for age vs CVD
ggplot(df, aes(x = AGE, y = as.numeric(CVD=="Yes"))) +
  geom_jitter(height=0.02, alpha=0.2) +
  geom_smooth(method="glm", method.args=list(family="binomial")) +
  ggtitle("Probability of CVD vs Age") +
  theme_minimal()
# Logistic smooth for systolic BP
ggplot(df, aes(x = SYSBP, y = as.numeric(CVD=="Yes"))) +
  geom_jitter(height=0.02, alpha=0.2) +
  geom_smooth(method="glm", method.args=list(family="binomial")) +
  ggtitle("Probability of CVD vs Systolic BP") +
  theme_minimal()

### Correlation analysis
# Choose continuous predictors
num_vars <- c("AGE","TOTCHOL","SYSBP","DIABP","BMI","HEARTRTE","GLUCOSE")
# Compute correlation matrix
cor_mat <- cor(df[num_vars], use="pairwise.complete.obs")
# Show numeric matrix
round(cor_mat, 2)
# Plot heatmap
library(corrplot)
corrplot(cor_mat, method="color", type="lower", tl.cex=0.7)

###Survival overview (Kaplan–Meier)
# Build survival object
library(survival)
library(survminer)
surv_obj <- Surv(time = df$TIMECVD, event = (df$CVD=="Yes"))
# KM curves by diabetes
fit_km <- survfit(surv_obj ~ DIABETES, data = df)
ggsurvplot(fit_km, data = df, conf.int = TRUE, pval = TRUE,
           risk.table = TRUE, ggtheme = theme_minimal(),
           legend.title="Diabetes status")
