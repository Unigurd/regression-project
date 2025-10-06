library(riskCommunicator)
library(dplyr)
library(ggplot2)
library(broom)
library(skimr)
library(GGally)
library(corrplot)

# Regression Project: Framingham Heart Study
setwd("C:/Users/SIMON/OneDrive/Dokumente/Uni/MQE/2025_26 WiSe (KU)/Regression/Assignment")
data <- read.csv("training_data.csv")

# Explanatory data analysis of the training data
# Check for:
# extreme observations / outliers
# missing values
# skewness / asymmetry in distributions

skim(data) 
# this reveals that all variables are numeric. 
# transform educ and PERIOD into factors
data <- data %>% 
  mutate(educ = as.factor(educ),
         PERIOD = as.factor(PERIOD))

# I would keep binary variables numeric and not transform into factors, but
# I think that, ultimately, it does not matter too much
# NAs:
# (relevant amount of) NAs in HDLC, LDLC, Glucose, Totchol, BPMEDS, educ.
# Maybe also cigpday, BMI are relevant NAs --> check whether there is any information about NAs in data description
# The completeness rate is still VERY high for these, so I would not be too worried about deleting them in doubt
# no NAs in outcome variable!

# check marginal distribution of outcome variable CVD 
table(data$CVD)
table(is.na(data$CVD)) 

table(data$PERIOD)
# create subsets of participants in each examination cycle to check to which extent they coincide
wave1 <- subset(data, PERIOD == 1)
wave2 <- subset(data, PERIOD == 2)
wave3 <- subset(data, PERIOD == 3)
length(intersect(wave2$RANDID, wave1$RANDID))
length(intersect(wave3$RANDID, wave1$RANDID)) 
length(intersect(wave3$RANDID, wave2$RANDID))
# can compare these numbers to total number of individuals
length(unique(data$RANDID)) # 4377

# can create a subset with only the first observation for each individual
first_obs <- data %>%
  group_by(RANDID) %>%
  arrange(PERIOD) %>%   # make sure periods are ordered
  slice(1) %>%          # take the first row per ID
  ungroup()

# want to know for how long a person was observed
hist(data$TIME)

# variable pre-screening
# Check correlation between possibly strongly correlated variabels to determine whether it is unwise to include both
cor(data$SYSBP, data$DIABP, use = "complete.obs") # not include both. cor = 0.71. Standard is to include Systolic
cor(data$CURSMOKE, data$CIGPDAY, use = "complete.obs")  # only include cigpday and not cursmoke (the information of which is contained in cigpday)
cor(data$HDLC, data$LDLC, use = "complete.obs") # low cor (-0.12). If we decide to use it, we should create a single DLC variable with 3 factors.
# Since we don't plan on using DLC (only available in period 3), I don't do that at this point


# list of plausible predictors
predictors <- c("SEX", "AGE", "educ", "TOTCHOL", "SYSBP", "CIGPDAY", "BMI", "DIABETES", "BPMEDS", "HEARTRTE", "GLUCOSE",
                "PREVCHD", "PREVAP", "PREVMI", "PREVHYP")
all_vars <- c("CVD", predictors)
# does not include: 
# indicators that show indication of diseases during follow up (these are no predictors!)
# Time variables
# DIABP because we have SYSBP, CURSMOKE because we have CIGPDAY, HDLC & LDLC because they only appear in period 3
# stroke because there is so little people with strokes

# check again for correlations between included predictors (multicollinearity)
# Make a numeric matrix for correlation (convert factors to numeric if necessary)
corr_data <- data %>%
  select(all_of(all_vars)) %>%
  mutate(across(everything(), ~ as.numeric(.)))  # ensures all columns are numeric

# Compute Spearman correlation
cp <- cor(
  data.matrix(corr_data), 
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

# insights from this:
# shows particularly strong correlations between PREVAP, PREVCHD, and PREVMI
cor(corr_data$PREVAP, corr_data$PREVCHD) # 0.85
cor(corr_data$PREVAP, corr_data$PREVMI) # 0.35
cor(corr_data$PREVMI, corr_data$PREVCHD) # 0.65
# also noticable correlation to CVD !
cor(corr_data$CVD, corr_data$PREVCHD)
# as well as between PREVHYP and SYSBP
cor(corr_data$PREVHYP, corr_data$SYSBP) # 0.67

table(data$PREVSTRK)

