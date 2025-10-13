###########################################################
#                 framingham data
###########################################################
#          - DEL 1: (EDA) EKSPLORATIV DATAANALYSE
# ___________________________________________________________
# Formål: Systematisk eksplorere datasættet før modellering
# Følger opgavens struktur:
# 1. Marginal distributions
# 2. Missing values
# 3. Selection of potential predictors
# 4. Pairwise associations
# 5. Marginal associations with response
# ___________________________________________________________


library(riskCommunicator) # Framingham data
library(dplyr)   # For at manipulere data (data_mamipulation)
library(ggplot2)   # Visualisering & plots.
library(tidyr)   # Data reshaping
library(Hmisc)  # Korrelation w. p-valus ()
library(corrplot)   # Korrelationsmatrix plots
library(psych)   # Deskriptiv statistik 
library(naniar)     # Missing data visualisering
library(GGally)     # Pairwise plots
library(RwR)    # For bike dataset
library(ggplot2)    # For plotting
library(broom)   # For tidy model output
library(dplyr)      # For data manipulation
library(MASS)
library(splines)    # For spline models
#install.packages("GGally")
# Load data
data(framingham, package = "riskCommunicator")

# Grundlæggende info
dim(framingham)            # Dimensioner: rækker x kolonner
str(framingham)            # Struktur: datatyper
head(framingham, 10)       # Første 10 rækker

# ============================================================================
# MARGINAL DISTRIBUTIONS OF VARIABLES INCLUDED
# ==================   ==========================================================

# --------------------------------------------------------------------------------------
#  Identificer variabeltyper
# ----------------------------------------------------------------------------
framingham
# Kontinuerte prædiktorer
continuous_vars <- c("AGE", "TOTCHOL", "SYSBP", "DIABP", "CIGPDAY", 
                     "BMI", "HEARTRTE", "GLUCOSE", "HDLC", "LDLC")

# Binære/kategoriske prædiktorer
binary_vars <- c("SEX", "CURSMOKE", "DIABETES", "BPMEDS", "CVD", "DEATH",
                 "PREVCHD", "PREVAP", "PREVMI", "PREVSTRK", "PREVHYP",
                 "ANGINA", "HOSPMI", "MI_FCHD", "ANYCHD", "STROKE", "HYPERTEN")

# Ordinal variabel
ordinal_vars <- c("educ")

# Study design variable
design_vars <- c("RANDID", "PERIOD", "TIME")

# Time-to-event variable
time_vars <- c("TIMECVD", "TIMEDTH", "TIMEMI", "TIMEMIFC", "TIMECHD", 
               "TIMESTRK", "TIMEAP", "TIMEHYP")

# ----------------------------------------------------------------------------
# Deskriptiv statistik - Kontinuerte variable
# ---------------------------------------------------------------------------------..-

# Detaljeret deskriptiv statistik
desc_continuous <- describe(framingham[continuous_vars])
print(desc_continuous)

# Summary statistik (base R)
summary(framingham[continuous_vars])

# Lav egen tabel med nøgletal (se om det giver bedre mening -> kunne være interessant, hvis vi så et mønstre og fokuserede på det)
desc_table <- framingham %>%
  select(all_of(continuous_vars)) %>%
  summarise(across(everything(), list(
    n = ~sum(!is.na(.)),                    # Antal ikke-missing
    mean = ~mean(., na.rm = TRUE),          # Gennemsnit
    sd = ~sd(., na.rm = TRUE),              # Standardafvigelse
    min = ~min(., na.rm = TRUE),            # Minimum
    q25 = ~quantile(., 0.25, na.rm = TRUE), # 1. kvartil
    median = ~median(., na.rm = TRUE),      # Median
    q75 = ~quantile(., 0.75, na.rm = TRUE), # 3. kvartil
    max = ~max(., na.rm = TRUE),            # Maximum
    missing = ~sum(is.na(.))                # Antal missing
  ))) %>%
  pivot_longer(everything()) %>%
  separate(name, into = c("Variable", "Statistic"), sep = "_(?=[^_]+$)") %>%
  pivot_wider(names_from = Statistic, values_from = value)

print(desc_table)

# ---------.- .----------------------------------------------------------__---_---------
#  Visualisering - Histogrammer for kontinuerte variable
# ---------  -------------------------------------------------------------------____----

# Histogrammer i facets
fig1_histograms <- framingham %>%
  select(all_of(continuous_vars)) %>%
  pivot_longer(everything(), names_to = "Variable", values_to = "Value") %>%
  ggplot(aes(x = Value)) +
  geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7, color = "black") +
  facet_wrap(~ Variable, scales = "free", ncol = 3) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    strip.text = element_text(size = 10, face = "bold")
  ) +
  labs(
    title = "Figur 1: Marginalfordelinger - Kontinuerte Variable",
    subtitle = "Alle observationer (n=11,627)",
    x = "Værdi",
    y = "Frekvens"
  )

print(fig1_histograms)


# ----------------------------------------------------------------------------
#  Q-Q plots for normalitetscheck
# ----------------------------------------------------------------------------

# Q-Q plots for at vurdere normalfordeling
par(mfrow = c(3, 4), mar = c(4, 4, 2, 1))
for (var in continuous_vars) {
  qqnorm(framingham[[var]], main = var, cex.main = 1)
  qqline(framingham[[var]], col = "red", lwd = 2)
}
par(mfrow = c(1, 1))

# ----------------------------------------------------------------------------
#  Deskriptiv statistik - Binære/kategoriske variable
# ----------------------------------------------------------------------------

# Frekvens---tabeller for alle binære variable
freq_tables <- lapply(framingham[binary_vars], function(x) {
  tab <- table(x, useNA = "ifany")
  prop <- prop.table(tab) * 100
  data.frame(
    Value = names(tab),
    N = as.numeric(tab),
    Percent = as.numeric(prop)
  )
})

#
print(freq_tables$SEX)

print(freq_tables$CVD)


print(freq_tables$DIABETES)

# Samlet tabel for binære variable ( ---- ----- 0,1,0,1,0k,1,10,1,) 

binary_summary <- framingham %>%
  select(all_of(binary_vars)) %>%
  summarise(across(everything(), list(
    n_total = ~sum(!is.na(.)),
    n_0 = ~sum(. == 0, na.rm = TRUE),
    n_1 = ~sum(. == 1, na.rm = TRUE),
    n_2 = ~sum(. == 2, na.rm = TRUE),  # For SEX
    pct_1 = ~mean(. == 1, na.rm = TRUE) * 100
  ))) %>%
  pivot_longer(everything()) %>%
  separate(name, into = c( ,  ), sep = ) %>%
  pivot_wider(names_from = Statistic, values_from = value)

print(binary_summary)

# ----------------------------------------------------------------------------
#  Ordinal variabel: educ
# ----------------------------------------------------------------------------

# Education level (1-4)
educ_freq <- table(framingham$educ, useNA = "ifany")
educ_prop <- prop.table(educ_freq) * 100

print(data.frame(
  Level = names(educ_freq),
  N = as.numeric(educ_freq),
  Percent = as.numeric(educ_prop)
))

# Barplot
barplot(educ_freq, 
        main = "Figur 2: Uddannelsesniveau (educ)",
        xlab = "Education Level",
        ylab = "Frekvens",
        col = "steelblue")

# ----------------------------------------------------------------------------
# Study design variable: PERIOD og TIME
# ----------------------------------------------------------------------------

# PERIOD fordeling
table(framingham$PERIOD)

# Antal målinger per person
measurements_per_person <- framingham %>%
  group_by(RANDID) %>%
  summarise(n_measurements = n()) %>%
  count(n_measurements, name = "n_persons")

print("Antal målinger per person:")
print(measurements_per_person)

# TIME statistik per PERIOD
time_by_period <- framingham %>%
  group_by(PERIOD) %>%
  summarise(
    n = n(),
    mean_TIME = mean(TIME),
    sd_TIME = sd(TIME),
    min_TIME = min(TIME),
    max_TIME = max(TIME),
    mean_AGE = mean(AGE, na.rm = TRUE)
  )

print("TIME statistik per PERIOD:")
print(time_by_period)

# ============================================================================
# . MISSING VALUES
# ============================================================================

# ----------------------------------------------------------------------------
# Samlet missing values analyse
# ----------------------------------------------------------------------------

# Beregn missing for ALLE variable
missing_summary <- framingham %>%
  summarise(across(everything(), ~sum(is.na(.) | . == "NA" | . == ""))) %>%
  pivot_longer(everything(), names_to = "Variable", values_to = "N_Missing") %>%
  mutate(
    N_Total = nrow(framingham),
    Pct_Missing = round(N_Missing / N_Total * 100, 2)
  ) %>%
  arrange(desc(Pct_Missing)) %>%
  filter(N_Missing > 0)  # Vis kun variable med missing

print(missing_summary, n = 30)

# ----------------------------------------------------------------------------
#  Missing patterns visualisering
# ----------------------------------------------------------------------------

# Vis missing per variabel (barplot)
fig3_missing_bar <- gg_miss_var(framingham) +
  theme_minimal() +
  labs(
    title = "Figur 3: Antal Missing Values per Variabel",
    x = "Variabel",
    y = "Antal Missing"
  )

print(fig3_missing_bar)
# ggsave("fig3_missing_bar.png", width = 10, height = 6, dpi = 300)

# Missing pattern heatmap (subset af variable)
key_vars <- c("CVD", "AGE", "SEX", "TOTCHOL", "SYSBP", "DIABP", "BMI", 
              "DIABETES", "GLUCOSE", "CURSMOKE", "HDLC", "LDLC", "educ")

fig4_missing_pattern <- vis_miss(framingham[key_vars]) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Figur 4: Missing Data Patterns - Key Variables")

print(fig4_missing_pattern)
# ggsave("fig4_missing_pattern.png", width = 10, height = 6, dpi = 300)

# ----------------------------------------------------------------------------
# Missing stratificeret på PERIOD
# ----------------------------------------------------------------------------

# Beregn missing per PERIOD for key variable
missing_by_period <- framingham %>%
  group_by(PERIOD) %>%
  summarise(
    n = n(),
    across(all_of(c("HDLC", "LDLC", "GLUCOSE", "educ", "BPMEDS", "TOTCHOL", "BMI")),
           ~sum(is.na(.)) / n() * 100,
           .names = "pct_missing_{.col}")
  ) %>%
  pivot_longer(
    cols = starts_with("pct_missing"),
    names_to = "Variable",
    values_to = "Pct_Missing",
    names_prefix = "pct_missing_"
  )

print(missing_by_period)

# Visualiser missing over PERIOD
fig5_missing_period <- ggplot(missing_by_period, 
                              aes(x = PERIOD, y = Pct_Missing, 
                                  color = Variable, group = Variable)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(
    title = "Figur 5: Missing Values per PERIOD",
    x = "PERIOD",
    y = "Procent Missing (%)",
    color = "Variabel"
  ) +
  theme(legend.position = "right")

print(fig5_missing_period)

# ----------------------------------------------------------------------------
#  Missing mechanism undersøgelse --> Hvis nødvendigt
# ----------------------------------------------------------------------------

# Er GLUCOSE missing relateret til andre variable?
# Sammenlign karakteristika: GLUCOSE missing vs ikke missing

glucose_missing_analysis <- framingham %>%
  filter(PERIOD == 1) %>%  # Kun baseline for nu
  mutate(GLUCOSE_missing = is.na(GLUCOSE)) %>%
  group_by(GLUCOSE_missing) %>%
  summarise(
    N = n(),
    Mean_Age = mean(AGE, na.rm = TRUE),
    Pct_Male = mean(SEX == 1, na.rm = TRUE) * 100,
    Pct_CVD = mean(CVD, na.rm = TRUE) * 100,
    Pct_Diabetes = mean(DIABETES, na.rm = TRUE) * 100,
    Mean_SYSBP = mean(SYSBP, na.rm = TRUE)
  )


print(glucose_missing_analysis)

# Chi-square test: Er GLUCOSE missing associeret med CVD?
baseline_data <- framingham %>% filter(PERIOD == 1)
chisq_glucose_cvd <- chisq.test(
  table(is.na(baseline_data$GLUCOSE), baseline_data$CVD)
)

print(chisq_glucose_cvd)

# T-test: Er GLUCOSE missing associeret med AGE?
t_test_glucose_age <- t.test(AGE ~ is.na(GLUCOSE), data = baseline_data)
print("T-test: Age difference GLUCOSE missing vs ikke-missing")
print(t_test_glucose_age)

# ============================================================================
#  SELECTION OF POTENTIAL PREDICTORS
# ============================================================================

# ----------------------------------------------------------------------------
# Klassificering baseret på klinisk relevans
# ----------------------------------------------------------------------------

# Opret prædiktor-tabel
predictors_classification <- data.frame(
  Variable = c("AGE", "SEX", "TOTCHOL", "SYSBP", "DIABP", "BMI", 
               "DIABETES", "CURSMOKE", "GLUCOSE", "HDLC", "LDLC",
               "BPMEDS", "educ", "HEARTRTE", "CIGPDAY",
               "PREVCHD", "PREVAP", "PREVMI", "PREVSTRK", "PREVHYP"),
  Type = c("Continuous", "Binary/Categorical", "Continuous", "Continuous", "Continuous",
           "Continuous", "Binary", "Binary", "Continuous", "Continuous", 
           "Continuous", "Binary", "Ordinal", "Continuous", "Continuous",
           "Binary", "Binary", "Binary", "Binary", "Binary"),
  Clinical_Relevance = c("High", "High", "High", "High", "High", "High",
                         "High", "High", "High", "High", "High",
                         "Medium", "Medium", "Low", "Medium",
                         "Medium", "Low", "Low", "Low", "Medium"),
  Known_CVD_Risk_Factor = c("Yes", "Yes", "Yes", "Yes", "Yes", "Yes",
                            "Yes", "Yes", "Yes", "Yes", "Yes",
                            "Unclear", "Social", "No", "Yes",
                            "Prevalent disease", "Prevalent disease", 
                            "Prevalent disease", "Prevalent disease", "Yes"),
  Notes = c("", "1=Male, 2=Female", "", "", "Correlated with SYSBP", "",
            "", "", "Missing 9%", "Missing 100% at baseline", "Missing 100% at baseline",
            "", "1-4 ordinal", "", "Redundant with CURSMOKE?",
            "Previous CHD - confounding?", "", "", "", "")
)

print("POTENTIAL PREDICTORS CLASSIFICATION:")
print(predictors_classification)

# Fokusér på primary kandidater (high clinical relevance, ikke previous disease)
primary_predictors <- c("AGE", "SEX", "TOTCHOL", "SYSBP", "BMI", 
                        "DIABETES", "CURSMOKE", "GLUCOSE")

print("PRIMARY PREDICTOR CANDIDATES:")
print(primary_predictors)

# ============================================================================
# . PAIRWISE ASSOCIATIONS OF PREDICTORS
# ============================================================================

# ----------------------------------------------------------------------------
#  Korrelationsmatrix - Kontinuerte variable
# ----------------------------------------------------------------------------

# Vælg kontinuerte prædiktorer (ekskl. HDLC/LDLC pga. missing)
cont_predictors <- c("AGE", "TOTCHOL", "SYSBP", "DIABP", "BMI", 
                     "HEARTRTE", "GLUCOSE", "CIGPDAY")


# ----------------------------------------------------------------------------
#  Visualisering - Korrelationsmatrix
# ---------------------------------------------------------------------------




#############################################################################################################################

#==================================================
# Alternativ udgangspunkt: Period 1 (aka. Baseline)
#==================================================

                  #### Struktur: ####

# 1. Marginal distributions
# 2. Missing values
# 3. Selection of potential predictors
# 4. Pairwise associations
# 5. Marginal associations with response

# Libraries
library(riskCommunicator)
library(dplyr)
library(ggplot2)
library(tidyr)
library(knitr)

data(framingham, package = "riskCommunicator")

#################### Data_structur ####################

# Hurtig gennemgang:
#______________________

dim(framingham) # nrow =  11,627. &  ncol = 39

# Purpose: Simplify analysis by using only first observation per person

framingham_baseline <- framingham %>% 
  filter(PERIOD == 1)     # Period 1 -> Can be examinated.
dim(framingham_baseline)   # 4434   39

framingham_Period_2 <- framingham %>%
  filter(PERIOD == 2)
dim(framingham_Period_2)  # 3930   39

framingham_Period_3 <- framingham %>%
  filter(PERIOD == 3) 
dim(framingham_Period_3) # 3263   39

# Check baseline dataset dimensions
dim(framingham_baseline)
# Result: Fewer rows now - only one per participant

#  SELECT KEY VARIABLES FOR INITIAL ANALYSIS ----
# Purpose: Focus on most clinically relevant predictors and outcomes
selected_vars <- framingham_baseline %>%
  select(
    # Demographic
    SEX, AGE,
    # Cholesterol
    TOTCHOL,
    # Blood Pressure
    SYSBP, DIABP, BPMEDS,
    # Lifestyle
    CURSMOKE, CIGPDAY, BMI,
    # Metabolic
    DIABETES, GLUCOSE,
    # Medical History
    PREVHYP,
    # Outcomes (will choose one as primary response later)
    CVD, ANYCHD, DEATH, HYPERTEN
  )

#  EXAMINE CONTINUOUS VARIABLES ----
# Purpose: Understand distribution, central tendency, spread, outliers

continuous_vars <- c("AGE", "TOTCHOL", "SYSBP", "DIABP", "BMI", "GLUCOSE")

#  Summary statistics for continuous variables
continuous_summary <- selected_vars %>%
  select(all_of(continuous_vars)) %>%
  summary()
print(continuous_summary)

#  Visualize distributions with histograms
par(mfrow = c(2, 3))  # Create 2x3 plot grid

for (var in continuous_vars) {
  hist(selected_vars[[var]], 
       main = paste("Distribution of", var),
       xlab = var,
       col = "lightblue",
       breaks = 20)
}

# Reset plot layout
par(mfrow = c(1, 1))

#  Check for outliers with boxplots
ggplot(gather(selected_vars %>% select(all_of(continuous_vars)), 
              key = "variable", value = "value"), 
       aes(x = variable, y = value)) +
  geom_boxplot(fill = "lightgreen") +
  theme_minimal() +
  labs(title = "Boxplots of Continuous Variables",
       x = "Variable",
       y = "Value") +
  coord_flip()  # Flip for better readability

# EXAMINE CATEGORICAL VARIABLES ----
# Purpose: Understand frequency distributions of categorical/binary variables

categorical_vars <- c("SEX", "CURSMOKE", "DIABETES", "BPMEDS", "PREVHYP")

# Frequency tables for categorical variables
cat_freq_tables <- list()

for (var in categorical_vars) {
  cat_freq_tables[[var]] <- table(selected_vars[[var]])
  print(paste("Frequency table for", var, ":"))
  print(cat_freq_tables[[var]])
}

#  Visualize categorical distributions with bar plots
par(mfrow = c(2, 3))

for (var in categorical_vars) {
  freq_table <- table(selected_vars[[var]])
  barplot(freq_table,
          main = paste("Distribution of", var),
          xlab = var,
          ylab = "Frequency",
          col = "lightcoral")
}

par(mfrow = c(1, 1))

# . EXAMINE POTENTIAL RESPONSE VARIABLES ----
# Purpose: Understand the outcomes we might want to predict

response_vars <- c("CVD", "ANYCHD", "DEATH", "HYPERTEN")

# 5.1 Frequency of outcomes in baseline data
response_freq <- list()

for (var in response_vars) {
  response_freq[[var]] <- table(selected_vars[[var]])
  print(paste("Prevalence of", var, "at baseline:"))
  print(round(prop.table(response_freq[[var]]) * 100, 1))

}

# . DATA QUALITY CHECK ----
# Purpose: Identify any obvious data quality issues

# Check for impossible values

range(selected_vars$AGE, na.rm = TRUE)
cat("SYSBP range:", range(selected_vars$SYSBP, na.rm = TRUE)
range(selected_vars$BMI, na.rm = TRUE)

# Check CIGPDAY when CURSMOKE = 0

selected_vars %>%
  filter(CURSMOKE == 0) %>%
  summarise(
    n = n(),
    cigpday_zero = sum(CIGPDAY == 0, na.rm = TRUE),
    cigpday_nonzero = sum(CIGPDAY > 0, na.rm = TRUE)
  )
# ==========================================================================



  

skim(data)
summary(data$CIGPDAY)
glimpse(data$CIGPDAY)
skim(data$CIGPDAY)
plot(data$CIGPDAY)
boxplot(data$CIGPDAY)
par(mfrow = c(1,2))
hist(data$CIGPDAY, breaks = 70)
hist(data$CIGPDAY, breaks = c(0,10, 20, 30, 40, 50, 60, 70))
plot( data$AGE, data$CIGPDAY)

table(training_data$CIGPDAY)
# Mange 0-værdier, få positive værdier med begrænset variation
data <- training_data
# Opret aldersgrupper
data_for_plot_CIGPDAY_AGE <- data %>%
  mutate(age_group = cut(AGE, 
                         breaks = c(30, 40, 50, 60, 70),
                         labels = c("30-39", "40-49", "50-59", "60-70")))

ggplot(data_for_plot_CIGPDAY_AGE, aes(x = age_group, y = CIGPDAY)) +
  geom_boxplot(fill = "lightblue", alpha = 0.7) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "red") +
  labs(
    title = "Cigarettes per Day by Age Group",
    subtitle = "Red diamond = mean, Box = IQR",
    x = "Age Group",
    y = "Cigarettes per Day"
  ) +
  theme_minimal()

ggplot(data_for_plot_CIGPDAY_AGE, aes(x = age_group, y = CIGPDAY)) +
  geom_violin(fill = "lightgreen", alpha = 0.6) +
  geom_boxplot(width = 0.2, fill = "white", alpha = 0.8) +
  labs(
    title = "Distribution of Cigarette Consumption by Age",
    x = "Age Group", 
    y = "Cigarettes per Day"
  ) +
  theme_minimal()

# Opret ryger kategorier
data_smoke_Categ <- data_for_plot_CIGPDAY_AGE %>%
  mutate(smoking_status = case_when(
    CIGPDAY == 0 ~ "Non-smoker (0)",
    CIGPDAY <= 10 ~ "Lightsmokers (1-10)",
    CIGPDAY <= 20 ~ "Moderate smoker (10-20)", 
    TRUE ~ "Heavysmokers > (20)"
  ))
summary(data_smoke_Categ)

# Mere detaljeret summary inklusive CVD
detailed_smoking_summary <- data_smoke_Categ %>%
  mutate(smoking_status = factor(smoking_status,
                                 levels = c("Non-smoker (0)",
                                            "Lightsmokers (1-10)",
                                            "Moderate smoker (10-20)",
                                            "Heavysmokers > (20)"))) %>% 
  group_by(smoking_status) %>%
  summarise(
    # Basis statistik
    n = n(),
    percentage = round(n / nrow(data_smoke_Categ) * 100, 1),
    
    # Alder statistik
    mean_age = round(mean(AGE), 1),
    median_age = median(AGE),
    age_range = paste(min(AGE), "-", max(AGE)),
    
    # CVD statistik
    n_cvd = sum(CVD == 1),
    pct_cvd = round(n_cvd / n * 100, 1),
    
    # Køn fordeling
    n_male = sum(SEX == 1),  # Antag SEX=1 er male
    pct_male = round(n_male / n * 100, 1),
    
    .groups = 'drop'
  ) %>% 
  arrange(smoking_status)

print(detailed_smoking_summary)

# Stacked bar plot
ggplot(data_smoke_Categ, aes(x = age_group, fill = smoking_status)) +
  geom_bar(position = "fill") +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Smoking Status Distribution by Age Group",
    subtitle = "Proportion of each smoking category",
    x = "Age Group",
    y = "Proportion",
    fill = "Smoking Status"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")


ggplot(data_smoke_Categ, aes(x = AGE, y = CIGPDAY)) +
  geom_jitter(alpha = 0.3, width = 0.5, height = 0.5, color = "blue") +
  geom_smooth(method = "loess", color = "red", se = TRUE) +
  labs(
    title = "Cigarette Consumption vs Age with Trend Line",
    x = "Age",
    y = "Cigarettes per Day"
  ) +
  theme_minimal()

# Facet plot med CVD status
ggplot(data_smoke_Categ , aes(x = age_group, y = CIGPDAY, fill = factor(CVD))) +
  geom_boxplot(alpha = 0.7) +
  facet_wrap(~ factor(CVD, labels = c("No CVD", "CVD")), ncol = 2) +
  scale_fill_manual(values = c("lightblue", "salmon")) +
  labs(
    title = "Cigarette Consumption by Age and CVD Status",
    x = "Age Group",
    y = "Cigarettes per Day",
    fill = "CVD Status"
  ) +
  theme_minimal() +
  theme(legend.position = "none")


# Fokus kun på rygere
smokers_data <- data_smoke_Categ %>% filter(CIGPDAY > 0)

ggplot(smokers_data, aes(x = AGE, y = CIGPDAY)) +
  geom_point(alpha = 0.5, color = "darkred") +
  geom_density_2d_filled(alpha = 0.7) +
  geom_smooth(method = "lm", color = "black", linetype = "dashed") +
  labs(
    title = "Density of Cigarette Consumption by Age (Smokers Only)",
    subtitle = "Contour lines show density of observations",
    x = "Age",
    y = "Cigarettes per Day"
  ) +
  theme_minimal()

# install.packages("plotly")
library(plotly)

p <- ggplot(data, aes(x = AGE, y = CIGPDAY, color = factor(CVD), 
                      text = paste("Age:", AGE, "<br>Cigarettes:", CIGPDAY))) +
  geom_jitter(alpha = 0.6) +
  scale_color_manual(values = c("blue", "red"), labels = c("No CVD", "CVD")) +
  labs(
    title = "Cigarette Consumption by Age and CVD Status",
    x = "Age",
    y = "Cigarettes per Day",
    color = "CVD Status"
  ) +
  theme_minimal()

ggplotly(p, tooltip = "text")




# Opret summary statistiker
smoking_summary <- data_smoke_Categ %>%
  group_by(age_group) %>%
  summarise(
    n = n(),
    n_smokers = sum(CIGPDAY > 0),
    pct_smokers = round(n_smokers / n * 100, 1),
    mean_cigarettes = round(mean(CIGPDAY[CIGPDAY > 0]), 1),
    median_cigarettes = median(CIGPDAY[CIGPDAY > 0]),
    .groups = 'drop'
  )

print(smoking_summary)


# 



