
library(tidyverse)
library(pROC)
library(car)
library(broom)


# Make sure thhe data corresponds to cleaned_version -> model_input   !!!!!!!!!

set.seed(8088)


train <- read_csv("model_input.csv") %>%
  select(-1) %>%
  mutate(
    CVD = factor(CVD, levels = c(0,1), labels = c("No","Yes")),
    SEX = factor(SEX, levels = c(1,2), labels = c("Male","Female")),
    DIABETES = factor(DIABETES, levels = c(0,1), labels = c("No","Yes")),
    PREVHYP = factor(PREVHYP, levels = c(0,1), labels = c("No","Yes")),
    educ = factor(educ, ordered = TRUE)
  )


nrow(train)
round(100*mean(train$CVD=="Yes"),1) # CVD-RATIO = yes in our data


# Models: 

# Model 1: Minimal
m1 <- glm(CVD ~ AGE + SEX, data = train, family = binomial)

# Model 2: Core
m2 <- glm(CVD ~ AGE + SEX + SYSBP + TOTCHOL + BMI + DIABETES, 
          data = train, family = binomial)

# Model 3: 
m3 <- glm(CVD ~ AGE + SEX + SYSBP + TOTCHOL + BMI + DIABETES + 
            CIGPDAY + educ + PREVHYP,
          data = train, family = binomial)

summary(m1)
summary(m2)
summary(m3)

ompare <- tibble(
  Model = c("Minimal","Core","m3"),
  n_params = c(length(coef(m1)), length(coef(m2)), length(coef(m3))),
  AIC = c(AIC(m1), AIC(m2), AIC(m3)),
  AUC = c(
    auc(roc(train$CVD, fitted(m1), quiet=T)),
    auc(roc(train$CVD, fitted(m2), quiet=T)),
    auc(roc(train$CVD, fitted(m3), quiet=T))
  )
)

compare

anova(m1,m2,m3, test = "LRT")

vif_vals <- vif(m3)
max(round(vif_vals, 2))
max(vif_vals)

# sætter m3 som primære model nu
primary <- m3


roc_obj <- roc(train$CVD, fitted(primary), quiet=TRUE)
roc_obj


# visueltmæssigt, lav kollonne for signifikant og se p-værdier i stjerne form
or_table <- tidy(primary, conf.int=TRUE, exponentiate=TRUE) %>%
  filter(term != "(Intercept)") %>%
  mutate(
    sig = ifelse(p.value < 0.001, "***",
                 ifelse(p.value < 0.01, "**",
                        ifelse(p.value < 0.05, "*", "")))
  ) %>%
  select(term, OR=estimate, Lower=conf.low, Upper=conf.high, p=p.value, sig)

print(or_table, digits=3)

# ===========================================

# lets plot some (se pdf for nogle flotte plot): 


par(mfrow=c(2,2), mar=c(4,4,2,1))

plot(roc_obj, main="ROC Curve", 
    col="steelblue", lwd=2, print.auc=TRUE, print.auc.y=0.4)
abline(0,1, lty=2, col="gray")


pred <- fitted(primary)
obs <- as.numeric(train$CVD=="Yes")
cal_data <- data.frame(pred, obs) %>%
  mutate(bin = cut(pred, breaks=10)) %>%
  group_by(bin) %>%
  summarise(pred_mean = mean(pred), obs_mean = mean(obs), .groups="drop")

plot(cal_data$pred_mean, cal_data$obs_mean, 
     xlim=c(0,1), ylim=c(0,1), pch=19, col="steelblue", cex=1.5,
     xlab="Predicted Probability", ylab="Observed Proportion",
     main="Calibration Plot")
abline(0,1, lty=2, col="red", lwd=2)


# 3. Forest Plot (Top 5)
or_data <- or_table %>% slice_head(n=7)
plot(or_data$OR, 1:nrow(or_data), pch=19, col="steelblue", cex=1.5,
     xlim=c(0.4, 3), yaxt="n", xlab="Odds Ratio (log scale)", 
     ylab="", main="Odds Ratios", log="x")
segments(or_data$Lower, 1:nrow(or_data), 
         or_data$Upper, 1:nrow(or_data), lwd=2, col="steelblue")
abline(v=1, lty=2, col="red", lwd=2)
axis(2, at=1:nrow(or_data), labels=or_data$term, las=1, cex.axis=0.8)


ages <- seq(30, 75, by=5)
risk_male <- predict(primary, 
                     newdata=data.frame(AGE=ages, SEX="Male", SYSBP=130, TOTCHOL=220, 
                                        BMI=26, DIABETES="No", CIGPDAY=0, educ=factor(2,ordered=T), 
                                        PREVHYP="No"),
                     type="response")

risk_female <- predict(primary,
                       newdata=data.frame(AGE=ages, SEX="Female", SYSBP=130, TOTCHOL=220,
                                          BMI=26, DIABETES="No", CIGPDAY=0, educ=factor(2,ordered=T),
                                          PREVHYP="No"),
                       type="response")

plot(ages, risk_male, type="l", col="steelblue", lwd=3,
     ylim=c(0,0.8), xlab="Age", ylab="Predicted CVD Risk",
     main="Risk by Age & Sex")
lines(ages, risk_female, col="salmon", lwd=3)
legend("topleft", c("Male","Female"), col=c("steelblue","salmon"), 
       lwd=3, bty="n")

par(mfrow=c(1,1))

# ============================================================================


# fast summary:
m1
m2
m3


# round(AIC(m1),1)
# round(AIC(m2),1)
# round(AIC(m3),1)

compare


round(auc(roc_obj),3)

round(max(vif_vals),2)




top1 <- or_table %>% arrange(p) %>% slice_head(n=3) # brug cat funktionen for bvisualisering af resultater - better output
for(i in 1) {
  cat("  ", i, ". ", top1$term[i], " (OR=", round(top1$OR[i],2), 
      ", p", format.pval(top1$p[i],digits=2), ")\n", sep="")
}
top3 <- or_table %>% arrange(p) %>% slice_head(n=3)
for(i in 1:3) {
  cat("  ", i, ". ", top3$term[i], " (OR=", round(top3$OR[i],2), 
      ", p", format.pval(top3$p[i],digits=2), ")\n", sep="")
} 
top5 <- or_table %>% arrange(p) %>% slice_head(n=5)
  for(i in 1:5) {
    cat("  ", i, ". ", top5$term[i], " (OR=", round(top5$OR[i],2), 
        ", p", format.pval(top5$p[i],digits=2), ")\n", sep="")
}
# 1. SEXFemale (OR=0.39, p<2e-16)
# 2. AGE (OR=1.05, p<2e-16)
# 3. DIABETESYes (OR=3.13, p4.5e-07)
# 4. SYSBP (OR=1.01, p3.6e-06)
# 5. TOTCHOL (OR=1, p5.3e-05)





w_2i <- exp(eta_2i)

# Round to 2 decimals
round(w_2i, 2)




" Roc_curve:
* X-akse: Specificity (1 - False Positive Rate)
* Y-akse: Sensitivity (True Positive Rate)
* Kurven: Viser trade-off mellem at fange CVD cases (sensitivity) og korrekt identificere non-CVD (specificity)
Fortolkning:
* AUC = 0.745 → "God" discrimination
* Kurven ligger langt over diagonalen (som er random guessing)
* Jo mere kurven "buer" mod øvre venstre hjørne, jo bedre model 
Can report this if  we will use it -> 
The ROC curve demonstrates good discriminative ability with an AUC of 0.745, indicating the model correctly ranks 74.5% of randomly selected CVD/non-CVD pairs. 
This performance exceeds the acceptable threshold (AUC > 0.7) for clinical risk prediction models."


"Calibration
* X-akse: Predicted Probability (model's forudsigelse)
* Y-akse: Observed Proportion (faktisk CVD rate)
* Rød stiplet linje: Perfekt calibration (45° linje)
* Blå points: Gennemsnit per decile
Fortolkning:
* Points ligger tæt på den røde linje → God calibration
* Betyder: Når model siger "20% risk", har ~20% faktisk CVD
* Lidt afvigelse ved høje probabilities, men overall acceptabelt

Can report this if  we will use it -> 
"The calibration plot shows good agreement between predicted probabilities and observed event rates. 
Points cluster closely around the 45° reference line, indicating the model provides well-calibrated risk estimates across the probability spectrum. 
Minor deviations at higher predicted probabilities suggest slight overestimation for highest-risk individuals."


"Odds-ratio
* X-akse: Odds Ratio (log scale)
* Y-akse: Prædiktorer
* Points: OR estimat
* Horisontale linjer: 95% Confidence Intervals
* Rød vertikal linje: OR = 1 (ingen effekt)
Fortolkning:
Højre side af rød linje (OR > 1) = Øget risk:
* DIABETES (længst til højre): Stærkeste effekt, OR ≈ 2.4
* CIGPDAY, BMI, TOTCHOL, SYSBP: Moderate positive effekter
Venstre side af rød linje (OR < 1) = Reduceret risk:
* SEX-Female: Beskyttende effekt (OR ≈ 0.55)
Tæt på rød linje:
* AGE: OR ≈ 1.08 (per år) - ser lille ud, men kraftig over tid" 

Can report this if  we will use it -> 
"Figure 3 displays odds ratios with 95% confidence intervals. Diabetes shows the strongest association (OR≈2.4), more than doubling CVD odds. 
Female sex demonstrates a protective effect (OR=0.55), with 45% lower odds compared to males. Age, though appearing close to unity on the log scale (OR=1.08 per year), compounds to substantial risk over decades. 
Traditional cardiovascular risk factors (SYSBP, TOTCHOL, BMI) show expected positive associations, while confidence intervals exclude unity for all significant predictors."

Risk by age
" * X-akse: Age (30-75 years)
* Y-akse: Predicted CVD Risk
* Blå linje: Males
* Rød linje: Females
Fortolkning:
* Clear age gradient: Risk stiger med alder for begge køn
* Sex difference persistent: Mænd har konsistent ~2× højere risk på tværs af alle aldre
* Acceleration: Risk stiger hurtigere efter ~55 år
* At age 70: Males ≈60% risk, Females ≈35% risk"

"Can report this if  we will use it -> 
Predicted CVD risk increases substantially with age for both sexes, with males consistently demonstrating approximately twice the risk of females across all age strata. 
At age 40, predicted risk is relatively low (males 15%, females 8%), but by age 70, males reach 60% predicted risk while females reach 35%. 
The parallel trajectories indicate the sex effect is multiplicative rather than additive, persisting throughout the age range. 
These predictions assume median values for other risk factors (no diabetes, no smoking, normal blood pressure)." 






"SAMLET! 
Predicted CVD risk increases substantially with age for both sexes, with males consistently demonstrating approximately twice the risk of females across all age strata. 
At age 40, predicted risk is relatively low (males 15%, females 8%), but by age 70, males reach 60% predicted risk while females reach 35%. 
The parallel trajectories indicate the sex effect is multiplicative rather than additive, persisting throughout the age range. 
These predictions assume median values for other risk factors (no diabetes, no smoking, normal blood pressure)."

