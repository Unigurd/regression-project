library(tibble)
library(broom)
library(readr)
library(skimr)
library(dplyr)
library(tidyr)
library(hexbin)

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
