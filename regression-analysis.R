library(splines)
library(Hmisc)
library(tibble)
library(broom)
library(readr)
library(skimr)
library(dplyr)
library(tidyr)
library(hexbin)
library(ggplot2)

library(riskCommunicator)
data(framingham, package="riskCommunicator")


input <- as_tibble(read.csv("training_data.csv"))

set.seed(8088)
n               <- nrow(framingham)
validation_size <- floor(n/10)
i               <- sample(n, validation_size)
validation_data <- framingham[i,]
training_data   <- framingham[-i,]

framingham_ordered_cols <- function(data) {
    data |>
        mutate(
            SEX      = ordered(SEX, levels=c(1,2), labels=c("M", "F")),
            CURSMOKE = as.ordered(CURSMOKE),
            DIABETES = as.ordered(DIABETES),
            BPMEDS   = as.ordered(BPMEDS),
            ## We don't actually know what education the educ levels stand for,
            ## but we assume they're ordered.
            educ     = as.ordered(educ),
            PREVCHD  = as.ordered(PREVCHD),
            PREVAP   = as.ordered(PREVAP),
            PREVMI   = as.ordered(PREVMI),
            PREVSTRK = as.ordered(PREVSTRK),
            PREVHYP  = as.ordered(PREVHYP),
            PERIOD   = as.ordered(PERIOD),
            DEATH    = as.ordered(DEATH),
            ANGINA   = as.ordered(ANGINA),
            HOSPMI   = as.ordered(HOSPMI),
            MI_FCHD  = as.ordered(MI_FCHD),
            ANYCHD   = as.ordered(ANYCHD),
            STROKE   = as.ordered(STROKE),
            CVD      = as.ordered(CVD),
            HYPERTEN = as.ordered(HYPERTEN),
            )
}

real_validation_data <- validation_data |> filter(PERIOD == 1) |> select(RANDID) |> distinct() |> left_join(framingham, by=join_by(RANDID)) |> framingham_ordered_cols()
framingham1          <- training_data   |> filter(PERIOD == 1) |> select(RANDID) |> distinct() |> left_join(framingham, by=join_by(RANDID)) |> framingham_ordered_cols()

## In order for us to use assumption A4 about the values of CVD being uncorrelated conditioned on the other variables, we need to get rid of all but one observation for each person. Otherwise there will be an obvious correlation between the values of CVD for the observations of the same person: They will be the same. We do this now so that when we start plotting the distributions of our variables we can better justify thinking of each observation as independent, making the plots easier to interpret. We cannot use the observations where the first checkup we have for a person is from period 2 or 3, since it will not be independent of any observations from earlier periods for that person that may be in the validation data. We therefore just us the observations from the baseline exam.
framingham2 <- framingham1 |> filter(PERIOD == 1)

## We also don't want to include people who dropped out of the study too early, as they had less time to get CVD, thus skewing our predictions for probability of getting CVD downwards. To do this we notice in the documentation of the data set that all the time-to-event columns like TIMECVD, TIMEDTH and so on register the real censoring time for each subject when the given event didn't occur. We have checked that these are all in agreement except for TIMEDTH which sometimes has the end-of-followup time even when the other time-to-event columns indicate the subject was lost to followup before that. The following piece of code calculates the observation duration under the assumption that TIME is always 0 in period 1 observations, which we have checked.
framingham3 <- framingham2 |>
    mutate(
        timeap   = ifelse(ANGINA   == 1, -1, TIMEAP),
        timemi   = ifelse(HOSPMI   == 1, -1, TIMEMI),
        timemifc = ifelse(MI_FCHD  == 1, -1, TIMEMIFC),
        timechd  = ifelse(ANYCHD   == 1, -1, TIMECHD),
        timestrk = ifelse(STROKE   == 1, -1, TIMESTRK),
        timecvd  = ifelse(CVD      == 1, -1, TIMECVD),
        timehyp  = ifelse(HYPERTEN == 1, -1, TIMEHYP),
        max      = pmax(timeap, timemi, timemifc, timechd, timestrk, timecvd, timehyp),
        ## if max == -1 then all of the time-to-event cols had the event happen.
        ## In that case we use TIMEDTH instead, since it contains the time of last contact
        ## with some person (or maybe a larger value as we saw above, but it's still our
        ## best guess) if the person didnt die, and also if the person did die.
        time_censor = ifelse(max == -1, TIMEDTH, max),
        .after=RANDID
    ) |>
    subset(select=-c(max, timeap, timemi, timemifc, timechd, timestrk, timecvd, timehyp))

## Now we can remove any observations from people who dropped out of the study before the halfway mark. The halfway mark was chosen pretty arbitrarily.
framingham4 <- framingham3 |> filter(time_censor > 12*365.25)

## A lot of the columns can already be ruled out as potential predictors. Either they're meaningless like RANDINT, or indicate a future event from the point of the baseline exam, or the time until a future event. This leaves us with the following potential predictors:

## We're now ready to examine the marginal distributions of the continuous variables:
dens1 <- ggplot(data=framingham4) + geom_density(fill=gray(0.7))
gridExtra::grid.arrange(
               dens1 + aes(TOTCHOL),  dens1 + aes(AGE),
               dens1 + aes(SYSBP),    dens1 + aes(DIABP),
               dens1 + aes(CIGPDAY),  dens1 + aes(BMI),
               dens1 + aes(HEARTRTE), dens1 + aes(GLUCOSE),
               ncol=4
           )

## We notice that GLUCOSE, BMI and mayyybe HEARTRTE, SYSBP and DIABP have skewed distributions, making them candidates for transformations, maybe by a logarithm. We do not have the domain knowledge to say whether this makes sense, but on an intuitive level it makes sense that if something can vary like glucose does, then it's not whether the glucose measurements differ by e.g. 10 that matters, but whether the measurements double or triple (This shouldn't be understood in a causative sense, naturally). So in this sense a log transform seems reasonable. We note that CIGPDAY is very skewed, but don't commit to any transformation of it yet.

## The log-transforms of the variables look much more reasonable, so we stick with that for now. The only one that doesn't improve much is GLUCOSE, and it keeps having a large tail regardless of how much we log-transform it. Maybe it's an extreme outlier causing the tail then.
gridExtra::grid.arrange(
               dens1 + aes(log(SYSBP)),   dens1 + aes(log(DIABP)),
               dens1 + aes(log(BMI)),     dens1 + aes(log(HEARTRTE)),
               dens1 + aes(log(GLUCOSE)), dens1 + aes(log(log(GLUCOSE))),
               dens1 + aes(log(log(log(GLUCOSE)))), ncol=4
           )

## All of our categorical variables except educ are binary. educ tapers off to the right so higher numbers likely corresponds to higher levels of education, which was our assumption as well. There are some missing values we wil need to handle, but the most striking thing is how skewed DIABETES, BPMEDS, PREVCHD, PREVAP, PREVMI and PREVSTRK are. The latter two in particular seem to have so few occurrences that they're candidates for being thrown out. If we don't have enough observations where MI or CHD, we cannot say anything about their effect with certainty anyway.
bar1 <- ggplot(data=framingham4) + geom_bar(fill=gray(0.7))
gridExtra::grid.arrange(
               bar1 + aes(SEX),      bar1 + aes(CURSMOKE),
               bar1 + aes(DIABETES), bar1 + aes(BPMEDS),
               bar1 + aes(educ),     bar1 + aes(PREVCHD),
               bar1 + aes(PREVAP),   bar1 + aes(PREVMI),
               bar1 + aes(PREVSTRK), bar1 + aes(PREVHYP),
               bar1 + aes(CVD),      ncol=4
           )


framingham4_log <- framingham4 |>
    mutate(log(GLUCOSE), log(DIABP), log(SYSBP), log(BMI)) |>
    select(-c(GLUCOSE, DIABP, SYSBP, BMI))


# Examining the correlation between our continuous variable, the strong colinearity of SYSBP and DIABP immediately catches the eye. One of them has to go, and we choose DIABP because it has higher correlation with BMI, the variables that both are second-most correlated to. No other correlations are high enough to remove variables at this point, but this plot does give us a convenient view of outliers. In particular, the tail of GLUCOSE is not caused by a few outliers, but rather it seems like a real trend in the data. Other than that, there are points that lie away from the other, but none that really necessitate removal.
ggally_hexbin <- function(data, mapping, ...) {
    ggplot(data, mapping) + stat_binhex(...)
}

GGally::ggpairs(
            framingham4_log,
            columns = c("AGE", "TOTCHOL", "log(SYSBP)", "log(DIABP)",
                        "CIGPDAY", "log(BMI)", "HEARTRTE", "log(GLUCOSE)"),
            lower   = list(continuous = GGally::wrap("hexbin", bins=20)),
            diag = list(continuous = "blankDiag")
        )


## Since all our categorical variables conveniently are ordered, we use the spearman correlation to get an idea of the correlation of all of our potentials predictors. There is naturally a strong correlation between CURSMOKE and CIGPDAY. If we believe that smoking any cigarettes at all is qualitatively different than not smoking we might choose to include both, but we don't have much reason to believe so we get rid of CURSMOKE.
cp <- cor(
    framingham4 |> select(AGE, TOTCHOL, log(SYSBP), log(DIABP),
                          CIGPDAY, log(BMI), HEARTRTE, log(GLUCOSE))
    |> data.matrix(),
    use    = "complete.obs",
    method = "spearman"
)
corrplot::corrplot(
              cp,
              diag    = FALSE,
              order   = "hclust",
              addrect = 14,
              tl.srt  = 45,
              tl.col  = "black",
              tl.cex  = 0.8
          )

## PREVCHD is very correlated to both PREVMI and PREVAP. Reading their descriptions we see that PREVAP is PREValent Angina Pectoris, PREVMI is PREValent Myocardial Infarction, and that PREVCHD is either of those or prevalent coronary insufficiency. Checking the data, we indeed see that PREVCHD is always 1 whenever PREVAP or PREVMI is, and that we only have two cases of PREVCHD that is not also PREVAP or PREVMI, that is, cases that we know must be prevalent coronary insufficiency. That is not information enough to sensibly impute prevalent coronary insufficiency for all the observations where either PREVAP or PREVMI is 1, so we just get rid of PREVCHD wholesale.
framingham4 |> filter(PREVCHD == 0, PREVMI == 1 | PREVAP == 1) |> nrow()
framingham4 |> filter(PREVCHD == 1, PREVMI == 0,  PREVAP == 0) |> nrow()

## The last strongly correlated group is SYSBP, DIABP and PREVHYP. We have already gotten rid of DIABP, and looking at a plot of smoothed CVD probability stratified by PREVHYP, we don't really see any reason to keep PREVHYP. The behavior seems exactly the same, just more fluctuating for PREVHYP=1 due to the fewer observations.
ggplot(framingham4, aes(log(SYSBP), as.numeric(CVD), color=PREVHYP, fill=PREVHYP)) + geom_smooth(alpha=0.5, method="loess")

## CVD didn't really seem strongly correlated with anything so we plot the density of the continuous variables stratified by CVD.
density_by_cvd_legend <- ggplot(data=framingham4, aes(fill=as.factor(CVD)))  + geom_density(alpha=0.4)
density_by_cvd        <- density_by_cvd_legend + theme(legend.position="none")
gridExtra::grid.arrange(
               density_by_cvd + aes(TOTCHOL),    density_by_cvd + aes(AGE),
               density_by_cvd + aes(log(SYSBP)), density_by_cvd_legend + aes(CIGPDAY),
               density_by_cvd + aes(log(BMI)),   density_by_cvd + aes(HEARTRTE),
               density_by_cvd + aes(log(GLUCOSE)), ncol=4
           )
## We see that some of the variables have their distribution translated slightly higher when CVD=1, but the biggest change is that risk of CVD clearly increases with age.

## Looking at barplots of the relative frequencies of CVD for the different values of the categorical variables, we see that a lot of the rare conditions indicate a great increase in the risk of CVD. PREVMI and PRESTRK in particular, which were the ones we considered removing because there were so few instances of them happening. Seeing this, we'll keep them in.
stacked_bars_legend <- ggplot(data=framingham4, aes(SEX, fill=CVD)) + geom_bar(position="fill", alpha=0.7) + ylab("")
stacked_bars        <- stacked_bars_legend + theme(legend.position="none")
gridExtra::grid.arrange(
               stacked_bars + aes(SEX),         stacked_bars + aes(CURSMOKE),
               stacked_bars + aes(DIABETES),    stacked_bars + aes(BPMEDS),
               stacked_bars_legend + aes(educ), stacked_bars + aes(PREVCHD),
               stacked_bars + aes(PREVAP),      stacked_bars + aes(PREVMI),
               stacked_bars + aes(PREVSTRK),    stacked_bars + aes(PREVHYP),
               ncol=5
           )


## Examining the missing values in our data, we see that GLUCOSE represents most of them. It makes sense that CIGPDAY and BMI are more likely to be missing when they're high, which would make them MNAR and thus not really possible to impute. We don't really have a reason to suspect any of the other columns are MNAR, so we will impute them
framingham4 |>
    skim() |>
    select(skim_variable, n_missing) |>
    filter(n_missing != 0, skim_variable != "HDLC", skim_variable != "LDLC") |>
    arrange(desc(n_missing))

## Luckily the missing CIGPDAY and BMI consist of less than 0.001 of the data, so if we impute the rest we
sum(is.na(framingham4$CIGPDAY) | is.na(framingham4$BMI))/nrow(framingham4)

table(rowSums(is.na(framingham4)))

cona4 <- framingham4 |>
    select(TOTCHOL, CIGPDAY, BPMEDS, GLUCOSE, educ)|>
    filter(rowSums(is.na(framingham4)) >= 2) |>
    as.matrix() |>
    is.na() |>
    crossprod()
round(sweep(cona4, 1, diag(cona4), "/"), 2)


period1 <- framingham1 |> filter(PERIOD == 1)
period2 <- framingham1 |> filter(PERIOD == 2) |>
    mutate(RANDID, CURSMOKE2=CURSMOKE, CIGPDAY2=CIGPDAY, BMI2=BMI,
           TOTCHOL2=TOTCHOL, BPMEDS2=BPMEDS, educ2=educ, GLUCOSE2=GLUCOSE,
           .keep="none"
           )
period3 <- framingham1 |> filter(PERIOD == 3) |>
    mutate(RANDID, CURSMOKE3=CURSMOKE, CIGPDAY3=CIGPDAY, BMI3=BMI,
           TOTCHOL3=TOTCHOL, BPMEDS3=BPMEDS, educ3=educ, GLUCOSE3=GLUCOSE,
           LDLC3=LDLC, HDLC3=HDLC,
           .keep="none"
)

## Add prevhyp and hyperten to impute BPMEDS.
## Don't bother with the other educs. They're all na if the first is na.

framingham5 <- period1 |>
    left_join(period2, by=join_by(RANDID)) |>
    left_join(period3, by=join_by(RANDID))

framingham7 <- framingham5
## Why are there more missing BMIs here than in framingham4? Oh
## because we've filtered rows with a short time off. We should
## probably because this on that, rather than framingham1.
## Not much to impute from though.
framingham7 |> filter(is.na(BMI)) |> mutate(B2=!is.na(BMI2), B3=!is.na(BMI3), .keep="none") |> table()
## Looks a bit better here. There might still be too many missing in
## GLUCOSE2 and GLUCOSE3 for it to be worth the time.
framingham7 |> filter(is.na(GLUCOSE)) |> mutate(B2=!is.na(GLUCOSE2), B3=!is.na(GLUCOSE3), .keep="none") |> table()
## Nothing to gain from educ.
framingham7 |> filter(is.na(educ)) |> mutate(B2=!is.na(educ), B3=!is.na(educ3), .keep="none") |> table()
## All nas in BPMEDS2
framingham7 |> filter(is.na(BPMEDS)) |> mutate(B2=!is.na(BPMEDS), B3=!is.na(BPMEDS3), .keep="none") |> table()
## All nas in TOTCHOL2
framingham7 |> filter(is.na(TOTCHOL)) |> mutate(B2=!is.na(TOTCHOL), B3=!is.na(TOTCHOL3), .keep="none") |> table()
## All nas in CIGPDAY2
framingham7 |> filter(is.na(CIGPDAY)) |> mutate(B2=!is.na(CIGPDAY), B3=!is.na(CIGPDAY3), .keep="none") |> table()

## We have removed PREVSTRK because it gave issues because bootstrap resamples would have too few unique values of it.
impute_form <- reformM(
    ~ CVD + SEX + TOTCHOL + AGE + SYSBP +
        CIGPDAY + DIABETES + HEARTRTE + PREVAP + PREVMI +
        + PREVHYP + BMI + BPMEDS + GLUCOSE + educ +
        CURSMOKE2 + CIGPDAY2 + BMI2 + TOTCHOL2 + BPMEDS2 + educ2 + GLUCOSE2 +
        CURSMOKE3 + CIGPDAY3 + BMI3 + TOTCHOL3 + BPMEDS3 + educ3 + GLUCOSE3 +
        LDLC3 + HDLC3,
    framingham5
)

imp <- aregImpute(
    impute_form,
    framingham7,
    # Set n.impute to proper value later
    n.impute=5,
    nk=4
)

plotting_data <- as_tibble(impute.transcan(
    imp,
    data=framingham7,
    imputation=1,
    list.out=TRUE,
    pr=FALSE
))

## THERE'S A PLOT METHOD FOR DISTRIBUTION OF IMPUTED VALUES
## MENTIONED IN  areg.impute docs


impute.all <- function(form, data, var, k, n.impute=5, nk=4) {
    n <- nrow(data)
    imputed <- list()
    for (i in seq(n.impute)) {
        imputed[i] <- list(rep(NA, n))
    }
    groups <- sample(rep(1:k, length.out=n))
    for (i in seq(k)) {
        data_na <- mutate(data, !!var:=ifelse(groups == i, NA, data[[var]]))
        imp <- aregImpute(
            form,
            data_na,
            n.impute=n.impute,
            nk=nk
        )
        for (j in seq(n.impute)) {
            imputed[[j]][groups==i] <- impute.transcan(
                imp,
                data=data_na,
                imputation=j,
                list.out=TRUE,
                pr=FALSE
            )[[var]][groups==i]
        }
    }
    imputed
}


imputed_glucose <- impute.all(impute_form, framingham7, "GLUCOSE", k=10)
write.csv(imputed_glucose, "imputed_glucose.csv")


imputed_educ <- impute.all(impute_form, framingham7, "educ", k=10, nk=0)

imputed_bpmeds <- impute.all(impute_form, framingham7, "BPMEDS", k=10)

## Plot against correlated variables

ggplot(mapping=aes(
           framingham7$GLUCOSE[!is.na(framingham7$GLUCOSE)],
           imputed_glucose[[2]][!is.na(framingham7$GLUCOSE)])
       ) +
    geom_abline(intercept=0, slope=1, color="blue") +
    geom_point()


## Maybe tables for the first two?
## Or overlapping bars or something.
gridExtra::grid.arrange(
               ggplot(mapping=aes(framingham7$educ[!is.na(framingham7$educ)])) +
               geom_bar(fill=gray(0.7)),
               ggplot(mapping=aes(imputed_educ[[1]][!is.na(framingham7$educ)])) +
               geom_bar(fill=gray(0.7)),
               ggplot(mapping=aes(imputed_educ[[2]][!is.na(framingham7$educ)])) +
               geom_bar(fill=gray(0.7)),
               ggplot(mapping=aes(imputed_educ[[3]][!is.na(framingham7$educ)])) +
               geom_bar(fill=gray(0.7)),
               ggplot(mapping=aes(imputed_educ[[4]][!is.na(framingham7$educ)])) +
               geom_bar(fill=gray(0.7)),
               ggplot(mapping=aes(imputed_educ[[5]][!is.na(framingham7$educ)])) +
               geom_bar(fill=gray(0.7)),
               ncol=3)



gridExtra::grid.arrange(
               ggplot(mapping=aes(framingham7$BPMEDS[!is.na(framingham7$BPMEDS)])) +
               geom_bar(fill=gray(0.7)),
               ggplot(mapping=aes(imputed_bpmeds[[1]][!is.na(framingham7$BPMEDS)])) +
               geom_bar(fill=gray(0.7)),
               ggplot(mapping=aes(imputed_bpmeds[[2]][!is.na(framingham7$BPMEDS)])) +
               geom_bar(fill=gray(0.7)),
               ggplot(mapping=aes(imputed_bpmeds[[3]][!is.na(framingham7$BPMEDS)])) +
               geom_bar(fill=gray(0.7)),
               ggplot(mapping=aes(imputed_bpmeds[[4]][!is.na(framingham7$BPMEDS)])) +
               geom_bar(fill=gray(0.7)),
               ggplot(mapping=aes(imputed_bpmeds[[5]][!is.na(framingham7$BPMEDS)])) +
               geom_bar(fill=gray(0.7)),
               ncol=3)


gridExtra::grid.arrange(
               ggplot(mapping=aes(framingham7$GLUCOSE[!is.na(framingham7$GLUCOSE)])) +
               geom_density(fill=gray(0.7)),
               ggplot(mapping=aes(imputed_glucose[[1]][!is.na(framingham7$GLUCOSE)])) +
               geom_density(fill=gray(0.7)),
               ggplot(mapping=aes(imputed_glucose[[2]][!is.na(framingham7$GLUCOSE)])) +
               geom_density(fill=gray(0.7)),
               ggplot(mapping=aes(imputed_glucose[[3]][!is.na(framingham7$GLUCOSE)])) +
               geom_density(fill=gray(0.7)),
               ggplot(mapping=aes(imputed_glucose[[4]][!is.na(framingham7$GLUCOSE)])) +
               geom_density(fill=gray(0.7)),
               ggplot(mapping=aes(imputed_glucose[[5]][!is.na(framingham7$GLUCOSE)])) +
               geom_density(fill=gray(0.7)),
               ncol=3)

table(framingham7$educ[!is.na(framingham7$educ)], imputed_educ[[1]][!is.na(framingham7$educ)])

table(framingham7$BPMEDS[!is.na(framingham7$BPMEDS)], imputed_bpmeds[[1]][!is.na(framingham7$BPMEDS)])

form1  <- CVD ~ SEX + AGE + SYSBP + TOTCHOL +
    PREVMI + DIABETES + PREVHYP + log(HEARTRTE) + CIGPDAY +
    log(BMI) + educ + log(GLUCOSE) + PREVAP + BPMEDS + PREVSTRK

model1 <- fit.mult.impute(
    form1,
    function(formula, data){glm(formula, data, family=binomial(link=logit))},
    imp,
    data=framingham7
)

# Which variables can explain CVD on their own?
# Not HEARTRTE. Kinda CIGPDAY.
### NOOOO GLM call directly!
null_model <- glm(CVD ~ 1, framingham7, family=binomial(link=logit))
add1(null_model, form1, test="LRT") |>
    ## p values get way smaller when I do this:
    ## as.data.frame() |>
    ## rownames_to_column("term") |>
    ## as.tibble() |>
    ## filter(term != "<none>") |>
    arrange(`Pr(>Chi)`, desc(LRT)) |>
    select(-AIC)

# Which variables have significance when dropped from the full model?
# Use the LRT when we don't need to estimate dispersion,
# like when doing logistic modeling.
drop1(model1, test="LRT") |>
    tidy() |>
    filter(term != "<none>") |>
    arrange(p.value, desc(LRT)) |>
    select(-AIC)


## 95% confidence intervals for parameters and their odds-ratios.
## Interpretation of OR:
## If the value of the variable were one higher,
## the odds of the person getting CVD would be OR times higher.
## Interpretation of odds:
## How many more times likely it is to get CVD than
## it is not to get CVD.
tidy(model1, conf.int=TRUE)[-1, c("term", "conf.low", "estimate", "conf.high")] |>
    arrange(estimate) |>
    mutate(
        `estimate 2.5%` = `conf.low`,
        `estimate 97.5%`= `conf.high`,
        `OR 2.5%`       = exp(`estimate 2.5%`),
        `OR`            = exp(estimate),
        `OR 97.5%`      = exp(`estimate 97.5%`),
        `conf.low`      = NULL,
        `conf.high`     = NULL,
        ) |>
    relocate(estimate, .after=`estimate 2.5%`) |>
    knitr::kable(digits=2)


## The OR for PREVMI in the above is wayyyy extreme, and so is the upper
## bound of the confidence interval of PREVCHD and the lower bound of the
## confidence interval of PREVAP. Notably the latter two includes no effect
## within the interval so the extreme ORs might just be due to few cases.
## Their p values are above 5% too.
## We see that there are only 71 cases of PREVMI, amounting to 2% of the data.
## But still the p value for PREVMI is about 1.33e-8, the third highest, so
## it has plenty of significance still. 97% of the PREVMI cases got CVD too,
## so that makes sense.
sum(as.numeric(framingham7$PREVMI)-1)
sum(as.numeric(framingham7$PREVMI)-1)/nrow(framingham7)
mean(as.numeric((framingham7 |> filter(as.numeric(PREVMI) == 2))$CVD)-1)



## Pearson residuals because their variance should be around 1.
model1_aug <- augment(model1, type.residuals = "pearson") |>
    mutate(.fitted_group = cut(.fitted, quantile(.fitted, probs = seq(5, 95, 1) / 100))) # |>



## Why are there nulls in .fitted_group?

## The estimate of the dispersion is close to 1 as it should be. Very nice.
## Based on (6.5)
model1_aug$.resid^2 |> sum() / df.residual(model1)

# Very nice .fitted around 0, bit of a weird lower tail
p1 <- model1_aug |>
    ggplot(aes(.fitted, .resid)) +
    geom_point(alpha=0.5) +
    geom_smooth() +
    xlab("fitted values") +
    ylab("Pearson residuals")

p2 <- model1_aug |>
    ggplot(aes(AGE, .resid)) +
    geom_point(alpha=0.5) +
    geom_smooth() +
    ylab("Pearson residuals")

## What is that grouping of residuals right above 0?
p3 <- model1_aug |>
    ggplot(aes(`log(GLUCOSE)`, .resid)) +
    geom_point(alpha=0.5) +
    geom_smooth(method = "lm") +
    ylab("Pearson residuals")

# Variance of pearson residuals are around 1 as they should be
p4 <- model1_aug |>
    group_by(.fitted_group) |>
    summarize(.fitted_local = mean(.fitted), .var_local = var(.resid)) |>
    na.omit() |>
    ggplot(aes(.fitted_local, .var_local)) +
    geom_point() +
    geom_hline(yintercept = 1, linetype = 2) +
    geom_smooth() +
    xlab("fitted values") +
    ylab("variance")

gridExtra::grid.arrange(p1, p2, p3, p4, ncol = 4)


## There are plenty of observations in the weird low tail
model1_aug |>
    filter(.fitted < -2) |>
    nrow()



B <- 8
sim_plots <- list()
sims <- simulate(model1, nsim=B)
for (b in 1:B) {
    fram_sim <- framingham7 |> mutate(CVD=sims[[b]])
    ## sim_glm_tmp <- glm(form1, data = fram_sim, family = binomial(link = "logit"))
    sim_glm_tmp <- fit.mult.impute(
        form1,
        function(formula, data){glm(formula, data, family=binomial(link=logit))},
        imp,
        data=fram_sim
    )
    sim_diag <- augment(sim_glm_tmp, type.residuals = "pearson") |>
        mutate(.fitted_group = cut(.fitted, quantile(.fitted, probs = seq(5, 95, 1) / 100))) # |>
    sim_plots[[b]] <- sim_diag |>
        group_by(.fitted_group) |>
        summarize(.fitted_local = mean(.fitted), .var_local = var(.resid)) |>
        na.omit() |>
        ggplot(aes(.fitted_local, .var_local)) +
        geom_point() +
        geom_hline(yintercept = 1, linetype = 2) +
        geom_smooth() +
        xlab("fitted values") +
        ylab("variance")
}

gridExtra::grid.arrange(sim_plots[[1]], sim_plots[[2]], sim_plots[[3]], sim_plots[[4]],
                        sim_plots[[5]], sim_plots[[6]], sim_plots[[7]], sim_plots[[8]],
                        ncol=4)


## Check that avoid0 is a numeric way to create the convention
## 0*log(0) == 0:
## avoid0(0)*log(avoid0(0))
avoid0 <- function(n, epsilon=0.000000000000001) {
    ifelse(n == 0,
           n + epsilon,
    ifelse(n == 1,
           n - epsilon,
           n
           )
    )
}


## Uses GA3
dev.res.sq <- function(data, muhat) {
    ## CVD is a factor...
    Y     <- as.numeric(data$CVD)-1
    m     <- rep(1, length(Y))
    term1 <- avoid0(Y)*log(avoid0(Y/(m*muhat)))
    term2 <- avoid0(m-Y)*log(avoid0((m-Y)/(m*(1-muhat))))
    dev   <- 2*(term1 + term2)
    dev
}

## Find the pearson loss too as it doesn't use GA3
cv <- function(form, data, k, m, family, loss) {
    n <- nrow(data)
    losses <- list()
    for (b in seq(m)) {
        muhat <- list()
        groups <- sample(rep(1:k, length.out=n))
        for (i in seq(k)) {
            ## TODO: Ensure all missings in train_data is present in val_data
            ## Preferably not in the same row.
            val_data   <- data[groups == i,]
            train_data <- data[groups != i,]
            train_imp <- aregImpute(
                impute_form,
                train_data,
                n.impute=5,
                nk=4
            )
            model <- fit.mult.impute(
                form,
                function(formula, data){glm(formula, data, family=binomial(link=logit))},
                train_imp,
                data=train_data
            )
            imputed_val_data <- impute.transcan(
                train_imp,
                data=train_data,
                newdata=val_data,
                imputation=(((b-1)*k+i)%%5)+1,
                list.out=TRUE,
                pr=FALSE
            )
            ## browser()
            ## muhat[groups==i] <- predict(model, newdata=imputed_val_data, type="response")
            muhat <- predict(model, newdata=imputed_val_data, type="response")
            losses[[(b-1)*k+i]] <- loss(imputed_val_data, muhat)
        }
        ## losses[[b]] <- loss(data, unlist(muhat))
    }
    losses |> unlist() |> mean()
}

## Need something to compare this against
## Do this with pearson residuals too
cv1 <- cv(
    form1,
    framingham7,
    5,
    2,
    ## 10,
    ## 20,
    family=binomial(link=logit),
    loss=dev.res.sq
)


## Which sample size does cv as done here correspond to?
## Plot cv results against PREVMI

# the total deviance is the sum of the squared deviance residuals.

## As defined in the side notes in the beginning of section 9.4
dev.loss <- function(data, mu_hat) {
    dev.res.sq(data, mu_hat)
}


unit.deviance <- function(data, mu) {
    y <- as.numeric(data$CVD)-1
    m <- 1
    2*(y*log(avoid0(y/mu)) + (m-y)*log(avoid0((m-y)/(m-mu))))
}

train.error <- function(model, loss) {
    sum(loss(model$data, fitted(model)))/nrow(model$data)
}

# Effective degrees of freedom for a GLM with the canonical link and dispersion 1.
eff.df <- function(model) {
    X      <- model.matrix(model)
    mu_hat <- fitted(model)
    V      <- mu_hat*(1-mu_hat)
    W      <- diag(V)
    S      <- X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
    sum(diag(S))
}


# In-sample error for a GLM with the canonical link and dispersion 1.
in_sample_error <- function(model) {
    train.error(model, unit.deviance) + 2*eff.df(model)/nrow(model$data)
}


ise1 <- in_sample_error(model1)


aic_risk <- function(model) {
    train.error(model, unit.deviance) + length(coef(model))/nrow(model$data)
}

aicr1 <- aic_risk(model1)

auc <- function(model) {
    eta0 <- predict(model, model$data |> filter(CVD == 0), type="response")
    eta1 <- predict(model, model$data |> filter(CVD == 1), type="response")
    mean(outer(eta0, eta1, `<`) + outer(eta0, eta1, `==`)/2)
}

auc1 <- auc(model1)

group_p <- function(y, p, n_probs = 100, filter=NULL) {
    if (!is.null(filter)) {
        y <- y[filter]
        p <- p[filter]
    }
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
    coord_equal() +
    scale_x_continuous(breaks=seq(0,1,0.1)) +
    scale_y_continuous(breaks=seq(0,1,0.1)) +
    xlab("predicted probability") +
    ylab("observed probability")


## Plot of the calibration. Looking good!
cal_plot + group_p(as.numeric(plotting_data$CVD)-1, fitted(model1))



## Calibration separated by various factors.
## The high probabilities are generally a bit too high.
gridExtra::grid.arrange(
               cal_plot + group_p(as.numeric(plotting_data$CVD)-1, fitted(model1),
                                  filter=plotting_data$SEX == "F"),
               cal_plot + group_p(as.numeric(plotting_data$CVD)-1, fitted(model1),
                                  filter=plotting_data$SEX == "M"),
               ncol=2
           )

gridExtra::grid.arrange(
               cal_plot + group_p(as.numeric(plotting_data$CVD)-1, fitted(model1),
                                  filter=plotting_data$educ == "1"),
               cal_plot + group_p(as.numeric(plotting_data$CVD)-1, fitted(model1),
                                  filter=plotting_data$educ == "2"),
               cal_plot + group_p(as.numeric(plotting_data$CVD)-1, fitted(model1),
                                  filter=plotting_data$educ == "3"),
               cal_plot + group_p(as.numeric(plotting_data$CVD)-1, fitted(model1),
                                  filter=plotting_data$educ == "4"),
               ncol=2
           )


gridExtra::grid.arrange(
               cal_plot + group_p(as.numeric(plotting_data$CVD)-1, fitted(model1),
                                  filter=plotting_data$CIGPDAY == 0),
               cal_plot + group_p(as.numeric(plotting_data$CVD)-1, fitted(model1),
                                  filter=(plotting_data$CIGPDAY > 0 & plotting_data$CIGPDAY <= 40)),
               cal_plot + group_p(as.numeric(plotting_data$CVD)-1, fitted(model1),
                                  filter=plotting_data$CIGPDAY >= 40),
               ncol=3
           )


model1_aug2 <- augment(model1, type.predict="response", type.residuals = "pearson")

## Looks good: SYSBP, AGE, TOTCHOL, HEARTRTE

# Looks kinda okay, fluctuates too much.
ggplot(data=model1_aug2, aes(x=`log(BMI)`)) +
    geom_smooth(aes(y=as.numeric(CVD)-1, color="Observed"), fill="blue", alpha=0.2) +
    geom_smooth(aes(y=.fitted, color="Fitted"),   fill="red",  alpha=0.2) +
    geom_point(aes(y=as.numeric(CVD)-1)) +
    scale_color_manual(values=c("Observed" = "blue", "Fitted"="red"))
    labs(color="Type")


## Bad tail behavior, but probably just randomness in data due to too few observations.
ggplot(data=model1_aug2, aes(x=`log(GLUCOSE)`)) +
    geom_smooth(aes(y=as.numeric(CVD)-1,     color="Observed"), fill="blue", alpha=0.2) +
    geom_smooth(aes(y=.fitted, color="Fitted"),   fill="red",  alpha=0.2) +
    geom_point(aes(y=as.numeric(CVD)-1)) +
    scale_color_manual(values=c("Observed" = "blue", "Fitted"="red"))
    labs(color="Type")

## Low tail is significant, high tail is probably noise.
framingham7 |> filter(log(GLUCOSE) < 4.25) |> nrow()
framingham7 |> filter(log(GLUCOSE) > 5.25) |> nrow()


## Something not captured here. Interactions with other risk factors?
ggplot(data=model1_aug2, aes_string(x="CIGPDAY")) +
    geom_smooth(aes(y=as.numeric(CVD)-1,     color="Observed"), fill="blue", alpha=0.2) +
    geom_smooth(aes(y=.fitted, color="Fitted"),   fill="red",  alpha=0.2) +
    geom_point(aes(y=as.numeric(CVD)-1), alpha=0.2) +
    scale_color_manual(values=c("Observed" = "blue", "Fitted"="red"))
    labs(color="Type")

framingham7 |> filter(CIGPDAY > 40) |> nrow()

## TODO: Fix the label
ggplot(model1_aug2, aes(x=DIABETES, fill=CVD)) +
    geom_bar(position="fill") +
    ## geom_boxplot(aes(y=.fitted, group=DIABETES), fill=NA) +
    stat_summary(aes(y=.fitted, group=DIABETES), fun=mean, geom="point", color="red") +
    ## geom_point(aes(y=.fitted), position=position_jitter(width=0.2), alpha=0.5) +
    labs(y="Proportion", fill="Response") +
    scale_fill_manual(values=(c("0" = "grey", "1"="steelblue")))



medians <- framingham7 |>
    ## filter(SEX == "F") |>
    ## mutate(
    ##     educ = ordered(educ, levels = c(1, 2, 3, 4), labels=c(1,2,3,4)),
    ##     SEX=NULL
    ## ) |>
    summarise(
        across(where(is.numeric), ~ median(.x, na.rm=1)),
        across(where(is.factor), ~ names(sort(table(.x), decreasing=TRUE))[1])
    )
    ## summarise(across(everything(), ~ median(.x)))


plot_marginal_predictions <- function(model, var, range) {
    grid          <- medians[rep(1,length(range)),]
    grid[["SEX"]] <- "F"
    grid[[var]]   <- range
    fittedF       <- predict(model, grid, type="response",)
    grid[["SEX"]] <- "M"
    fittedM       <- predict(model, grid, type="response")
    max_obs       <- max(plotting_data[[var]])
    min_obs       <- min(plotting_data[[var]])
    alpha         <- 0.002
    ggplot(mapping=aes(grid[[var]])) +
        geom_vline(xintercept=max_obs, linetype="dashed") +
        geom_vline(xintercept=min_obs, linetype="dashed") +
        geom_rect(aes(xmin=-Inf, xmax=min_obs, ymin=-Inf, ymax=Inf), fill="red", alpha=alpha) +
        geom_rect(aes(xmin=max_obs, xmax=Inf, ymin=-Inf, ymax=Inf), fill="red", alpha=alpha) +
        geom_rect(aes(xmin=min_obs, xmax=max_obs, ymin=-Inf, ymax=Inf), fill="green", alpha=alpha) +
        geom_line(mapping=aes(y=fittedM, color="Male")) +
        geom_line(mapping=aes(y=fittedF, color="Female")) +
        scale_color_manual(values=c("Male"="purple", "Female" = "yellow"),
                           breaks=c("Male", "Female")) +
        labs(x=var, y="Predicted probability", color="Gender")
}

plot_marginal_predictions(model1, "HEARTRTE", seq(0,200))

plot_marginal_predictions(model1, "TOTCHOL", seq(0,600))

plot_marginal_predictions(model1, "BMI", seq(0,75))

plot_marginal_predictions(model1, "GLUCOSE", seq(1,500))

plot_marginal_predictions(model1, "SYSBP", seq(1,350))



## Looks like low BMI means increased CVD risk. Not very certain though.
ggplot(plotting_data, aes(log(BMI), as.numeric(CVD)-1)) + geom_smooth(method="loess")
## There are enough low observations that it might be real.
plotting_data |> filter(log(BMI) < 3) |> nrow()

## Oooh low glucose def looks promising
ggplot(plotting_data, aes(log(GLUCOSE), as.numeric(CVD)-1)) + geom_smooth(method="loess")
## Pleenty of observations!
plotting_data |> filter(log(BMI) < 4.25) |> nrow()

ggplot(plotting_data, aes(SYSBP, as.numeric(CVD)-1)) +
    geom_smooth(method="loess")
plotting_data |> filter(SYSBP < 100)


knots_sysbp    <- as.vector(quantile(framingham7$SYSBP, probs=c(0.33,0.67), na.rm=1))
knots_log_bmi  <- as.vector(quantile(log(framingham7$BMI), probs=c(0.33,0.67), na.rm=1))
knots_log_gluc <- as.vector(quantile(log(framingham7$GLUCOSE), probs=c(0.33,0.67), na.rm=1))
bnds_sysbp     <- range(framingham7$SYSBP, na.rm=1)
bnds_log_bmi   <- range(log(framingham7$BMI), na.rm=1)
bnds_log_gluc  <- range(log(framingham7$GLUCOSE), na.rm=1)
form2  <- CVD ~ SEX + AGE +
    ns(SYSBP, knots=knots_sysbp, Boundary.knots=bnds_sysbp) +
    TOTCHOL +
    PREVMI + DIABETES + PREVHYP +
    log(HEARTRTE) +
    CIGPDAY +
    ns(log(BMI), knots=knots_log_bmi, Boundary.knots=bnds_log_bmi) +
    educ +
    ns(log(GLUCOSE), knots=knots_log_gluc, Boundary.knots=bnds_log_gluc) +
    PREVSTRK + PREVAP + BPMEDS
model2 <- fit.mult.impute(
    form2,
    function(formula, data){glm(formula, data, family=binomial(link=logit))},
    imp,
    data=framingham7
)


anova(model1, model2, test="LRT")

# Which variables can explain CVD on their own?
# Not HEARTRTE. Kinda CIGPDAY.
add1(null_model, form2, test="LRT") |>
    ## p values get way smaller when I do this:
    ## as.data.frame() |>
    ## rownames_to_column("term") |>
    ## as.tibble() |>
    ## filter(term != "<none>") |>
    arrange(`Pr(>Chi)`, desc(LRT)) |>
    select(-AIC)

## splines on glucose is definitely good. splines on
## bmi is more so-so
drop1(model2, test="LRT") |>
    tidy() |>
    filter(term != "<none>") |>
    arrange(p.value, desc(LRT)) |>
    select(-AIC)

## Glucose has an extreme OR
tidy(model2, conf.int=TRUE)[-1, c("term", "conf.low", "estimate", "conf.high")] |>
    arrange(estimate) |>
    mutate(
        `estimate 2.5%` = `conf.low`,
        `estimate 97.5%`= `conf.high`,
        `OR 2.5%`       = exp(`estimate 2.5%`),
        `OR`            = exp(estimate),
        `OR 97.5%`      = exp(`estimate 97.5%`),
        `conf.low`      = NULL,
        `conf.high`     = NULL,
        ) |>
    relocate(estimate, .after=`estimate 2.5%`) |>
    knitr::kable(digits=2)

model2_aug <- augment(model2, type.residuals = "pearson") |>
    mutate(.fitted_group = cut(.fitted, quantile(.fitted, probs = seq(5, 95, 1) / 100)))

## The estimate of the dispersion is close to 1 as it should be. Very nice.
## Based on (6.5)
model2_aug$.resid^2 |> sum() / df.residual(model2)

# Very nice .fitted around 0, bit of a weird lower tail
p5 <- model2_aug |>
    ggplot(aes(.fitted, .resid)) +
    geom_point(alpha=0.5) +
    geom_smooth() +
    xlab("fitted values") +
    ylab("Pearson residuals")

p6 <- model2_aug |>
    ggplot(aes(TOTCHOL, .resid)) +
    geom_point(alpha=0.5) +
    geom_smooth() +
    ylab("Pearson residuals")

## What is that grouping of residuals right above 0?
p7 <- model2_aug |>
    ggplot(aes(log(plotting_data$GLUCOSE), .resid)) +
    geom_point(alpha=0.5) +
    geom_smooth(method = "lm") +
    ylab("Pearson residuals")

# Variance of pearson residuals are around 1 as they should be
p8 <- model2_aug |>
    group_by(.fitted_group) |>
    summarize(.fitted_local = mean(.fitted), .var_local = var(.resid)) |>
    na.omit() |>
    ggplot(aes(.fitted_local, .var_local)) +
    geom_point() +
    geom_hline(yintercept = 1, linetype = 2) +
    geom_smooth() +
    xlab("fitted values") +
    ylab("variance")

gridExtra::grid.arrange(
               p5,
               p6, p7,
               p8, ncol = 4)


cv2 <- cv(
    form2,
    framingham7,
    5,
    2,
    ## 10,
    ## 20,
    family=binomial(link=logit),
    loss=dev.res.sq
)

ise2 <- in_sample_error(model2)

aicr2 <- aic_risk(model2)

auc2 <- auc(model2)


## Plot of the calibration. Looking good!
cal_plot + group_p(as.numeric(plotting_data$CVD)-1, fitted(model2))


gridExtra::grid.arrange(
               cal_plot + group_p(as.numeric(plotting_data$CVD)-1, fitted(model2),
                                  filter=plotting_data$SEX == "F"),
               cal_plot + group_p(as.numeric(plotting_data$CVD)-1, fitted(model2),
                                  filter=plotting_data$SEX == "M"),
               ncol=2
           )

gridExtra::grid.arrange(
               cal_plot + group_p(as.numeric(plotting_data$CVD)-1, fitted(model2),
                                  filter=plotting_data$educ == "1"),
               cal_plot + group_p(as.numeric(plotting_data$CVD)-1, fitted(model2),
                                  filter=plotting_data$educ == "2"),
               cal_plot + group_p(as.numeric(plotting_data$CVD)-1, fitted(model2),
                                  filter=plotting_data$educ == "3"),
               cal_plot + group_p(as.numeric(plotting_data$CVD)-1, fitted(model2),
                                  filter=plotting_data$educ == "4"),
               ncol=2
           )

gridExtra::grid.arrange(
               cal_plot + group_p(as.numeric(plotting_data$CVD)-1, fitted(model2),
                                  filter=plotting_data$CIGPDAY == 0),
               cal_plot + group_p(as.numeric(plotting_data$CVD)-1, fitted(model2),
                                  filter=(plotting_data$CIGPDAY > 0 & plotting_data$CIGPDAY <= 40)),
               cal_plot + group_p(as.numeric(plotting_data$CVD)-1, fitted(model2),
                                  filter=plotting_data$CIGPDAY >= 40),
               ncol=2
           )



model2_aug2 <- augment(model2, type.predict="response", type.residuals = "pearson")

ggplot(data=model2_aug2, aes(x=log(plotting_data$BMI))) +
    geom_smooth(aes(y=as.numeric(CVD)-1,     color="Observed"), fill="blue", alpha=0.2, method="loess") +
    geom_smooth(aes(y=.fitted, color="Fitted"),   fill="red",  alpha=0.2, method="loess") +
    geom_point(aes(y=as.numeric(CVD)-1)) +
    scale_color_manual(values=c("Observed" = "blue", "Fitted"="red"))
    labs(color="Type")

ggplot(data=model2_aug2, aes(x=log(plotting_data$GLUCOSE))) +
    geom_smooth(aes(y=as.numeric(CVD)-1,     color="Observed"), fill="blue", alpha=0.2, method="loess") +
    geom_smooth(aes(y=.fitted, color="Fitted"),   fill="red",  alpha=0.2, method="loess") +
    geom_point(aes(y=as.numeric(CVD)-1)) +
    scale_color_manual(values=c("Observed" = "blue", "Fitted"="red"))
    labs(color="Type")

ggplot(data=model2_aug2, aes(x=plotting_data$SYSBP)) +
    geom_smooth(aes(y=as.numeric(CVD)-1,     color="Observed"), fill="blue", alpha=0.2, method="loess") +
    geom_smooth(aes(y=.fitted, color="Fitted"),   fill="red",  alpha=0.2, method="loess") +
    geom_point(aes(y=as.numeric(CVD)-1)) +
    scale_color_manual(values=c("Observed" = "blue", "Fitted"="red"))
    labs(color="Type")


## Does this type of graph make sense for CIGPDAY?
ggplot(data=model2_aug2, aes_string(x="CIGPDAY")) +
    geom_smooth(aes(y=as.numeric(CVD)-1,     color="Observed"), fill="blue", alpha=0.2) +
    geom_smooth(aes(y=.fitted, color="Fitted"),   fill="red",  alpha=0.2) +
    geom_point(aes(y=as.numeric(CVD)-1)) +
    scale_color_manual(values=c("Observed" = "blue", "Fitted"="red"))
    labs(color="Type")


## TODO: Fix the label
ggplot(model2_aug2, aes(x=DIABETES, fill=CVD)) +
    geom_bar(position="fill") +
    ## geom_boxplot(aes(y=.fitted, group=DIABETES), fill=NA) +
    stat_summary(aes(y=.fitted, group=DIABETES), fun=mean, geom="point", color="red") +
    ## geom_point(aes(y=.fitted), position=position_jitter(width=0.2), alpha=0.5) +
    labs(y="Proportion", fill="Response") +
    scale_fill_manual(values=(c("0" = "grey", "1"="steelblue")))


plot_marginal_predictions(model2, "HEARTRTE", seq(1,200))

plot_marginal_predictions(model2, "TOTCHOL", seq(0,600))

plot_marginal_predictions(model2, "BMI", seq(1,75))

plot_marginal_predictions(model2, "GLUCOSE", seq(1,500))

## There aren't that many GLUCOSE observations above 150, so it's probably part of the
## noise discussed above. Model generalizes poorly.
framingham7 |> filter(GLUCOSE > 150)

plot_marginal_predictions(model2, "SYSBP", seq(0,350))

## Too few observations up in the high end. Model generalizes poorly.
framingham7 |> filter(SYSBP > 200)


quartiles <- function(vec) {
    cut(
        vec,
        breaks=quantile(vec, probs=seq(0,1,0.25)),
        labels=FALSE,
        include.lowest=TRUE,
        )
}

splits <- function(vec) {
    cut(
        vec,
        breaks=seq(min(vec), max(vec), (max(vec) - min(vec))/4),
        labels=FALSE,
        include.lowest=TRUE,
        )
}


## DIABETES * BMI * GLUCOSE * CIGPDAY
ggplot(plotting_data, aes(log(BMI), as.numeric(CVD)-1)) +
    facet_grid(. ~ DIABETES, label= label_both) +
    geom_smooth(method="loess")

ggplot(plotting_data, aes(GLUCOSE, as.numeric(CVD)-1)) +
    facet_grid(. ~ DIABETES , label= label_both) +
    geom_smooth(method="loess")


ggplot(plotting_data, aes(CIGPDAY, as.numeric(CVD)-1)) +
    facet_grid(. ~ quartiles(BMI) , label= label_both) +
    geom_smooth(method="loess")


## BPMEDS * PREVHYP
## Not included because bpmeds only when prevhyp
ggplot(plotting_data |> group_by(PREVHYP, BPMEDS) |> summarise(prop=mean(as.numeric(CVD)-1)),
       aes(as.numeric(PREVHYP), prop, fill=as.numeric(BPMEDS))) +
    geom_bar(stat="identity")


ggplot(plotting_data, aes(educ, as.numeric(CVD)-1#, color=SEX
                          )) +
    facet_grid(. ~ SEX, label= label_both) +
    geom_bar(stat="identity")
    ## geom_smooth(method="lm")

ggplot(plotting_data, aes(AGE, as.numeric(CVD)-1#, color=SEX
                          )) +
    facet_grid(. ~ educ, label= label_both) +
    ## geom_bar(stat="identity")
    ## geom_smooth(method="lm")
    geom_smooth(method="lm")


form3  <- CVD ~ DIABETES * log(BMI) * log(GLUCOSE)  * CIGPDAY +
    AGE + SEX * educ +
    BPMEDS + PREVHYP +
    SYSBP + TOTCHOL +
    PREVMI + log(HEARTRTE) +
    PREVSTRK + PREVAP
model3 <- fit.mult.impute(
    form3,
    function(formula, data){glm(formula, data, family=binomial(link=logit))},
    imp,
    data=framingham7
)



anova(model1, model3, test="LRT")

# Which variables can explain CVD on their own?
# Not HEARTRTE. Kinda CIGPDAY.
add1(null_model, form3, test="LRT") |>
    ## p values get way smaller when I do this:
    ## as.data.frame() |>
    ## rownames_to_column("term") |>
    ## as.tibble() |>
    ## filter(term != "<none>") |>
    arrange(`Pr(>Chi)`, desc(LRT)) |>
    select(-AIC)

## splines on glucose is definitely good. splines on
## bmi is more so-so
drop1(model3, test="LRT") |>
    tidy() |>
    filter(term != "<none>") |>
    arrange(p.value, desc(LRT)) |>
    select(-AIC)

## Glucose has an extreme OR
tidy(model3, conf.int=TRUE)[-1, c("term", "conf.low", "estimate", "conf.high")] |>
    arrange(estimate) |>
    mutate(
        `estimate 2.5%` = `conf.low`,
        `estimate 97.5%`= `conf.high`,
        `OR 2.5%`       = exp(`estimate 2.5%`),
        `OR`            = exp(estimate),
        `OR 97.5%`      = exp(`estimate 97.5%`),
        `conf.low`      = NULL,
        `conf.high`     = NULL,
        ) |>
    relocate(estimate, .after=`estimate 2.5%`) |>
    knitr::kable(digits=2)

model3_aug <- augment(model3, type.residuals = "pearson") |>
    mutate(.fitted_group = cut(.fitted, quantile(.fitted, probs = seq(5, 95, 1) / 100)))

## The estimate of the dispersion is close to 1 as it should be. Very nice.
## Based on (6.5)
model3_aug$.resid^2 |> sum() / df.residual(model3)

p9 <- model3_aug |>
    ggplot(aes(.fitted, .resid)) +
    geom_point(alpha=0.5) +
    geom_smooth() +
    xlab("fitted values") +
    ylab("Pearson residuals")

p10 <- model3_aug |>
    ggplot(aes(TOTCHOL, .resid)) +
    geom_point(alpha=0.5) +
    geom_smooth() +
    ylab("Pearson residuals")

p11 <- model3_aug |>
    ggplot(aes(log(plotting_data$GLUCOSE), .resid)) +
    geom_point(alpha=0.5) +
    geom_smooth(method = "lm") +
    ylab("Pearson residuals")

p12 <- model3_aug |>
    group_by(.fitted_group) |>
    summarize(.fitted_local = mean(.fitted), .var_local = var(.resid)) |>
    na.omit() |>
    ggplot(aes(.fitted_local, .var_local)) +
    geom_point() +
    geom_hline(yintercept = 1, linetype = 2) +
    geom_smooth() +
    xlab("fitted values") +
    ylab("variance")

gridExtra::grid.arrange(
               p9,
               p10, p11,
               p12, ncol = 4)


cv3 <- cv(
    form3,
    framingham7,
    5,
    2,
    ## 10,
    ## 20,
    family=binomial(link=logit),
    loss=dev.res.sq
)

ise3 <- in_sample_error(model3)

aicr3 <- aic_risk(model3)

auc3 <- auc(model3)

cal_plot + group_p(as.numeric(plotting_data$CVD)-1, fitted(model3))


gridExtra::grid.arrange(
               cal_plot + group_p(as.numeric(plotting_data$CVD)-1, fitted(model3),
                                  filter=plotting_data$SEX == "F"),
               cal_plot + group_p(as.numeric(plotting_data$CVD)-1, fitted(model3),
                                  filter=plotting_data$SEX == "M"),
               ncol=2
           )

gridExtra::grid.arrange(
               cal_plot + group_p(as.numeric(plotting_data$CVD)-1, fitted(model3),
                                  filter=plotting_data$educ == "1"),
               cal_plot + group_p(as.numeric(plotting_data$CVD)-1, fitted(model3),
                                  filter=plotting_data$educ == "2"),
               cal_plot + group_p(as.numeric(plotting_data$CVD)-1, fitted(model3),
                                  filter=plotting_data$educ == "3"),
               cal_plot + group_p(as.numeric(plotting_data$CVD)-1, fitted(model3),
                                  filter=plotting_data$educ == "4"),
               ncol=2
           )

gridExtra::grid.arrange(
               cal_plot + group_p(as.numeric(plotting_data$CVD)-1, fitted(model3),
                                  filter=plotting_data$CIGPDAY == 0),
               cal_plot + group_p(as.numeric(plotting_data$CVD)-1, fitted(model3),
                                  filter=(plotting_data$CIGPDAY > 0 & plotting_data$CIGPDAY <= 40)),
               cal_plot + group_p(as.numeric(plotting_data$CVD)-1, fitted(model3),
                                  filter=plotting_data$CIGPDAY >= 40),
               ncol=3
           )


model3_aug2 <- augment(model3, type.predict="response", type.residuals = "pearson")

## Looks good: SYSBP, AGE, TOTCHOL, HEARTRTE

ggplot(data=model3_aug2, aes(x=log(model3_aug2$`log(BMI)`))) +
    geom_smooth(aes(y=as.numeric(CVD)-1,     color="Observed"), fill="blue", alpha=0.2, method="loess") +
    geom_smooth(aes(y=.fitted, color="Fitted"),   fill="red",  alpha=0.2, method="loess") +
    geom_point(aes(y=as.numeric(CVD)-1)) +
    scale_color_manual(values=c("Observed" = "blue", "Fitted"="red"))
    labs(color="Type")

ggplot(data=model3_aug2, aes(x=log(plotting_data$GLUCOSE))) +
    geom_smooth(aes(y=as.numeric(CVD)-1,     color="Observed"), fill="blue", alpha=0.2, method="loess") +
    geom_smooth(aes(y=.fitted, color="Fitted"),   fill="red",  alpha=0.2, method="loess") +
    geom_point(aes(y=as.numeric(CVD)-1)) +
    scale_color_manual(values=c("Observed" = "blue", "Fitted"="red"))
    labs(color="Type")

ggplot(data=model3_aug2, aes(x=plotting_data$SYSBP)) +
    geom_smooth(aes(y=as.numeric(CVD)-1,     color="Observed"), fill="blue", alpha=0.2, method="loess") +
    geom_smooth(aes(y=.fitted, color="Fitted"),   fill="red",  alpha=0.2, method="loess") +
    geom_point(aes(y=as.numeric(CVD)-1)) +
    scale_color_manual(values=c("Observed" = "blue", "Fitted"="red"))
    labs(color="Type")


## Does this type of graph make sense for CIGPDAY?
ggplot(data=model3_aug2, aes_string(x="CIGPDAY")) +
    geom_smooth(aes(y=as.numeric(CVD)-1,     color="Observed"), fill="blue", alpha=0.2, method="loess") +
    geom_smooth(aes(y=.fitted, color="Fitted"),   fill="red",  alpha=0.2, method="loess") +
    geom_point(aes(y=as.numeric(CVD)-1)) +
    scale_color_manual(values=c("Observed" = "blue", "Fitted"="red"))
    labs(color="Type")


## TODO: Fix the label
ggplot(model3_aug2, aes(x=DIABETES, fill=CVD)) +
    geom_bar(position="fill") +
    ## geom_boxplot(aes(y=.fitted, group=DIABETES), fill=NA) +
    stat_summary(aes(y=.fitted, group=DIABETES), fun=mean, geom="point", color="red") +
    ## geom_point(aes(y=.fitted), position=position_jitter(width=0.2), alpha=0.5) +
    labs(y="Proportion", fill="Response") +
    scale_fill_manual(values=(c("0" = "grey", "1"="steelblue")))


plot_marginal_predictions(model3, "HEARTRTE", seq(1,200))

plot_marginal_predictions(model3, "TOTCHOL", seq(0,600))

plot_marginal_predictions(model3, "BMI", seq(1,75))

plot_marginal_predictions(model3, "GLUCOSE", seq(1,500))

plot_marginal_predictions(model3, "SYSBP", seq(0,350))


## left out: PREVMI + DIABETES + PREVHYP + log(HEARTRTE) +
## CIGPDAY + educ + log(GLUCOSE) + PREVAP + BPMEDS + PREVSTRK
form4  <- CVD ~ AGE + SEX +
    ns(SYSBP, knots=knots_sysbp, Boundary.knots=bnds_sysbp) +
    TOTCHOL + CIGPDAY * log(BMI)
model4 <- fit.mult.impute(
    form4,
    function(formula, data){glm(formula, data, family=binomial(link=logit))},
    imp,
    data=framingham7
)

anova(model4, model1, test="LRT")

drop1(model4, test="LRT") |>
    tidy() |>
    filter(term != "<none>") |>
    arrange(p.value, desc(LRT)) |>
    select(-AIC)

## Glucose has an extreme OR
tidy(model4, conf.int=TRUE)[-1, c("term", "conf.low", "estimate", "conf.high")] |>
    arrange(estimate) |>
    mutate(
        `estimate 2.5%` = `conf.low`,
        `estimate 97.5%`= `conf.high`,
        `OR 2.5%`       = exp(`estimate 2.5%`),
        `OR`            = exp(estimate),
        `OR 97.5%`      = exp(`estimate 97.5%`),
        `conf.low`      = NULL,
        `conf.high`     = NULL,
        ) |>
    relocate(estimate, .after=`estimate 2.5%`) |>
    knitr::kable(digits=2)

model4_aug <- augment(model4, type.residuals = "pearson") |>
    mutate(.fitted_group = cut(.fitted, quantile(.fitted, probs = seq(5, 95, 1) / 100)))

## The estimate of the dispersion is close to 1 as it should be. Very nice.
## Based on (6.5)
model4_aug$.resid^2 |> sum() / df.residual(model4)

p13 <- model4_aug |>
    ggplot(aes(.fitted, .resid)) +
    geom_point(alpha=0.5) +
    geom_smooth() +
    xlab("fitted values") +
    ylab("Pearson residuals")

p14 <- model4_aug |>
    ggplot(aes(TOTCHOL, .resid)) +
    geom_point(alpha=0.5) +
    geom_smooth() +
    ylab("Pearson residuals")

p15 <- model4_aug |>
    ggplot(aes(log(plotting_data$AGE), .resid)) +
    geom_point(alpha=0.5) +
    geom_smooth(method = "lm") +
    ylab("Pearson residuals")

p16 <- model4_aug |>
    group_by(.fitted_group) |>
    summarize(.fitted_local = mean(.fitted), .var_local = var(.resid)) |>
    na.omit() |>
    ggplot(aes(.fitted_local, .var_local)) +
    geom_point() +
    geom_hline(yintercept = 1, linetype = 2) +
    geom_smooth() +
    xlab("fitted values") +
    ylab("variance")

gridExtra::grid.arrange(p13, p14, p15, p16, ncol = 4)


cv4 <- cv(
    form4,
    framingham7,
    5,
    2,
    ## 10,
    ## 20,
    family=binomial(link=logit),
    loss=dev.res.sq
)

ise4 <- in_sample_error(model4)

aicr4 <- aic_risk(model4)

auc4 <- auc(model4)

cal_plot + group_p(as.numeric(plotting_data$CVD)-1, fitted(model4))


gridExtra::grid.arrange(
               cal_plot + group_p(as.numeric(plotting_data$CVD)-1, fitted(model4),
                                  filter=plotting_data$SEX == "F"),
               cal_plot + group_p(as.numeric(plotting_data$CVD)-1, fitted(model4),
                                  filter=plotting_data$SEX == "M"),
               ncol=2
           )

gridExtra::grid.arrange(
               cal_plot + group_p(as.numeric(plotting_data$CVD)-1, fitted(model4),
                                  filter=plotting_data$educ == "1"),
               cal_plot + group_p(as.numeric(plotting_data$CVD)-1, fitted(model4),
                                  filter=plotting_data$educ == "2"),
               cal_plot + group_p(as.numeric(plotting_data$CVD)-1, fitted(model4),
                                  filter=plotting_data$educ == "3"),
               cal_plot + group_p(as.numeric(plotting_data$CVD)-1, fitted(model4),
                                  filter=plotting_data$educ == "4"),
               ncol=2
           )

gridExtra::grid.arrange(
               cal_plot + group_p(as.numeric(plotting_data$CVD)-1, fitted(model4),
                                  filter=plotting_data$CIGPDAY == 0),
               cal_plot + group_p(as.numeric(plotting_data$CVD)-1, fitted(model4),
                                  filter=(plotting_data$CIGPDAY > 0 & plotting_data$CIGPDAY <= 40)),
               cal_plot + group_p(as.numeric(plotting_data$CVD)-1, fitted(model4),
                                  filter=plotting_data$CIGPDAY >= 40),
               ncol=3
           )


model4_aug2 <- augment(model4, type.predict="response", type.residuals = "pearson")

## Looks good: SYSBP, AGE, TOTCHOL, HEARTRTE

ggplot(data=model4_aug2, aes(x=log(plotting_data$BMI))) +
    geom_smooth(aes(y=as.numeric(CVD)-1,     color="Observed"), fill="blue", alpha=0.2, method="loess") +
    geom_smooth(aes(y=.fitted, color="Fitted"),   fill="red",  alpha=0.2, method="loess") +
    geom_point(aes(y=as.numeric(CVD)-1)) +
    scale_color_manual(values=c("Observed" = "blue", "Fitted"="red"))
    labs(color="Type")

ggplot(data=model4_aug2, aes(x=log(plotting_data$GLUCOSE))) +
    geom_smooth(aes(y=as.numeric(CVD)-1,     color="Observed"), fill="blue", alpha=0.2, method="loess") +
    geom_smooth(aes(y=.fitted, color="Fitted"),   fill="red",  alpha=0.2, method="loess") +
    geom_point(aes(y=as.numeric(CVD)-1)) +
    scale_color_manual(values=c("Observed" = "blue", "Fitted"="red"))
    labs(color="Type")

ggplot(data=model4_aug2, aes(x=plotting_data$SYSBP)) +
    geom_smooth(aes(y=as.numeric(CVD)-1,     color="Observed"), fill="blue", alpha=0.2, method="loess") +
    geom_smooth(aes(y=.fitted, color="Fitted"),   fill="red",  alpha=0.2, method="loess") +
    geom_point(aes(y=as.numeric(CVD)-1)) +
    scale_color_manual(values=c("Observed" = "blue", "Fitted"="red"))
    labs(color="Type")


## Does this type of graph make sense for CIGPDAY?
ggplot(data=model4_aug2, aes_string(x="CIGPDAY")) +
    geom_smooth(aes(y=as.numeric(CVD)-1,     color="Observed"), fill="blue", alpha=0.2, method="loess") +
    geom_smooth(aes(y=.fitted, color="Fitted"),   fill="red",  alpha=0.2, method="loess") +
    geom_point(aes(y=as.numeric(CVD)-1)) +
    scale_color_manual(values=c("Observed" = "blue", "Fitted"="red"))
    labs(color="Type")

## Plots highlighting how this model doesn't predict well on the discarded variables.

plot_marginal_predictions(model4, "TOTCHOL", seq(0,600))

plot_marginal_predictions(model4, "BMI", seq(1,75))

plot_marginal_predictions(model4, "SYSBP", seq(0,350))


dens <- ggplot(model4_aug2, aes(.fitted)) + geom_density(alpha=0.4)
gridExtra::grid.arrange(
               dens + aes(fill=plotting_data$PREVMI),
               dens + aes(fill=plotting_data$PREVSTRK),
               ncol=2
           )


dens <- ggplot(model3_aug2, aes(.fitted)) + geom_density(alpha=0.4)
gridExtra::grid.arrange(
               dens + aes(fill=plotting_data$PREVMI),
               dens + aes(fill=plotting_data$PREVSTRK),
               ncol=2
           )
