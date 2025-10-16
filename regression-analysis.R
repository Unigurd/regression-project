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

set.seed(8088)
framingham1 <- as_tibble(read.csv("training_data.csv"))


framingham2 <- framingham1 |>
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
        time_to_cvd    = TIMECVD-TIME,
        time_censor    = ifelse(max == -1, TIMEDTH, max),
        time_to_censor = time_censor-TIME,
        .after=RANDID
    ) |>
    subset(select=-c(max, timeap, timemi, timemifc, timechd, timestrk, timecvd, timehyp))

framingham3 <- framingham2 |> filter(PERIOD == min(PERIOD), .by=RANDID)

## All the people are represented
all(unique(framingham2$RANDID) == unique(framingham3$RANDID))
# exactly once.
n_distinct(framingham2$RANDID) == nrow(framingham3)

framingham4 <- framingham3 |> filter(CVD == 1 | time_to_censor > 12*365.25)
framingham5 <- framingham4 |> select(-c(HDLC, LDLC))

framingham5 |>
    skim() |>
    select(skim_variable, n_missing) |>
    filter(n_missing != 0) |>
    arrange(desc(n_missing))

table(rowSums(is.na(framingham5)))

framingham5 |>
    select(TOTCHOL, CIGPDAY, BPMEDS, GLUCOSE, educ)|>
    filter(rowSums(is.na(framingham5)) >= 2) |>
    print(n=64)


framingham6 <- framingham5 |>
    filter(!is.na(CIGPDAY), !is.na(BMI))

framingham6 |>
    skim() |>
    select(skim_variable, n_missing) |>
    filter(n_missing != 0) |>
    arrange(desc(n_missing))

table(rowSums(is.na(framingham6)))

framingham6 |>
    select(TOTCHOL, BPMEDS, GLUCOSE, educ)|>
    filter(rowSums(is.na(framingham6)) >= 2) |>
    print(n=54)


## TODO: remove outliers
framingham7 <- framingham6 |>
    mutate(
        AGE,
        TOTCHOL,
        SYSBP,
        CIGPDAY,
        BMI,
        HEARTRTE,
        GLUCOSE,
        CVD      = as.ordered(CVD),
        SEX      = ordered(SEX, levels=c(1,2), labels=c("M", "F")),
        educ     = as.ordered(educ),
        DIABETES = as.ordered(DIABETES),
        BPMEDS   = as.ordered(BPMEDS),
        PREVAP   = as.ordered(PREVAP),
        PREVMI   = as.ordered(PREVMI),
        PREVSTRK = as.ordered(PREVSTRK),
        PREVHYP  = as.ordered(PREVHYP),
        .keep="none"
    )
## select(CVD, SEX, AGE, educ, TOTCHOL, SYSBP, CIGPDAY, BMI, DIABETES,
##        BPMEDS, HEARTRTE, GLUCOSE, PREVAP, PREVMI, PREVSTRK, PREVHYP)



impute_form <- ~ CVD + SEX + TOTCHOL + AGE + SYSBP +
    CIGPDAY + DIABETES + HEARTRTE + PREVAP + PREVMI +
    PREVSTRK + PREVHYP + BMI + BPMEDS + GLUCOSE + educ

## Set n.impute accordingly later
reformM(impute_form, framingham7)

imp <- aregImpute(
    impute_form,
    framingham7,
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
    ## The use of `confint.lm` instead of `confint` is because the latter is slower.
    ## We just want approximate intervals right now.
    ## cbind(confint.lm(model1)[-1,]) |>
    ## cbind(confint(model1)[-1,]) |>
    mutate(
        `estimate 2.5%` = `conf.low`,
        `estimate 97.5%`= `conf.high`,
        `OR 2.5%`       = exp(`estimate 2.5%`),
        `OR`            = exp(estimate),
        `OR 97.5%`      = exp(`estimate 97.5%`),
        `conf.low`      = NULL,
        `conf.high`     = NULL,
        ) |>
    ## mutate(
    ##     `estimate 2.5%` = `2.5 %`,
    ##     `estimate 97.5%`= `97.5 %`,
    ##     `OR 2.5%`       = exp(`2.5 %`),
    ##     `OR`            = exp(estimate),
    ##     `OR 97.5%`      = exp(`97.5 %`),
    ##     `2.5 %`         = NULL,
    ##     `97.5 %`        = NULL,
    ##     ) |>
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


knots_sysbp <- as.vector(quantile(framingham7$SYSBP, probs=c(0.33,0.67), na.rm=1))
knots_log_bmi   <- as.vector(quantile(log(framingham7$BMI), probs=c(0.33,0.67), na.rm=1))
knots_log_gluc  <- as.vector(quantile(log(framingham7$GLUCOSE), probs=c(0.33,0.67), na.rm=1))
bnds_sysbp <- range(framingham7$SYSBP, na.rm=1)
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
    ## The use of `confint.lm` instead of `confint` is because the latter is slower.
    ## We just want approximate intervals right now.
    ## cbind(confint.lm(model2)[-1,]) |>
    ## cbind(confint(model2)[-1,]) |>
    mutate(
        `estimate 2.5%` = `conf.low`,
        `estimate 97.5%`= `conf.high`,
        `OR 2.5%`       = exp(`estimate 2.5%`),
        `OR`            = exp(estimate),
        `OR 97.5%`      = exp(`estimate 97.5%`),
        `conf.low`      = NULL,
        `conf.high`     = NULL,
        ) |>
    ## mutate(
    ##     `estimate 2.5%` = `2.5 %`,
    ##     `estimate 97.5%`= `97.5 %`,
    ##     `OR 2.5%`       = exp(`2.5 %`),
    ##     `OR`            = exp(estimate),
    ##     `OR 97.5%`      = exp(`97.5 %`),
    ##     `2.5 %`         = NULL,
    ##     `97.5 %`        = NULL,
    ##     ) |>
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
    ## The use of `confint.lm` instead of `confint` is because the latter is slower.
    ## We just want approximate intervals right now.
    ## cbind(confint.lm(model3)[-1,]) |>
    ## cbind(confint(model3)[-1,]) |>
    mutate(
        `estimate 2.5%` = `conf.low`,
        `estimate 97.5%`= `conf.high`,
        `OR 2.5%`       = exp(`estimate 2.5%`),
        `OR`            = exp(estimate),
        `OR 97.5%`      = exp(`estimate 97.5%`),
        `conf.low`      = NULL,
        `conf.high`     = NULL,
        ) |>
    ## mutate(
    ##     `estimate 2.5%` = `2.5 %`,
    ##     `estimate 97.5%`= `97.5 %`,
    ##     `OR 2.5%`       = exp(`2.5 %`),
    ##     `OR`            = exp(estimate),
    ##     `OR 97.5%`      = exp(`97.5 %`),
    ##     `2.5 %`         = NULL,
    ##     `97.5 %`        = NULL,
    ##     ) |>
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



form4  <- CVD ~ AGE + SEX + SYSBP + TOTCHOL + CIGPDAY + log(BMI)
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
    ## The use of `confint.lm` instead of `confint` is because the latter is slower.
    ## We just want approximate intervals right now.
    ## cbind(confint.lm(model4)[-1,]) |>
    ## cbind(confint(model4)[-1,]) |>
    mutate(
        `estimate 2.5%` = `conf.low`,
        `estimate 97.5%`= `conf.high`,
        `OR 2.5%`       = exp(`estimate 2.5%`),
        `OR`            = exp(estimate),
        `OR 97.5%`      = exp(`estimate 97.5%`),
        `conf.low`      = NULL,
        `conf.high`     = NULL,
        ) |>
    ## mutate(
    ##     `estimate 2.5%` = `2.5 %`,
    ##     `estimate 97.5%`= `97.5 %`,
    ##     `OR 2.5%`       = exp(`2.5 %`),
    ##     `OR`            = exp(estimate),
    ##     `OR 97.5%`      = exp(`97.5 %`),
    ##     `2.5 %`         = NULL,
    ##     `97.5 %`        = NULL,
    ##     ) |>
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


## cal_plot + group_p(framingham$CVD, fitted(model4))

## gridExtra::grid.arrange(
##                cal_plot + group_p(framingham$CVD, fitted(model4), filter=framingham$SEX == "F"),
##                cal_plot + group_p(framingham$CVD, fitted(model4), filter=framingham$SEX == "M"),
##                ncol=2
##            )

## gridExtra::grid.arrange(
##                cal_plot + group_p(framingham$CVD, fitted(model4), filter=framingham$educ == "1"),
##                cal_plot + group_p(framingham$CVD, fitted(model4), filter=framingham$educ == "2"),
##                cal_plot + group_p(framingham$CVD, fitted(model4), filter=framingham$educ == "3"),
##                cal_plot + group_p(framingham$CVD, fitted(model4), filter=framingham$educ == "4"),
##                ncol=2
##            )

## gridExtra::grid.arrange(
##                cal_plot + group_p(framingham$CVD, fitted(model4), filter=framingham$CIGPDAY == 0),
##                cal_plot + group_p(framingham$CVD, fitted(model4), filter=framingham$CIGPDAY != 0),
##                ncol=2
##            )


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
