library(ggplot2)
library(tibble)
library(broom)
library(splines)
library(readr)
library(skimr)
library(dplyr)
library(rlang)
library(hexbin)
library(DHARMa)

framingham <- as_tibble(read.csv("model_input.csv")) |>
    mutate(
        X         = NULL,
        ## GLUCOSE   = log(GLUCOSE),
        ## BMI       = log(BMI)
        ## HEARTRE   = log(HEARTRTE)
        SEX       = factor(SEX, levels = c(1, 2), labels=c("M", "F")),
        educ      = factor(educ, levels = c(1, 2, 3, 4)),
    )

fill = gray(0.7)

## Used to check what to log-transform
ggplot(data=framingham, aes(GLUCOSE)) + geom_density(fill=fill)
ggplot(data=framingham, aes(log(GLUCOSE))) + geom_density(fill=fill)

ggplot(data=framingham, aes(BMI)) + geom_density(fill=fill)
ggplot(data=framingham, aes(log(BMI))) + geom_density(fill=fill)

# NO
ggplot(data=framingham, aes(TOTCHOL)) + geom_density(fill=fill)
ggplot(data=framingham, aes(log(TOTCHOL))) + geom_density(fill=fill)

ggplot(data=framingham, aes(HEARTRTE)) + geom_density(fill=fill)
ggplot(data=framingham, aes(log(HEARTRTE))) + geom_density(fill=fill)

# Nah
ggplot(data=framingham, aes(SYSBP)) + geom_density(fill=fill)
ggplot(data=framingham, aes(log(SYSBP))) + geom_density(fill=fill)


## Variables:
## - CVD
## - SEX
## - TOTCHOL
## - AGE
## - SYSBP
## - CIGPDAY
## - BMI
## - DIABETES
## - BPMEDS
## - HEARTRTE
## - GLUCOSE
## - educ
## - PREVCHD
## - PREVAP
## - PREVMI
## - PREVHYP

## Assumptions:
## - GA1. The conditional expectation of Y given X is E(Y|X) = mu(eta)
##        with mu: R -> (0, infty) the mean value function.
##        The set J \subseteq R denotes the range of the mean
##        value function.
## - GA2. The conditional variance of Y given X is V(Y|X) = psi*V(mu(eta)),
##        with V: J -> (0, infty) the variance function and psi > 0 the
##        dispersion parameter.
## - GA3. The conditional distribution of Y given X is the
##        (theta(eta), vau_psi)-exponential dispersion distribution
##        Y|X ~ E(theta(eta), vau_psi),
##        (E is not the expectation but some fancy E)

## For logistic regression:
## g(mu)      = ln(mu/(1-mu))
## mu(eta)    = 1/(exp(-eta)+1)
## V(mu(eta)) = mu*(1-mu)
## mu'        = V(mu(eta))  -- Corollary 5.1

## Pearson residuals:  GA1 and GA2
## Deviance residuals: GA1, GA2 and GA3
## In-sample error:    A5

form1  <- CVD ~ SEX + AGE + SYSBP + TOTCHOL +
    PREVMI + DIABETES + PREVHYP + log(HEARTRTE) + CIGPDAY +
    log(BMI) + educ + log(GLUCOSE) + PREVCHD + PREVAP + BPMEDS
model1 <- glm(form1, framingham, family=binomial(link=logit))

# Which variables can explain CVD on their own?
# Not HEARTRTE. Kinda CIGPDAY.
null_model <- glm(CVD ~ 1, framingham, family=binomial(link=logit))
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
sum(framingham$PREVMI)                             # 71
sum(framingham$PREVMI)/nrow(framingham)            # 2%
(framingham |> filter(PREVMI == 1))$CVD |> mean()  # 97%

## Pearson residuals because their variance should be around 1.
model1_aug <- augment(model1, type.residuals = "pearson") |>
    mutate(.fitted_group = cut(.fitted, quantile(.fitted, probs = seq(5, 95, 1) / 100)))

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
    ggplot(aes(TOTCHOL, .resid)) +
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

gridExtra::grid.arrange(
               p1,
               p2, p3,
               p4, ncol = 4)

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
    Y     <- data$CVD
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
            val_data   <- data[groups == i,]
            train_data <- data[groups != i,]
            model <- glm(form, data=train_data, family=family)
            muhat[groups==i] <- predict(model, newdata=val_data, type="response")
        }

        losses[[b]] <- loss(data, unlist(muhat))
    }
    losses |> unlist() |> mean()
}

## Need something to compare this against
## Do this with pearson residuals too
cv1 <- cv(
    form1,
    framingham,
    10,
    20,
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
    y <- data$CVD
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

auc <- function(model) {
    eta0 <- predict(model, model$data |> filter(CVD == 0), type="response")
    eta1 <- predict(model, model$data |> filter(CVD == 1), type="response")
    mean(outer(eta0, eta1, `<`) + outer(eta0, eta1, `==`)/2)
}

auc1 <- auc(model1)

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
cal_plot + group_p(framingham$CVD, fitted(model1))

## Calibration separated by various factors.
## The high probabilities are generally a bit too high.
gridExtra::grid.arrange(
               cal_plot + group_p(framingham$CVD, fitted(model1), filter=framingham$SEX == "F"),
               cal_plot + group_p(framingham$CVD, fitted(model1), filter=framingham$SEX == "M"),
               ncol=2
           )

gridExtra::grid.arrange(
               cal_plot + group_p(framingham$CVD, fitted(model1), filter=framingham$educ == "1"),
               cal_plot + group_p(framingham$CVD, fitted(model1), filter=framingham$educ == "2"),
               cal_plot + group_p(framingham$CVD, fitted(model1), filter=framingham$educ == "3"),
               cal_plot + group_p(framingham$CVD, fitted(model1), filter=framingham$educ == "4"),
               ncol=2
           )

gridExtra::grid.arrange(
               cal_plot + group_p(framingham$CVD, fitted(model1), filter=framingham$CIGPDAY == 0),
               cal_plot + group_p(framingham$CVD, fitted(model1), filter=framingham$CIGPDAY != 0),
               ncol=2
           )


model1_aug2 <- augment(model1, type.predict="response", type.residuals = "pearson")

## Looks good: SYSBP, AGE, TOTCHOL, HEARTRTE

# Looks kinda okay, fluctuates too much.
ggplot(data=model1_aug2, aes(x=`log(BMI)`)) +
    geom_smooth(aes(y=CVD,     color="Observed"), fill="blue", alpha=0.2) +
    geom_smooth(aes(y=.fitted, color="Fitted"),   fill="red",  alpha=0.2) +
    geom_point(aes(y=CVD)) +
    scale_color_manual(values=c("Observed" = "blue", "Fitted"="red"))
    labs(color="Type")


## Bad tail behavior, but probably just randomness in data due to too few observations.
ggplot(data=model1_aug2, aes(x=`log(GLUCOSE)`)) +
    geom_smooth(aes(y=CVD,     color="Observed"), fill="blue", alpha=0.2) +
    geom_smooth(aes(y=.fitted, color="Fitted"),   fill="red",  alpha=0.2) +
    geom_point(aes(y=CVD)) +
    scale_color_manual(values=c("Observed" = "blue", "Fitted"="red"))
    labs(color="Type")


## Does this type of graph make sense for CIGPDAY?
ggplot(data=model1_aug2, aes_string(x="CIGPDAY")) +
    geom_smooth(aes(y=CVD,     color="Observed"), fill="blue", alpha=0.2) +
    geom_smooth(aes(y=.fitted, color="Fitted"),   fill="red",  alpha=0.2) +
    geom_point(aes(y=CVD)) +
    scale_color_manual(values=c("Observed" = "blue", "Fitted"="red"))
    labs(color="Type")


## TODO: Fix the label
ggplot(model1_aug2, aes(x=DIABETES, fill=factor(CVD))) +
    geom_bar(position="fill") +
    ## geom_boxplot(aes(y=.fitted, group=DIABETES), fill=NA) +
    stat_summary(aes(y=.fitted, group=DIABETES), fun=mean, geom="point", color="red") +
    ## geom_point(aes(y=.fitted), position=position_jitter(width=0.2), alpha=0.5) +
    labs(y="Proportion", fill="Response") +
    scale_fill_manual(values=(c("0" = "grey", "1"="steelblue")))


medians <- framingham |>
    ## filter(SEX == "F") |>
    ## mutate(
    ##     educ = ordered(educ, levels = c(1, 2, 3, 4), labels=c(1,2,3,4)),
    ##     SEX=NULL
    ## ) |>
    summarise(
        across(where(is.numeric), ~ median(.x)),
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
    max_obs       <- max(framingham[[var]])
    min_obs       <- min(framingham[[var]])
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


## Looks like low BMI means increased CVD risk. Not very certain though.
ggplot(framingham, aes(log(BMI), CVD)) + geom_smooth(method="loess")
## There are enough low observations that it might be real.
framingham |> filter(log(BMI) < 3) |> nrow()

## Oooh low glucose def looks promising
ggplot(framingham, aes(log(GLUCOSE), CVD)) + geom_smooth(method="loess")
## Pleenty of observations!
framingham |> filter(log(BMI) < 4.25) |> nrow()

ggplot(framingham, aes(SYSBP, CVD)) +
    geom_smooth(method="loess")
framingham |> filter(SYSBP < 100)


knots_sysbp <- c(100, 125) ## as.vector(quantile(framingham$SYSBP, probs=c(0.33,0.67)))
knots_log_bmi   <- as.vector(quantile(log(framingham$BMI), probs=c(0.33,0.67)))
knots_log_gluc  <- as.vector(quantile(log(framingham$GLUCOSE), probs=c(0.33,0.67)))
bnds_sysbp <- range(framingham$SYSBP)
bnds_log_bmi   <- range(log(framingham$BMI))
bnds_log_gluc  <- range(log(framingham$GLUCOSE))
form2  <- CVD ~ SEX + AGE +
    ns(SYSBP, knots=knots_sysbp, Boundary.knots=bnds_sysbp) +
    TOTCHOL +
    PREVMI + DIABETES + PREVHYP +
    log(HEARTRTE) +
    CIGPDAY +
    ns(log(BMI), knots=knots_log_bmi, Boundary.knots=bnds_log_bmi) +
    educ +
    ns(log(GLUCOSE), knots=knots_log_gluc, Boundary.knots=bnds_log_gluc) +
    PREVCHD + PREVAP + BPMEDS
environment(form2) <- environment()
model2 <- glm(form2, framingham, family=binomial(link=logit))

# Hmm, this is only a real change for the ~176 observations with BMI < 20,
# so it would be difficult to get low p value regardless how much better
# the change describes those points.
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
    ggplot(aes(log(framingham$GLUCOSE), .resid)) +
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
    framingham,
    10,
    20,
    family=binomial(link=logit),
    loss=dev.res.sq
)

ise2 <- in_sample_error(model2)

aicr2 <- aic_risk(model2)

auc2 <- auc(model2)

cal_plot + group_p(framingham$CVD, fitted(model2))

gridExtra::grid.arrange(
               cal_plot + group_p(framingham$CVD, fitted(model2), filter=framingham$SEX == "F"),
               cal_plot + group_p(framingham$CVD, fitted(model2), filter=framingham$SEX == "M"),
               ncol=2
           )

gridExtra::grid.arrange(
               cal_plot + group_p(framingham$CVD, fitted(model2), filter=framingham$educ == "1"),
               cal_plot + group_p(framingham$CVD, fitted(model2), filter=framingham$educ == "2"),
               cal_plot + group_p(framingham$CVD, fitted(model2), filter=framingham$educ == "3"),
               cal_plot + group_p(framingham$CVD, fitted(model2), filter=framingham$educ == "4"),
               ncol=2
           )

gridExtra::grid.arrange(
               cal_plot + group_p(framingham$CVD, fitted(model2), filter=framingham$CIGPDAY == 0),
               cal_plot + group_p(framingham$CVD, fitted(model2), filter=framingham$CIGPDAY != 0),
               ncol=2
           )


model2_aug2 <- augment(model2, type.predict="response", type.residuals = "pearson")

## Looks good: SYSBP, AGE, TOTCHOL, HEARTRTE

ggplot(data=model2_aug2, aes(x=log(framingham$BMI))) +
    geom_smooth(aes(y=CVD,     color="Observed"), fill="blue", alpha=0.2, method="loess") +
    geom_smooth(aes(y=.fitted, color="Fitted"),   fill="red",  alpha=0.2, method="loess") +
    geom_point(aes(y=CVD)) +
    scale_color_manual(values=c("Observed" = "blue", "Fitted"="red"))
    labs(color="Type")

ggplot(data=model2_aug2, aes(x=log(framingham$GLUCOSE))) +
    geom_smooth(aes(y=CVD,     color="Observed"), fill="blue", alpha=0.2, method="loess") +
    geom_smooth(aes(y=.fitted, color="Fitted"),   fill="red",  alpha=0.2, method="loess") +
    geom_point(aes(y=CVD)) +
    scale_color_manual(values=c("Observed" = "blue", "Fitted"="red"))
    labs(color="Type")

ggplot(data=model2_aug2, aes(x=framingham$SYSBP)) +
    geom_smooth(aes(y=CVD,     color="Observed"), fill="blue", alpha=0.2, method="loess") +
    geom_smooth(aes(y=.fitted, color="Fitted"),   fill="red",  alpha=0.2, method="loess") +
    geom_point(aes(y=CVD)) +
    scale_color_manual(values=c("Observed" = "blue", "Fitted"="red"))
    labs(color="Type")


## Does this type of graph make sense for CIGPDAY?
ggplot(data=model2_aug2, aes_string(x="CIGPDAY")) +
    geom_smooth(aes(y=CVD,     color="Observed"), fill="blue", alpha=0.2) +
    geom_smooth(aes(y=.fitted, color="Fitted"),   fill="red",  alpha=0.2) +
    geom_point(aes(y=CVD)) +
    scale_color_manual(values=c("Observed" = "blue", "Fitted"="red"))
    labs(color="Type")


## TODO: Fix the label
ggplot(model2_aug2, aes(x=DIABETES, fill=factor(CVD))) +
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

plot_marginal_predictions(model2, "SYSBP", seq(0,350))

## DIABETES * BMI * GLUCOSE
ggplot(framingham, aes(log(BMI), CVD)) +
    facet_grid(. ~ DIABETES, label= label_both) +
    geom_smooth(method="loess")

ggplot(framingham, aes(GLUCOSE, CVD)) +
    facet_grid(. ~ DIABETES , label= label_both) +
    geom_smooth(method="loess")

## BPMEDS * PREVHYP
ggplot(framingham |> group_by(PREVHYP, BPMEDS) |> summarise(prop=mean(CVD)),
       aes(PREVHYP, prop, fill=BPMEDS)) +
    geom_bar(stat="identity")


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

ggplot(framingham, aes(AGE, CVD, color=SEX)) +
    facet_grid(. ~ quartiles(SYSBP), label= label_both) +
    geom_smooth(method="lm")


form3  <- CVD ~ DIABETES * log(BMI) * log(GLUCOSE) +
    AGE * SYSBP +
    BPMEDS + PREVHYP +
    SEX + TOTCHOL +
    PREVMI + log(HEARTRTE) + CIGPDAY +
    educ + PREVCHD + PREVAP
model3 <- glm(form3, framingham, family=binomial(link=logit))


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
    ggplot(aes(log(framingham$GLUCOSE), .resid)) +
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
    framingham,
    10,
    20,
    family=binomial(link=logit),
    loss=dev.res.sq
)

ise3 <- in_sample_error(model3)

aicr3 <- aic_risk(model3)

auc3 <- auc(model3)

cal_plot + group_p(framingham$CVD, fitted(model3))

gridExtra::grid.arrange(
               cal_plot + group_p(framingham$CVD, fitted(model3), filter=framingham$SEX == "F"),
               cal_plot + group_p(framingham$CVD, fitted(model3), filter=framingham$SEX == "M"),
               ncol=2
           )

gridExtra::grid.arrange(
               cal_plot + group_p(framingham$CVD, fitted(model3), filter=framingham$educ == "1"),
               cal_plot + group_p(framingham$CVD, fitted(model3), filter=framingham$educ == "2"),
               cal_plot + group_p(framingham$CVD, fitted(model3), filter=framingham$educ == "3"),
               cal_plot + group_p(framingham$CVD, fitted(model3), filter=framingham$educ == "4"),
               ncol=2
           )

gridExtra::grid.arrange(
               cal_plot + group_p(framingham$CVD, fitted(model3), filter=framingham$CIGPDAY == 0),
               cal_plot + group_p(framingham$CVD, fitted(model3), filter=framingham$CIGPDAY != 0),
               ncol=2
           )


model3_aug2 <- augment(model3, type.predict="response", type.residuals = "pearson")

## Looks good: SYSBP, AGE, TOTCHOL, HEARTRTE

ggplot(data=model3_aug2, aes(x=log(framingham$BMI))) +
    geom_smooth(aes(y=CVD,     color="Observed"), fill="blue", alpha=0.2, method="loess") +
    geom_smooth(aes(y=.fitted, color="Fitted"),   fill="red",  alpha=0.2, method="loess") +
    geom_point(aes(y=CVD)) +
    scale_color_manual(values=c("Observed" = "blue", "Fitted"="red"))
    labs(color="Type")

ggplot(data=model3_aug2, aes(x=log(framingham$GLUCOSE))) +
    geom_smooth(aes(y=CVD,     color="Observed"), fill="blue", alpha=0.2, method="loess") +
    geom_smooth(aes(y=.fitted, color="Fitted"),   fill="red",  alpha=0.2, method="loess") +
    geom_point(aes(y=CVD)) +
    scale_color_manual(values=c("Observed" = "blue", "Fitted"="red"))
    labs(color="Type")

ggplot(data=model3_aug2, aes(x=framingham$SYSBP)) +
    geom_smooth(aes(y=CVD,     color="Observed"), fill="blue", alpha=0.2, method="loess") +
    geom_smooth(aes(y=.fitted, color="Fitted"),   fill="red",  alpha=0.2, method="loess") +
    geom_point(aes(y=CVD)) +
    scale_color_manual(values=c("Observed" = "blue", "Fitted"="red"))
    labs(color="Type")


## Does this type of graph make sense for CIGPDAY?
ggplot(data=model3_aug2, aes_string(x="CIGPDAY")) +
    geom_smooth(aes(y=CVD,     color="Observed"), fill="blue", alpha=0.2, method="loess") +
    geom_smooth(aes(y=.fitted, color="Fitted"),   fill="red",  alpha=0.2, method="loess") +
    geom_point(aes(y=CVD)) +
    scale_color_manual(values=c("Observed" = "blue", "Fitted"="red"))
    labs(color="Type")


## TODO: Fix the label
ggplot(model3_aug2, aes(x=DIABETES, fill=factor(CVD))) +
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
