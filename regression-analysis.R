## START import
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
## END import


library(riskCommunicator)
data(framingham, package="riskCommunicator")
dir.create("resources")



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

validation_data2 <- framingham_ordered_cols(validation_data)
framingham1      <- framingham_ordered_cols(training_data)

## In order for us to use assumption A4 about the values of CVD being
## uncorrelated conditioned on the other variables, we need to get rid
## of all but one observation for each person. Otherwise there will be
## an obvious correlation between the values of CVD for the
## observations of the same person: They will be the same. We do this
## now so that when we start plotting the distributions of our
## variables we can better justify thinking of each observation as
## independent, making the plots easier to interpret. We cannot use
## the observations where the first checkup we have for a person is
## from period 2 or 3, since it will not be independent of any
## observations from earlier periods for that person that may be in
## the validation data. We therefore just us the observations from the
## baseline exam.

## We also don't want to include people who dropped out of the study
## too early, as they had less time to get CVD, thus skewing our
## predictions for probability of getting CVD downwards. To do this we
## notice in the documentation of the data set that all the
## time-to-event columns like TIMECVD, TIMEDTH and so on register the
## real censoring time for each subject when the given event didn't
## occur. We have checked that these are all in agreement except for
## TIMEDTH which sometimes has the end-of-followup time even when the
## other time-to-event columns indicate the subject was lost to
## followup before that. The following piece of code calculates the
## observation duration under the assumption that TIME is always 0 in
## period 1 observations, which we have checked.

## Now we can remove any observations from people who dropped out of
## the study before the halfway mark. The halfway mark was chosen
## pretty arbitrarily.
framingham2 <- framingham1 |> filter(PERIOD == 1)
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
framingham4 <- framingham3 |> filter(time_censor > 12*365.25)

## Rename to framingham5
## Remove observations where TIMECVD == 0?
## Understand the high correlation between PREVCHD and PREVAP again
framingham4 <- framingham4 |> filter(PREVSTRK == 0 | is.na(PREVSTRK), PREVMI == 0 | is.na(PREVSTRK), BMI < 55 | is.na(BMI))



## We're now ready to examine the marginal distributions of the
## continuous variables:
dens1 <- ggplot(data=framingham4) + geom_density(fill=gray(0.7))
p1 <- gridExtra::grid.arrange(
               dens1 + aes(TOTCHOL),  dens1 + aes(AGE),
               dens1 + aes(SYSBP),    dens1 + aes(DIABP),
               dens1 + aes(CIGPDAY),  dens1 + aes(BMI),
               dens1 + aes(HEARTRTE), dens1 + aes(GLUCOSE),
               ncol=4
           )

ggsave("resources/marg-cont.pdf", p1, width=5000, height=2000, units="px")

## We notice that GLUCOSE, BMI and mayyybe HEARTRTE, SYSBP and DIABP
## have skewed distributions, making them candidates for
## transformations, maybe by a logarithm. We do not have the domain
## knowledge to say whether this makes sense, but on an intuitive
## level it makes sense that if something can vary like glucose does,
## then it's not whether the glucose measurements differ by e.g. 10
## that matters, but whether the measurements double or triple (This
## shouldn't be understood in a causative sense, naturally). So in
## this sense a log transform seems reasonable. We note that CIGPDAY
## is very skewed, but don't commit to any transformation of it yet.

## The log-transforms of the variables look much more reasonable, so
## we stick with that for now. The only one that doesn't improve much
## is GLUCOSE, and it keeps having a large tail regardless of how much
## we log-transform it. Maybe it's an extreme outlier causing the tail
## then.
p2 <- gridExtra::grid.arrange(
               dens1 + aes(log(SYSBP)),   dens1 + aes(log(DIABP)),
               dens1 + aes(log(BMI)),     dens1 + aes(log(HEARTRTE)),
               dens1 + aes(log(GLUCOSE)), dens1 + aes(log(log(GLUCOSE))),
               dens1 + aes(log(log(log(GLUCOSE)))), ncol=4
           )

ggsave("resources/marg-log-cont.pdf", p2, width=5000, height=2000, units="px")

## All of our categorical variables except educ are binary. educ
## tapers off to the right so higher numbers likely corresponds to
## higher levels of education, which was our assumption as well. There
## are some missing values we wil need to handle, but the most
## striking thing is how skewed DIABETES, BPMEDS, PREVCHD, PREVAP,
## PREVMI and PREVSTRK are. The latter two in particular seem to have
## so few occurrences that they're candidates for being thrown out. If
## we don't have enough observations where MI or CHD, we cannot say
## anything about their effect with certainty anyway.
bar1 <- ggplot(data=framingham4) + geom_bar(fill=gray(0.7))
p3 <- gridExtra::grid.arrange(
               bar1 + aes(SEX),      bar1 + aes(CURSMOKE),
               bar1 + aes(DIABETES), bar1 + aes(BPMEDS),
               bar1 + aes(educ),     bar1 + aes(PREVCHD),
               bar1 + aes(PREVAP),   bar1 + aes(PREVMI),
               bar1 + aes(PREVSTRK), bar1 + aes(PREVHYP),
               bar1 + aes(CVD),      ncol=4
           )

ggsave("resources/marg-factor.pdf", p3, width=5000, height=2000, units="px")


## Examining the correlation between our continuous variable, the
## strong colinearity of SYSBP and DIABP immediately catches the eye.
## One of them has to go, and we choose DIABP because it has higher
## correlation with BMI, the variables that both are second-most
## correlated to. No other correlations are high enough to remove
## variables at this point, but this plot does give us a convenient
## view of outliers. In particular, the tail of GLUCOSE is not caused
## by a few outliers, but rather it seems like a real trend in the
## data. Other than that, there are points that lie away from the
## other, but none that really necessitate removal.
framingham4_log <- framingham4 |>
    mutate(log(GLUCOSE), log(DIABP), log(SYSBP), log(BMI)) |>
    select(-c(GLUCOSE, DIABP, SYSBP, BMI))

ggally_hexbin <- function(data, mapping, ...) {ggplot(data, mapping) + stat_binhex(...)}
p4 <- GGally::ggpairs(
            framingham4_log,
            columns = c("AGE", "TOTCHOL", "log(SYSBP)", "log(DIABP)",
                        "CIGPDAY", "log(BMI)", "HEARTRTE", "log(GLUCOSE)"),
            lower = list(continuous = GGally::wrap("hexbin", bins=20)),
            diag  = list(continuous = "blankDiag")
        )

ggsave("resources/corr-cont.pdf", p4, width=2000, height=2000, units="px")

## Since all our categorical variables conveniently are ordered, we
## use the spearman correlation to get an idea of the correlation of
## all of our potentials predictors. There is naturally a strong
## correlation between CURSMOKE and CIGPDAY. If we believe that
## smoking any cigarettes at all is qualitatively different than not
## smoking we might choose to include both, but we don't have much
## reason to believe so we get rid of CURSMOKE.

pdf("resources/corr-all.pdf")


cp <- cor(
    framingham4_log |> select(AGE, TOTCHOL, `log(SYSBP)`, `log(DIABP)`,
                              CIGPDAY, `log(BMI)`, HEARTRTE, `log(GLUCOSE)`,
                              SEX, CURSMOKE, DIABETES, BPMEDS, educ,
                              PREVCHD,
                              PREVAP,
                              ## PREVMI, PREVSTRK,
                              PREVHYP , CVD
                              )
    |> data.matrix(),
    use    = "complete.obs",
    method = "spearman"
)
p5 <- corrplot::corrplot(cp, diag=FALSE, order="hclust",
                         addrect=14, tl.srt=45, tl.col="black", tl.cex=0.8)

dev.off()

## PREVCHD is very correlated to both PREVMI and PREVAP. Reading their
## descriptions we see that PREVAP is PREValent Angina Pectoris,
## PREVMI is PREValent Myocardial Infarction, and that PREVCHD is
## either of those or prevalent coronary insufficiency. Checking the
## data, we indeed see that PREVCHD is always 1 whenever PREVAP or
## PREVMI is, and that we only have two cases of PREVCHD that is not
## also PREVAP or PREVMI, that is, cases that we know must be
## prevalent coronary insufficiency. That is not information enough to
## sensibly impute prevalent coronary insufficiency for all the
## observations where either PREVAP or PREVMI is 1, so we just get rid
## of PREVCHD wholesale.
framingham4 |> filter(PREVCHD == 0, PREVMI == 1 | PREVAP == 1) |> nrow()
framingham4 |> filter(PREVCHD == 1, PREVMI == 0,  PREVAP == 0) |> nrow()

## The last strongly correlated group is SYSBP, DIABP and PREVHYP. We
## have already gotten rid of DIABP, and looking at a plot of smoothed
## CVD probability stratified by PREVHYP, we don't really see any
## reason to keep PREVHYP. The behavior seems exactly the same, just
## more fluctuating for PREVHYP=1 due to the fewer observations.
p6 <- ggplot(framingham4, aes(log(SYSBP), as.numeric(CVD), color=PREVHYP, fill=PREVHYP)) +
    geom_smooth(alpha=0.5, method="loess") ## +


ggsave("resources/sysbp-prevhyp.pdf", p6, width=5000, height=2000, units="px")


## CVD didn't really seem strongly correlated with anything so we plot
## the density of the continuous variables stratified by CVD.
density_by_cvd_legend <- ggplot(data=framingham4, aes(fill=as.factor(CVD)))  + geom_density(alpha=0.4)
density_by_cvd        <- density_by_cvd_legend + theme(legend.position="none")
p7 <- gridExtra::grid.arrange(
               density_by_cvd + aes(TOTCHOL),    density_by_cvd + aes(AGE),
               density_by_cvd + aes(log(SYSBP)), density_by_cvd_legend + aes(CIGPDAY),
               density_by_cvd + aes(log(BMI)),   density_by_cvd + aes(HEARTRTE),
               density_by_cvd + aes(log(GLUCOSE)), ncol=4
           )

ggsave("resources/resp-pred-cont.pdf", p7, width=5000, height=2000, units="px")

## We see that some of the variables have their distribution
## translated slightly higher when CVD=1, but the biggest change is
## that risk of CVD clearly increases with age.


## Looking at barplots of the relative frequencies of CVD for the
## different values of the categorical variables, we see that a lot of
## the rare conditions indicate a great increase in the risk of CVD.
## PREVMI and PREVSTRK in particular, which were the ones we
## considered removing because there were so few instances of them
## happening. Seeing this, we'll keep them in.
stacked_bars_legend <- ggplot(data=framingham4, aes(SEX, fill=CVD)) + geom_bar(position="fill", alpha=0.7) + ylab("")
stacked_bars        <- stacked_bars_legend + theme(legend.position="none")
p8 <- gridExtra::grid.arrange(
               stacked_bars + aes(SEX),         stacked_bars + aes(CURSMOKE),
               stacked_bars + aes(DIABETES),    stacked_bars + aes(BPMEDS),
               stacked_bars_legend + aes(educ), stacked_bars + aes(PREVCHD),
               stacked_bars + aes(PREVAP),      stacked_bars + aes(PREVMI),
               stacked_bars + aes(PREVSTRK),    stacked_bars + aes(PREVHYP),
               ncol=5
           )

ggsave("resources/resp-pred-factor.pdf", p8, width=5000, height=2000, units="px")

## Examining the missing values in our data, we see that GLUCOSE
## represents most of them. It makes sense that CIGPDAY and BMI are
## more likely to be missing when they're high, which would make them
## MNAR and thus not really possible to impute. But dropping those
## observations is not a good choice either. We tried examining
## whether it was possible to impute the missing values based on the
## later observations by adding them as separate columns, e.g.
## GLUCOSE2 and GLUCOSE3, but that didn't seem to improve the
## accuracy. The later measurements would often be missing too for
## those variables, and they would be missing for the people who were
## no longer part of the study, making imputation really slow. We
## therefore dropped that approach.
table1 <- framingham4 |>
    skim() |>
    select(skim_variable, n_missing) |>
    filter(n_missing != 0, skim_variable != "HDLC", skim_variable != "LDLC") |>
    arrange(desc(n_missing)) |>
    knitr::kable(digits=2)

cat(table1, file="resources/table1.txt", sep="\n")

## Luckily the missing CIGPDAY and BMI consist of less than 0.001 of
## the data, so it doesn't matter much what we do with them. Let's
## just impute them. Why not.
sum(is.na(framingham4$CIGPDAY) | is.na(framingham4$BMI))/nrow(framingham4)

table(rowSums(is.na(framingham4)))

cona4 <- framingham4 |>
    select(TOTCHOL, CIGPDAY, BPMEDS, GLUCOSE, educ)|>
    filter(rowSums(is.na(framingham4)) >= 2) |>
    as.matrix() |>
    is.na() |>
    crossprod()
round(sweep(cona4, 1, diag(cona4), "/"), 2)

## We impute our data using aregImpute from the Hmisc package as
## suggested in the chapter on missing data in the Harrell book. It
## uses bayesian bootstrap, predictive mean matching and multiple
## imputation to capture as much randomness as possible. The multiple
## imputation means it imputes multiple complete data sets and we're
## then to impute our model to all of those data sets using
## fit.mult.impute from the same package.
## START impute
impute_form <- reformM(
    ~ CVD + SEX + TOTCHOL + AGE + SYSBP + CIGPDAY + DIABETES + HEARTRTE + PREVAP +
        ## PREVMI + PREVSTRK +
        PREVHYP + BMI + BPMEDS + GLUCOSE + educ, framingham4
)
imp <- aregImpute(impute_form, framingham4, n.impute=14, nk=4)
## END impute

## We pick the first imputed data set arbitrarily to use for plotting.
plotting_data <- as_tibble(impute.transcan(
    imp,
    data=framingham4,
    imputation=1,
    list.out=TRUE,
    pr=FALSE
))


## We want to compare our imputed values against the non-missing
## values. We do this by repeatedly removing a fraction of the values
## and then imputing them, finally gathering these fractions into a
## fully imputed vector of values. Quite like K-fold cross-validation.
## START impute.all
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
imputed_glucose <- impute.all(impute_form, framingham4, "GLUCOSE", k=10)
## END impute.all

## To check our imputed values we plotted the marginal densities of
## the imputed values, the collinearities between imputed values and
## various covariates, as well as imputed values against the real
## values for the observations where the variable wasn't missing. The
## show plots here is an example of these plots. The two first look
## fine but the last doesn't really show that the imputed values hit
## the real values with any accuracy. It looks like they could just as
## well have been imputed from their marginal distributiones. We don't
## believe this to be a problem, however, it's just an indication that
## there's no real predictive signal for aregImpute to use.
p9 <- gridExtra::grid.arrange(
               ggplot(mapping=aes(imputed_glucose[[1]])) + geom_density(fill=gray(0.7)),
               ggally_hexbin(
                   cbind(GLUCOSE=log(imputed_glucose[[1]]), HEARTRTE=framingham4$HEARTRTE),
                   aes(GLUCOSE, HEARTRTE)),
               ggplot(mapping=aes(
                          framingham4$GLUCOSE[!is.na(framingham4$GLUCOSE)],
                          imputed_glucose[[1]][!is.na(framingham4$GLUCOSE)])
                      ) + geom_abline(intercept=0, slope=1, color="blue") + geom_point(),
               ncol=3)

ggsave("resources/imputed-glucose.pdf", p9, width=5000, height=2000, units="px")

## Our first model will be an additive logistic one with all the
## predictors we didn't find any reason to remove.
## START model1
form1  <- CVD ~ SEX + AGE + SYSBP + TOTCHOL +
    # PREVMI +
    DIABETES + PREVHYP + log(HEARTRTE) + CIGPDAY +
    log(BMI) + educ + log(GLUCOSE) + BPMEDS + PREVAP # + PREVSTRK
model1 <- fit.mult.impute(
    form1, function(formula, data){glm(formula, data, family=binomial(link=logit))},
    imp, data=framingham4)
## END model1


## To assess the model fit we do the usual residual plots. We use the
## pearson residuals since their variance should be around 1 if the
## variance structure is correctly specified. Note that since the
## logistic regression has a fixed mean/variance relationship,
## checking this is another check of the mean structure, not an
## independent check of the variance structure. All the plots look as
## they should and we see nothing to contradict GA1 and GA2. We do not
## have a check of GA3, but we'll allow ourselves to make use of this
## assumption later on as needed.
model1_diag <- augment(model1, type.residuals = "pearson") |>
    mutate(.fitted_group = cut(.fitted, quantile(.fitted, probs = seq(5, 95, 1) / 100))) # |>

p10 <- model1_diag |>
    ggplot(aes(.fitted, .resid)) +
    geom_point(alpha=0.5) + geom_smooth() +
    xlab("fitted values") + ylab("Pearson residuals")

p11 <- model1_diag |>
    ggplot(aes(AGE, .resid)) +
    geom_point(alpha=0.5) + geom_smooth() +
    ylab("Pearson residuals")

p12 <- model1_diag |>
    ggplot(aes(`log(GLUCOSE)`, .resid)) +
    geom_point(alpha=0.5) + geom_smooth(method = "lm") +
    ylab("Pearson residuals")

p13 <- model1_diag |>
    group_by(.fitted_group) |>
    summarize(.fitted_local = mean(.fitted), .var_local = var(.resid)) |>
    na.omit() |>
    ggplot(aes(.fitted_local, .var_local)) +
    geom_point() +
    geom_hline(yintercept = 1, linetype = 2) +
    geom_smooth() +
    xlab("fitted values") +
    ylab("variance")

p14 <- gridExtra::grid.arrange(p10, p11, p12, p13, ncol = 4)

ggsave("resources/resid1.pdf", p14, width=5000, height=2000, units="px")

## Before we started imputing values we consistently saw that the left
## tails of the variance plot was a bit too high, around 1.25. We do
## not see this in the above plot, so we simulate from our model and
## look at the variance plots.

## START sim
B <- 8
sim_plots <- list()
sims <- simulate(model1, nsim=B)
for (b in 1:B) {
    fram_sim <- framingham4 |> mutate(CVD=sims[[b]])
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
## END sim

p15 <- gridExtra::grid.arrange(
    sim_plots[[1]], sim_plots[[2]], sim_plots[[3]], sim_plots[[4]],
    sim_plots[[5]], sim_plots[[6]], sim_plots[[7]], sim_plots[[8]],
    ncol=4)

ggsave("resources/simul.pdf", p15, width=5000, height=2000, units="px")

## Overall they look alright. The left tail does veer around more than
## the rest, but that can be explained by the fact that each locally
## computed point of variance varies more for some reason. The left
## tail just doesn't seem to have converged to 1 quite as well as the
## rest yet. The high left tail in our previous models would just be
## an artifact of this, then.

## Due to the fixed mean/variance relationship in the logistic model
## we don't have to estimate the dispersion, but we do it to check
## whether our model is over/under-dispersed. It should be close to 1,
## and indeed it is. This means we can use the deviance test statistic
## rather than the F-test statistic.
## START disp1
model1_diag$.resid^2 |> sum() / df.residual(model1)
## END disp1

## We check the significance of our predictors using likelihood-ratio
## tests assuming GA3. A lot of our predictors are very significant.
## In particular PREVMI and PREVSTRK have a really high significance
## given how infrequent their occurrence is.
table2 <- drop1(model1, test="LRT") |>
    tidy() |>
    filter(term != "<none>") |>
    arrange(p.value, desc(LRT)) |>
    select(-AIC) |>
    mutate(p.value=as.character(signif(p.value, 3))) |>
    knitr::kable(digits=2)

cat(table2, file="resources/table2.txt", sep="\n")

## We now show the MLEs of the parameters and their odds-ratios along
## with 95% confidence intervals of them. The odds-ratios of PREVMI
## and PREVSTRK are striking! According to our model, the odds of
## getting CVD for a person who've already had a stroke is 72436 times
## higher than the odds for a similar person who's never had a stroke
## before. But the confidence interval is so large that the upper
## bound could not be computed, nor could even the upper bound of the
## parameter estimate. The lower bound is below 1, which means that
## the confidence interval doesn't rule out the possibility that
## having had a stroke before actually decreases the odds of getting
## CVD. There is a tension here between the low p-value of PREVSTRK in
## our previous table and the inclusion of no effect in our confidence
## interval. PREVMI has a large OR too, as would be expected, but
## nowhere near as dramatic and it does not include 1 in the
## confidence interval.
## START OR1
table3 <- tidy(model1, conf.int=TRUE)[-1, c("term", "conf.low", "estimate", "conf.high")] |>
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
## END OR1

cat(table3, file="resources/table3.txt", sep="\n")

## Our model is overall very well calibrated although we overpredict
## in the right tail. Considering how relatively few observations
## there are there, we are not surprised or concerned. We use code
## from the book to create this plot.
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


p16 <- cal_plot + group_p(as.numeric(plotting_data$CVD)-1, fitted(model1))

ggsave("resources/cal1.pdf", p16, width=5000, height=2000, units="px")

## We also looked at our CVD predictions and the smoothed observed CVD
## means against various variables. The plot that stood out to us was
## the log(GLUCOSE) plot. The left tail of the observed means goes up
## again! There are enough observations down there that we cannot just
## rule this behavior out as noise, so maybe the risk of CVD actually
## does increase when your glucose levels get too low. That does make
## sense, both too low measurements as well as too high measurements
## generally indicates that something is wrong. SYSBP shows slight
## signs of this too with the left tail going up a bit. We want to
## examine this.

estimate_plot <- function(model_diag, var, xlab) {
    ggplot(data=model_diag, aes(!!var)) +
        geom_smooth(aes(y=as.numeric(CVD)-1, color="Observed"), fill="blue", alpha=0.2, method="loess") +
        geom_smooth(aes(y=.fitted, color="Fitted"),   fill="red",  alpha=0.2, method="loess") +
        geom_point(aes(y=as.numeric(CVD)-1), alpha=0.2) +
        scale_color_manual(values=c("Observed" = "blue", "Fitted"="red")) +
        xlab(xlab) +
        labs(color="Type")
}

model1_diag2 <- augment(model1, type.predict="response", type.residuals = "pearson")
p17 <- gridExtra::grid.arrange(
    estimate_plot(model1_diag2, model1_diag2$`log(BMI)`, xlab="log(BMI)"),
    estimate_plot(model1_diag2, model1_diag2$`log(GLUCOSE)`, xlab="log(GLUCOSE)"),
    estimate_plot(model1_diag2, model1_diag2$SYSBP, xlab="SYSBP"),
    estimate_plot(model1_diag2, model1_diag2$`log(HEARTRTE)`, xlab="log(HEARTRTE)"),
    ncol=4)

ggsave("resources/estim1.pdf", p17, width=5000, height=2000, units="px")

## To do this we make plots of how our model predicts on the median
## person where we vary just one single parameter (except we always
## plot both sexes). The vertical dotted lines are the minimum and
## maximum observed value for that variable. Anything outside them is
## pure extrapolation. We show here the five continuous variables that
## represent actual measurements. It is only these where it makes
## sense that both too high and too low values are unhealthy.
medians <- framingham4 |>
    summarise(
        across(where(is.numeric), ~ median(.x, na.rm=1)),
        across(where(is.factor),  ~ names(sort(table(.x), decreasing=TRUE))[1])
    )

plot_marginal_predictions <- function(model, var, range, legend=FALSE) {
    grid          <- medians[rep(1,length(range)),]
    grid[["SEX"]] <- "F"
    grid[[var]]   <- range
    fittedF       <- predict(model, grid, type="response",)
    grid[["SEX"]] <- "M"
    fittedM       <- predict(model, grid, type="response")
    max_obs       <- max(plotting_data[[var]])
    min_obs       <- min(plotting_data[[var]])
    alpha         <- 0.002
    plot          <- ggplot(mapping=aes(grid[[var]])) +
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
    if (legend) {
        plot
    } else {
        plot + theme(legend.position="none")
    }
}

p18 <- gridExtra::grid.arrange(
    plot_marginal_predictions(model1, "HEARTRTE", seq(0,200)),
    plot_marginal_predictions(model1, "TOTCHOL", seq(0,600)),
    plot_marginal_predictions(model1, "BMI", seq(0,75)),
    plot_marginal_predictions(model1, "GLUCOSE", seq(1,500)),
    plot_marginal_predictions(model1, "SYSBP", seq(1,350), legend=TRUE),
    nrow=1
)

ggsave("resources/median1.pdf", p18, width=5000, height=2000, units="px")

## We decided to add splines to GLUCOSE, BMI and SYSBP. We didn't add
## splines to TOTCHOL because we there are too few values below 120,
## which is what would be considered unhealthily low according to
## first page search results. We do not remember why we did not
## include HEARTRTE. This might have been an oversight.



## One well-placed knot might be enough to allow the model to express
## this non-linearity we're exploring but we didn't want to fiddle
## with that, so instead we placed two knots: one at the 33th
## percentile and one at the 67th percentile. This is probably enough
## freedom to allow the bend to bend, if it exists at all.
## START model2
knots_sysbp    <- as.vector(quantile(framingham4$SYSBP,        probs=c(0.33,0.67), na.rm=1))
knots_log_bmi  <- as.vector(quantile(log(framingham4$BMI),     probs=c(0.33,0.67), na.rm=1))
knots_log_gluc <- as.vector(quantile(log(framingham4$GLUCOSE), probs=c(0.33,0.67), na.rm=1))
bnds_sysbp     <- range(framingham4$SYSBP,        na.rm=1)
bnds_log_bmi   <- range(log(framingham4$BMI),     na.rm=1)
bnds_log_gluc  <- range(log(framingham4$GLUCOSE), na.rm=1)
form2  <- CVD ~ SEX + AGE +
    ns(SYSBP, knots=knots_sysbp, Boundary.knots=bnds_sysbp) +
    TOTCHOL + # PREVMI +
    DIABETES + PREVHYP + log(HEARTRTE) + CIGPDAY +
    ns(log(BMI), knots=knots_log_bmi, Boundary.knots=bnds_log_bmi) +
    ns(log(GLUCOSE), knots=knots_log_gluc, Boundary.knots=bnds_log_gluc) +
    educ + #PREVSTRK +
    PREVAP + BPMEDS
model2 <- fit.mult.impute(
    form2, function(formula, data){glm(formula, data, family=binomial(link=logit))},
    imp, data=framingham4)
## END model2


## The estimate of the dispersion is still close to 1, and the
## residual plots look pretty indistinguishable from the first ones.
model2_diag <- augment(model2, type.residuals = "pearson") |>
    mutate(.fitted_group = cut(.fitted, quantile(.fitted, probs = seq(5, 95, 1) / 100)))

model2_diag$.resid^2 |> sum() / df.residual(model2)

p19 <- model2_diag |>
    ggplot(aes(.fitted, .resid)) +
    geom_point(alpha=0.5) +
    geom_smooth() +
    xlab("fitted values") +
    ylab("Pearson residuals")

p20 <- model2_diag |>
    ggplot(aes(TOTCHOL, .resid)) +
    geom_point(alpha=0.5) +
    geom_smooth() +
    ylab("Pearson residuals")

p21 <- model2_diag |>
    ggplot(aes(log(plotting_data$GLUCOSE), .resid)) +
    geom_point(alpha=0.5) +
    geom_smooth(method = "lm") +
    ylab("Pearson residuals")

p22 <- model2_diag |>
    group_by(.fitted_group) |>
    summarize(.fitted_local = mean(.fitted), .var_local = var(.resid)) |>
    na.omit() |>
    ggplot(aes(.fitted_local, .var_local)) +
    geom_point() +
    geom_hline(yintercept = 1, linetype = 2) +
    geom_smooth() +
    xlab("fitted values") +
    ylab("variance")

p23 <- gridExtra::grid.arrange(p19, p20, p21, p22, ncol = 4)

ggsave("resources/resid2.pdf", p23, width=5000, height=2000, units="px")

## We see that we get a p-value of 8%. It's not particularly low,
## certainly not enough to reject the additive model. Interestingly
## the p-values for our spline models were much lower before we
## started imputing missing data. It could be that the observations
## with missing values we dropped caused the data to be worse
## described by an additive model somehow, but that explanation
## doesn't really feel plausible based on gut feeling. It could be
## that the imputation doesn't work that well and it creates data that
## both models are worse at explaining.
table4 <- anova(model1, model2, test="LRT") |> knitr::kable(digits=2)

cat(table4, file="resources/table4.txt", sep="\n")

## Examining the splined variables with an LRT test, we see that the
## p-value for both BMI and SYSBP got worse, while GLUCOSE got better
## with a factor of about 10. This does make sense, the glucose plot
## was the only one that really showed evidence of this non-linear
## behavior we hypothesized.
table5 <- drop1(model2, test="LRT") |>
    tidy() |>
    filter(term != "<none>") |>
    arrange(p.value, desc(LRT)) |>
    select(-AIC) |>
    mutate(term=gsub("(,[^)]+)", "", term)) |>
    mutate(p.value=as.character(signif(p.value, 3))) |>
    filter(grepl("^ns", term)) |>
    knitr::kable(digits=2)

cat(table5, file="resources/table5.txt", sep="\n")

table5.5 <- tidy(model2, conf.int=TRUE)[-1, c("term", "conf.low", "estimate", "conf.high")] |>
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
    mutate(term=gsub("(,[^)]+)", "", term)) |>
    ## filter(grepl("^ns", term)) |>
    knitr::kable(digits=2)


cat(table5.5, file="resources/table5.5.txt", sep="\n")

## The calibration plots of the new model looks similar to the plots
## for the additive model and so we don't show them. But now we have
## two models we can compare their predictive power!


## START cv
## This function is used to make sure we don't take the logarithm of 0.
avoid0 <- function(n, epsilon=0.000000000000001) {
    ifelse(n == 0, n + epsilon, ifelse(n == 1, n - epsilon, n))}

## Squared deviance residuals
## Uses GA3
dev.res.sq <- function(data, muhat) {
    ## CVD is a factor so we need to convert back to the right numerics.
    Y     <- as.numeric(data$CVD)-1
    m     <- rep(1, length(Y))
    term1 <- avoid0(Y)*log(avoid0(Y/(m*muhat)))
    term2 <- avoid0(m-Y)*log(avoid0((m-Y)/(m*(1-muhat))))
    dev   <- 2*(term1 + term2)
    dev
}

cv <- function(form, data, k, m, family, loss) {
    n <- nrow(data)
    losses <- list()
    for (b in seq(m)) {
        muhat <- list()
        groups <- sample(rep(1:k, length.out=n))
        for (i in seq(k)) {
            val_data   <- data[groups == i,]
            train_data <- data[groups != i,]
            ## There needs to be a case of PREVSTRK
            ## Every time we train the data so we add one.
            train_data[1,"PREVSTRK"] <- "1"
            train_imp <- aregImpute(impute_form, train_data, n.impute=5, nk=4)
            model <- fit.mult.impute(
                form, function(formula, data){glm(formula, data, family=binomial(link=logit))},
                train_imp, data=train_data)
            imputed_val_data <- impute.transcan(
                train_imp, data=train_data, newdata=val_data,
                imputation=(((b-1)*k+i)%%5)+1, list.out=TRUE, pr=FALSE)
            muhat <- predict(model, newdata=imputed_val_data, type="response")
            losses[[(b-1)*k+i]] <- loss(imputed_val_data, muhat)
        }
    }
    losses |> unlist() |> mean()
}

cv1 <- cv(form1, framingham4, 10, 10, family=binomial(link=logit), loss=dev.res.sq)
cv2 <- cv(form2, framingham4, 10, 10, family=binomial(link=logit), loss=dev.res.sq)
## END cv


## TODO: write this based on actual computations
## They're pretty similar


## We also look at the estimated in-sample error and AUC for both
## models. The estimate of the generalization error based on AIC is
## the same as the estimate of the in-sample error but using the spent
## degrees of freedom rather than the effective degrees of freedom, so
## we omit that.

## START pred-error
train.error <- function(model, loss) {
    sum(loss(model$data, fitted(model)))/nrow(model$data)
}

eff.df <- function(model) {
    X      <- model.matrix(model)
    mu_hat <- fitted(model)
    V      <- mu_hat*(1-mu_hat)
    W      <- diag(V)
    S      <- X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
    sum(diag(S))
}

in_sample_error <- function(model) {
    train.error(model, dev.res.sq) + 2*eff.df(model)/nrow(model$data)
}

auc <- function(model) {
    eta0 <- predict(model, model$data |> filter(CVD == 0), type="response")
    eta1 <- predict(model, model$data |> filter(CVD == 1), type="response")
    mean(outer(eta0, eta1, `<`) + outer(eta0, eta1, `==`)/2)
}

ise1 <- in_sample_error(model1)
ise2 <- in_sample_error(model2)
auc1 <- auc(model1)
auc2 <- auc(model2)
## END pred-error

## The in-sample error estimates are remarkably similar, suggesting
## that either model generalizes equally well. The auc is slightly
## higher for the spline model, suggesting that it would be a bit
## easier to set thresholds for the spline model's predictions.
model2_diag2 <- augment(model2, type.predict="response", type.residuals = "pearson")

p24 <- gridExtra::grid.arrange(
    estimate_plot(model2_diag2, log(plotting_data$BMI), xlab="log(BMI)"),
    estimate_plot(model2_diag2, log(plotting_data$GLUCOSE), xlab="log(GLUCOSE)"),
    estimate_plot(model2_diag2, plotting_data$SYSBP, xlab="SYSBP"),
    estimate_plot(model2_diag2, log(plotting_data$HEARTRTE), xlab="log(HEARTRTE)"),
    ncol=4)

ggsave("resources/estim2.pdf", p24, width=5000, height=2000, units="px")

## But we do see that GLUCOSE in particular, but also BMI, exhibit the
## behavior we hypothesized. At a glance the spline model seem to
## behave reasonably outside the observed range of the variables too.
p25 <- gridExtra::grid.arrange(
    plot_marginal_predictions(model2, "BMI", seq(1,75)),
    plot_marginal_predictions(model2, "GLUCOSE", seq(1,500)),
    plot_marginal_predictions(model2, "SYSBP", seq(1,350), legend=TRUE),
    nrow=1)

ggsave("resources/median2.pdf", p25, width=5000, height=2000, units="px")

## But looking back at the plots exploring the predictions of GLUCOSE
## in the additive model, we do see a bit of an upwards slope in the
## left tail. This suggests that other variables have values that make
## the predictions higher, even though GLUCOSE is low. Maybe there are
## interactions we can find that makes the predictions follow the
## observations on the GLUCOSE plot.

quartiles <- function(vec) {
    cut(
        vec, breaks=quantile(vec, probs=seq(0,1,0.25)),
        labels=FALSE, include.lowest=TRUE,
        )
}

p26 <- gridExtra::grid.arrange(
    ggplot(plotting_data, aes(log(BMI), as.numeric(CVD)-1, color=DIABETES)) + geom_smooth(method="loess"),
    ggplot(plotting_data, aes(GLUCOSE, as.numeric(CVD)-1, color=DIABETES)) + geom_smooth(method="loess"),
    ggplot(plotting_data, aes(CIGPDAY, as.numeric(CVD)-1, color=as.ordered(quartiles(BMI)))) + geom_smooth(method="loess"),
    ggplot(plotting_data |>
           group_by(educ, SEX) |>
           summarise(meancvd=mean(as.numeric(CVD)-1), .groups="drop"),
           aes(educ, meancvd, color=SEX, fill=SEX)) + geom_bar(stat="identity", position="dodge"),
    nrow=1)

ggsave("resources/int-expl.pdf", p26, width=5000, height=2000, units="px")

## These are the best candidates for interactions we found. We also
## throw in an interaction between educ and SEX just because we
## noticed how the mean of CVD goes up again for men in educ groups 3
## and 4 but not for women in those groups. Let's get modeling!
## START model3
form3  <- CVD ~ DIABETES * log(BMI) * log(GLUCOSE)  * CIGPDAY +
    AGE + SEX * educ +
    BPMEDS + PREVHYP +
    SYSBP + TOTCHOL +
    # PREVMI +
    log(HEARTRTE) +
    # PREVSTRK +
    PREVAP
model3 <- fit.mult.impute(
    form3, function(formula, data){glm(formula, data, family=binomial(link=logit))},
    imp, data=framingham4)
## END model3


## The model diagnostics look just fine again, so we examine p-values.
## We see on the p-value from the LRT test that the interaction model
## would explain the data this much better only 1% of the time given
## that the additive model is the true one. This is an improvement
## from the spline model. The p-value for the interaction centered
## around GLUCOSE is pretty high, certainly much higher than the
## p-value for just GLUCOSE in the additive model. The p-value for the
## SEX*educ interaction is somewhat middling, but it is an interaction
## between two variables that had a very low and a pretty high p-value
## in the additive model.
table6 <- anova(model1, model3, test="LRT") |> knitr::kable(digits=2)
table7 <- drop1(model3, test="LRT") |>
    tidy() |>
    filter(term != "<none>") |>
    arrange(p.value, desc(LRT)) |>
    select(-AIC)


cat(table6, file="resources/table6.txt", sep="\n")
cat(table7, file="resources/table7.txt", sep="\n")


## But looking at the odds-ratios, we now see an extreme OR for plain
## DIABETES of 5.6e+39 with a 95% confidence interval that goes right
## down to 0.00 after rounding. This does feel quite wonky. A lot of
## the other variables produced by the interaction does counteract it,
## but their confidence intervals also includes the possibity that
## they don't have any effect. This is likely just all because there
## are too few observations of DIABETES to include it in an
## interaction like this.
table8 <- tidy(model3, conf.int=TRUE)[-1, c("term", "conf.low", "estimate", "conf.high")] |>
    arrange(estimate) |>
    mutate(
        `estim. 2.5%`   = `conf.low`,
        `estim.`        = `estimate`,
        `estim. 97.5%`  = `conf.high`,
        `OR 2.5%`       = exp(`estim. 2.5%`),
        `OR`            = exp(estim.),
        `OR 97.5%`      = exp(`estim. 97.5%`),
        `conf.low`      = NULL,
        `conf.high`     = NULL,
        `estimate`      = NULL,
        ) |>
    ## relocate(estim., .after=`estim. 2.5%`) |>
    mutate(term=gsub("\\.L", "", term)) |>
    mutate(term=gsub("(log\\(([^)]*)\\))", "\\2", term)) |>
    filter(grepl("DIABETES|BMI|GLUCOSE|CIGPDAY", term)) |>
    knitr::kable(digits=2)

cat(table8, file="resources/table8.txt", sep="\n")

## We do see the left tail of GLUCOSE going up a bit, but not as much
## as the observed mean or it did in model 2. It goes up too much for
## BMI.
model3_diag2 <- augment(model3, type.predict="response", type.residuals = "pearson")

p27 <- gridExtra::grid.arrange(
    estimate_plot(model3_diag2, log(plotting_data$BMI), xlab="log(BMI)"),
    estimate_plot(model3_diag2, log(plotting_data$GLUCOSE), xlab="log(GLUCOSE)"),
    nrow=1)

ggsave("resources/estim3.pdf", p27, width=5000, height=2000, units="px")

## This time we see that the predictions for BMI and GLUCOSE by
## themselves increase monotonically. So this model claims that low
## GLUCOSE and BMI don't increase the risk of CVD themselves, but
## GLUCOSE interacts with other variables in a way such that in
## practice you do have higher risk of CVD if you have low glucose.
p28 <- gridExtra::grid.arrange(
    plot_marginal_predictions(model3, "BMI", seq(1,75)),
    plot_marginal_predictions(model3, "GLUCOSE", seq(1,500)),
    nrow=1)

ggsave("resources/median3.pdf", p28, width=5000, height=2000, units="px")

## And the calibration is even better than the two previous models.
p29 <- cal_plot + group_p(as.numeric(plotting_data$CVD)-1, fitted(model3))

ggsave("resources/cal3.pdf", p29, width=5000, height=2000, units="px")

## So far we've only been discussing models with all the variables we
## deemed to be possible predictors. The variables have been chosen
## ad-hoc in the modelling process. Roughly they have been chosen
## because they are mostly actual measurements as would be done at a
## health exam, and we had a feeling they would have an effect early
## on in the process. We've have included most of the interaction from
## the previous model, omitting DIABETES since we believe it made the
## odds-ratios crazy.
## START model4
form4  <- CVD ~ AGE + SEX + SYSBP + TOTCHOL + CIGPDAY * log(BMI) * log(GLUCOSE)
model4 <- fit.mult.impute(
    form4, function(formula, data){glm(formula, data, family=binomial(link=logit))},
    imp, data=framingham4)
## END model4

## Unsurprisingly the additive model is very significant. Apart from
## this the model behaves pretty sane. Perfectly fine residuals, much
## more reasonable odds-ratios. We'll skip straight to the
## generalization error estimates for this as well as the full
## interaction model.
table9 <- anova(model4, model1, test="LRT") |> knitr::kable(digits=2)

cat(table9, file="resources/table9.txt", sep="\n")


cv3 <- cv(form3, framingham4, 10, 10, family=binomial(link=logit), loss=dev.res.sq)
cv4 <- cv(form4, framingham4, 10, 10, family=binomial(link=logit), loss=dev.res.sq)
ise3 <- in_sample_error(model3)
ise4 <- in_sample_error(model4)
auc3 <- auc(model3)
auc4 <- auc(model4)

errors_df <- data.frame(CV=c(cv1,cv2,cv3,cv4), ISE=c(ise1, ise2, ise3, ise4), AUC=c(auc1, auc2, auc3, auc4))
rownames(errors_df) <- c("Additive", "Splines", "Interactions", "Fewer variables")
table10 <- errors_df |> knitr::kable(digits=3)

cat(table10, file="resources/table10.txt", sep="\n")



## We see the interaction model seems to predict the best and the
## model with fewer variables the worst. The latter does provide a
## nice sense of how big a change matter. Before it was pretty
## difficult to tell if a change from 0.92 to 0.91 matters at all, but
## with this reduced model we get the sense that while it's not a wild
## change, it does indeed indicate something real.


## This model captures the left tail behavior of GLUCOSE just about as
## well as the full interaction model did, and the left tail of BMI
## goes up too high as well. On the variables it does model it
## performs about just as well.
model4_diag2 <- augment(model4, type.predict="response", type.residuals = "pearson")
p30 <- gridExtra::grid.arrange(
    estimate_plot(model4_diag2, log(plotting_data$BMI)),
    estimate_plot(model4_diag2, log(plotting_data$GLUCOSE)),
    nrow=1)

ggsave("resources/estim4.pdf", p30, width=5000, height=2000, units="px")

## But we also see how it predicts markedly differently for important
## indicators we didn't include, like PREVMI and PREVSTRK.
dens3 <- ggplot(model3_diag2, aes(.fitted)) + geom_density(alpha=0.4)
dens4 <- ggplot(model4_diag2, aes(.fitted)) + geom_density(alpha=0.4)
p31 <- gridExtra::grid.arrange(
    dens4 + aes(fill=plotting_data$PREVMI),
    dens4 + aes(fill=plotting_data$PREVSTRK),
    dens3 + aes(fill=plotting_data$PREVMI),
    dens3 + aes(fill=plotting_data$PREVSTRK),
    nrow=1
)

ggsave("resources/median4.pdf", p31, width=5000, height=2000, units="px")
