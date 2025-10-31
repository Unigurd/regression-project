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
library(ggforce)
library(riskCommunicator)
## END import

dir.create("resources")


## START intro
data(framingham, package="riskCommunicator")
set.seed(8088)
n               <- nrow(framingham)
validation_size <- floor(n/10)
i               <- sample(n, validation_size)
validation_data <- framingham[i,]
training_data   <- framingham[-i,]

fix_framingham <- function(data) {
    data |> mutate(
                ## Factorize
                SEX  = ordered(SEX, levels=c(1,2), labels=c("M", "F")),
                CVDF = as.ordered(CVD),
                across(c(CURSMOKE, BPMEDS, PREVCHD, PREVMI, PREVHYP, DEATH, HOSPMI, ANYCHD, educ), as.ordered),
                cig=cut(CIGPDAY, breaks=c(-1,0,10, 20, Inf),
                        labels=c("0", "1-10", "11-20", "21+")),
                ## Calculate censor time
                timeap   = ifelse(ANGINA   == 1, -1, TIMEAP),   timemi   = ifelse(HOSPMI   == 1, -1, TIMEMI),
                timemifc = ifelse(MI_FCHD  == 1, -1, TIMEMIFC), timechd  = ifelse(ANYCHD   == 1, -1, TIMECHD),
                timestrk = ifelse(STROKE   == 1, -1, TIMESTRK), timecvd  = ifelse(CVD      == 1, -1, TIMECVD),
                timehyp  = ifelse(HYPERTEN == 1, -1, TIMEHYP),
                max      = pmax(timeap, timemi, timemifc, timechd, timestrk, timecvd, timehyp),
                ## if max == -1 then all of the time-to-event cols had the event happen. In that
                ## case we use TIMEDTH instead, since it always contains the time of last contact.
                time_censor = ifelse(max == -1, TIMEDTH, max),
                ) |>
        subset(select=-c(max, timeap, timemi, timemifc, timechd, timestrk, timecvd, timehyp)) |>
        filter (PERIOD == 1, time_censor > 12*365.25,
                PREVSTRK == 0 | is.na(PREVSTRK), PREVMI == 0 | is.na(PREVSTRK),
                BMI < 55 | is.na(BMI))
}

framingham4 <- fix_framingham(training_data)
## END intro

dens1 <- ggplot(data=framingham4) + geom_density(fill=gray(0.7))
p1 <- gridExtra::grid.arrange(
               dens1 + aes(TOTCHOL),  dens1 + aes(AGE),
               dens1 + aes(SYSBP),    dens1 + aes(DIABP),
               dens1 + aes(CIGPDAY),  dens1 + aes(BMI),
               dens1 + aes(HEARTRTE), dens1 + aes(GLUCOSE),
               ncol=4
           )

ggsave("resources/marg-cont.pdf", p1, width=5000, height=2000, units="px")


bar1 <- ggplot(data=framingham4) + geom_bar(fill=gray(0.7))
p3 <- gridExtra::grid.arrange(
               bar1 + aes(SEX),      bar1 + aes(CURSMOKE),
               bar1 + aes(DIABETES), bar1 + aes(BPMEDS),
               bar1 + aes(educ),     bar1 + aes(PREVAP),
               bar1 + aes(CVDF) + xlab("CVD"), nrow=2
           )

ggsave("resources/marg-factor.pdf", p3, width=5000, height=2000, units="px")


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


pdf("resources/corr-all.pdf")
cp <- cor(
    framingham4_log |> select(AGE, TOTCHOL, `log(SYSBP)`, `log(DIABP)`,
                              CIGPDAY, `log(BMI)`, HEARTRTE, `log(GLUCOSE)`,
                              SEX, CURSMOKE, DIABETES, BPMEDS, educ,
                              PREVAP, CVD)
    |> data.matrix(),
    use    = "complete.obs",
    method = "spearman"
)
p5 <- corrplot::corrplot(cp, diag=FALSE, order="hclust",
                         addrect=10, tl.srt=45, tl.col="black", tl.cex=0.8)
dev.off()

density_by_cvd_legend <- ggplot(data=framingham4, aes(fill=CVDF)) + labs(fill="CVD") + geom_density(alpha=0.4)
density_by_cvd        <- density_by_cvd_legend + theme(legend.position="none")
p7 <- gridExtra::grid.arrange(
               density_by_cvd + aes(TOTCHOL),    density_by_cvd + aes(AGE),
               density_by_cvd + aes(log(SYSBP)), density_by_cvd_legend + aes(CIGPDAY),
               density_by_cvd + aes(log(BMI)),   density_by_cvd + aes(HEARTRTE),
               density_by_cvd + aes(log(GLUCOSE)), ncol=4
           )

ggsave("resources/resp-pred-cont.pdf", p7, width=5000, height=2000, units="px")


stacked_bars_legend <- ggplot(data=framingham4, aes(SEX, fill=CVDF)) + labs(fill="CVD") + geom_bar(position="fill", alpha=0.7) + ylab("")
stacked_bars        <- stacked_bars_legend + theme(legend.position="none")
p8 <- gridExtra::grid.arrange(
               stacked_bars + aes(SEX),         stacked_bars + aes(CURSMOKE),
               stacked_bars_legend + aes(DIABETES),    stacked_bars + aes(BPMEDS),
               stacked_bars_legend + aes(educ), stacked_bars + aes(PREVAP),
               nrow=2
           )

ggsave("resources/resp-pred-factor.pdf", p8, width=5000, height=2000, units="px")

table1 <- framingham4 |>
    skim() |>
    select(skim_variable, n_missing) |>
    filter(n_missing != 0, skim_variable != "HDLC", skim_variable != "LDLC") |>
    arrange(desc(n_missing)) |>
    knitr::kable(digits=2)

cat(table1, file="resources/table1.txt", sep="\n")

sum(is.na(framingham4$CIGPDAY) | is.na(framingham4$BMI))/nrow(framingham4)

table(rowSums(is.na(framingham4)))

cona4 <- framingham4 |>
    select(TOTCHOL, CIGPDAY, BPMEDS, GLUCOSE, educ)|>
    filter(rowSums(is.na(framingham4)) >= 2) |>
    as.matrix() |>
    is.na() |>
    crossprod()
round(sweep(cona4, 1, diag(cona4), "/"), 2)

## START impute
impute_form <- reformM(
    ~ CVDF + SEX + TOTCHOL + AGE + SYSBP + cig + DIABETES + HEARTRTE + PREVAP +
         BMI + BPMEDS + GLUCOSE + educ, framingham4
)
imp <- aregImpute(impute_form, framingham4, n.impute=14, nk=4)
## END impute

## START impute.all
impute.all <- function(form, data, var, k, n.impute=14, nk=4) {
    imputed <- list()
    for (i in seq(n.impute)) {imputed[i] <- list(rep(NA, nrow(data)))}
    groups <- sample(rep(1:k, length.out=nrow(data)))
    for (i in seq(k)) {
        data_na <- mutate(data, !!var:=ifelse(groups == i, NA, data[[var]]))
        imp <- aregImpute(form, data_na, n.impute=n.impute, nk=nk)
        for (j in seq(n.impute)) {
            imputed[[j]][groups==i] <- impute.transcan(
                imp, data=data_na, imputation=j, list.out=TRUE, pr=FALSE)[[var]][groups==i]
        }
    }
    imputed
}
imputed_glucose <- impute.all(impute_form, framingham4, "GLUCOSE", k=10)
## END impute.all

p9 <- gridExtra::grid.arrange(
                     ggplot(mapping=aes(imputed_glucose[[1]])) + geom_density(fill=gray(0.7))
                     + xlab("imputed GLUCOSE"),
                     ggally_hexbin(
                         cbind(GLUCOSE=log(imputed_glucose[[1]]), HEARTRTE=framingham4$HEARTRTE),
                         aes(GLUCOSE, HEARTRTE)),
                     ## ggally_hexbin(
                     ##     cbind(`imputed GLUCOSE`=imputed_glucose[[1]][!is.na(framingham4$GLUCOSE)],
                     ##           GLUCOSE=framingham$GLUCOSE[!is.na(framingham4$GLUCOSE)]),
                     ##     aes(`imputed GLUCOSE`, GLUCOSE)),
                     ggplot(cbind(GLUCOSE=framingham4$GLUCOSE[!is.na(framingham4$GLUCOSE)],
                                  `imputed GLUCOSE`=imputed_glucose[[1]][!is.na(framingham4$GLUCOSE)]),
                            mapping=aes(GLUCOSE, `imputed GLUCOSE`)
                            ) + geom_abline(intercept=0, slope=1, color="blue") + geom_point(alpha=0.5),
                     ncol=3)
ggsave("resources/imputed-glucose.pdf", p9, width=5000, height=2000, units="px")

## START model1
form1  <- CVDF ~ SEX + AGE + SYSBP + TOTCHOL + DIABETES +  cig +
    log(BMI) + educ + log(GLUCOSE) + BPMEDS + PREVAP + log(HEARTRTE)
model1 <- fit.mult.impute(form1, xtrans=imp, data=framingham4,
    fitter=function(formula, data){glm(formula, data, family=binomial(link=logit))})

model.diag <- function(model) {
    augment(model, type.residuals = "pearson") |>
        mutate(.fitted_group = cut(.fitted, quantile(.fitted, probs = seq(5, 95, 1) / 100)))
}

disp1 <- model.diag(model1)$.resid^2 |> sum() / df.residual(model1)

plot.resid <- function(model) {
    ggplot(model.diag(model), aes(y=.resid)) +
        geom_point(alpha=0.5) + geom_smooth(method="loess", span=1.0) + ylab("Pearson residuals")
}

plot.resid.var <- function(model) {
    model.diag(model) |> group_by(.fitted_group) |>
        summarize(.fitted_local = mean(.fitted), .var_local = var(.resid)) |>
        na.omit() |> ggplot(aes(.fitted_local, .var_local)) +
        geom_point() + geom_smooth(method="loess") +
        geom_hline(yintercept = 1, linetype = 2) +
        xlab("fitted values") + ylab("variance")
}

p14 <- gridExtra::grid.arrange(plot.resid(model1) + aes(.fitted),        plot.resid(model1) + aes(log(SYSBP)),
                               plot.resid(model1) + aes(`log(GLUCOSE)`), plot.resid.var(model1), nrow=1)
## END model1
ggsave("resources/resid1.pdf", p14, width=5000, height=2000, units="px")

## START sim
B <- 8
sim_plots <- list()
sims <- simulate(model1, nsim=B)
for (b in 1:B) {
    fram_sim <- framingham4 |> mutate(CVDF=sims[[b]])
    sim_glm_tmp <- fit.mult.impute(form1, xtrans=imp, data=fram_sim,
        fitter=function(formula, data){glm(formula, data, family=binomial(link=logit))})
    sim_plots[[b]] <- plot.resid.var(sim_glm_tmp)
}
## END sim

p15 <- gridExtra::grid.arrange(
    sim_plots[[1]], sim_plots[[2]], sim_plots[[3]], sim_plots[[4]],
    sim_plots[[5]], sim_plots[[6]], sim_plots[[7]], sim_plots[[8]],
    ncol=4)

ggsave("resources/simul.pdf", p15, width=5000, height=2000, units="px")


## START plot_pred
eta <- 0
plot.pred <- function(model) {
    data <- model$data |> mutate(.fitted = predict(model))
    ggplot(data |> filter(.fitted > eta)) +
        geom_point(aes(y=CVD), alpha=0.2) +
        geom_smooth(aes(y=CVD), method="loess", color="blue", se=FALSE, span=1.0) +
        geom_smooth(data=data |> filter(.fitted <= eta), mapping=aes(y=CVD), method="loess", color="green", se=FALSE, span=1.0) +
        geom_smooth(aes(y=1/(1+exp(-.fitted))), method="loess", color="red", se=FALSE, span=1.0) +
        geom_hline(yintercept=1/(1+exp(-eta)), linetype="dotted") #+
}
p16 <- gridExtra::grid.arrange(
               plot.pred(model1) + aes(TOTCHOL),  plot.pred(model1) + aes(AGE),
               plot.pred(model1) + aes(SYSBP),    plot.pred(model1) + aes(CIGPDAY),
               plot.pred(model1) + aes(BMI),      plot.pred(model1) + aes(HEARTRTE),
               plot.pred(model1) + aes(GLUCOSE),  nrow=2
)
## END plot_pred
ggsave("resources/plot_pred.pdf", p16, width=4000, height=2000, units="px")


## START int-expl
tmpdens <- ggplot(framingham4 |> mutate(.fitted = predict(model1, type="response")) |> filter(!is.na(GLUCOSE)),
                  aes(y=CVD)) +
    geom_point() + geom_smooth(method="loess", span=1.0, na.rm=TRUE)
p17 <- gridExtra::grid.arrange(tmpdens + aes(TOTCHOL, color=SYSBP<150),
                               tmpdens + aes(log(HEARTRTE), color=GLUCOSE<100), nrow=1)
## END int-expl
ggsave("resources/int-expl.pdf", p17, width=4000, height=2000, units="px")



## START model2
form2  <- CVDF ~ log(GLUCOSE) * log(HEARTRTE) + TOTCHOL * SYSBP +
    SEX + AGE + log(BMI) + SYSBP + DIABETES + cig + educ + BPMEDS + PREVAP
model2 <- fit.mult.impute(form2, xtrans=imp, data=framingham4,
    fitter=function(formula, data){glm(formula, data, family=binomial(link=logit))})
disp2 <- model.diag(model2)$.resid^2 |> sum() / df.residual(model2)
p18 <- gridExtra::grid.arrange(plot.resid(model2) + aes(.fitted), plot.resid(model2) + aes(log(SYSBP)),
                               plot.resid(model2) + aes(`log(GLUCOSE)`), plot.resid.var(model2), nrow=1)
## END model2
ggsave("resources/model2.pdf", p18, width=4000, height=2000, units="px")

## START tables2
p.value.table <- function(model) {
    drop1(model, test="LRT") |> tidy() |> filter(term != "<none>") |> arrange(p.value, desc(LRT)) |>
        select(-c(AIC, deviance, LRT)) |> mutate(p.value=as.character(signif(p.value, 3)))
}

conf.int.table <- function(model, threshold=-Inf) {
    tidy(model, conf.int=TRUE)[-1, c("term", "conf.low", "estimate", "conf.high")] |>
        mutate(`est.` = estimate, `est. 2.5%` = `conf.low`, `est. 97.5%`= `conf.high`,
               `OR   2.5%` = exp(`est. 2.5%`), `OR` = exp(estimate), `OR   97.5%`= exp(`est. 97.5%`),
               term=term, `conf.low` = NULL, `conf.high` = NULL, estimate=NULL) |>
        arrange(est.) |> relocate(est., .after=`est. 2.5%`) |> filter(abs(est.) > threshold)
}

table2p <- p.value.table(model2)
table2c <- conf.int.table(model2, threshold=1)
## END tables2
cat(table2p |> mutate(term=gsub("(,[^)]+)", "", term)) |> knitr::kable(digits=4) , file="resources/table2p.txt", sep="\n")
cat(table2c |>
    mutate(term=gsub("(,[^)]+)", "", term)) |>
    pivot_longer(-term) |>
    pivot_wider(names_from=term) |>
    mutate(HEART=`log(HEARTRTE)`, GLUC=`log(GLUCOSE)`,
           `GLUC:HEART`=`log(GLUCOSE):log(HEARTRTE)`,
           name=name,
           .keep="none") |>
    knitr::kable(digits=4),
    file="resources/table2c.txt", sep="\n")



## START model3
knots_log_gluc <- as.vector(quantile(log(framingham4$GLUCOSE), probs=c(0.33,0.67), na.rm=1))
bnds_log_gluc  <- range(log(framingham4$GLUCOSE), na.rm=1)
form3  <- CVDF ~ ns(log(GLUCOSE), knots=knots_log_gluc, Boundary.knots=bnds_log_gluc) * log(HEARTRTE) +
    TOTCHOL * SYSBP + SEX + AGE + log(BMI) + SYSBP + DIABETES + cig + educ + BPMEDS + PREVAP
model3 <- fit.mult.impute(form3, xtrans=imp, data=framingham4,
    fitter=function(formula, data){glm(formula, data, family=binomial(link=logit))})
p19 <- gridExtra::grid.arrange(plot.resid(model3) + aes(.fitted), plot.resid(model3) + aes(log(SYSBP)),
         plot.resid(model3)+aes(log(framingham4$GLUCOSE))+xlab("log(GLUCOSE)"), plot.resid.var(model3), nrow=1)
## END model3
ggsave("resources/model3.pdf", p19, width=4000, height=2000, units="px")


## START tables3
table3p <- p.value.table(model3)
table3c <- conf.int.table(model3, threshold=1)
## END tables3
cat(table3c |> mutate(term=gsub("(,[^)]+)", "", term)) |> knitr::kable(digits=4), file="resources/table3c.txt", sep="\n")



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



## START cal
p20 <- gridExtra::grid.arrange(cal_plot + group_p(framingham4$CVD, fitted(model1)),
  cal_plot + group_p(framingham4$CVD, fitted(model2)), cal_plot + group_p(framingham4$CVD, fitted(model3)), nrow=1)
## END cal
ggsave("resources/cal.pdf", p20, width=5000, height=2000, units="px")

## START pred-error
## This function is used to make sure we don't take the logarithm of 0.
avoid0 <- function(n, epsilon=0.000000000000001)
    {ifelse(n == 0, n + epsilon, ifelse(n == 1, n - epsilon, n))}

## Squared deviance residuals
## Uses GA3
dev.res.sq <- function(data, muhat) {
    Y     <- as.numeric(data$CVDF)-1
    term1 <- avoid0(Y)   * log(avoid0(   Y  / (muhat)))
    term2 <- avoid0(1-Y) * log(avoid0((1-Y) / (1*(1-muhat))))
    2*(term1 + term2)
}

in_sample_error <- function(model) {
    X           <- model.matrix(model)
    mu_hat      <- fitted(model)
    W           <- diag(mu_hat*(1-mu_hat))
    eff_df      <- sum(diag(X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W))
    train_error <- sum(dev.res.sq(model$data, fitted(model)))/nrow(model$data)
    train.error + 2*eff_df/nrow(model$data)}

auc <- function(model) {
    eta0 <- predict(model, model$data |> filter(CVD == 0), type="response")
    eta1 <- predict(model, model$data |> filter(CVD == 1), type="response")
    mean(outer(eta0, eta1, `<`) + outer(eta0, eta1, `==`)/2)
}

cv <- function(form, k, m) {
    losses <- list()
    for (b in seq(m)) {
        muhat <- list()
        groups <- sample(rep(1:k, length.out=nrow(framingham4)))
        for (i in seq(k)) {
            val_data   <- framingham4[groups == i,]
            train_data <- framingham4[groups != i,]
            train_imp <- aregImpute(impute_form, train_data, n.impute=14, nk=4)
            model <- fit.mult.impute(
                form, function(formula, data){glm(formula, data, family=binomial(link=logit))},
                train_imp, data=train_data)
            imputed_val_data <- impute.transcan(
                train_imp, data=train_data, newdata=val_data,
                imputation=(((b-1)*k+i)%%5)+1, list.out=TRUE, pr=FALSE)
            muhat <- predict(model, newdata=imputed_val_data, type="response")
            losses[[(b-1)*k+i]] <- dev.res.sq(imputed_val_data, muhat)
        }
    }
    losses |> unlist() |> mean()
}

errors_df <- data.frame(CV=c(cv(form1, 10, 10), cv(form2, 10, 10), cv(form3, 10, 10)),
                        ISE=c(in_sample_error(model1), in_sample_error(model2), in_sample_error(model3)),
                        AUC=c(auc(model1), auc(model2), auc(model3)))
rownames(errors_df) <- c("model1", "model2", "model3")
## END pred-error
pred_table <- errors_df |> knitr::kable(digits=3)
cat(pred_table, file="resources/pred-error.txt", sep="\n")


## cv1 <- cv(form1, framingham4, 3, 3, family=binomial(link=logit), loss=dev.res.sq)
## cv2 <- cv(form2, framingham4, 3, 3, family=binomial(link=logit), loss=dev.res.sq)
## cv3 <- cv(form3, framingham4, 3, 3, family=binomial(link=logit), loss=dev.res.sq)

## START validation
validation_data2 <- fix_framingham(validation_data)
validation_data3 <- impute.transcan(
    imp, data=framingham4, newdata=validation_data2,
    imputation=1, list.out=TRUE, pr=FALSE)
gen_error1 <- mean(dev.res.sq(validation_data3, predict(model1, newdata=validation_data3, type="response")))
gen_error2 <- mean(dev.res.sq(validation_data3, predict(model2, newdata=validation_data3, type="response")))
gen_error3 <- mean(dev.res.sq(validation_data3, predict(model3, newdata=validation_data3, type="response")))
## END validation

errors_df2 <- data.frame(`Generalization error estimated
on validation data`=c(gen_error1, gen_error2, gen_error3))
rownames(errors_df2) <- c("model1", "model2", "model3")
colnames(errors_df2) <- c("Validation error")
## pred_table2 <- mutate(errors_df, Val=c(gen_error1, gen_error2, gen_error3))
cat(errors_df2 |> knitr::kable(digits=3), file="resources/pred-error2.txt", sep="\n")
