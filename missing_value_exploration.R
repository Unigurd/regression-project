library(tibble)
library(broom)
library(readr)
library(skimr)
library(dplyr)
library(tidyr)
library(hexbin)
library(Hmisc)

set.seed(8088)
framingham_raw <- read.csv("training_data.csv") 
framingham <- framingham_raw |>
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
    ) 
## |>
    ## drop_na()

## ~ SEX + CIGPDAY + DIABETES + AGE + CVD + BMI + BPMEDS + GLUCOSE + educ

fram_imp <- framingham |>
    mutate_at(c("CVD", "SEX", "DIABETES", "BPMEDS", "educ",
                "PREVCHD", "PREVAP", "PREVMI", "PREVHYP"),
              as.ordered)
impute_form <- ~ CVD + SEX + TOTCHOL + AGE + SYSBP + CIGPDAY + DIABETES + HEARTRTE + PREVCHD + PREVAP + PREVMI + PREVHYP + BMI + BPMEDS + GLUCOSE + educ 
imp <- aregImpute(
    impute_form,
    fram_imp,
    n.impute=5,
    nk=4
)

n_distinct(framingham$PREVCHD)

Vectorize(function(col){n_distinct(framingham[[col]])})(c("CVD", "SEX", "TOTCHOL", "AGE", "SYSBP", "CIGPDAY", "BMI", "DIABETES", "BPMEDS", "HEARTRTE", "GLUCOSE", "educ", "PREVCHD", "PREVAP", "PREVMI", "PREVHYP"))

fill = gray(0.7)

# Plot the marginal distribution of a column
marginal_plot <- function(data, col) {
    if (is.factor(data[[col]]) | col == "CVD" | col == "CIGPDAY") {
        plot <- ggplot(data=data) +
            aes_string(col) +
            xlab(col) +
            ylab("") +
            geom_bar(fill=fill)
    } else {
        plot <- ggplot(data=data) +
            aes_string(col) +
            xlab(col) +
            ylab("") +
            geom_density(fill=fill)
    }
    plot
}

marginal_plots <- function(data, cols=NULL) {
    if (is.null(cols)) {
        cols <- names(data)
    }
    vec.marg.plot <- Vectorize(marginal_plot, vectorize.args="col")
    tmpplots <- vec.marg.plot(data, cols)
    do.call(gridExtra::grid.arrange, tmpplots)}

# SEX 2 has more missing GLUCOSEs
gridExtra::grid.arrange(
    ggplot(data=framingham_raw |> filter(!is.na(GLUCOSE)), aes(SEX)) +
    geom_bar(fill=fill),
    ggplot(data=framingham_raw |> filter(is.na(GLUCOSE)), aes(SEX)) +
    geom_bar(fill=fill),
    ncol=2
)

# For SEX 1 with missing GLUCOSEs, they're slightly more frequently
# smokers, die slightly less, get slightly less ANYCHD and get
# slightly less CVD. But maybe there's just too little data there and
# it's random fluctuations?

# For SEX 2 with missing GLUCOSE, they're slightly more frequently
# smokers. The rest of the changes for SEX 1 are not really present.

fram <- framingham_raw |>
    filter(
        ## PERIOD == 1,
        ## SEX == 2
           ) |>
    mutate_at(
        c("SEX", "CURSMOKE", "DIABETES", "BPMEDS", "educ",
          "PREVCHD", "PREVAP", "PREVMI", "PREVSTRK",
          "PREVHYP", "PERIOD", "DEATH", "ANGINA", "HOSPMI",
          "MI_FCHD", "ANYCHD", "STROKE", "HYPERTEN",
          "CVD"),
        as.ordered
    ) |>
    mutate(X=NULL, RANDID=NULL, HDLC=NULL, LDLC=NULL, TIME=NULL, PERIOD=NULL)

cp <- cor(
    data.matrix(fram),
    use    = "complete.obs",
    method = "spearman"
)
corrplot::corrplot(
              cp,
              diag    = FALSE,
              order   = "hclust",
              addrect = 13,
              tl.srt  = 45,
              tl.col  = "black",
              tl.cex  = 0.8
          )

# Diabetes and age are most correlated with GLUCOSE.
gluc_cps <- tibble(vars=names(cp[,"GLUCOSE"]), values=as.numeric(cp[,"GLUCOSE"])) |>
    arrange(values)
print(gluc_cps, n=40)

# SEX, CURSMOKE, DIABETES, AGE, CVD, BMI


gp(marginal_plots(fram |> filter(!is.na(GLUCOSE))))
gp(marginal_plots(fram |> filter(is.na(GLUCOSE))))

tmp <- lapply(names(framingham), marginal_plot)
do.call(gridExtra::grid.arrange, c(tmp, ncol=5))


ggplot(framingham, aes(BMI, CIGPDAY)) +
    geom_point(data=subset(framingham, !is.na(GLUCOSE)), color="black") + 
    geom_point(data=subset(framingham, is.na(GLUCOSE)), color="red")


table(rowSums(is.na(framingham |> mutate(GLUCOSE=NULL, educ=NULL))))

table(rowSums(is.na(framingham )))
framingham |> filter(rowSums(is.na(framingham)) == 2)

framingham |> filter(rowSums(is.na(framingham)) == 2) |> select(CIGPDAY, BMI, BPMEDS, GLUCOSE, educ) |> naclus()

naplot(framingham)

# Taken from the book chapter 3
ggally_hexbin <- function(data, mapping, ...) {
    ggplot(data, mapping) + stat_binhex(...) 
    ## +
    ## geom_point(data=subset(data, is.na(GLUCOSE)), color="red", alpha=0.2)
}

# We see that SYSBP and DIABP are the most correlated with a Corr score of 0.722.
# We see a single outlier for HEARTRTE, four for TOTCHOL, maybe some for SYSBP?
# So we should maybe pick DIABP over it?
# CIGPDAY and educ are very skewed, as are most of the binary factors.
# Should we make groupings for CIGPDAY?
# Should we make groupings for time_to_censor? four of them?
ggally <- GGally::ggpairs(
            framingham |> filter(is.na(GLUCOSE)),
            columns = c("AGE", "TOTCHOL", "SYSBP", "CIGPDAY", "BMI", "HEARTRTE"## , "GLUCOSE"
                        )
           ,
            lower   = list(continuous = GGally::wrap("hexbin", bins=20)),
            diag = list(continuous = "blankDiag")
        )
gp(ggally)


ggplot()framingham |> filter(is.na(GLUCOSE))


cp <- cor(
    data.matrix(fram),
    use    = "complete.obs",
    method = "spearman"
)
corrplot::corrplot(
              cp,
              diag    = FALSE,
              order   = "hclust",
              addrect = 13,
              tl.srt  = 45,
              tl.col  = "black",
              tl.cex  = 0.8
          )

cpplot

gp(cpplot)

gp <- function(plot) {
    x11()
    print(plot)
}

cor(data.matrix(mutate(framingham, x=log(GLUCOSE), y=SYSBP, .keep="none")),
    use="complete.obs",
    method="spearman")


gp(ggplot(framingham, aes(log(GLUCOSE), SYSBP)) +
   geom_point())

skim(framingham)

# cigpday: Use the median 0? That goes against the
# suspicion that they're missing for high values

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
