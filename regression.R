library(ggplot2)
library(tibble)
library(broom)
library(splines)
library(readr)
library(skimr)
library(dplyr)
library(rlang)

set.seed(8088)
framingham_raw <- as_tibble(read.csv("training_data.csv"))


# Framingham is ordered so that lower RANDIDs come first, and when there are
# multiple rows for the same RANDID, the earlier checkups come first.
# This is evidenced by the below check.
all(arrange(framingham_raw[,c("RANDID","TIME")]) == framingham_raw[,c("RANDID","TIME")])

# Thankfully noone attended more than three examinations.
max(table(framingham_raw[,"RANDID"]))


# We use this knowledge to add the column days_till_censoring which shows
# how many years of observation are left until the censoring date. This is
# based on some loose estimate of the censoring date being after 24 years.
# We also drop X since it just counts rows and reencode a bunch of columns
# as factors. Not CVD though since it's our response and so it has to be
# numeric if we want to fit regular lms on it.
framingham_all <- mutate(
    framingham_raw,
    SEX      = factor(SEX,      ordered=TRUE),
    CURSMOKE = factor(CURSMOKE, ordered=TRUE),
    DIABETES = factor(DIABETES, ordered=TRUE),
    BPMEDS   = factor(BPMEDS,   ordered=TRUE),
    educ     = factor(educ,     ordered=TRUE),
    PREVCHD  = factor(PREVCHD,  ordered=TRUE),
    PREVAP   = factor(PREVAP,   ordered=TRUE),
    PREVMI   = factor(PREVMI,   ordered=TRUE),
    PREVSTRK = factor(PREVSTRK, ordered=TRUE),
    PREVHYP  = factor(PREVHYP,  ordered=TRUE),
    PERIOD   = factor(PERIOD,   ordered=TRUE),
    DEATH    = factor(DEATH,    ordered=TRUE),
    ANGINA   = factor(ANGINA,   ordered=TRUE),
    HOSPMI   = factor(HOSPMI,   ordered=TRUE),
    MI_FCHD  = factor(MI_FCHD,  ordered=TRUE),
    ANYCHD   = factor(ANYCHD,   ordered=TRUE),
    STROKE   = factor(STROKE,   ordered=TRUE),
    HYPERTEN = factor(HYPERTEN, ordered=TRUE),
    # We can see that the censoring date was after 8766 days because that's
    # the maximum value in all the TIME* columns, and I think people who didn't
    # have a specific event happen were registered as having it occur after the
    # whole examination period.
    days_till_censoring=8766-TIME,
    X=NULL,
)

framingham <- select(
    framingham_all,
    CVD,
    AGE,
    SEX,
    educ,
    days_till_censoring,
    TOTCHOL,
    SYSBP,
    DIABP,
    BPMEDS,
    CURSMOKE,
    CIGPDAY,
    BMI,
    HEARTRTE,
    GLUCOSE,
    DIABETES,
    PREVCHD,
    PREVAP,
    PREVMI,
    PREVSTRK,
    PREVHYP
)

# The following expression tells us that over half of the observations
# are in the second and third checkups, so it would be a shame to have to
# throw all that away.
table(framingham_all$PERIOD)

fill = gray(0.7)

# Plot the marginal distribution of a column
marginal_plot <- function(col) {
    if (is.factor(framingham[[col]]) | col == "CVD") {
        plot <- ggplot(data=framingham) +
            aes_string(col) +
            xlab(col) +
            ylab("") +
            geom_bar(fill=fill)
    } else {
        if (is.null(timecols[[col]])) {
            data <- framingham
        } else {
            data <- timefilter(col)
        }
        plot <- ggplot(data=data) +
            aes_string(col) +
            xlab(col) +
            ylab("") +
            geom_density(fill=fill)
    }
    plot
}

tmp <- lapply(names(framingham), marginal_plot)
do.call(gridExtra::grid.arrange, c(tmp, ncol=5))


# Taken from the book chapter 3
ggally_hexbin <- function(data, mapping, ...) {
    ggplot(data, mapping) + stat_binhex(...)
}


GGally::ggpairs(
            framingham,
            columns = c("AGE", "days_till_censoring", "TOTCHOL", "SYSBP", "DIABP", "CIGPDAY", "BMI", "HEARTRTE", "GLUCOSE"),
            lower   = list(continuous = GGally::wrap("hexbin", bins=20)),
            diag = list(continuous = "blankDiag")
        )

skim(framingham)

