library(ggplot2)
library(tibble)
library(broom)
library(splines)
library(readr)
library(skimr)
library(dplyr)
library(rlang)
library(hexbin)

set.seed(8088)
framingham_raw <- read.csv("training_data.csv") |>
    as_tibble() |>
    mutate_at(
        c("SEX", "CURSMOKE", "DIABETES", "BPMEDS", "educ", "PREVCHD", "PREVAP", "PREVMI", "PREVSTRK",
          "PREVHYP", "PERIOD", "DEATH", "ANGINA", "HOSPMI", "MI_FCHD", "ANYCHD", "STROKE", "HYPERTEN"),
        as.ordered
    ) |>
    mutate(X=NULL)


# Framingham is ordered so that lower RANDIDs come first, and when there are
# multiple rows for the same RANDID, the earlier checkups come first.
# This is evidenced by the below check.
all(arrange(framingham_raw[,c("RANDID","TIME")]) == framingham_raw[,c("RANDID","TIME")])

# Thankfully noone attended more than three examinations.
max(table(framingham_raw[,"RANDID"]))

# This shows that if a CVD event happened to a person, then the CVD col is 1 for
# all rows corresponding to that person, so we have to do some work to figure out
# when the event happened relative to the various checkups. An observation after an event
# happened cannot be used to predict that event.
all((framingham_raw |> select(RANDID, CVD) |> mutate(eq=max(CVD) == min(CVD), .by=RANDID))$eq)

# ctc = check_time_consistency
ctc <- function(max, col) {
    (col == max | col == -1)
}

# I want to figure out for each checkup how long that person was observed for
# afterwards. The study ended after 8766 days, so it can be that at most. But for
# the events that didn't happen to a person, the corresponding time-to-event columns
# contain the time at which the study lost contact to that person, i.e. how long that
# person was observed for. I'd like to use that to determine a more accurate censoring
# date for each person.
# First we create a tibble to use for some sanity checks:
time_checking <- framingham_raw |>
    mutate(
        # We make copies of the time columns where
        # they're now set to the invalid time -1 when
        # the event happened.
        timeap   = ifelse(ANGINA   == 1, -1, TIMEAP),
        timemi   = ifelse(HOSPMI   == 1, -1, TIMEMI),
        timemifc = ifelse(MI_FCHD  == 1, -1, TIMEMIFC),
        timechd  = ifelse(ANYCHD   == 1, -1, TIMECHD),
        timestrk = ifelse(STROKE   == 1, -1, TIMESTRK),
        timecvd  = ifelse(CVD      == 1, -1, TIMECVD),
        timedth  = ifelse(DEATH    == 1, -1, TIMEDTH),
        timehyp  = ifelse(HYPERTEN == 1, -1, TIMEHYP),
        # We use these to find the maximum time recorded
        # in any of those columns in each row. The timedth column
        # is added separately because, as we will see, it doesn't
        # behave well.
        max = pmax(
            timeap,
            timemi,
            timemifc,
            timechd,
            timestrk,
            timecvd,
            timehyp),
        maxd = pmax(max, timedth),
        # We now check that all of those columns are either equal to the max in that row,
        # or equal to -1 (which means that the event happened and we want to disregard
        # that row in that column). If any column anywhere has a different value than
        # either of those, something is wrong with our assumption that time-to-event columns
        # don't record the time until the study lost contact with that person when their event
        # didn't happen.
        check = ctc(max, timeap) &
            ctc(max, timemi) &
            ctc(max, timemifc) &
            ctc(max, timechd) &
            ctc(max, timestrk) &
            ctc(max, timecvd) &
            ctc(max, timehyp),
        checkd = ctc(max, timeap) &
            ctc(max, timemi) &
            ctc(max, timemifc) &
            ctc(max, timechd) &
            ctc(max, timestrk) &
            ctc(max, timecvd) &
            ctc(max, timedth) &
            ctc(max, timehyp)
    ) |>
    select(RANDID, DEATH,
           timeap, timemi, timemifc, timechd, timestrk, timecvd, timedth, timehyp,
           max, maxd, check, checkd
    )

# We see that there are 68 rows where our assumption didn't hold.
time_checking |> filter(!checkd) |> nrow()

# But examining this printout we see that it's always because TIMEDTH was set to
# 8766, the end of the study. That's why we handled that separately above.
print(time_checking |> filter(!checkd) |> distinct(RANDID, .keep_all=TRUE), n=100)

# If we don't include TIMEDTH we see that 0 rows break our assumption.
time_checking |> filter(!check) |> nrow()

# We also check that for each person, the time until lost contact recorded in the
# time-to-event columns are the same for the different observations.
(time_checking |> mutate(id_check=n_distinct(max) == 1, .by=RANDID))$id_check |> all()

# Using the above approach we add columns containing the accurate time until censor
# as well as the time until a CVD event happens
framingham_all <- framingham_raw |>
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


# We can now remove all observations that occur after the CVD event happened and all observations
# where the person was never seen again. We also remove all columns that we don't want as predictors.
framingham <- framingham_all |>
    ## Should this be < or <=? There are rows with equality but all except two are on day 0.
    filter(time_to_cvd >= 0, time_to_censor > 0) |>
    subset(select=-c(
                       ## RANDID shouldn't be a predictor.
                       RANDID,
                       ## These have too many missing values because they were only
                       ## tested for at the last checkup.
                       HDLC, LDLC,
                       ## These have played their role and should not be used to predict
                       time_censor, time_to_cvd, TIME, PERIOD,
                       # We cannot use future events to predict
                       ANGINA, HOSPMI, MI_FCHD,  ANYCHD,  STROKE,            DEATH,   HYPERTEN,
                       TIMEAP, TIMEMI, TIMEMIFC, TIMECHD, TIMESTRK, TIMECVD, TIMEDTH, TIMEHYP)
           )


# The following expression tells us that over half of the observations
# are in the second and third checkups, so it would be a shame to have to
# throw all that away.
table(framingham_all$PERIOD)


# Issues to look out for:
# - Extreme observations and potential outliers
# - Missing values
# - Skewness or asymmetry of marginal distributions
# - collinearity of predictors

fill = gray(0.7)

# Plot the marginal distribution of a column
marginal_plot <- function(col) {
    if (is.factor(framingham[[col]]) | col == "CVD" | col == "CIGPDAY") {
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


# We see that SYSBP and DIABP are the most correlated with a Corr score of 0.722.
# We see a single outlier for HEARTRTE, four for TOTCHOL, maybe some for SYSBP?
# So we should maybe pick DIABP over it?
# CIGPDAY and educ are very skewed, as are most of the binary factors.
# Should we make groupings for CIGPDAY?
# Should we make groupings for time_to_censor? four of them?
GGally::ggpairs(
            framingham,
            columns = c("AGE", "time_to_censor", "TOTCHOL", "SYSBP",
                        "DIABP", "CIGPDAY", "BMI", "HEARTRTE", "GLUCOSE"),
            lower   = list(continuous = GGally::wrap("hexbin", bins=20)),
            diag = list(continuous = "blankDiag")
        )


cp <- cor(
    data.matrix(framingham),
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


# Plot the density of variables when cvd == 0 and when cvd == 1
plot_density_cvd <- function(col, data=framingham) {
    ggplot(data) +
        facet_grid(. ~ CVD, label=label_both) +
        aes_string(col) +
        geom_density(fill=fill)
}

numeric_names <- names(framingham)[!sapply(framingham, is.factor) &
                                   names(framingham) != "CVD"]

# Interestingly, the time_to_censor is more spread out
# for the cases where the person did get CVD. People
# who get CVD are more likely to go no-contact with the
# study? Do they die?
do.call(gridExtra::grid.arrange, c(lapply(numeric_names, plot_density_cvd), ncol=3))

