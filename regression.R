library(ggplot2)
library(tibble)
library(broom)
library(splines)
library(readr)
library(skimr)
library(dplyr)

framingham <- as_tibble(read.csv("training_data.csv"))


# Framingham is ordered so that lower RANDIDs come first, and when there are
# multiple rows for the same RANDID, the lower ages come first.
# This is evidenced by the below check.
all(arrange(framingham[,c("RANDID","AGE")]) == framingham[,c("RANDID","AGE")])

# We use this knowledge to add the following columns:
# - age_diff shows how much each person has aged since their first checkup.
# - checkup shows which checkup number the row is.
# - years_till_censoring shows how many years of observation are left until the censoring date.
# We also drop X since it just counts rows.
fram_checkup <- mutate(
    framingham,
    age_diff=AGE-AGE[1],
    # Choices for cutoffs are argued below.
    checkup=ifelse(age_diff == 0, 1, ifelse(age_diff < 10, 2, 3)),
    years_till_censoring=24-age_diff,
    X=NULL,
    .by=RANDID,
    .after=RANDID
)

# The following expression shows how many examinations were done roughly how
# many years after each person's first checkup. Note that this is not the
# same as calendar years since the first observation.
# The checkups are approimately 6 years apart, so it seems reasonable that
# the one observation after 8 years would belong to the second checkup round,
# while the three checkups after ten years would belong to the third. This way, each
# checkup round runs over a consecutive number of years (as measured by age_diff).
table(fram_checkup$age_diff)

# Thankfully noone attended more than three examinations.
max(fram_checkup$checkup)

# The following expression tells use that over half of the observations
# are in the second and third checkups, so it would be a shame to have to
# throw all that away.
table(fram_checkup$checkup)

