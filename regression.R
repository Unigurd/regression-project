library(ggplot2)
library(tibble)
library(broom)
library(splines)
library(readr)
library(skimr)
library(dplyr)

framingham <- as_tibble(read.csv("training_data.csv"))


