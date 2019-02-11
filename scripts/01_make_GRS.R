#install.packages( "http://www.well.ox.ac.uk/~gav/resources/rbgen_v1.1.4.tgz", repos = NULL, type = "source" )
setwd("/home/jm2294/repos/never_lupus")
lapply(c("dplyr", "ggplot2", "rbgen"), require, character.only = T)

ranges = data.frame(
  chromosome = "01",
  start = 100000,
  end = 120000
)
data <- bgen.load("SOMAS_round_all_SOMAS_round_all_impute_1_interval_filtered2_reworked.bgen", ranges)