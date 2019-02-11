#install.packages( "http://www.well.ox.ac.uk/~gav/resources/rbgen_v1.1.4.tgz", repos = NULL, type = "source" )
lapply(c("dplyr", "ggplot2", "rbgen"), require, character.only = T)
