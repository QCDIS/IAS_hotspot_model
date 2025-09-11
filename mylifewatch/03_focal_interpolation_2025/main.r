
lapply(list.files("functions", full.names = TRUE),source)
lapply(list.files("/wrp/functions", full.names = TRUE, pattern = "\\.R$"), source)

args = args_parse(commandArgs(trailingOnly = TRUE))



source("code/focal_interpolation_2025.R")