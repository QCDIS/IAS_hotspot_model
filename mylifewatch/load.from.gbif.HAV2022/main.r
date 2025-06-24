
lapply(list.files("functions", full.names = TRUE),source)
lapply(list.files("/wrp/functions", full.names = TRUE, pattern = "\\.R$"), source)

args = args_parse(commandArgs(trailingOnly = TRUE))



source("code/load.from.gbif.HAV2022.r")