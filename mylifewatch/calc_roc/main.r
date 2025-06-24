# Load the functions/SEanalytics.functions.r

source("functions/SEanalytics.functions.r")
lapply(list.files("/wrp/functions", full.names = TRUE, pattern = "\\.R$"), source)
args = args_parse(commandArgs(trailingOnly = TRUE))

# Call the read.and.extract function

Speciespath <- "/mnt/inputs/species_to_use"
Outpath <- "/mnt/outputs/RF.indata"

Stack <- file.path(Stackpath, "globalStack.rda") # "globalStack.rda" or "europeStack.rda"
Iterations.path <- "/mnt/outputs/iterations" # Path to the iterations zip file

RF_results <- "/mnt/outputs/RF.results.cv"

all.ROC <- lapply(all.rep, function(curr.rep)
 calc.ROC(rf.output.cv[[curr.rep]],.true.class)
)