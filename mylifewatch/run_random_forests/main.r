# Load the functions/SEanalytics.functions.r

source("functions/SEanalytics.functions.r")
lapply(list.files("/wrp/functions", full.names = TRUE, pattern = "\\.R$"), source)
args = args_parse(commandArgs(trailingOnly = TRUE))

# Call the read.and.extract function

Speciespath <- "/mnt/inputs/species_to_use"
Outpath <- "/mnt/outputs/RF.indata"

Species <- read.species(Speciespath) # Read the species list from the input path
Stack <- file.path(Stackpath, "globalStack.rda") # "globalStack.rda" or "europeStack.rda"
Iterations.path <- "/mnt/outputs/iterations" # Path to the iterations zip file

RF_results <- "/mnt/outputs/RF.results.cv"

rf.output.list <- run.random.forests(species = Species,
            selvar = "all",
            indata.path = Outpath,
            iterations.path= Iterations.path)

rf_path = paste(RF_results,"RF.results.CV",sep="/")
rf.output.cv <- rf.output.list[[rf_path]]