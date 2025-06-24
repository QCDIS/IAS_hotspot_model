# Load the functions/SEanalytics.functions.r

source("functions/SEanalytics.functions.r")

# Call the read.and.extract function

Speciespath <- "/mnt/inputs/species_to_use"
Outpath <- "/mnt/outputs/RF.indata"

Species <- read.species(Speciespath) # Read the species list from the input path
Stack <- file.path(Stackpath, "globalStack.rda") # "globalStack.rda" or "europeStack.rda"
Iterations.path <- "/mnt/outputs/iterations" # Path to the iterations zip file



  split.data(species = Species,
             indata.path = Outpath,
             iterations.path= Iterations.path)