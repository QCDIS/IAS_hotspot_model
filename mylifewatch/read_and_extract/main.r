# Load the functions/SEanalytics.functions.r

source("functions/SEanalytics.functions.r")

# Call the read.and.extract function



data_table <- "/mnt/inputs/data.table.csv"
Speciespath <- "/mnt/inputs/species_to_use"
Stackpath <- "/mnt/inputs/rasterstacks"
Plotpath_input <- "/mnt/inputs/plots/"
Outpath <- "/mnt/outputs//RF.indata"

Stack <- file.path(Stackpath, "globalStack.rda") # "globalStack.rda" or "europeStack.rda"

 species.data <-  read.and.extract(data.table = Data.table,
                                   species = Species,
                                   stack = Stack,
                                   speciespath = Speciespath,
                                   stackpath = Stackpath,
                                   plotpath = Plotpath_input,
                                   outpath = Outpath)