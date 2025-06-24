# Load the functions/SEanalytics.functions.r

source("functions/SEanalytics.functions.r")

# Call the read.and.extract function

Speciespath <- "/mnt/inputs/species_to_use"
Outpath <- "/mnt/outputs/RF.indata"
data_table = "/mnt/inputs/data_table.csv"


Stack <- file.path(Stackpath, "globalStack.rda") # "globalStack.rda" or "europeStack.rda"
Iterations.path <- "/mnt/outputs/iterations" # Path to the iterations zip file

RF_results <- "/mnt/outputs/RF.results.cv"


Data.table <- read.csv2(data_table,header=TRUE, sep=";",stringsAsFactors = F)


for(Species in Data.table$species){
    rf.output.list <- run.random.forests(species = Species,
                selvar = "all",
                indata.path = Outpath,
                iterations.path= Iterations.path)

    rf_path = paste(RF_results,"RF.results.CV",sep="/")
    rf.output.cv <- rf.output.list[[rf_path]]
}