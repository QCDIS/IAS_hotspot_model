# Load the functions/SEanalytics.functions.r

source("functions/SEanalytics.functions.r")

# Call the read.and.extract function

Speciespath <- "/mnt/inputs/species_to_use"
Outpath <- "/mnt/outputs/RF.indata"
data_table = "/mnt/inputs/data_table.csv"
RF_results <- "/mnt/outputs/RF.results.cv"

Stack <- file.path(Stackpath, "globalStack.rda") # "globalStack.rda" or "europeStack.rda"
Iterations.path <- "/mnt/outputs/iterations" # Path to the iterations zip file

Data.table <- read.csv2(data_table,header=TRUE, sep=";",stringsAsFactors = F)


for(Species in Data.table$species){
    print(Species)
    species.data <- split.data(species = Species, indata.path = Outpath, iterations.path= Iterations.path)
    write.csv2(species.data, file = paste(Outpath,"/",Species,"indata.csv", sep=""),row.names=F)
}