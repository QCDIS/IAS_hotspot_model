# Load the functions/SEanalytics.functions.r

source("functions/SEanalytics.functions.r")
lapply(list.files("/wrp/functions", full.names = TRUE, pattern = "\\.R$"), source)

# Call the read.and.extract function

Speciespath <- "/mnt/inputs/species_to_use"
Outpath <- "/mnt/outputs/RF.indata"
Modelpath <- "/mnt/outputs/models"
Plotpath_output <- "/mnt/outputs/plots/"
mappath <- "/mnt/outputs/maps"

brk <- c(seq(0, 1,by=0.1),1.05)
for(Species in Data.table$species){
 map <-  plot.maps(species = Species,
            indata.path = Outpath,
            modelpath = Modelpath,
            plotpath = Plotpath_output,
            colors = c(colorsBr, "lightgrey"),
            brk = brk
            )
 save(map, file= paste(mappath,"/predicted.map.",Species,".rda",sep=""))
}