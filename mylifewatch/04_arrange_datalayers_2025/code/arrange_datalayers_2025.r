# Required libraries
require(raster)
library(sf)
require(fBasics)

# Define directories
# --- Define file paths and output folders ---
inputs_path = "/mnt/inputs/"
outputs_path <- "/mnt/outputs/"

biooracle_dir <- paste0(inputs_path, "biooracle","/rasterstacks/")
rasterstacks_outputs <- outputs_path

# Get scenario folders
dataset_scenarios <- list.dirs(
  biooracle_dir,
  full.names = FALSE,
  recursive = FALSE
)



# --- Scenario loop ---
for (scenario in dataset_scenarios) {
    print(paste("Processing scenario:", scenario))
    #  outdir <- paste(dir,datalayer_dir,scenario,"/",sep="")

    rasterpath <- paste(biooracle_dir,scenario,"/",sep="")
    stackpath <- paste(rasterstacks_outputs,"/rasterstacks/",scenario,"/",sep="")
    if (!dir.exists(stackpath)) dir.create(stackpath, recursive = TRUE)

    print(paste("Raster path:", rasterpath))
    lista.ras<- Sys.glob(paste(rasterpath,"/*",".tif",sep=""))
    print(paste("Number of raster files found:", length(lista.ras)))
    mystack <- stack(lista.ras)
    for(i in lista.ras){
        print(i)
        gc()
    }
    filename <- paste(stackpath,"/","Biooracle.global.tif", sep="")
    writeRaster(mystack, filename, format="GTiff",overwrite=TRUE)

    e <- extent(-25, 45, 30, 72) #xmin, xmax,ymin,ymax
    rasterstack.Europe.2025 <- crop(mystack, e)
    plot(rasterstack.Europe.2025)
    png(paste(stackpath,"/","Biooracle.Europe2025.png", sep=""))

    filename2 <- paste(stackpath,"/","Biooracle.Europe2025.tif", sep="")
    writeRaster(rasterstack.Europe.2025, filename2, format="GTiff",overwrite=TRUE)


    layernames <- names(mystack)
    save(layernames, file = paste(stackpath,"/layernames.rda" ,sep =""))

}