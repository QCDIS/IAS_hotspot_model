library(terra)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(raster)
library(ncdf4)
library(RNetCDF)

biooracle_dir <- "/mnt/inputs/biooracle/"
rasterstacks_dir <- "/mnt/outputs/rasterstacks"
if (!dir.exists(rasterstacks_dir)) dir.create(rasterstacks_dir, recursive = TRUE)
dataset_scenarios <- c("baseline", "ssp119", "ssp126", "ssp245", "ssp370", "ssp460", "ssp585")
dec.vec <- c("", "dec50", "dec100")

# --- Outer scenario loop ---
for (sel.sen in c(3:7)) {
  scenario <- dataset_scenarios[[sel.sen]]
  # --- Decadal loop ---
  for (dec in dec.vec[c(2, 3)]) {
    print(paste(scenario, dec))

    stackpath <- paste(rasterstacks_dir, scenario, dec, "/", sep = "")
    rasterpath <- paste(biooracle_dir, "datalayer.nc/", scenario, dec, sep = "")


    print(paste("stackpath: ", stackpath))
    print(paste("rasterpath: ", rasterpath))

    lista.ras <- Sys.glob(paste(rasterpath, "/*.nc", sep = ""))
    layers <- rast(lista.ras)
#
#     template <- layers[[1]]
#
#     land <- ne_countries(scale = "medium", returnclass = "sf")
#     land <- st_transform(land, crs(template))
#     land_vect <- vect(land)
#     plot(land_vect)
#
#     land_mask <- rasterize(land_vect, template, field = 1, background = NA)
#     land_mask[!is.na(land_mask)] <- 1
#     print(unique(values(land_mask)))
#
#     masked_list <- lapply(1:nlyr(layers), function(i) {
#       this_layer <- layers[[i]]
#       mask(this_layer, land_mask, maskvalue = 1)
#     })
#     masked_layers <- rast(masked_list)
#
#     filled_list <- lapply(1:nlyr(masked_layers), function(i) {
#       this_layer <- masked_layers[[i]]
#       focal(this_layer, w = 3, fun = function(x) {
#         if (is.na(x[5])) {
#           mean(x, na.rm = TRUE)
#         } else {
#           x[5]
#         }
#       })
#     })
#     filled_layers <- rast(filled_list)
#
#     print("write raster")
#     summary(values(template))
#     summary(values(filled_layers[[1]]))
#     writeRaster(filled_layers, paste(stackpath, "filled_layers_new.tif", sep = ""), overwrite = TRUE, filetype = "GTiff")
#
#     mystack <- stack(paste(stackpath, "filled_layers_new.tif", sep = ""))
#     rm(masked_layers)
#     rm(land_mask)
#     rm(this_layer)
#     gc()
#     rm(filled_layers)
#     rm(land_mask)
#     gc()
#
#     e <- extent(-25, 45, 30, 72)
#     rasterstack.filled.layers.Europe.2025 <- crop(mystack, e)
#
#     filename <- paste(stackpath, "/", "Biooracle.filled.layers.global2025.tif", sep = "")
#     print("save filled layers)")
#     writeRaster(mystack, filename, format = "GTiff", overwrite = TRUE)
#     print(filename)
#
#     filename2 <- paste(stackpath, "/", "Biooracle.filled.layers.Europe2025.tif", sep = "")
#     writeRaster(rasterstack.filled.layers.Europe.2025, filename2, format = "GTiff", overwrite = TRUE)
#
#     print("save layernames")
#     layernames <- names(mystack)
#     save(layernames, file = paste(stackpath, "/layernames.filled.rda", sep = ""))
#     rm(mystack)
#     rm(rasterstack.filled.layers.Europe.2025)
#     gc()
  } # --- End decadal loop ---
} # --- End outer scenario loop ---

# get full layernames
# load(paste(stackpath, "/layernames.rda", sep = ""))
# names(mystack) <- layernames