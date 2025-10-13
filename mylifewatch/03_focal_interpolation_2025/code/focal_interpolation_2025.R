library(terra)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

library(raster)
library(ncdf4)
library(RNetCDF)

inputs_path = "/mnt/inputs/"
outputs_path <- "/mnt/outputs/"

input_dir <- paste0(inputs_path, "biooracle/")
output_dir <- paste0(outputs_path, "biooracle/")
plots_path <- paste0(outputs_path, "plots/")
if (!dir.exists(plots_path)) dir.create(plots_path, recursive = TRUE)

dataset_scenarios <- list.dirs(
  paste(input_dir, "datalayer.nc", sep = ""),
  full.names = FALSE,
  recursive = FALSE
)

for (scenario in dataset_scenarios) {
  rasterpath <- paste(input_dir, "datalayer.nc/", scenario, sep = "")
  stackpath <- paste(output_dir, "rasterstacks/", scenario, "/", sep = "")
  if (!dir.exists(stackpath)) dir.create(stackpath, recursive = TRUE)
  print(paste("Crated stackpath:", stackpath))

  lista.ras <- Sys.glob(paste(rasterpath, "/*.nc", sep = ""))
  layers <- rast(lista.ras)
  template <- layers[[1]]

  # Load and reproject the land polygons to the template CRS:
  land <- ne_countries(scale = "medium", returnclass = "sf")
  land <- st_transform(land, crs(template))

  # Convert to a SpatVector:
  land_vect <- vect(land)
  plot(land_vect)
  # Save plot to file
  png(filename = paste(plots_path, "land_polygons.png", sep = ""), width = 800, height = 600)
  dev.off()

  # Rasterize the land polygons onto the template.
  land_mask <- rasterize(land_vect, template, field = 1, background = NA)
  land_mask[!is.na(land_mask)] <- 1

  # Check the unique values:
  print(unique(values(land_mask)))

  masked_list <- lapply(1:nlyr(layers), function(i) {
    this_layer <- layers[[i]]
    mask(this_layer, land_mask, maskvalue = 1)
  })
  masked_layers <- rast(masked_list)

  filled_list <- lapply(1:nlyr(masked_layers), function(i) {
    this_layer <- masked_layers[[i]]
    focal(this_layer, w = 3, fun = function(x) {
      if (is.na(x[5])) {
        mean(x, na.rm = TRUE)
      } else {
        x[5]
      }
    })
  })
  filled_layers <- rast(filled_list)

  print("write raster")
  summary(values(template))
  summary(values(filled_layers[[1]]))
  writeRaster(filled_layers, paste(stackpath, "filled_layers_new.tif", sep = ""), overwrite = TRUE, filetype = "GTiff")

  mystack <- stack(paste(stackpath, "filled_layers_new.tif", sep = ""))
  rm(masked_layers)
  rm(land_mask)
  rm(this_layer)
  gc()
  rm(filled_layers)
  rm(land_mask)
  gc()

  e <- extent(-25, 45, 30, 72) #xmin, xmax, ymin, ymax
  rasterstack.filled.layers.Europe.2025 <- crop(mystack, e)

  filename <- paste(stackpath, "/", "Biooracle.filled.layers.global2025.tif", sep = "")
  print("save filled layers)")
  writeRaster(mystack, filename, format = "GTiff", overwrite = TRUE)
  print(filename)

  filename2 <- paste(stackpath, "/", "Biooracle.filled.layers.Europe2025.tif", sep = "")
  writeRaster(rasterstack.filled.layers.Europe.2025, filename2, format = "GTiff", overwrite = TRUE)

  print("save layernames")
  layernames <- names(mystack)
  save(layernames, file = paste(stackpath, "/layernames.filled.rda", sep = ""))
  rm(mystack)
  rm(rasterstack.filled.layers.Europe.2025)
  gc()
}