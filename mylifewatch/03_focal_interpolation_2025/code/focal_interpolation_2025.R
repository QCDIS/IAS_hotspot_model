library(terra)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
#library(rgdal)

library(raster)
library(ncdf4)
library(RNetCDF)

input_dir <- "/mnt/inputs/biooracle/"
output_dir <- "/mnt/outputs/biooracle/"
plots_path <- "/mnt/outputs/plots/"
if (!dir.exists(plots_path)) dir.create(plots_path, recursive = TRUE)

dataset_scenarios <- list.dirs(
  paste(input_dir, "datalayer.nc", sep = ""),
  full.names = FALSE,
  recursive = FALSE
)


#for(sel.sen in c(2:7)){
 # sel.sen <- 1
 # scenario <- dataset_scenarios[[sel.sen]]
  #  outdir <- paste(dir,"/datalayer.tiff/",scenario,"/",sep="")

for (scenario in dataset_scenarios) {

rasterpath <- paste(input_dir,"datalayer.nc/",scenario,sep="")
stackpath <- paste(output_dir,"rasterstacks/",scenario,"/",sep="")
if (!dir.exists(stackpath)) dir.create(stackpath, recursive = TRUE)
print(paste("Crated stackpath:", stackpath))
#paste(path, c("data/Biooracle_datalager/datalager.tiff"), sep="/")

lista.ras<- Sys.glob(paste(rasterpath,"/*.nc",sep=""))
layers <- rast(lista.ras)

# Use the first layer from your NetCDF stack as the template:
#layers()
#mystacknc <- stack(lista.ras[1])
#mystacknc <- ncdf4::nc_open(lista.ras[1])
#
template <- layers[[1]]

# Load and reproject the land polygons to the template CRS:
land <- ne_countries(scale = "medium", returnclass = "sf")
land <- st_transform(land, crs(template))

# Convert to a SpatVector:
land_vect <- vect(land)
plot(land_vect)
# Save plot to file
png(filename = paste(plots_path, "land_polygons.png", sep = ""), width = 800, height = 600)

#plot(template)
# Rasterize the land polygons onto the template.
# Use background = NA so that cells not covered by a polygon remain NA.
land_mask <- rasterize(land_vect, template, field = 1, background = NA)

# Now force any cell that got a value (i.e. where land is present) to 1:
land_mask[!is.na(land_mask)] <- 1

# Check the unique values:
print(unique(values(land_mask)))
# Expect something like: [1]  1  NA

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

# PLOTS TO CHECK RESULTS
print("plot land mask")
plot(land_mask, main = "Binary Land Mask (1 = Land, NA = Water)")
png(filename = paste(plots_path, "land_mask.png", sep = ""), width = 800, height = 600)
dev.off()

# # Overlay the original land polygons for reference
print("plot land vect")
plot(land_vect, add = TRUE, border = "red")
png(filename = paste(plots_path, "land_vect.png", sep = ""), width = 800, height = 600)
dev.off()

print("plot masked layer")
plot(masked_layers[[1]], main = "Masked Layer (Land = NA)")
png(filename = paste(plots_path, "masked_layer.png", sep = ""), width = 800, height = 600)
dev.off()

# Overlay the land polygons
print("plot land vect blue")
plot(land_vect, add = TRUE, border = "blue")
png(filename = paste(plots_path, "land_vect.png", sep = ""), width = 800, height = 600)
dev.off()

print("plot filled layer")
plot(filled_layers[[1]], main = "Filled Layer (Coastal Interpolated)")
png(filename = paste(plots_path, "filled_layer.png", sep = ""), width = 800, height = 600)
dev.off()

# Overlay the land polygons to see the boundary
print("plot land vect green")
plot(land_vect, add = TRUE, border = "green")
png(filename = paste(plots_path, "land_vect_green.png", sep = ""), width = 800, height = 600)
dev.off()

print("plot 3 panel")
par(mfrow = c(1, 3))
plot(layers[[1]], main = "Original Layer")
png(filename = paste(plots_path, "original_layer.png", sep = ""), width = 800, height = 600)
dev.off()

print("plot masked layer 3 panel")
plot(masked_layers[[1]], main = "Masked Layer")
png(filename = paste(plots_path, "masked_layer_3panel.png", sep = ""), width = 800, height = 600)
dev.off()

print("plot filled layer 3 panel")
plot(filled_layers[[1]], main = "Filled Layer")
png(filename = paste(plots_path, "filled_layer_3panel.png", sep = ""), width = 800, height = 600)
dev.off()

#par(mfrow = c(1, 1))
print("write raster")

summary(values(template))
summary(values(filled_layers[[1]]))
writeRaster(filled_layers, paste(stackpath ,"filled_layers_new.tif",sep=""), overwrite = TRUE, filetype = "GTiff")
#saveRDS(filled_layers, paste(stackpath ,"filled_layers.rds",sep=""))
#
# ########## load and create stacks for use
# #filled_layers <- readRDS(paste(stackpath ,"filled_layers_new.tif",sep=""))
mystack <- stack(paste(stackpath ,"filled_layers_new.tif",sep=""))
rm(masked_layers)
rm(land_mask)
rm(this_layer)
gc()
#mystack <- filled_layers
rm(filled_layers)
rm(land_mask)
gc()

print("plot mystack")
plot(mystack)
png(filename = paste(plots_path, "mystack.png", sep = ""), width = 800, height = 600)
dev.off()

e <- extent(-25, 45, 30, 72) #xmin, xmax,ymin,ymax
rasterstack.filled.layers.Europe.2025 <- crop(mystack, e)
#plot(rasterstack.filled.layers.Europe.2025)

filename <- paste(stackpath,"/","Biooracle.filled.layers.global2025.tif", sep="")
print("save filled layers)")

writeRaster(mystack, filename, format="GTiff",overwrite=TRUE)
print(filename)
filename2 <- paste(stackpath,"/","Biooracle.filled.layers.Europe2025.tif", sep="")
writeRaster(rasterstack.filled.layers.Europe.2025, filename2, format="GTiff",overwrite=TRUE)

print("save layernames")
layernames <- names(mystack)
save(layernames, file = paste(stackpath,"/layernames.filled.rda" ,sep =""))
rm(mystack)
rm(rasterstack.filled.layers.Europe.2025)
gc()
}
# # get full layernames
# load(paste(stackpath,"/layernames.rda" ,sep =""))
#
# names(mystack) <- layernames