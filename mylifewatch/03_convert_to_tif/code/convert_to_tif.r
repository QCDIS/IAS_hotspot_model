library(raster)


inputs_path = "/mnt/inputs/"
outputs_path <- "/mnt/outputs/"


converted_biooracle_dir <- paste0(outputs_path,"biooracle/datalayer_tif/")

biooracle_datalayer_path <- paste0(inputs_path,"biooracle/datalayer_nc/")

dataset_scenarios <- list.dirs(
  paste(biooracle_datalayer_path, sep = ""),
  full.names = FALSE,
  recursive = FALSE
)

for (scenario in dataset_scenarios) {
    nc_dir <- paste(biooracle_datalayer_path, scenario, "/", sep = "")
    tif_dir <- paste(converted_biooracle_dir, scenario, "/", sep = "")
    if (!dir.exists(tif_dir)) dir.create(tif_dir, recursive = TRUE)
    nc_files <- list.files(nc_dir, pattern = "\\.nc$", full.names = TRUE)
    for (nc_file in nc_files) {
      r <- brick(nc_file)
      tif_file <- paste0(
        tif_dir,
        tools::file_path_sans_ext(basename(nc_file)),
        ".tif"
      )
      writeRaster(r, tif_file, format = "GTiff", overwrite = TRUE)
      print(paste("Converted", nc_file, "to", tif_file))
  }
}