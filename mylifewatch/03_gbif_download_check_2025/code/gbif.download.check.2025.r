## Load required libraries
require(CoordinateCleaner)
require(speciesgeocodeR)
require(sf)
require(raster)
require(biooracler)

# --- Define file paths and output folders ---
inputs_path = "/mnt/inputs/"
outputs_path <- "/mnt/outputs/"

species_file <- paste0(inputs_path,"0010903-240202131308920.csv")
download_zip_data_if_not_present_and_unzip(
    data_path = species_file,
    data_url = args$species_file_url,
    dest_path = inputs_path
    )
nis_list_path <- paste0(inputs_path,"NIS_list_combined_Mar2025_v2.csv")
download_zip_data_if_not_present_and_unzip(
    data_path = nis_list_path,
    data_url = args$nis_list_url,
    dest_path = inputs_path
    )

specie_splots_dir <- paste0(outputs_path,"speciesplots/")
if (!dir.exists(specie_splots_dir)) dir.create(specie_splots_dir, recursive = TRUE)

test_plot <- "test_plot.jpg"
cleanput_plot <- paste(specie_splots_dir, "cleanput_plot.jpg", sep = "")
filled_layers_plot <- paste(specie_splots_dir, "filled_layers_plot.jpg", sep = "")

testplot_cleanput_mars_2025 <- "testplot_cleanput_mars_2025.jpg"
test_plot_cleanput_2024 <- "testplot_cleanput_mar2024.jpg"
coordinates_output = paste0(outputs_path,"cleaned_coordinates/")
if (!dir.exists(coordinates_output)) dir.create(coordinates_output, recursive = TRUE)
filtered_clean_coordinates_output_mars_2025 <- paste(coordinates_output, "filtered_clean_coordinates_output_mars_2025.rda", sep = "")
filtered_clean_coordinates_output_mars_2025_2 <- paste(coordinates_output, "filtered_clean_coordinates_output_mars_2025_2.rda", sep = "")
filtered_clean_marine_coordinates_output_mars_2025 <- paste(coordinates_output, "filtered_clean_marine_coordinates_output_mars_2025.rda", sep = "")
filtered_clean_coordinates_mars_2024_output <- paste(coordinates_output, "filtered_clean_coordinates_mars_2024_output.rda", sep = "")
filled_layers_file <- paste0(inputs_path,"biooracle/rasterstacks/baselinedec50/filled_layers_new.tif")

speciespath <- paste0(outputs_path,"species/")
if (!dir.exists(speciespath)) dir.create(speciespath, recursive = TRUE)

cleanput_cache <- paste(inputs_path, "cleaned_coordinates.rds", sep = "")
shapefile_zip <- paste(inputs_path , "shapefile.shp.zip", sep = "")

# --- Download and unzip shapefile if needed ---
download_zip_data_if_not_present_and_unzip(
    data_path = shapefile_zip,
    data_url = args$shapefile_link,
    dest_path = inputs_path
    )
shapefile <- list.files(inputs_path, pattern = "\\.shp$", full.names = TRUE)[1]
print(paste("Using shapefile:", shapefile))


# --- Read and clean input data ---
input <- read.delim(species_file, header = TRUE, sep = "\t", na = c("", "NA"))
to.use <- complete.cases(input[, c("decimalLatitude", "decimalLongitude")])
input <- input[to.use, ]

# --- Plot initial data ---
print("Plotting initial data")
factor <- as.numeric(as.factor(input$species))
fac1 <- ceiling(factor / 10)
fac2 <- factor - 10 * (fac1 - 1)
world1 <- st_read(shapefile, layer = "CNTR_RG_01M_2020_4326")
plot(world1$geometry)

par(mfrow = c(1, 1))
xlim <- c(-100, 20)
ylim <- c(8, 68)
test_plot_file_name <- paste(specie_splots_dir, test_plot)
jpeg(test_plot_file_name, width = 10 * (xlim[2] - xlim[1]), height = 10 * (ylim[2] - ylim[1]), pointsize = 4)
plot(world1$geometry, xlim = xlim, ylim = ylim, col = "light grey")
points(input$decimalLongitude, input$decimalLatitude, col = fac2, pch = fac1)
dev.off()
print(paste("plot written to", test_plot))

# --- Clean coordinates ---
print("Clean coordinates")
if (file.exists(cleanput_cache)) {
  print(paste("Filtered clean coordinates file already exists:", cleanput_cache))
  cleanput <- readRDS(cleanput_cache)
} else {
  print(paste("Filtered clean coordinates file does not exist, will create:", cleanput_cache))
  cleanput <- input[-c(1), ]
  cleanput[, 23] <- as.numeric(cleanput[, 23])
  names(cleanput)[c(22:23)] <- c("decimalLatitude", "decimalLongitude")
  cleanput[, 22] <- as.numeric(cleanput[, 22])
  cleanput <- clean_coordinates(x = cleanput)
  jpeg(cleanput_plot, width = 1000, height = 1000, pointsize = 10)
  plot(cleanput)
  dev.off()
  print(paste("plot written to", cleanput_plot))
  save(cleanput, file = filtered_clean_coordinates_output_mars_2025)
  saveRDS(cleanput, cleanput_cache)
}

# --- Plot cleaned data ---
print("plot on world map again")
factor <- as.numeric(as.factor(cleanput$species))
fac1 <- ceiling(factor / 10)
fac2 <- factor - 10 * (fac1 - 1)
testplot_cleanput_mars_2025_filename <- paste(specie_splots_dir, testplot_cleanput_mars_2025, sep = "")
jpeg(testplot_cleanput_mars_2025_filename, width = 90 * (xlim[2] - xlim[1]), height = 90 * (ylim[2] - ylim[1]), pointsize = 10)
plot(world1$geometry, xlim = xlim, ylim = ylim, col = "light grey")
points(cleanput$decimalLongitude, cleanput$decimalLatitude, col = fac2, pch = fac1)
dev.off()

# --- Plot each species separately ---
print("Plot each species separately")
all.species <- unique(cleanput$species)
xlim <- c(-180, 180)
ylim <- c(-40, 80)
for (s in all.species) {
  print(paste("Plotting species", s))
  plot_filename <- paste(specie_splots_dir, s, ".jpg", sep = "")
  if (file.exists(plot_filename)) {
    print(paste("Plot already exists, skipping:", plot_filename))
    next
  }
  jpeg(plot_filename, width = 18 * (xlim[2] - xlim[1]), height = 18 * (ylim[2] - ylim[1]), pointsize = 4)
  plot(world1$geometry, xlim = xlim, ylim = ylim, col = "light grey")
  points(cleanput$decimalLongitude[cleanput$species == s], cleanput$decimalLatitude[cleanput$species == s],
         col = ifelse(cleanput$occurrenceStatus == "PRESENT", 2, 4),
         pch = ifelse(cleanput$occurrenceStatus == "PRESENT", "+", "o"), cex = 4)
  text(-170, -10, paste(s, "number of findings=", length(which(cleanput$species == s))), cex = 20)
  p <- length(which(cleanput$species == s & cleanput$occurrenceStatus == "PRESENT"))
  a <- length(which(cleanput$species == s & cleanput$occurrenceStatus == "ABSENT"))
  text(-170, -20, paste("n presence =", p, ": n absence =", a), cex = 20)
  title(main = s, cex = 20)
  dev.off()
}

# --- Check if filter worked ---
print("check if filter worked")
for (u in unique(cleanput$basisOfRecord)) {
  n <- length(which(cleanput$basisOfRecord == u))
  print(paste("BasisOfRecord", u, ":", n, "cases"))
}
exclude.Basis <- c("MACHINE_OBSERVATION", "PRESERVED_SPECIMEN", "FOSSIL_SPECIMEN",
                   "MATERIAL_CITATION", "MATERIAL_SAMPLE", "LIVING_SPECIMEN")
exclude.Basis.pattern <- paste(exclude.Basis, collapse = "|")
exlude.index <- grep(exclude.Basis.pattern, cleanput$basisOfRecord)
print("Excluding entries")
filtered.cleanput <- cleanput[-exlude.index, ]
for (u in unique(filtered.cleanput$basisOfRecord)) {
  n <- length(which(filtered.cleanput$basisOfRecord == u))
  print(paste("BasisOfRecord", u, ":", n, "cases"))
}
save(filtered.cleanput, file = filtered_clean_coordinates_output_mars_2025_2)
save(filtered.cleanput, file = filtered_clean_coordinates_mars_2024_output)

# --- Filter out data outside marine environment ---
print("Filter out data outside marine environment")
mystack <- raster(filled_layers_file)
mask <- mean(mystack)
jpeg(filled_layers_plot, width = 1000, height = 1000, pointsize = 10)
plot(mask)
print("mask plotted")


if (!file.exists(filtered_clean_marine_coordinates_output_mars_2025)) {
    observations <- unique(filtered.cleanput[, c("gbifID", "decimalLongitude", "decimalLatitude", "occurrenceStatus")])
    names(observations) <- c("ID", "Lon", "Lat", "occurrenceStatus")
    observations$Lon <- as.numeric(observations$Lon)
    observations$Lat <- as.numeric(observations$Lat)
    coord <- as.data.frame(observations[, c("Lon", "Lat")])
    names(coord) <- c("Lon", "Lat")
    observations2 <- as.data.frame(observations[, -c(2, 3)])
    points <- SpatialPointsDataFrame(coord, observations2, proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
    points2 <- extract(mask, points, sp = TRUE)
    plot(mask)
    plot(world1$geometry)
    points(points, col = ifelse(is.na(points2$layer), 2, 3))
    print(paste("Number of observations in marine environment:", length(which(!is.na(points2$layer)))))
    is.marine <- sapply(filtered.cleanput$gbifID, function(i) points2$layer[which(points2$ID == i)])
    filtered.cleanput.marine <- filtered.cleanput[!is.na(is.marine), ]
    save(filtered.cleanput.marine, file = filtered_clean_marine_coordinates_output_mars_2025)
} else {
    print(paste("Loading filtered clean marine coordinates file:", filtered_clean_marine_coordinates_output_mars_2025))
    load(filtered_clean_marine_coordinates_output_mars_2025)
}


# --- Prepare input files for species distribution modelling ---
all.species <- unique(filtered.cleanput.marine$species)
all.species <- c(all.species, "Asparagopsis armata")
print("Prepare input files for species distribution modelling")
selected.species <- read.csv2(nis_list_path, sep = ",")
cathegories <- unique(selected.species$category)

tab_colnames <- c("species", "present.data", "absence.data", "pseudoabsence.data", "n.present", "n.absent", "comment on Taxon")
tab <- data.frame(matrix(ncol = length(tab_colnames), nrow = 0))
colnames(tab) <- tab_colnames


print(paste("Number of species to process:", length(all.species)))
print(paste("Categories:", cathegories))
print(paste("selected.species: ", head(selected.species)))

for (s in all.species) {
    comment <- "no match"
    id <- which(!is.na(match(selected.species$Taxon.name, s)))
    if (length(id) < 1) {
        spaces_name <- strsplit(s, " ")[[1]][1]
        id <- grep(spaces_name, selected.species$Taxon.name)
        comment <- selected.species$Taxon.name[id]
    } else {
        comment <- "names match"
    }
    cat <- selected.species$category[id]
    if (length(cat) < 1) cat <- "no_match"
    if (length(comment) < 1) comment <- "no_match"
    if (length(which(filtered.cleanput.marine$species == s)) < 1) {
        print(paste("No records for species", s, "skipping to next species"))
        next
    }
    temp <- filtered.cleanput.marine[filtered.cleanput.marine$species == s,
                                   c("gbifID", "occurrenceID", "species", "occurrenceStatus", "decimalLongitude", "decimalLatitude",
                                     "coordinateUncertaintyInMeters", "depth", "depthAccuracy", "eventDate")]
    names(temp) <- c("gbifID", "occurrenceID", "species", "occurrenceStatus", "decimalLongitude", "decimalLatitude",
                   "coordinateUncertaintyInMeters", "depth", "depthAccuracy", "eventDate")

    write.csv(temp, file = paste(speciespath, s, ".csv", sep = ""), row.names = FALSE)
    my.filename <- paste(s, ".csv", sep = "")
    my.pseudoname <- paste("pseudoabsences.marine.excludebox", cat, ".csv", sep = "")
    tab <- rbind(tab, data.frame(
      species = s,
      present.data = my.filename,
      absence.data = my.filename,
      pseudoabsence.data = my.pseudoname,
      n.present = length(which(temp$occurrenceStatus == "PRESENT")),
      n.absent = length(which(temp$occurrenceStatus == "ABSENT")),
      `comment on Taxon` = comment,
      stringsAsFactors = FALSE
    ))

}


data_table_path = paste(speciespath, "data_table.csv", sep = "")
write.csv(tab, file = data_table_path, row.names = FALSE)
print(paste("Wrote species data table to", paste(speciespath, "data_table.csv", sep = "")))

# --- Define a box to exclude pseudoabsences ---
filtered.cleanput <- filtered.cleanput.marine
print("define a box to exclude pseudoabsences")
boxxlim <- c(12, 30)
boxylim <- c(50, 66)
xlim <- c(-180, 180)
ylim <- c(-60, 84)
plot(world1$geometry, xlim = xlim, ylim = ylim, col = "light grey")
lines(boxxlim[c(1, 2, 2, 1, 1)], boxylim[c(1, 1, 2, 2, 1)])

excludebox <- intersect(
  which(filtered.cleanput$decimalLongitude > boxxlim[1] & filtered.cleanput$decimalLongitude < boxxlim[2]),
  which(filtered.cleanput$decimalLatitude > boxylim[1] & filtered.cleanput$decimalLatitude < boxylim[2])
)
points(filtered.cleanput$decimalLongitude[excludebox], filtered.cleanput$decimalLatitude[excludebox], col = 2)
points(filtered.cleanput$decimalLongitude[-excludebox], filtered.cleanput$decimalLatitude[-excludebox], col = 3)
filtered.cleanput.unbox <- filtered.cleanput[-excludebox, ]

# --- Define groups of species and generate pseudoabsences ---
selected.species <- read.csv2(nis_list_path, sep = ",")
cathegories <- unique(selected.species$category)
for (my.cathegory in cathegories) {
    my.species.list <- selected.species$Taxon.name[selected.species$category == my.cathegory]
    filtered.cleanput.subset <- filtered.cleanput.unbox[!is.na(match(filtered.cleanput.unbox$species, my.species.list)), ]
    n_samples <- min(1000, length(filtered.cleanput.subset$gbifID))
    if (n_samples <= 0) {
        print(paste("No samples for category", my.cathegory, "skipping to next category"))
        next
    }
    locationsamples <- sample(1:length(filtered.cleanput.subset$gbifID), n_samples, replace =F)
    pseudoabsences <- filtered.cleanput.subset[locationsamples,
                                             c("gbifID", "occurrenceID", "species", "occurrenceStatus", "decimalLongitude", "decimalLatitude",
                                               "coordinateUncertaintyInMeters", "depth", "depthAccuracy", "eventDate")]
    pseudoabsences$gbifID <- paste("pseudo", seq(1:1000), sep = "")
    pseudoabsences$species <- NA
    pseudoabsences$occurrenceStatus <- "ABSENT"
    pseudoabsences$coordinateUncertaintyInMeters <- NA
    pseudoabsences$depthAccuracy <- NA
    pseudoabsences$eventDate <- NA
    names(pseudoabsences) <- c("gbifID", "occurrenceID", "species", "occurrenceStatus", "decimalLongitude", "decimalLatitude",
                             "coordinateUncertaintyInMeters", "depth", "depthAccuracy", "eventDate")
    xlim <- c(-180, 180)
    ylim <- c(-60, 84)
    test_plot_pseudoabsences = paste(specie_splots_dir, "testplot.pseudoabsences.", my.cathegory, ".jpg", sep = "")
    print(paste("Plotting test_plot_pseudoabsences: ", test_plot_pseudoabsences))
    jpeg(test_plot_pseudoabsences,
       width = 18 * (xlim[2] - xlim[1]), height = 18 * (ylim[2] - ylim[1]), pointsize = 4)
    plot(world1$geometry, xlim = xlim, ylim = ylim, col = "light grey")
    lines(boxxlim[c(1, 2, 2, 1, 1)], boxylim[c(1, 1, 2, 2, 1)])
    points(pseudoabsences$decimalLongitude, pseudoabsences$decimalLatitude, col = "red", pch = "*", cex = 5)
    dev.off()
    pseudoabsences_marine_excludebox = paste(speciespath, "pseudoabsences.marine.excludebox", my.cathegory, ".csv", sep = "")
    print(paste("pseudoabsences_marine_excludebox: ", pseudoabsences_marine_excludebox))
    write.csv(pseudoabsences, file = pseudoabsences_marine_excludebox, row.names = FALSE)
} # --- End category loop ---


# --- Final plot of cleaned data ---
filtered.cleanput <- filtered.cleanput[, c("gbifID", "occurrenceID", "species", "occurrenceStatus", "decimalLongitude", "decimalLatitude",
                                           "coordinateUncertaintyInMeters", "depth", "depthAccuracy", "eventDate")]
world1 <- st_read(shapefile, layer = "CNTR_RG_01M_2020_4326")
factor <- as.numeric(as.factor(filtered.cleanput$species))
fac1 <- ceiling(factor / 10)
fac2 <- factor - 10 * (fac1 - 1)
xlim <- c(-110, 40)
ylim <- c(-2, 68)
testplot_file = paste(specie_splots_dir, "testplot.CLEANPUT.SMALLER.world.mar2025.jpg", sep = "")
print(paste("Plotting testplot_file: ", testplot_file))
jpeg(testplot_file,
     width = 10 * (xlim[2] - xlim[1]), height = 10 * (ylim[2] - ylim[1]), pointsize = 4)
plot(world1, xlim = xlim, ylim = ylim, col = "light grey")
points(filtered.cleanput$decimalLongitude, filtered.cleanput$decimalLatitude, col = fac2, pch = fac1)
dev.off()