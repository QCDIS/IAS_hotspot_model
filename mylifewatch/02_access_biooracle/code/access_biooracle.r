
library(biomod2)
library(terra)
library(raster)
library(devtools)

library(biooracler)
library(rgbif)
library(sf)
library(dplyr)
library(ggplot2)
library(maps)
library(glue)
library(readxl)# to rename files

# Explore datasets in the package
#list_datasets()

# Explore layers in a dataset
#list_layers() 



output_dir <- "/mnt/outputs/biooracle/"
if(!dir.exists(output_dir)){
  dir.create(output_dir, recursive = TRUE)
}

och_variables_url <- "https://github.com/QCDIS/IAS_hotspot_model/raw/refs/heads/master/Modelling_trial_2025/input%20data/dataset%20och%20variabler%202025.xlsx"
destfile = "/mnt/inputs/cache/dataset_och variabler_2025.xlsx"


cahche_dir <- "/mnt/inputs/cache/"
if(!dir.exists(cahche_dir)){
  dir.create(cahche_dir, recursive = TRUE)
}

#
### BIO-ORACLE LAYERS
# Explore environmental variables from Bio-Oracle
layers.bio2 <- as.data.frame(list_layers(simplify=F) )#simplify=F if need more info

#remove irrelevant
exclude <- c(grep("terrain",layers.bio2$dataset_id),grep("kdpar",layers.bio2$dataset_id))
layers.bio2<- layers.bio2[-exclude,]

#Find all variables
print("Find all variables")
dataset_variables <- unique(lapply(layers.bio2$dataset_id, function(i)
  strsplit(i,"_")[[1]][1]
))
dataset_scenarios <- unique(lapply(layers.bio2$dataset_id, function(i)
  strsplit(i,"_")[[1]][2]
))

print(paste("dataset_scenarios: ",length(dataset_scenarios)))
dataset_scenarios<- dataset_scenarios[-8]#remove mean
dataset_scenarios_titles <- unique(layers.bio2$title)
print(paste("dataset_scenarios_titles: ",length(dataset_scenarios_titles)))

#get variables from par
print("et variables from par")
info_layer("par_mean_baseline_2000_2020_depthsurf")

#layers.bio2$dataset_id[grep("_mean",layers.bio2$dataset_id)]

#get variables from layer
#layers.bio2[1,]
#layers.bio2[1,c("title")]

#testlayer <- layers.bio2$dataset_id[1]
#testlayerspellout <- layers.bio2$title[1]


###
#scenario <- dataset_scenarios[1]
#selected.datasets.nr <-grep(scenario,layers.bio2$dataset_id )
#layers.bio2$dataset_id[selected.datasets.nr ]
#for(i in selected.datasets.nr){
#  print(paste(layers.bio2$title[i], layers.bio2$dataset_id[i], sep=";"))
 # print(paste(unlist(info_layer(layers.bio2$dataset_id[i])[["variables"]]["variable_name"]),collapse=";"))
#}

if(file.exists(destfile)){
    print("File already exists, no need to download")
}else{
    print("Downloading file")
    #download xlsx file with variable selection from och_variables_url
    download.file(
      url = och_variables_url,
      destfile = destfile,
      mode = "wb"
    )
}

variable.selection <- read_excel(destfile)
print(paste("variable.selection len: ",length(variable.selection)))
#remove not used
print("remove not used")
variable.selection <- variable.selection[!is.na(variable.selection$Final),]
variable.selection <-variable.selection[-3,]# chl does not exist at depthmean ssp119
print(paste("variable.selection len: ",length(variable.selection)))

dec.vec <- c("", "dec50", "dec100")
# --- Outer scenario loop ---
for (sel.sen in c(4:7)) {
    scenario <- dataset_scenarios[[sel.sen]][1]
    # --- Decadal loop ---
    for (dec in dec.vec[c(2, 3)]) {

        ssp_scenarios <- list()
        constraints.base <- list(
            #latitude = c(25, 80),
            #longitude = c(-15, 40),
            time = c("2010-01-01T00:00:00Z", "2010-01-01T00:00:00Z")
        )
        constraints.ssp <- list(
            #latitude = c(25, 80),
            #longitude = c(-15, 40),
            time = c("2020-01-01T00:00:00Z", "2020-01-01T00:00:00Z")
        )
        constraints.dec50 <- list(
            #latitude = c(25, 80),
            #longitude = c(-15, 40),
            time = c("2040-01-01T00:00:00Z", "2040-01-01T00:00:00Z")
        )
        constraints.dec100 <- list(
            #latitude = c(25, 80),
            #longitude = c(-15, 40),
            time = c("2090-01-01T00:00:00Z", "2090-01-01T00:00:00Z")
        )
        i <- 1
        if (scenario == "baseline") {
            print("Scenario is baseline")
            constraints <- constraints.base
        } else {
            constraints <- constraints.ssp
            if (dec == "dec50") constraints <- constraints.dec50
            if (dec == "dec100") constraints <- constraints.dec100
        }

        scenario.info <- list()
        scenario.info[[scenario]] <- list()
        scenario.info[[scenario]][["datasets"]] <- list()
        number <- ifelse(scenario == "baseline", length(variable.selection$Spellout), length(variable.selection$Spellout) - 1)

        # --- Variable selection loop ---
        for (i in 1:number) {
            scenario.info[[scenario]][["datasets"]][[i]] <- list(
                dataset_id = variable.selection$dataset[i],
                variables = gsub(";", ",", variable.selection$Final[i]),
                constraints = constraints
            )
        } # --- End variable selection loop ---

        # --- Dataset loop ---
        for (dataset.nr in 1:length(scenario.info[[scenario]][["datasets"]])) {
            dataset <- scenario.info[[scenario]][["datasets"]][[dataset.nr]]
            dataset_id <- dataset$dataset_id
            if (!(scenario == "baseline")) {
                scenario.pos <- which(!is.na(match(strsplit(dataset_id, "_")[[1]], "baseline")))
                variable.vector <- strsplit(dataset_id, "_")[[1]]
                variable.vector[scenario.pos] <- scenario
                variable.vector[scenario.pos + 1] <- 2020
                variable.vector[scenario.pos + 2] <- 2100
                dataset_id <- paste(variable.vector, collapse = "_")
            }
            strsplit(dataset_id, "_")[[1]]
            variables <- dataset$variables
            constraints <- dataset$constraints
            # --- Variable loop ---
            for (variable in strsplit(variables, ",")[[1]]) {
#                 info_layer(dataset_id) # We get Error: 429 Too Many Requests
                filename1_nc <- glue("{output_dir}{dataset_id}_{variables}.nc")
                outdir <- paste(output_dir, "datalayer.nc/", scenario, dec, "/", sep = "")
#                 outdir <- paste(output_dir, "datalayer.tif/", scenario, dec, "/", sep = "")
                if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
                depth <- strsplit(dataset_id, "_")[[1]][length(strsplit(dataset_id, "_")[[1]])]
#                 filename_tif <- paste(outdir, variable, "_", depth, ".tif", sep = "")
                filename_nc <- paste(outdir, variable, "_", depth, ".nc", sep = "")
                if (file.exists(filename_nc)) {
                    next
                }
                print(paste("file: ", filename_nc," does not exist, downloading..."))
                a <- download_layers(dataset_id, variables = variable, constraints = constraints, directory = output_dir)
                filename_with_ext <- basename(terra::sources(a))[1]
                file.rename(
                    from = glue("{output_dir}{filename_with_ext}"),
                    to = glue("{output_dir}{dataset_id}_{variables}.nc")
                )
                b <- brick(glue("{output_dir}{dataset_id}_{variables}.nc"))
#                 print(paste("Renaming ", filename1_nc, " to ", filename_tif))
#                 file.rename(filename1_nc, filename_tif)
                print(paste("Renaming ", filename1_nc, " to ", filename_nc))
                file.rename(filename1_nc, filename_nc)
            } # --- End variable loop ---
        } # --- End dataset loop ---
    } # --- End decadal loop ---
} # --- End outer scenario loop ---

# ##################################################################
# ######## copy the missing layer to the new folders
# ##################################################################
dec.vec <- c("", "dec50", "dec100")
# --- Scenario loop for copying ---
for (sel.sen in c(2:7)) {
    scenario <- dataset_scenarios[[sel.sen]][1]
    # clear whitespace
    scenario <- gsub(" ", "", scenario)
    # --- Decadal loop for copying ---
    for (dec in dec.vec[c(2, 3)]) {
        indir <- paste(output_dir, "datalayer.nc/", dataset_scenarios[[1]][1], "/", sep = "")
        indir <- gsub(" ", "", indir)
        outdir <- paste(output_dir, "datalayer.nc/", scenario, dec, "/", sep = "")

        f_name_dest = paste(outdir, "par_mean_mean_depthsurf.nc", sep = "")
        f_name_source = paste(indir, "par_mean_mean_depthsurf.nc", sep = "")
#         print(paste(outdir, f_name_dest," exists", file.exists(f_name_dest)))
        file.copy(f_name_source, f_name_dest, overwrite = TRUE)
#         print(paste(f_name_dest, file.exists(f_name_dest)))
    } # --- End decadal loop for copying ---
} # --- End scenario loop for copying ---

# Rename datalayer.nc/ to datalayer.tiff/

# ##################################################################
# ####################################################
# ## rename files
# # --- Scenario loop for renaming ---
# for (sel.sen in c(4:7)) {
#     scenario <- dataset_scenarios[[sel.sen]][1]
#     print(paste("Renaming files for scenario:", scenario))
#     indir <- paste(output_dir, "datalayer.nc/", scenario, dec, "/", sep = "")
#     outdir <- paste(output_dir, "datalayer.tiff/", scenario, dec, "/", sep = "")
#     if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
#     print(paste("outdir:", outdir))
#     a <- list.files(indir, pattern = "*_1.tif")
#     # --- File renaming loop ---
#     for (i in a) {
#         print(paste("Renaming file:", i))
#         i <- paste(outdir, i, sep = "")
#         print(paste("Renaming file:", i))
#         file.rename(i, gsub("_1.tif", ".tif", i))
#     } # --- End file renaming loop ---
# } # --- End scenario loop for renaming ---