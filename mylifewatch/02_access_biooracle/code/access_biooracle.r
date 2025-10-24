
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



inputs_path = "/mnt/inputs/"
outputs_path <- "/mnt/outputs/"


biooracle_datalayer_path <- paste0(outputs_path,"biooracle/datalayer.nc")
if(!dir.exists(biooracle_datalayer_path)){
  dir.create(biooracle_datalayer_path, recursive = TRUE)
}


och_variables_url=args$och_variables_url
constraints_base_time_from = args$constraints_base_time_from
constraints_base_time_to = args$constraints_base_time_to
constraints_ssp_time_from = args$constraints_ssp_time_from
constraints_ssp_time_to = args$constraints_ssp_time_to
constraints_dec50_time_from = args$constraints_dec50_time_from
constraints_dec50_time_to = args$constraints_dec50_time_to
constraints_dec100_time_from = args$constraints_dec100_time_from
constraints_dec100_time_to = args$constraints_dec100_time_to

dataset_och_variabler = paste0(inputs_path,"dataset_och_variabler_2025.xlsx")

download_zip_data_if_not_present_and_unzip(
    data_path = dataset_och_variabler,
    data_url = och_variables_url,
    dest_path = inputs_path
    )

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

print(paste("dataset_scenarios: ",dataset_scenarios))
dataset_scenarios<- dataset_scenarios[-8]#remove mean
dataset_scenarios_titles <- unique(layers.bio2$title)

#get variables from par
print("et variables from par")
info_layer("par_mean_baseline_2000_2020_depthsurf")

variable.selection <- read_excel(dataset_och_variabler)
print(paste("variable.selection len: ",length(variable.selection)))
#remove not used
print("remove not used")
variable.selection <- variable.selection[!is.na(variable.selection$Final),]
variable.selection <-variable.selection[-3,]# chl does not exist at depthmean ssp119
print(paste("variable.selection len: ",length(variable.selection)))

dec.vec <- c("baseline", "dec50", "dec100")
# --- Outer scenario loop ---
for (scenario in dataset_scenarios){
    print(paste("Scenario:", scenario))
    # --- Decadal loop ---
    for (dec in dec.vec[c(2, 3)]) {
            constraints.base <- list(
            #latitude = c(25, 80),
            #longitude = c(-15, 40),
            time = c(constraints_base_time_from, constraints_base_time_to)
        )
        constraints.ssp <- list(
            #latitude = c(25, 80),
            #longitude = c(-15, 40),
            time = c(constraints_ssp_time_from, constraints_ssp_time_to)
        )
        constraints.dec50 <- list(
            #latitude = c(25, 80),
            #longitude = c(-15, 40),
            time = c(constraints_dec50_time_from, constraints_dec50_time_to)
        )
        constraints.dec100 <- list(
            #latitude = c(25, 80),
            #longitude = c(-15, 40),
            time = c(constraints_dec100_time_from, constraints_dec100_time_to)
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
                filename1_nc <- glue("{biooracle_datalayer_path}{dataset_id}_{variables}.nc")
                outdir_datalayer_nc <- paste(biooracle_datalayer_path, scenario, dec, "/", sep = "")
                if (!dir.exists(outdir_datalayer_nc)) dir.create(outdir_datalayer_nc, recursive = TRUE)
                depth <- strsplit(dataset_id, "_")[[1]][length(strsplit(dataset_id, "_")[[1]])]
                filename_nc <- paste(outdir_datalayer_nc, variable, "_", depth, ".nc", sep = "")
                if (file.exists(filename_nc)) {
                    next
                }
                print(paste("file: ", filename_nc," does not exist, downloading..."))
                a <- download_layers(dataset_id, variables = variable, constraints = constraints, directory = biooracle_datalayer_path)
                filename_with_ext <- basename(terra::sources(a))[1]
                file.rename(
                    from = glue("{biooracle_datalayer_path}{filename_with_ext}"),
                    to = glue("{biooracle_datalayer_path}{dataset_id}_{variables}.nc")
                )
                print(paste("Renaming ", filename1_nc, " to ", filename_nc))
                file.rename(filename1_nc, filename_nc)
                break
            } # --- End variable loop ---
            break
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
        indir <- paste(biooracle_datalayer_path, dataset_scenarios[[1]][1], "/", sep = "")
        indir <- gsub(" ", "", indir)
        outdir <- paste(biooracle_datalayer_path, scenario, dec, "/", sep = "")

        f_name_dest = paste(outdir, "par_mean_mean_depthsurf.nc", sep = "")
        f_name_source = paste(indir, "par_mean_mean_depthsurf.nc", sep = "")
#         print(paste(outdir, f_name_dest," exists", file.exists(f_name_dest)))
        file.copy(f_name_source, f_name_dest, overwrite = TRUE)
#         print(paste(f_name_dest, file.exists(f_name_dest)))
    } # --- End decadal loop for copying ---
} # --- End scenario loop for copying ---



# ##################################################################
# ####################################################
# tRANSOM NC TO TIF



# ## rename files
# # --- Scenario loop for renaming ---
# for (sel.sen in c(4:7)) {
#     scenario <- dataset_scenarios[[sel.sen]][1]
#     print(paste("Renaming files for scenario:", scenario))
#     indir <- paste(biooracle_datalayer_path, "datalayer.nc/", scenario, dec, "/", sep = "")
#     outdir <- paste(biooracle_datalayer_path, "datalayer.tiff/", scenario, dec, "/", sep = "")
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