##require(rgdal)
require(CoordinateCleaner)
require(speciesgeocodeR)# installed via zip-tar frpm CRAN
require(sf)
require(raster)
require(biooracler)


species_file <- "/mnt/inputs/0010903-240202131308920.csv"# Arter från Matthias Mars 2024
file_url <- "/mnt/inputs/NIS_list_combined_Mar2025_v2.csv"
plots_folder <- "/mnt/outputs/speciesplots/"
if (!dir.exists(plots_folder)){
  dir.create(plots_folder, recursive = TRUE)
}
test_plot <- "test_plot.jpg"
cleanput_plot = "/mnt/outputs/cleanput_plot.jpg"
filled_layers_plot = "/mnt/outputs/filled_layers_plot.jpg"
testplot_cleanput_mars_2025 = "testplot_cleanput_mars_2025.jpg"
test_plot_cleanput_2024 = "testplot_cleanput_mar2024.jpg"
filtered_clean_coordinates_output_mars_2025 = "/mnt/outputs/filtered_clean_coordinates_output_mars_2025.rda"
filtered_clean_coordinates_output_mars_2025_2 = "/mnt/outputs/filtered_clean_coordinates_output_mars_2025_2.rda"
filtered_clean_marine_coordinates_output_mars_2025 = "/mnt/outputs/filtered_clean_marine_coordinates_output_mars_2025.rda"
filled_layers_file <- "/mnt/inputs/filled_layers_new.tif"

speciespath <- "/mnt/outputs/species/"

if(!dir.exists(speciespath)){
  dir.create(speciespath, recursive = TRUE)
}

cahche_dir <- "/mnt/inputs/cache/"
if(!dir.exists(cahche_dir)){
  dir.create(cahche_dir, recursive = TRUE)
}
cleanput_cache = paste(cahche_dir, "cleaned_coordinates.rds", sep = "")
shapefile_zip = paste(cahche_dir, "shapefile.shp.zip", sep = "")

# If shapefile_zip exits skip download
if(file.exists(shapefile_zip)){
    print(paste("Shapefile zip file already exists:", shapefile_zip))
    } else {
    print(paste("Shapefile zip file does not exist, will download to:", shapefile_zip))
    shapefile_link=args$shapefile_link
    # Download file from shapefile_link
    download.file(shapefile_link, destfile = shapefile_zip, mode = "wb")
}
# Unzip the downloaded file
unzip(shapefile_zip, exdir = ".")
shapefile <- list.files(".", pattern = "\\.shp$", full.names = TRUE)[1]

input <- read.delim(species_file,header =T,sep="\t",
                    na = c("", "NA"))#[ 1:1000,]
to.use <- complete.cases(input[,c("decimalLatitude" ,"decimalLongitude")])

input <- input[to.use,]
names(input)
#selcol <- c(1,3,12,16,62,68,70,71,71,73,76,77,84,108,108,109,115,116,122,137,138,139,175,189)
#head(input[,selcol])


################################# ################################# #################################
################################# plot  ################################# #################################
print("Plotting initial data")
unique(input$species)
factor <- as.numeric(as.factor(input$species))
fac1 <- ceiling(factor/10)
fac2 <- factor - 10*(fac1-1)

unique(fac1)
#names(input)[selcol]

#unique(input[,84])# occurance status
#which(input[,84]== "ABSENT")
#which(input[,84]== "PRESENT")

#unique(input[,112])# habitat

#unique(input[,115])# sampling effort
#unique(input[,137])# locationremarks


 #<- readOGR(dsn=paste(path,"data/ref-countries-2020-01m.shp/CNTR_RG_01M_2020_4326.shp",sep=""),layer="CNTR_RG_01M_2020_4326")

world1 <- st_read(shapefile,layer="CNTR_RG_01M_2020_4326")
print("Read shapefile")
#If you only want the outlines, i.e. only the geometry, you need to plot the geometry part of the sf object. Try plot(st_geometry(NJ_Map_Road)) or plot(NJ_Map_Road$geometry), both should work.
plot(world1$geometry)


# #xlim <- c(-180,180)
# #ylim <- c(-60, 84)
par(mfrow=c(1,1))
xlim <- c(-100,20)
ylim <- c(8, 68)
#jpeg("testplot2.jpg", width = 90*( xlim[2] - xlim[1]),height = 90*( ylim[2] - ylim[1]), pointsize = 10)
test_plot_file_name =  paste(plots_folder,test_plot)
jpeg(test_plot_file_name, width = 10*( xlim[2] - xlim[1]),height = 10*( ylim[2] - ylim[1]), pointsize = 4)
plot(world1$geometry, xlim = xlim, ylim = ylim, col = "light grey")
points(input$"decimalLongitude",input$"decimalLatitude", col = fac2, pch = fac1)
dev.off()
print(paste("plot written to", test_plot))

################################ ################################# #################################
####################### clean coordinates #################################
################################ ################################# #################################
print("Clean coordinates")
# if filtered_clean_coordinates_output_mars_2025 exists, skip cleaning
if(file.exists(cleanput_cache)){
    print(paste("Filtered clean coordinates file already exists:", cleanput_cache))
    cleanput <- readRDS(cleanput_cache)
    } else {
    print(paste("Filtered clean coordinates file does not exist, will create:", cleanput_cache))
    summary(input)
    names(input)[c(22:23)]
    cleanput <- input[-c(1),]
    cleanput[,23] <- as.numeric(cleanput[,23])

    names(cleanput)[c(22:23)] <- c("decimalLatitude" , "decimalLongitude")
    hist(cleanput[,22])
    hist(cleanput[,23])
    which(!is.numeric(cleanput[,22]))
    cleanput[,22]<- as.numeric(cleanput[,22])

    print("Running clean_coordinates")
    cleanput <- clean_coordinates(x = cleanput)
    print("Done clean_coordinates")

    jpeg(cleanput_plot, width = 1000,height=1000, pointsize = 10)
    plot(cleanput)
    dev.off()
    print(paste("plot written to", cleanput_plot))

    save(cleanput, file = filtered_clean_coordinates_output_mars_2025)
    saveRDS(cleanput, cleanput_cache)
}


# lines below used to merge datasets
#cleanput2 <- cleanput
#load(filtered_clean_coordinates_output_mars_2025)

# #Merge two downloads
# print("Merge two downloads")
#cleanput <- rbind(cleanput,cleanput2)
# #################################### plotta på världskarta igen ################################################
print("plot on world map again")
factor <- as.numeric(as.factor(cleanput$species))
fac1 <- ceiling(factor/10)
fac2 <- factor - 10*(fac1-1)

unique(fac1)
testplot_cleanput_mars_2025_filename = paste(plots_folder,testplot_cleanput_mars_2025, sep = "")
jpeg(testplot_cleanput_mars_2025_filename, width = 90*( xlim[2] - xlim[1]),height = 90*( ylim[2] - ylim[1]), pointsize = 10)
plot(world1$geometry, xlim = xlim, ylim = ylim, col = "light grey")
points(cleanput$"decimalLongitude",cleanput$"decimalLatitude", col = fac2, pch = fac1)
dev.off()

#
# # ####################################artvisa plottar ################################################
print("Plot each species separately")
all.species <- unique(cleanput$species)
xlim <- c(-180,180)
ylim <- c(-40, 80)

#for(s in all.species[c(8,9)]){
for(s in all.species){
    print(paste("Plotting species", s))
    plot_filename = paste(plots_folder, s,".jpg", sep = "")
    jpeg(plot_filename, width = 18*( xlim[2] - xlim[1]),height = 18*( ylim[2] - ylim[1]), pointsize = 4)
    plot(world1$geometry, xlim = xlim, ylim = ylim, col = "light grey")
    points(cleanput$"decimalLongitude"[which(cleanput$species == s)],cleanput$"decimalLatitude"[which(cleanput$species == s)],
           col = ifelse(cleanput$occurrenceStatus == "PRESENT",2,4),
           pch = ifelse(cleanput$occurrenceStatus == "PRESENT","+","o"), cex=4)
    text(-170, -10, paste(s, "number of findings=", length(which(cleanput$species == s))), cex=20)
    p <-length(which(cleanput$species == s & cleanput$occurrenceStatus == "PRESENT"))
    a <- length(which(cleanput$species == s & cleanput$occurrenceStatus == "ABSENT"))
    text(-170, -20, paste("n presence =", p, ": n absence =", a), cex=20)

    title(main = s, cex=20)
    dev.off()
}

############################## check if filter worked
print("check if filter worked")

for( u in unique(cleanput$basisOfRecord)){
  n <- length(which(cleanput$basisOfRecord == u))
  print(paste("BasisOfRecord", u,":", n,"cases"))
}
#[1] "HUMAN_OBSERVATION"   "PRESERVED_SPECIMEN"  "FOSSIL_SPECIMEN"     "OCCURRENCE"          "MATERIAL_CITATION"   "MATERIAL_SAMPLE"
#[7] "OBSERVATION"         "MACHINE_OBSERVATION" "LIVING_SPECIMEN"

exclude.Basis <- c("MACHINE_OBSERVATION","PRESERVED_SPECIMEN","FOSSIL_SPECIMEN",
             "MATERIAL_CITATION"  , "MATERIAL_SAMPLE" ,"LIVING_SPECIMEN")
exclude.Basis.pattern <- paste(exclude.Basis, collapse = "|")
exlude.index <- grep (exclude.Basis.pattern, cleanput$basisOfRecord)
print("Excluding entries")
filtered.cleanput <- cleanput[-exlude.index,]
unique(filtered.cleanput$basisOfRecord)
for( u in unique(filtered.cleanput$basisOfRecord)){
  n <- length(which(filtered.cleanput$basisOfRecord == u))
  print(paste("BasisOfRecord", u,":", n,"cases"))
}
save(filtered.cleanput, file = filtered_clean_coordinates_output_mars_2025_2)

save(filtered.cleanput, file = "filtered.clean.coordinates. mars 2024.output.rda")
##################### speciesgeododeR
#load( filtered_clean_coordinates_output_mars_2025_2)
#####load( "filtered.clean.coordinates. download feb 15.output.rda")


########## filter out data outside marine environment
print("Filter out data outside marine environment")
# stackpath <- paste(path, c("data/Biooracle.download/rasterstacks/baseline"), sep="/")
# filled_layers_file <- paste(stackpath,"/","filled_layers_new.tif", sep="")
mystack <- raster(filled_layers_file)
mask <- mean(mystack)
jpeg(filled_layers_plot,width=1000,height=1000, pointsize=10)
plot(mask)
print("mask plotted")

observations <- unique(filtered.cleanput[,c("gbifID","decimalLongitude", "decimalLatitude", "occurrenceStatus")])
names(observations) <- c("ID","Lon","Lat","occurrenceStatus")

#coordinatsystem. WGS 84 decimal-latitude. epsg:4326
observations$Lon <- as.numeric(observations$Lon)
observations$Lat <- as.numeric(observations$Lat)
coord<-as.data.frame(observations[,c("Lon","Lat")])
names(coord) <- c("Lon","Lat")

observations1 <- observations
#coord <- coord[c(1:5),]
#coordinates(coord) <- ~Lon+Lat

# head(observations)
print(paste("Number of unique observations to check for marine occurrences:", length(observations$ID)))
observations2 <- as.data.frame(observations[,-c(2,3)])
#points<-SpatialPointsDataFrame(coord,
#                              observations2, proj4string=CRS("+init=epsg:4326")) # expression decapriated
points<-SpatialPointsDataFrame(coord,
                               observations2, proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
xlim <- c(extent(points)[1],extent(points)[2]) + c(-0.01, 0.01)
ylim <- c(extent(points)[3],extent(points)[4]) + c(-0.01, 0.01)

extent(mask)
# #
# # xlim2 <- extent(mask)[1:2]
# # ylim2 <- extent(mask)[3:4]
# #
# # points2<-extract(mask, points, sp=TRUE)#t
# #
# # head(points2)
# # plot(mask)
# # plot(world1$geometry)
# # points(points, col=ifelse(is.na(points2$layer),2,3))
# #
# # head(points2)
# # names(points2)
# # print(paste("Number of observations in marine environment:", length(which(!is.na(points2$layer)))))
# #
# # is.marine <- sapply(filtered.cleanput$gbifID, function(i)
# #   points2$layer[which( points2$ID == i)]
# # )
# # filtered.cleanput.marine <- filtered.cleanput[!is.na(is.marine),]
# # save(filtered.cleanput.marine, file = filtered_clean_marine_coordinates_output_mars_2025)
# # ##############################################################
# # ##################################################################
# # all.species <- unique(filtered.cleanput.marine$species)
# #
# # #for(s in all.species){
# # #   s <- all.species[18]
# # #  temp <- filtered.cleanput[ which(filtered.cleanput$species == s & filtered.cleanput$occurrenceStatus == "PRESENT"),
# # #                    c("occurrenceID", "decimallongitude","decimallatitude")]
# # #  geocode1 <- SpeciesGeoCoder(temp, world1, areanames = "CNTR_NAME")
# #
# # head(filtered.cleanput.marine)
# # print("Prepare input files for species distribution modelling")
# # #unique(filtered.cleanput$occurrenceStatus)
# # # s <- all.species[18]
# # #prepare input file
# # print(paste("selected species from", file_url))
# # selected.species <- read.csv2(file_url, sep=",")
# # cathegories <- unique(selected.species$category)
# # #my.cathegory <- cathegories[1]
# # tab <- c()
# # print(paste("Number of species to process:", length(all.species)))
# # print(paste("Categories:", cathegories))
# # print(paste("selected.species: ", head(selected.species)))
# #
# # for(s in all.species){
# #   s <- "Asparagopsis armata"
# #   print(paste("Working on species", s))
# #   comment <- "no match"
# #   id <- which(!is.na(match(selected.species$Taxon.name,s)))
# #   print(paste("species id:", id))
# #   # If we didn't find a match, try to match only the first part of the name
# #   if(length(id)<1){
# #     print("id len <1")
# #     spaces_name <- strsplit(s," ")[[1]][1]
# #     print(paste("Try to find with spaces_name:", spaces_name))
# #     print(paste("selected.species$Taxon.name len", length(selected.species$Taxon.name)))
# #     id <- grep(spaces_name,selected.species$Taxon.name)
# #     print(paste("new species id:", id))
# #     comment <- selected.species$Taxon.name[id]
# #   }else{
# #     print("id len ok")
# #     comment <- "names match"
# #   }
# #   print(paste("comment:", comment))
# #   cat <- selected.species$category[id]
# #   if (length(cat)<1){
# #     cat <- "no_match"
# #   }
# #   print(paste("category:", cat))
# #
# #   temp <- filtered.cleanput.marine[ which(filtered.cleanput$species == s ),
# #                              c("gbifID","occurrenceID","species", "occurrenceStatus", "decimalLongitude","decimalLatitude","coordinateUncertaintyInMeters",
# #                                "depth", "depthAccuracy","eventDate")]
# #
# #
# #   #fix names of deciomal coordinates. capital L
# #   print("Fix names of decimal coordinates")
# #   names(temp)<- c("gbifID","occurrenceID","species", "occurrenceStatus", "decimalLongitude","decimalLatitude","coordinateUncertaintyInMeters",
# #                   "depth", "depthAccuracy","eventDate")
# #   #head(temp)
# #   write.csv(temp,file = paste(speciespath,s,".csv", sep =""), row.names = F)
# #   print(paste("wrote: ", speciespath,s,".csv", sep =""))
# #   print(paste("no positives: ",length(which(temp$occurrenceStatus == "PRESENT")), sep =""))
# #   print(paste("no negatives: ",length(which(temp$occurrenceStatus == "ABSENT")), sep =""))
# #   my.filename <- paste(s,".csv",sep="")
# #   print(paste("my.filename:", my.filename))
# #   my.pseudoname <- paste("pseudoabsences.marine.excludebox",cat,".csv",sep="")
# #   tab <- rbind(tab, c(s, my.filename,my.filename,my.pseudoname,length(which(temp$occurrenceStatus == "PRESENT")),length(which(temp$occurrenceStatus == "ABSENT")),comment))
# #   print(paste("my.pseudoname:", my.pseudoname))
# #   print(paste("tab in loop: ", tab))
# #   print(paste("tab len: ", length(tab)))
# # }
# # print(paste("tab : ", tab))
# # # print(paste("tab len: ", length(tab)))
# # # print(paste("tab summary: ", summary(tab)))
# # colnames(tab)<- c("species","present.data","absence.data","pseudoabsence.data","n.present","n.absent","comment on Taxon")
# # write.csv2(tab, file=paste(speciespath, "data_table_mars2025.csv",sep=""))
# # print(paste("Wrote species data table to", paste(speciespath, "data_table_mars2025.csv",sep="")))
# #
# #
# # # # ###################################################################
# # # # ###################################################################
# # # filtered.cleanput<-filtered.cleanput.marine
# # # #define a box to exclude pseudoabsences
# # # print("define a box to exclude pseudoabsences")
# # # boxxlim=c(12,30)
# # # boxylim = c(50,66)
# # # xlim <- c(-180,180)
# # # ylim <- c(-60, 84)
# # # plot(world1$geometry, xlim = xlim, ylim = ylim, col = "light grey")
# # # lines( boxxlim[c(1,2,2,1,1)], boxylim[c(1,1,2,2,1)])
# # #
# # # excludebox <- intersect(
# # #   which(filtered.cleanput$decimalLongitude > boxxlim[1] & filtered.cleanput$decimalLongitude < boxxlim[2]),
# # #   which(filtered.cleanput$decimalLatitude > boxylim[1] & filtered.cleanput$decimalLatitude < boxylim[2])
# # # )
# # # #plot excluded points
# # #
# # # points(filtered.cleanput$decimalLongitude[excludebox], filtered.cleanput$decimalLatitude[excludebox], col=2)
# # # #plot included points
# # # points(filtered.cleanput$decimalLongitude[-excludebox], filtered.cleanput$decimalLatitude[-excludebox], col=3)
# # #
# # # filtered.cleanput.unbox <- filtered.cleanput[-excludebox,]
# # #
# # # ##############################################
# # # ######### define groups of species      ######
# # #
# # # selected.species <- read.csv2(file_url, sep=",")
# # # cathegories <- unique(selected.species$category)
# # #
# # # #my.cathegory <- cathegories[1]
# # # for (my.cathegory in cathegories[-5]){
# # #
# # # my.species.list <- selected.species$Taxon.name[which(selected.species$category == my.cathegory)]
# # #
# # # #get entries matching
# # # filtered.cleanput.subset <- filtered.cleanput.unbox[which(!is.na(match(filtered.cleanput.unbox$species,my.species.list))),]
# # #
# # # locationsamples <- sample(1:length(filtered.cleanput.subset$gbifID), 1000, replace =F)
# # #
# # # #locationsamples <- sample(1:length(filtered.cleanput.unbox$gbifID), 10000, replace =F)
# # # pseudoabsences <- filtered.cleanput.subset[ locationsamples,
# # #                                            c("gbifID","occurrenceID","species", "occurrenceStatus", "decimalLongitude","decimalLatitude","coordinateUncertaintyInMeters",
# # #                                              "depth", "depthAccuracy","eventDate")]
# # # #pseudoabsences <- filtered.cleanput.unbox[ locationsamples,
# # # #                                     c("gbifID","occurrenceID","species", "occurrenceStatus", "decimalLongitude","decimalLatitude","coordinateUncertaintyInMeters",
# # # #                                       "depth", "depthAccuracy","eventDate")]
# # # pseudoabsences$gbifID <- paste("pseudo",seq(1:1000), sep ="")
# # # pseudoabsences$species <- NA
# # # pseudoabsences$occurrenceStatus <- "ABSENT"
# # # pseudoabsences$coordinateUncertaintyInMeters <- NA
# # # pseudoabsences$depthAccuracy <- NA
# # # pseudoabsences$eventDate <- NA
# # # head(pseudoabsences)
# # # names(pseudoabsences)<- c("gbifID","occurrenceID","species", "occurrenceStatus", "decimalLongitude","decimalLatitude","coordinateUncertaintyInMeters",
# # #                 "depth", "depthAccuracy","eventDate")
# # #
# # # xlim <- c(-180,180)
# # # ylim <- c(-60, 84)
# # #
# # # #boxxlim=c(12,30)
# # # #boxylim = c(50,66)
# # # jpeg(paste(path,"/speciesplots/","testplot.pseudoabsences.",my.cathegory,".jpg",sep=""), width = 18*( xlim[2] - xlim[1]),height = 18*( ylim[2] - ylim[1]), pointsize = 4)
# # # plot(world1$geometry, xlim = xlim, ylim = ylim, col = "light grey")
# # # lines( boxxlim[c(1,2,2,1,1)], boxylim[c(1,1,2,2,1)])
# # # points(pseudoabsences$"decimalLongitude",pseudoabsences$"decimalLatitude", col = "red", pch = "*", cex = 5)
# # # dev.off()
# # #
# # # write.csv(pseudoabsences,file = paste(speciespath,"pseudoabsences.marine.excludebox",my.cathegory,".csv", sep =""), row.names = F)
# # # print(paste("wrote: ", speciespath,"/pseudoabsences.marine.excludebox",my.cathegory,".csv", sep =""))
# # #
# # #
# # # }#end cathegory
# # # ##################xlim <- c(-100,20)
# # # filtered.cleanput <- filtered.cleanput[ ,
# # #                            c("gbifID","occurrenceID","species", "occurrenceStatus", "decimalLongitude","decimalLatitude","coordinateUncertaintyInMeters",
# # #                              "depth", "depthAccuracy","eventDate")]
# # # #fix names of deciomal coordinates. capital L
# # # #names(filtered.cleanput)<- c("gbifID","occurrenceID","species", "occurrenceStatus", "decimalLongitude","decimalLatitude","coordinateUncertaintyInMeters",
# # # #                "depth", "depthAccuracy","eventDate")
# # #
# # # world1 <- st_read(shapefile,layer="CNTR_RG_01M_2020_4326")
# # #
# # # factor <- as.numeric(as.factor(filtered.cleanput$species))
# # # fac1 <- ceiling(factor/10)
# # # fac2 <- factor - 10*(fac1-1)
# # # xlim <- c(-110,40)
# # # ylim <- c(-2, 68)#jpeg("testplot2.jpg", width = 90*( xlim[2] - xlim[1]),height = 90*( ylim[2] - ylim[1]), pointsize = 10)
# # # jpeg(paste(path,"speciesplots/","testplot.CLEANPUT.SMALLER.world.mar2025.jpg", sep=""), width = 10*( xlim[2] - xlim[1]),height = 10*( ylim[2] - ylim[1]), pointsize = 4)
# # # plot(world1, xlim = xlim, ylim = ylim, col = "light grey")
# # # points(filtered.cleanput$"decimalLongitude",filtered.cleanput$"decimalLatitude", col = fac2, pch = fac1)
# # # dev.off()
# #
# # ########################
