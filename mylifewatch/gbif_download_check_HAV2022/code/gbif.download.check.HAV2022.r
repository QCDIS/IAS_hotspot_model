require(rgdal)
require(CoordinateCleaner)
require(speciesgeocodeR)# installed via zip-tar frpm CRAN
#path_home
gbifs <- "/mnt/inputs/GBIF" # The path to the root of your project
jpgs <- "/mnt/outputs/jpg" # Path to save outputs
rdas <- "/mnt/outputs/rda" # Path to save outputs
speciespath = "/mnt/inputs/speciesIndata2022" # Path to save species data
ogr_path <- "/mnt/inputs/CNTR_RG_01M_2020_4326.shp" # Path to save ogr data

#Set file with GBIF data. Path relative to project dir.
file <- "GBIF alla arter 0429558-210914110416597.csv"
#file <- "Data2022/GBIF download.2.dec 13/0213464-220831081235567.csv" # alternative file

# Read data and select complete cases
input <- read.delim(paste(gbifs,file,sep=""),header =T,sep="\t",
                    na = c("", "NA"))#[ 1:1000,]
to.use <- complete.cases(input[,c("decimalLatitude" ,"decimalLongitude")])

input <- input[to.use,]

#names(input)


################################# ################################# #################################
################################# plot raw data ################################# #################################
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


world1 <- readOGR(dsn=ogr_path,layer="CNTR_RG_01M_2020_4326")
#xlim <- c(-180,180)
#ylim <- c(-60, 84)
xlim <- c(-100,20)
ylim <- c(8, 68)
# plot at high resolution to inspect data or low resolution for illustration

#jpeg("testplot.jpg", width = 90*( xlim[2] - xlim[1]),height = 90*( ylim[2] - ylim[1]), pointsize = 10)

jpeg("testplot.SMALLER.jpg", width = 10*( xlim[2] - xlim[1]),height = 10*( ylim[2] - ylim[1]), pointsize = 4)
plot(world1, xlim = xlim, ylim = ylim, col = "light grey")
points(input$"decimalLongitude",input$"decimalLatitude", col = fac2, pch = fac1)
dev.off()

################################# ################################# #################################
######################## clean coordinates #################################
################################# ################################# #################################
summary(input)
names(input)[c(22:23)]
cleanput <- input[-c(1),]
cleanput[,23] <- as.numeric(cleanput[,23])
names(cleanput)[c(22:23)] <- c("decimallatitude" , "decimallongitude")
hist(cleanput[,22])
hist(cleanput[,23])
which(!is.numeric(cleanput[,23]))
head(cleanput)
cleanput <- clean_coordinates(x = cleanput)
summary(cleanput)

#plot to visualize
res = paste(jpgs, "plot.cleanput.result.jpg", sep = "/")
jpeg(res, width = 1000,height=1000, pointsize = 10)
plot(cleanput)
dev.off()

# Save processed data as intermediate step
coordinates_path = paste(rdas, "clean.coordinates.output.filename.rda", sep = "/")
save(cleanput, file = coordinates_path)
#load again
load(coordinates_path)
#################################### plot worldmap igen with cleaned data ################################################
factor <- as.numeric(as.factor(cleanput$species))
fac1 <- ceiling(factor/10)
fac2 <- factor - 10*(fac1-1)

unique(fac1)
testplot = paste(jpgs, "testplot.cleanput.filename.jpg", sep = "/")
jpeg(testplot, width = 90*( xlim[2] - xlim[1]),height = 90*( ylim[2] - ylim[1]), pointsize = 10)
plot(world1, xlim = xlim, ylim = ylim, col = "light grey")
points(cleanput$"decimallongitude",cleanput$"decimallatitude", col = fac2, pch = fac1)
dev.off()


#################################### specieswise plots  ################################################
all.species <- unique(cleanput$species)
for(s in all.species){

testplot2 =  paste(jpgs, "testplot.cleanput.filename3.jpg", sep = "/")
jpeg(testplot2, width = 90*( xlim[2] - xlim[1]),height = 90*( ylim[2] - ylim[1]), pointsize = 10)
plot(world1, xlim = xlim, ylim = ylim, col = "light grey")
points(cleanput$"decimallongitude"[which(cleanput$species == s)],cleanput$"decimallatitude"[which(cleanput$species == s)],
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
head(cleanput)
unique(cleanput$basisOfRecord)
for( u in unique(cleanput$basisOfRecord)){
  n <- length(which(cleanput$basisOfRecord == u))
  print(paste("BasisOfRecord", u,":", n,"cases"))
}
#}
#[1] "HUMAN_OBSERVATION"   "PRESERVED_SPECIMEN"  "FOSSIL_SPECIMEN"     "OCCURRENCE"          "MATERIAL_CITATION"   "MATERIAL_SAMPLE"
#[7] "OBSERVATION"         "MACHINE_OBSERVATION" "LIVING_SPECIMEN"

# Exclude observations with irrelevant types of observation
exclude.Basis <- c("MACHINE_OBSERVATION","PRESERVED_SPECIMEN","FOSSIL_SPECIMEN",
             "MATERIAL_CITATION"  , "MATERIAL_SAMPLE" ,"LIVING_SPECIMEN")
exclude.Basis.pattern <- paste(exclude.Basis, collapse = "|")
exlude.index <- grep (exclude.Basis.pattern, cleanput$basisOfRecord)

filtered.cleanput <- cleanput[-exlude.index,]
unique(filtered.cleanput$basisOfRecord)
for( u in unique(filtered.cleanput$basisOfRecord)){
  n <- length(which(filtered.cleanput$basisOfRecord == u))
  print(paste("BasisOfRecord", u,":", n,"cases"))
}
# Save intermediate work
filtered_path = paste(rdas,"filtered.clean.coordinates.output.filename.rda", sep = "/")
save(filtered.cleanput, file = filtered_path)
###################### speciesgeododeR
# Reload data if taking up work
#load(  "filtered.clean.coordinates.output.filename.rda")

all.species <- unique(filtered.cleanput$species)

head(filtered.cleanput)
#unique(filtered.cleanput$occurrenceStatus)
  # s <- all.species[18]

for(s in all.species){

 temp <- filtered.cleanput[ which(filtered.cleanput$species == s ),
                   c("gbifID","occurrenceID","species", "occurrenceStatus", "decimallongitude","decimallatitude","coordinateUncertaintyInMeters",
                     "depth", "depthAccuracy","eventDate")]
 #fix names of deciomal coordinates. capital L
 names(temp)<- c("gbifID","occurrenceID","species", "occurrenceStatus", "decimalLongitude","decimalLatitude","coordinateUncertaintyInMeters",
                 "depth", "depthAccuracy","eventDate")
#head(temp)
write.csv(temp,file = paste(speciespath,s,".csv", sep =""), row.names = F)
print(paste("wrote: ", speciespath,s,".csv", sep =""))
print(paste("no positives: ",length(which(temp$occurrenceStatus == "PRESENT")), sep =""))
print(paste("no negatives: ",length(which(temp$occurrenceStatus == "ABSENT")), sep =""))

}


# Sample a rest of observations to be used as pseudoabsences
locationsamples <- sample(1:length(filtered.cleanput$gbifID), 10000, replace =F)
pseudoabsences <- filtered.cleanput[ locationsamples,
                                     c("gbifID","occurrenceID","species", "occurrenceStatus", "decimallongitude","decimallatitude","coordinateUncertaintyInMeters",
                                       "depth", "depthAccuracy","eventDate")]
pseudoabsences$gbifID <- paste("pseudo",seq(1:10000), sep ="")
pseudoabsences$species <- NA
pseudoabsences$occurrenceStatus <- "ABSENT"
pseudoabsences$coordinateUncertaintyInMeters <- NA
pseudoabsences$depthAccuracy <- NA
pseudoabsences$eventDate <- NA
head(pseudoabsences)
names(pseudoabsences)<- c("gbifID","occurrenceID","species", "occurrenceStatus", "decimalLongitude","decimalLatitude","coordinateUncertaintyInMeters",
                "depth", "depthAccuracy","eventDate")

xlim <- c(-180,180)
ylim <- c(-60, 84)
plot(world1, xlim = xlim, ylim = ylim, col = "light grey")
points(pseudoabsences$"decimallongitude",pseudoabsences$"decimallatitude", col = "red", pch = "*", cex = 0.2)

# Write file with pseudoabsence data
write.csv(pseudoabsences,file = paste(speciespath,"/pseudoabsences",".csv", sep =""), row.names = F)
print(paste("wrote: ", speciespath,"/pseudoabsences",".csv", sep =""))

################## make plots to illustrate the distribution of filtered data ########

filtered.cleanput <- filtered.cleanput[ ,
                           c("gbifID","occurrenceID","species", "occurrenceStatus", "decimallongitude","decimallatitude","coordinateUncertaintyInMeters",
                             "depth", "depthAccuracy","eventDate")]
#fix names of deciomal coordinates. capital L
names(filtered.cleanput)<- c("gbifID","occurrenceID","species", "occurrenceStatus", "decimalLongitude","decimalLatitude","coordinateUncertaintyInMeters",
                "depth", "depthAccuracy","eventDate")


world1 <- readOGR(dsn=ogr_path,layer="CNTR_RG_01M_2020_4326")

factor <- as.numeric(as.factor(filtered.cleanput$species))
fac1 <- ceiling(factor/10)
fac2 <- factor - 10*(fac1-1)
xlim <- c(-110,40)
ylim <- c(-2, 68)
# High res for inspection
#jpeg("testplot.jpg", width = 90*( xlim[2] - xlim[1]),height = 90*( ylim[2] - ylim[1]), pointsize = 10)
# Moderate res for visualization
SMALLER_world = paste(jpgs, "testplot.CLEANPUT.SMALLER.world.jpg", sep = "/")
jpeg(SMALLER_world, width = 10*( xlim[2] - xlim[1]),height = 10*( ylim[2] - ylim[1]), pointsize = 4)
plot(world1, xlim = xlim, ylim = ylim, col = "light grey")
points(filtered.cleanput$"decimalLongitude",filtered.cleanput$"decimalLatitude", col = fac2, pch = fac1)
dev.off()