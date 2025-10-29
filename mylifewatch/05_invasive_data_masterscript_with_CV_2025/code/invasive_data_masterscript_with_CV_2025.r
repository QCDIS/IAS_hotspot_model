require(raster)
#require(rgdal)
require(zip)
require(sf)
require(sp)
require(C50)
require(unix)
library(tools)
library(C50)
library(raster)

########################################################################################
### Ecological niche modelling for invasive species
### Gunnar Andersson, National Veterinary Institute, SVA, Sweden  gunnar.andersson@sva.se
### Prepared fro SE-analytics Mars 2020
########################################################################################

suffixes <- c(".marine.new.base")#c(".marine")
suffix <- suffixes[1]

inputs_path = "/mnt/inputs/"

# Paths to folders where shapefiels are stored describing ICES areas. Used as overlays for plots only.
ICESpath <- paste(inputs_path,"ICES_areas", sep="")
download_zip_data_if_not_present_and_unzip(
    data_path = ICESpath,
    data_url = args$ICESpath_url,
    dest_path = inputs_path
    )
ICESecopath <- paste(inputs_path,"ICES_ecoregions", sep="")
download_zip_data_if_not_present_and_unzip(
    data_path = ICESecopath,
    data_url = args$ICESecopath_url,
    dest_path = inputs_path
    )

#Folder where the rasterstacks are stored
baseline_path <- "baselinedec100"
filled_rasterstacks <- paste(inputs_path,"filled_rasterstacks", sep="")
arranged_rasterstacks <- paste(inputs_path,"/arranged_rasterstacks", sep="") # For the layernames.rda"
arranged_baseline_path <- paste(arranged_rasterstacks,"/",baseline_path, sep="")
filtered_baseline_path <- paste(filled_rasterstacks,"/",baseline_path, sep="")

species_path = paste(inputs_path, "species", sep="")
download_zip_data_if_not_present_and_unzip(
    data_path = species_path,
    data_url = args$species_url,
    dest_path = inputs_path
    )

data_table_path <- paste(species_path,"/data_table.csv", sep="")



traffic_path <- paste(inputs_path, "traffic_layers", sep="")
download_zip_data_if_not_present_and_unzip(
    data_path = traffic_path,
    data_url = args$ais_shipping_density_url,
    dest_path = traffic_path
    )

stringsAsFactors = F
#Folder with the original rasterdata

outputs_path = "/mnt/outputs/"
invasive_species_path <- paste(outputs_path, "invasive_species/", sep="")
if (!dir.exists(invasive_species_path)) dir.create(invasive_species_path, recursive = TRUE)

#This folder contains files describing the cross validation scheme for each species

iterations_path <- paste(outputs_path, "Iterations/", sep="")
if (!dir.exists(iterations_path)) dir.create(iterations_path, recursive = TRUE)

#ROC_path <- paste(path,"Data2022/ROC.no.chlora", sep ="/")
ROC_path <- paste(outputs_path, "ROC/", sep ="/")
if (!dir.exists(ROC_path)) dir.create(ROC_path, recursive = TRUE)

# where plots should be stored
Plotpath <- paste(outputs_path, "Plots/", sep="")
if (!dir.exists(Plotpath)) dir.create(Plotpath, recursive = TRUE)

#The "outpath" contains the .csv files with species data and extracted environmental variables
Outpath <- paste( outputs_path, "RF_indata/", sep="")
if (!dir.exists(Outpath)) dir.create(Outpath, recursive = TRUE)

# the modelpath is where the models from random forest is saved. The files will also contain predictions from the cross validation,
# used as input when ROC curves are prpares
Modelpath <- paste(outputs_path, "Models/", sep="")
if (!dir.exists(Modelpath)) dir.create(Modelpath, recursive = TRUE)

resultpath <- paste(outputs_path, "Results/", sep="")
if (!dir.exists(resultpath)) dir.create(resultpath, recursive = TRUE)


# This contains the predicted species distribuition maps as rda files. These are used input for the plots.
mappath <- paste(outputs_path, "Maps/", sep="")
if (!dir.exists(mappath)) dir.create(mappath, recursive = TRUE)

# stores the prediction resutls as geotiff rasters. Both data for individual species and for the different measures
rastermappath <- paste(outputs_path, "Rastermaps/" , sep="")
if (!dir.exists(rastermappath)) dir.create(rastermappath, recursive = TRUE)


exclude_species = paste(outputs_path,"exclude_species",suffix,".rda",sep="")
stats_csv = paste(outputs_path,"exclude_species",suffix,".csv",sep="")

all_AUC_path = paste(outputs_path,"all.AUC",suffix,".csv", sep="")

# read in the region map used in the plots, and transfor to WGS84 - World Geodetic System 1984
# actually the map will not be needed until later. Readin git could be omitted here to save memory.

# Read the
Data.table <- read.csv(data_table_path,header=TRUE, sep=",",stringsAsFactors = F)



which(duplicated(Data.table))

names(Data.table)[1] <- "species" # just to check

# select one species to try out the code. Not used in the loop.
species = Data.table$species[2] # Ficopomatus enigmaticus  "Neogobius melanostomus"
# Define which stack to used when extracting environmental data-
# not using alternative rasterstacks

biooracle_filled_layers = paste(filtered_baseline_path,"/Biooracle.filled.layers.global2025",".tif", sep="")
Stack <- stack(biooracle_filled_layers) #"globalStack.rda"#"globalStack.rda" or "europeStack.rda"
print(paste("Loading raster stack from:", biooracle_filled_layers))

layernames_path <- paste(inputs_path,"/layernames",".rda" ,sep ="")
download_zip_data_if_not_present_and_unzip
(
    data_path = layernames_path,
    data_url = args$layernames_url,
    dest_path = inputs_path
    )

load(layernames_path)
print(paste("Loaded layer names from:", layernames_path))
names(Stack) <- layernames

#Prevent that strings in csv files are red in as factors variables
stringsAsFactors= FALSE

stats <- c()
species.stats = paste(outputs_path,"species.stats.rda", sep="")
print(paste("Data.table$species: ",Data.table$species))
for(i in seq_along(Data.table$species[-1])) {
    species <- Data.table$species[-1][i]
    print(paste("species: ", species))
    species_indata_file =  paste(Outpath,species,"_indata.csv", sep="")
    if (file.exists(species_indata_file)) {
        print(paste("File already exists:", species_indata_file))
        next
    }
    # Read in present absent and pseudoabsent points,
    #convert these points to spatial coordinates and extract environmental variables from rasterstack
    if (!file.exists(species_path)) {
        print(paste("Species path does not exist:", species_path))
    }
    if (length(list.files(species_path, pattern = paste0("^", species, ".*\\.csv$"))) == 0) {
        print(paste("No CSV file found for species:", species))
    }
    if (!dir.exists(Plotpath)){
        print(paste("Plotpath does not exist:", Plotpath))
    }
    species.data.list <-  read.and.extract(data.table = Data.table,
                                   species = species,
                                   stack = Stack,
                                   speciespath = species_path,
                                   plotpath = Plotpath,
                                   outpath = Outpath)
    species.data <- species.data.list[["complete.points"]]
    stats.temp <- c(species,species.data.list[["stats"]])
    stats <- rbind(stats, stats.temp)
    ## Check that there are no erroneous entries in the files.
    # the most common synonyms for "present" and "absent" are identified and entries harmonized.
    #      For entires where status cannot be determiend the line is removed
    # If data was checked in the earleir stage these lines should not find any mistakes
    my.data <- species.data
    print(paste(species,paste(unique(my.data$occurrenceStatus) ))  )
    my.data$occurrenceStatus <- as.character(my.data$occurrenceStatus)
    present.synonyms <-  which(my.data$occurrenceStatus == "present"|my.data$occurrenceStatus == "Present"| my.data$occurrenceStatus == "established"| my.data$occurrenceStatus == "Established")
    my.data$occurrenceStatus[present.synonyms] <- "present"
    absent.synonyms <-  which(my.data$occurrenceStatus == "Absent"| my.data$occurrenceStatus == "")

    my.data$occurrenceStatus[absent.synonyms] <- "absent"
    remove <- which(!my.data$occurrenceStatus == "absent"  & !my.data$occurrenceStatus == "present")
    print(paste("remove: ",remove))# promt if lines are removed
    if(length(remove > 0)){
        my.data <- my.data[-remove,]
    }
    species.data <- my.data
    # Save the species data with environmental variables in selected folder
    print(paste(species, "npos=", length(which(species.data$occurrenceStatus == "present")),
                                          "nneg=", length(which(species.data$occurrenceStatus == "absent") )))
    print(paste("Writing to file:", species_indata_file))
    write.csv(species.data, file = species_indata_file,row.names=F)
    # save the stats variable
    save(stats, file = species.stats)

}

if (file.exists(species.stats)) load(species.stats)
colnames(stats)<- c("species","npos","nabs", "npos.unique", "nabs.unique", "npos.unique.complete", "nabs.unique.complete")
write.csv2(as.data.frame(stats),file =  paste(outputs_path,"species.stats",suffix,".csv", sep=""))

##
stats <- as.data.frame(stats)
exclude <-  union(which(as.numeric(stats$npos.unique.complete) <6),
           which(as.numeric(stats$nabs.unique.complete) <6))
#krasch <- c(11,17,34)
#exclude <- union(exclude, krasch)
#stats[exclude,]
#stats <- read.csv2(file = "stats.mar 2025.csv", sep=",")
write.csv2(stats, file = stats_csv)
# save(exclude, file = exclude_species)


Data.table$species[exclude]
print(paste("Excluding species:", paste(Data.table$species[exclude], collapse = ", ")))
#load(file = "excludemar2025.csv")
# load(file= exclude_species)
################################################################
## prepare iterations
#################################################################
#
#
#Call the split data function to generate a list object that describes which entries are
# used as training data in each repetition and CV-fold in the cross-validation
# the resutls is stored as a rda file in the folder defined in "iterations path"

#The settings are made in teh "split data" function
#  n.repeat <- 5 this defines the number of times the cross validation process is performed
#  CV.level <- 5 this defienes the level of the cross validation
# to speed up thing it is possible to make  2-fold cross validation and just one repeat.
# when the data is split into ICES statistical rectangles based on their latitude and longitide.
#If data from other parts of the globe are used this may have to be modified
#
exclude2 <- c()

for(species in Data.table$species){
    lista.csv<- Sys.glob(paste(Outpath,"*.csv",sep="/"))
    files_to_read <- lista.csv[grep(species, lista.csv)]
    if (length(files_to_read) == 0) {
        next
    }
    nsites <- split.data(species = species,
                 indata.path = Outpath,
                 iterations.path= iterations_path)
    if(nsites[1]<5 | nsites[2]<5){
        exclude2<- rbind(exclude2, c(species, nsites[1],nsites[2]))
    }
}
exclude2.index <- match(exclude2[,1],Data.table$species)
exclude <- union(exclude, exclude2.index)
##################################################################
## Run MCMC  the function will prepare data and execute.
##################################################################
# note to self exclude3 <- c(8,20,23,27,42,44,46,47,50, 53, 57, 58,62, 64, 65,68, 69,70,71,72,73,74,75,76,77,78,79 ,80) #Mnemiopsis leidyi,Halothrix lumbricalis,Haloa japonicaCelleporaria brunnea,Apionsoma misakianum,
# Pleurosira laevis,Fenestrulina delicia Gonionemus vertens Smittoidea prolifica Corambe obscura Haminella solitaria Sinelobus vanhaareni Cephalothrix simula
#Torquigener flavimaculosus Boccardia proboscidea Paracerceis sculpta Oithona davisae Pseudodiaptomus marinus Polydora websteri Xenostrobus securis Aurelia solida,
#Evadne anonyx Parathalestris harpactoides Stylochus ellipticus Marenzelleria neglecta Marenzelleria arctia Aporrectodea caliginosa Polycerella emertoni
#intersect(exclude, exclude3) #exclude 3 == exclude...
#exclude3 <- c(8,20,23,27,42,44,46,47,50, 53, 57, 58,62, 64, 65,68, 69,70,71,72,73,74,75,76,77,78,79 ,80)
#exclude <- union(exclude, exclude3)

mcmc <- T
if(mcmc){
  require("rmcfs")
  #for j48 trees
  require("RWeka")
  require("parallel")

  print(paste("Loading layer names from:", layernames_path))
  lista.rda<- Sys.glob(paste(iterations_path,"*.rda",sep="/"))

    for(species in Data.table$species[-exclude]){
        selected_vars = paste(resultpath,"/selected.vars.",species,".rda",sep="")
        if (file.exists(selected_vars)) {
#             print(paste("Selected vars file already exists:", selected_vars))
            next
        }

    indata.path = Outpath
#     print(paste("Loading data from:", indata.path))
    lista.csv<- Sys.glob(paste(indata.path,"*.csv",sep="/"))
    if (length(lista.csv[grep(species,lista.csv)]) == 0) {
        next
    }
    my.data <- read.csv(lista.csv[grep(species,lista.csv)],header=T)
    print("Loaded my.data")
    my.data$RANDOMVAR <- runif(length(my.data$ID),0,1)
    my.data$RANDOMVAR2 <- runif(length(my.data$ID),0,1)
    my.data$RANDOMVAR3 <- runif(length(my.data$ID),0,1)

    print(paste(species, "n.abs=", length(which(my.data$occurrenceStatus == "absent")),
    ":n.pos=" ,length(which(my.data$occurrenceStatus == "present"))))
    unique(my.data$occurrenceStatus )
##    We are using occurrence iters in this project, not NUTS iters
     use.site.inters <- TRUE
    if(use.site.inters){
      load(lista.rda[
        intersect(grep("occurance.iters_",lista.rda),grep(species,lista.rda))
      ])
    }
     # get info from iterations data structure
     nrep <- length(all.occurance.iters)
     CV.level <- length(all.occurance.iters[[1]])
     process <- seq(1,nrep*CV.level)
     rep <- sort(rep(seq(1,nrep),CV.level))
     iter <- rep(seq(1,CV.level),nrep)
     process.plan <- cbind(process,rep,iter)
     print("Prepared process plan")
    # initiate MCMC
     ptm2 <- proc.time()
     no_cores <- detectCores() -1
     no_cores <- min(no_cores, length(iter))
     cl <- makeCluster(no_cores)
     #clusterExport(cl,"dummy.process")
     #result <- parLapply(cl,1:60,function(i) dummy.process(i))
     clusterExport(cl,c("MCMC.process", "my.data", "all.occurance.iters","process.plan"))
     #clusterEvalQ(cl,library(randomForest))
     clusterEvalQ(cl,library(rmcfs))

     # start parallell execution
     MCMCresult <- try(parLapply(cl,1:length(process.plan[,1]),
                             function(i) MCMC.process(.mydata =my.data,
                                                      .iters =all.occurance.iters,
                                                      .process.plan = process.plan,
                                                      .process.ID = i)))
     stopCluster(cl)
     time.to.complete2 <- proc.time()-ptm2
     print(time.to.complete2)

     if(length(grep("Error",MCMCresult))>0){
       print("Error - fallback")
       myreturn <- list()
       for(prid in 1:25){
           print(paste(" my.data len: ", length(my.data$ID)))
           print(paste(" all.occurance.iters len: ", length(all.occurance.iters)))
           print(paste(" process.plan len: ", length(process.plan[,1])))
           myreturn[[prid]] <- MCMC.process(.mydata = my.data,
                                            .iters = all.occurance.iters,
                                            .process.plan = process.plan,
                                            .process.ID = prid)
     }
     MCMCresult <- myreturn
     }
     # Check if MCMCresult is NULL or empty
        if (is.null(MCMCresult) || length(MCMCresult) == 0) {
#             print(paste("MCMCresult is NULL or empty for species:", species))
            next
        }

     selected_vars = paste(resultpath,"/selected.vars.",species,".rda",sep="")
     save(MCMCresult, file= selected_vars)
     rm(MCMCresult)
     gc()
  }
}

##################################################################
## plot selected variables
##################################################################
#Plotpath
#resultpath
#species <-Data.table$species[1]){
lista.rda<- Sys.glob(paste(iterations_path,"*.rda",sep="/"))
lista.selection <- Sys.glob(paste(resultpath,"*.rda",sep="/"))
indata.path = Outpath
lista.csv<- Sys.glob(paste(indata.path,"*.csv",sep="/"))

for(species in Data.table$species[-exclude]){
    #  for(species in Data.table$species[-exclude][-c(4,7,9,10,12,13)]){
    if (length(lista.selection[intersect(grep(species,lista.selection),
                                       grep("selected.vars",lista.selection))]) == 0) {
        print(paste("No selection rda file found for species:", species))
        print("Skipping to next species")
        next
    }
    print(paste("Plotting selected variables for species:", species))
    my.data <- read.csv(lista.csv[grep(species,lista.csv)],header=T, sep=",")
    absent_length = length(which(my.data$occurrenceStatus == "absent"))
    present_length = length(which(my.data$occurrenceStatus == "present"))


    selection_indices <- intersect(
      grep(species, lista.selection),
      grep("selected.vars", lista.selection)
    )
    selection_file <- lista.selection[selection_indices]
    load(selection_file)

    all.attributes <- sort(MCMCresult[[1]]$RI$attribute.name)

    cutoffs <- mean(sapply(1:length(MCMCresult), function(x) MCMCresult[[x]]$cutoff))
    RI_cutoff <- mean( sapply(1:length(MCMCresult), function(x) MCMCresult[[x]]$RI$RI[  as.numeric(MCMCresult[[x]]$cutoff)] ) )

  selection <- as.data.frame(all.attributes)
  selection<- cbind(selection,
                    sapply(1:length(MCMCresult),function(i)
                      as.numeric(MCMCresult[[i]]$RI[ order(MCMCresult[[i]]$RI$attribute.name),c("RI")]
                      )))

  selection <- as.data.frame(selection)
  selection[,-1]<- sapply(2:(length(MCMCresult)+1),function(i)
    as.numeric(selection[,i]))
  RI.mean <- rowMeans(selection[,-1])

  RI.min <- sapply(1:length(selection[,1]),function(i)
    min(selection[i,-1]))
  RI.max <- sapply(1:length(selection[,1]),function(i)
    max(selection[i,-1]))

  summary <- as.data.frame(all.attributes)
  summary <- as.data.frame(cbind(summary,RI.mean))
  summary <- cbind(summary,RI.min)
  summary <- cbind(summary,RI.max)
  summary.order <- summary[rev(order(RI.mean)),]

  raw.labels=summary.order$all.attributes#[1:24]

  nrow <- ceiling(length(raw.labels)/6)

  #plot pdf selection
  plotname <- paste(Plotpath,"/","_",species,"_seleced_variables",".png",sep="")

  png(plotname)  #to make file
  par(mar=c(14,2,2,2))
  #plot(seq(1:7),summary.order$RI.mean[seq(1:7)],ylim=c(0,max(summary.order$RI.max)),
   #    axes=F,  xlab = "")
  plot(seq(1:length(summary.order$RI.mean)),summary.order$RI.mean,ylim=c(-0.1*(max(summary.order$RI.max)),max(summary.order$RI.max)),
           axes=F,  xlab = "")
  box()
  axis(side=1, at=seq(1:length(summary.order$RI.mean)),
       labels=summary.order$all.attributes[1:length(summary.order$RI.mean)],las=2,cex.axis=0.7)
  axis(side = 2)
  for(i in 1:length(summary.order$RI.mean)){
    lines(c(i,i),c(summary.order$RI.min[i],summary.order$RI.max[i]))
  }
  lines( c(0,length(summary.order$RI.mean)),rep(as.numeric(RI_cutoff),2), col=2)
  title( species)
  dev.off()# close plotfile

  # plot pdf var distribution
  pdf(paste(Plotpath,"/var_distribution",species,".pdf", sep=""), paper = "a4r")# height= 11.7 , width= 16.6)#paper = "a4r"

  random <- grep("RANDOM",raw.labels )
  if(length(random) >0){ raw.labels <- raw.labels[-random]}
  par(mfrow= c(nrow,3))
  par(mar=c(4,2,2,0))
  for(l in 1:length(raw.labels)){
    my.var <- as.character(raw.labels[l])

    predictor <- my.data[,my.var]
    # my.data[,"bio_1"][which(is.na(my.data[,"bio_1"]))]
    response <- as.factor(my.data$occurrenceStatus)

    #hist(predictor)
    unique(response)
   #temptab <- cbind(predictor,response)# old version, gives wrong order of classes
    temptab <- cbind(predictor,ifelse(response== "present",1,2))# new puts "present" as class 1
    colnames(temptab)[2]<- "response"# add colname
    temptab[,"response"]

    temptab <- as.data.frame(temptab)
    temptab[,"response"] <- gsub(2,0,  temptab[,"response"])
    #head(temptab)

    temptab <- temptab[order(temptab$predictor),]
    nsplit <- 5
    chunks <- cut(seq(1:length(temptab[,1])),nsplit, labels = F)
    summary(chunks)
    table <- c()
    for(i in 1: nsplit){
      meanp <- mean(as.numeric(temptab$predictor[which(chunks == i)]))
      maxp <- max(as.numeric(temptab$predictor[which(chunks == i)]))
      minp <- min(as.numeric(temptab$predictor[which(chunks == i)]))
      meanr <- mean(as.numeric(temptab$response[which(chunks == i)]) )
      table <- rbind(table,c(minp, maxp,meanp,meanr))
    }
    colnames(table)<- c("minx","maxx","meanx", "meany")
    ymax <- ceiling(20* max(table[,"meany"]))*5

    xbreaks <- c(table[,"minx"], max(table[,"maxx"])  )
    betweenbreaks <- unlist( sapply(1:(length(xbreaks)-1), function(b) mean(xbreaks[b:(b+1)])) )
    breaklabels <-   unlist( sapply(1:(length(xbreaks)-1), function(b) paste(round(xbreaks[b:(b+1)],1), collapse = "-"   )  ) )

    plot(table[,"meanx"],table[,"meany"]*100, xlim= c(min(table[,"minx"]), max(table[,"maxx"]) ), ylim=c(0,ymax),
         xlab= "", ylab="presence (%)", axes = F)
    title(main= paste(species,my.var), cex.main = 0.7)
    box()
    axis(side = 1, at = xbreaks, labels = FALSE, las=0, cex=0.7)
    axis(side = 2, at = seq(0,ymax, length.out = 6), labels = seq(0,ymax, length.out = 6))
    # axis(side = 2)
    #  text(x=xbreaks, y=par()$usr[3]-0.05*(par()$usr[4]-par()$usr[3]),
    #      labels=xbreaks, srt=90, adj=1, xpd=TRUE)
    for(i in 1:length(xbreaks)){
      lines(rep(xbreaks[i],2),c(0, ymax), lty=2, col="grey")
    }
    text(x=betweenbreaks, y=par()$usr[3]-0.02*(par()$usr[4]-par()$usr[3]),
         labels=breaklabels, srt=45, adj=1, xpd=TRUE, cex=0.8)
    #text(x=c(min(table[,"minx"]), max(table[,"maxx"]) ), y=par()$usr[3]-0.02*(par()$usr[4]-par()$usr[3]),
    #     labels=c(min(table[,"minx"]), max(table[,"maxx"]) ), srt=45, adj=1, xpd=TRUE, cex=0.7)

    lines(table[,"meanx"],table[,"meany"]*100, lty=2)
  }
  dev.off()
}

print("Completed plotting selected variables")
##################################################################
  ## Run random forests. the function will prepare data and execute.
##################################################################
rm(Stack)
rm(map)
rm(my.map)
rm(my.maps)
rm(mystack)
rm(my.data)
gc()
#rlimit_all()
#rlimit_as(1e12)  #increases to ~12GB
#rlimit_stack(100000000)

# for(species in Data.table$species){
 for(species in Data.table$species[-exclude]){
    random_forest_model = paste(Modelpath,"/RF.model.and.predictions.eur.wt.",species,".rda",sep="")
    cross_validation_file = paste(Modelpath,"/RF.model.and.predictions.CV.eur.wt.",species,".rda",sep="")
    if (file.exists(random_forest_model) & file.exists(cross_validation_file)) {
        next
    }
    file_to_load <- paste(indata.path, species, "_indata.csv", sep="")
    if (!file.exists(file_to_load)) {
        print(paste("file_to_load does not exist:", file_to_load))
        next
    }
    rf.output.list <- run.random.forests(species = species ,
                selvar = "all",
                indata.path = Outpath,
                iterations.path = iterations_path)

    # run a garbage collection to free memory
    gc()
    # the rf.output.cv object contains the results from the cross validations analysis.
    #Predictions for each entrie are stored for each repeat and cv-fold. The random forest model is not saved, to save memory
    rf.output.cv <- rf.output.list[["RF.results.CV"]]
    # this object contains the random forests model obtainied when all indata is used as trainingset.
    # the predicted probability of being present is stored for each data poing
    rf.output <- rf.output.list[["RF.results.alldata"]]

    #print results to see progression
    names(rf.output)
    print(rf.output$"RF.selected")
    # print(rf.output$"response.sel")

    #save cross validation results for later calculation of ROC

    save(rf.output.cv, file= cross_validation_file)
    #save the random forest model for late prediction of map
    save(rf.output, file= random_forest_model)
    #remove the large object generated to make space for the next species
    rm(rf.output.list)
    rm(rf.output.cv)
    rm(rf.output)
    gc()
}
# ####################
# try C50 rules
print("Generating C5.0 rules")
#for(species in Data.table$species){
for(species in Data.table$species[-exclude]){
    lista.csv<- Sys.glob(paste(species_path,"*.csv",sep="/"))
    if (length(lista.csv[grep(species,lista.csv)]) == 0) {
        print(paste("No CSV file found for species:", species))
        next
    }
    csv_file <- lista.csv[grep(species,lista.csv)]
    if (!file.exists(csv_file))  {
        print(paste("File does not exist:", csv_file))
        next
    }
    print(paste("Processing species:", species, "from file:", csv_file))
    my.data <- read.csv(csv_file,header=T)
    print(paste("Loaded my.data for species:", species))
    column_names = names(my.data)
    # Validate predictor columns and skip when none or all-NA
    vars_max_index <- min(length(column_names), 23)
    vars_min_index <- 5
    vars <- names(my.data)[vars_min_index:vars_max_index]

    # keep only non-NA names that actually exist in my.data
    valid_vars <- vars[!is.na(vars) & vars %in% names(my.data)]

    if (length(valid_vars) == 0) {
      warning(paste("No valid predictor columns (5:23) found for species:", species, "- skipping C5.0"))
      next
    }

    # subset safely and preserve data.frame structure
    x_var <- my.data[, valid_vars, drop = FALSE]

    # skip if all predictor values are NA
    if (all(sapply(x_var, function(col) all(is.na(col))))) {
      warning(paste("Predictor columns are all NA for species:", species, "- skipping C5.0"))
      next
    }

    # proceed to model (example)
    rule_mod <- tryCatch(
      C50::C5.0(x = x_var, y = as.factor(my.data$occurrenceStatus), rules = TRUE),
      error = function(e) { warning(paste("C5.0 failed for", species, ":", conditionMessage(e))); return(NULL) }
    )
    if (is.null(rule_mod)) next

    print(paste("Writing C5.0 rules to file:", c50_rules_file))
    sink(c50_rules_file)
    print(species)
    print("/n")
    print(summary(rule_mod))
    sink()
}
# ##################################################################
# # calculate ROC curves and plot
# ##################################################################
all.AUC <- c()
print("Calculating and plotting ROC curves")
for(species in Data.table$species[-exclude]){

    plotname <- paste(ROC_path,"/","plotROC_",species,"allvars",".png",sep="")
    if (file.exists(plotname)) {
#         print(paste("ROC plot already exists for species:", species))
        next
    }

    cross_validation_file = paste(Modelpath,"/RF.model.and.predictions.CV.eur.wt.",species,".rda",sep="")
    if (!file.exists(cross_validation_file)) {
#       print(paste("No cross validation file:", cross_validation_file))
      next
    }
    print(paste("Processing species:", species))
    load(cross_validation_file)
    load(paste(iterations_path,"/occurance.iters_",species,".rda",sep=""))

    lista.csv<- Sys.glob(paste(Outpath,"*.csv",sep="/"))
    my.data <- read.csv(lista.csv[grep(species,lista.csv)],header=T)

    .true.class <-my.data$occurrenceStatus

    #######recreate the experiment setup
    nrep <- length(all.occurance.iters)
    CV.level <- length(all.occurance.iters[[1]])
    process <- seq(1,nrep*CV.level)
    rep <- sort(rep(seq(1,nrep),CV.level))
    iter <- rep(seq(1,CV.level),nrep)
    process.plan <- cbind(process,rep,iter)# may be obsolete when not running parallell script

    ###
    all.rep <- unique(process.plan[,2])
    curr.rep <- 2
    i <-1
    method <- "allvars"
    #### calculate one ROC curve per repetitio of the Cross validation
    all.ROC <- lapply(all.rep, function(curr.rep)
     calc.ROC(   rf.output.cv[[curr.rep]],.true.class)
    )

    # prepared a average ROC curve based on all repetitions togeather
    mean.ROC <- calc.ROC(lapply(1:length(process.plan[,1]), function(i)
      rf.output.cv[[   process.plan[i,2]  ]][[ process.plan[i,3]  ]] ) ,.true.class)

    # .RF.result <-lapply(1:length(process.plan[,1]), function(i)  rf.output.cv[[   process.plan[i,2]  ]][[ process.plan[i,3]  ]] )
    #      dim( rf.output.cv[[curr.rep]][[i]]$predict.selected )
    #rf.output.cv[[curr.rep]][[i]]$predict.selected[1,]

    #Plot the ROC cuve in the dedicated folder
    mean.AUC <-  plot.ROC(.ROC.path = ROC_path,
                        .my.species = species,
                        .mean.ROC = mean.ROC,
                        .all.ROC = all.ROC)
    #remove objects t make list for next run
    rm(all.occurance.iters)
    rm(rf.output.cv)
    rm(all.ROC)
    rm(mean.ROC)
    gc()
    all.AUC <- rbind(all.AUC, c(species,mean.AUC))
}
write.csv2(as.data.frame(all.AUC),file = all_AUC_path )
##################################################################
## predict maps
##################################################################

dataset_scenarios <- list.dirs(
  paste(rasterstacks_path, sep = ""),
  full.names = FALSE,
  recursive = FALSE
)

# --- Scenario loop ---
for (sel.sen in 1:length(dataset_scenarios)) {

    scenario <- dataset_scenarios[[sel.sen]]
        rastermappath_scenario <- paste(rastermappath,scenario,sep ="")
        if (!dir.exists(rastermappath_scenario)) dir.create(rastermappath_scenario)

        Biooracle_scenario_path <- paste(rasterstacks_path,"/",scenario,"/",sep="")
        Biooracle.filled.layers.global = paste(Biooracle_scenario_path,"Biooracle.filled.layers.global2025",".tif", sep="")
        print(paste("Loading raster stack from:", Biooracle.filled.layers.global))
        if (!file.exists(Biooracle.filled.layers.global)) {
            print(paste("No raster stack file found:", Biooracle.filled.layers.global))
            next
        }
        Stack <- stack(Biooracle.filled.layers.global)
        rasterstacks_path.base <- paste(arranged_rasterstacks,"/baselinedec50", sep="/")
        layernames_path <- paste(rasterstacks_path.base,"/layernames",".rda" ,sep ="")
        if (!file.exists(layernames_path)) {
            print(paste("No layer names file found:", layernames_path))
            next
        }
        load(layernames_path)
        # Check lengths before assignment
        if (length(layernames) == nlayers(Stack)) {
          names(Stack) <- layernames
          print(length(layernames))
          print(nlayers(Stack))
        } else {
          print("Length of layernames does not match number of layers in Stack. Names not assigned.")
          print(length(layernames))
          print(nlayers(Stack))
          next
        }
        bar <- Stack
        #call the function to predict the species distribution at raster level.
        #The function return the prediction as a raster object but also make plots as .png
        # the lines between png() and dev.off() may be removed/inactivated if the png plots are not wanted.
        exclude_species = Data.table$species[-exclude]
        print(paste("exclude_species: ", exclude_species))
        print(paste("Number of species to process:", length(Data.table$species[-exclude])))
        if (length(Data.table$species[-exclude]) == 0) {
            print("No species to process for map prediction.")
            next
        }
        for(Species in Data.table$species[-exclude]){
            newfile <- paste(rastermappath_scenario,"/linear.prob.global.",Species,'.tif',sep="")
            if (file.exists(newfile)) {
                print(paste("Prediction map already exists:", newfile))
                next
            }
            model_file = paste(Modelpath,"RF.model.and.predictions.eur.wt.",species,".rda",sep="")
            print(paste("Loading model from:", model_file))
            if (!file.exists(model_file)) {
                print(paste("No random forest model file found for species:", Species))
                next
            }
            print(paste("Predicting map for species:", Species))
            map <-  predict.maps(species = Species,
                    modelpath = model_file
                    )
            if (is.null(map)) {
                print(paste("Prediction map is NULL for species:", Species))
                next
            }
            # save predicted probabilities as rasterfiles.
            print(paste("Writing raster to:", newfile))
            writeRaster(map, filename= newfile, format = "GTiff", overwrite=TRUE)
        }
}#end scenario

######################################################################################################################################
#####################################################################################################################################
## ###################################################################plot maps
##################################################################
######################################################################################################################################

#restore path
#rastermappath <- paste(path,"/data/Rastermaps", suffix,sep ="")



#colorsBrBG <- rev(divPalette(n=12, name = c( "BrBG") ) )
colorsBrBG2 <- rev(divPalette(n=19, name = c( "BrBG") ) )

plot(seq(1:25), seq(1:25), col=colorsBrBG2, pch=16)
colorsBr <- c("white",colorsBrBG2[11:19])
colorsBlue <- seqPalette(n=,12, name = c( "Blues") )
#shape.ecoregions <- readOGR(dsn=ICESecopath,
#                            layer="ICES_ecoregions_20171207_erase_ESRI",encoding="UTF-8")
shape.ecoregions <- st_read(paste(ICESecopath,"/ICES_ecoregions_20171207_erase_ESRI.shp",sep=""), )
print("Loaded ecoregions shapefile")
#plot(shape.ecoregions$geometry)

brk <- c(seq(0, 1,by=0.1),1.05)
# temporary solution
shape2 <-shape.ecoregions
xlim=c(0,30)# for swe
ylim = c(50,70)#for swe
#xlim=c(-180,180)
#ylim=c(-90,90)

if(!dir.exists(Plotpath)){dir.create(Plotpath)}

probability.all.species.file = paste(Plotpath,"mean.probability.all.species.png",sep="")
print(paste("Looking Mean probability plot:", probability.all.species.file))
if (file.exists(probability.all.species.file)) {
    print(paste("Mean probability plot already exists:", probability.all.species.file))
} else {
    rastermappath_scenario <- paste(rastermappath,dataset_scenarios[2],sep ="")
    lista.ras <- Sys.glob(paste(rastermappath_scenario,"*linear.prob.global*",sep="/"))
    print(paste("Found", length(lista.ras), "raster maps in", rastermappath_scenario))

    predvar<-stack(lista.ras)
    png(filename = probability.all.species.file, width = 800, height = 600)
    plot(mean(predvar))
    dev.off()
    rm(predvar)
    gc()
}

for(species in Data.table$species[-exclude]){
    plot_file_Swe = paste(Plotpath,"/",Species, ".lin.trainpoints.Swe.png",sep="")
    plot_file_world = paste(Plotpath,"/",Species, ".lin.trainpoints.world.png",sep="")
    if (file.exists(plot_file_Swe) & file.exists(plot_file_world)) {
        next
    }
    print(paste("Plotting map for species:", species))
    rastermap <- lista.ras[grep(species,lista.ras)]
    plot.maps(.Species=Species,
            indata.path = Outpath,
            .rastermap=rastermap,
            .plotpath=Plotpath,
            .colors=colorsBr,
            .brk=brk,
            .shape=shape2,
            xlim=xlim,
            ylim= ylim
            )
}
######################################################################################################################################
###################################################plot futire scenarios#########################################################
Plotpath2 <- gsub( ".marine","ssp",Plotpath)
if(!dir.exists(Plotpath2)){dir.create(Plotpath2)}

listatiff <- c()
# scenarios <- c( "baseline" ,"ssp119" ,  "ssp126"  , "ssp245"   ,"ssp370" ,  "ssp460"  , "ssp585"  )
# dec.vec <- c("", "dec50", "dec100")
for (sel.sen in 1:length(dataset_scenarios)) {
    scenario <- dataset_scenarios[[sel.sen]]
    rastermappath_scenario <- paste(rastermappath,scenario,sep ="")
    if(dir.exists(rastermappath_scenario)){
        lista.tif.temp<- Sys.glob(paste(rastermappath_scenario,"/","*.tif",sep=""))
        listatiff <- c(listatiff,lista.tif.temp)
    }else{
        print(paste("no dir",rastermappath_scenario))
    }
}

for (sel.sen in 1:length(dataset_scenarios)) {
    scenario <- dataset_scenarios[[sel.sen]]
    for(species in Data.table$species[-exclude]){
        plot_file = paste(Plotpath2,species, ".",scenario,".png",sep="")
        if (file.exists(plot_file)) {
            next
        }
        #get files for scenario
        cases <- union(grep("baseline",listatiff),grep(scenario,listatiff))
        print(paste("Found", length(cases), "files for scenario:", scenario))
        my.rasterfiles<- listatiff[ intersect(grep(species,listatiff), cases)  ]

        my.stack <- try(stack(my.rasterfiles))
        r <- 1

        my.layernames <- sapply(1:length(my.rasterfiles), function(r)
        strsplit(strsplit( my.rasterfiles[r],"Rastermaps")[[1]][2],"/")[[1]][1])
        names(my.stack)<- my.layernames
        if (nlayers(my.stack) >= 3) {
            my.change <- my.stack[[3]] - my.stack[[1]]
        }else if (nlayers(my.stack) >= 2) {
            my.change <- my.stack[[2]] - my.stack[[1]]
        } else {
            warning("Not enough layers in my.stack to compute change.")
            print("Skipping to next species.")
            next
        }
        if(length(grep("Error",my.stack[1]))<1 ){
          print(paste(names(my.stack), "stack success"))
          par(mfrow=c(3,1))
        }else{
          print(paste(names(my.stack), "stack fail"))
          break
        }
        colorsBrBG2 <- rev(divPalette(n=19, name = c( "BrBG") ) )
        colorsBrBG3 <- rev(divPalette(n=10, name = c( "BrBG") ) )


        colorsBr <- c("white",colorsBrBG2[11:19])
        colorsBlue <- seqPalette(n=,12, name = c( "Blues") )
        brk <- c(seq(0, 1,by=0.1),1.05)
        brk2 <- c(seq(-1, 1,by=0.1),1.05)

        png(plot_file,  width = 180, height = 180, units = "mm", res=1200)

        par(mfrow=c(2,2))
        par(mar=c(2,2,2,6))

        for(i in 1:nlayers(my.stack)){
          map<-my.stack[[i]]
          # plot(map)

          map2 <- map
          map2[is.na(map2)]<- 1.005
          my.colors <- c(colorsBr,"lightgrey")

          plot(map2, col=my.colors, breaks = brk)
          if(i == 1){
            title(main = species)}else{
          title(main = names(my.stack)[i])
            }
        }
      map<-my.change
      my.colors <- c("white",colorsBrBG2,"lightgrey")
      map2 <- map
      map2[is.na(map2)]<- 1.005
      plot(map2, col=my.colors, breaks = brk2)
      title(main = "dec100 - base")
      dev.off()
  }
}# end sel.sen species
print("Finished plotting future scenario maps")
######################################################################################################################################
### ###################################################################Plot logmaps
#####################################################################################################################################
######################################################################################################################################
#species  <- Data.table$species[1]

colorsBrBG <- c(rev(divPalette(n=12, name = c( "BrBG") ) ),"lightgrey","lightgrey","lightgrey")
colorsBlue <- c(seqPalette(n=,12, name = c( "Blues") ),"lightgrey","lightgrey","lightgrey" )
colorsBrBG2 <- rev(divPalette(n=22, name = c( "BrBG") ) )

#plot(seq(1:25), seq(1:23), col=colorsBrBG2, pch=16)
colorsBr <- c("white",colorsBrBG2[12:22], "lightgrey")
brk.log <- c(seq(-3, 0,by=0.25),0.25)
length(brk.log)
# define the area to plot. In this case the coordinates fro teh swedish map
xlim=c(0,30)
ylim = c(50,70)


lista.ras <- Sys.glob(paste(rastermappath,"*linear.prob.*",sep="/"))
predvar<-stack(lista.ras)
print("Calculating mean probability of presence for all species")
for(species in Data.table$species){
    plot_file_eu = paste(Plotpath,"/",species, ".logprob.Eur.png",sep="")
    plot_file_sw = paste(Plotpath,"/",species, ".logprob.Swe.png",sep="")
    if (file.exists(plot_file_eu) & file.exists(plot_file_sw)) {
        next
    }
    rastermap <- lista.ras[grep(species,lista.ras)]
    if (length(rastermap) == 0) {
        next
    }
    map<-raster(rastermap)
    lista.csv<- Sys.glob(paste(Outpath,"*.csv",sep="/"))
    # Print all the filenames
    print(lista.csv[grep(species,lista.csv)])
    if (length(lista.csv[grep(species,lista.csv)]) == 0) {
        next
    }
    my.data <- read.csv(lista.csv[grep(species,lista.csv)],header=T)
    print(paste("Plotting log map for species:", species))

    # prepare a rasterlayer with the 10-log of the predicted probability of presence.
    # An arbitrary small number is added, as log10(0) is not defined
    logmap <- log10(map+0.001) #log10(max(map , 0.001, na.rm=T))

    # an arbitrary large number is inserted at locations where presence is undefined, in other words land.
    logmap[is.na(mean(bar))]<- 0.24

    # plot maps at European and Swedish scale this presence points added.
    # Plot without indicating xlim and ylim. This gives a plot area defined by the raster extent
    png(plot_file_eu,  width = 180, height = 180, units = "mm", res=1200)
    plot(logmap,col=colorsBr, breaks=brk.log)
    plot(shape2$geometry,add = TRUE, xlim=xlim, ylim=ylim, border = 1, lwd = 0.1)
    points(my.data$Lon[which(my.data$occurrenceStatus == "absent")],
         my.data$Lat[which(my.data$occurrenceStatus == "absent")], col=4, pch=1, cex=0.3, lwd=0.2)
    points(my.data$Lon[which(my.data$occurrenceStatus == "present")],
         my.data$Lat[which(my.data$occurrenceStatus == "present")], col=2, pch=1, cex=0.3, lwd=0.2)
    dev.off()
  ### same for Sweden
  # make the plot using the xlim and ylim defined earlier.
  png(plot_file_sw,  width = 180, height = 180, units = "mm", res=1200)
  plot(logmap,col=colorsBr, breaks=brk.log , xlim=xlim, ylim=ylim )
  plot(shape2$geometry,add = TRUE, xlim=xlim, ylim=ylim, border = 1, lwd = 0.1)
  points(my.data$Lon[which(my.data$occurrenceStatus == "absent")],
         my.data$Lat[which(my.data$occurrenceStatus == "absent")], col=4, pch=1, cex=1, lwd=0.25)
  points(my.data$Lon[which(my.data$occurrenceStatus == "present")],
         my.data$Lat[which(my.data$occurrenceStatus == "present")], col=2, pch=1, cex=1, lwd=0.25)
    dev.off()

    writelograster <- F
    if(writelograster){
        newfile <- paste(rastermappath,"/log.prob.",species,sep="")
        writeRaster(logmap, filename= newfile, format = "GTiff", suffix='.tif', overwrite=TRUE)
    }
  rm(logmap)
  rm(map)
  gc()
}
print("Finished plotting log maps for all species")
######################################################################################################################################

###################################################################
### #################################################Plot combined maps
##################################################################
######################################################################################################################################
# define color palettes. White and lightgrey to incidate values outside range.
colorsBrBG <- c(rev(divPalette(n=12, name = c( "BrBG") ) ),"lightgrey","lightgrey","lightgrey")
colorsBlue <- c(seqPalette(n=,12, name = c( "Blues") ),"lightgrey","lightgrey","lightgrey" )
colorsBrBG2 <- rev(divPalette(n=22, name = c( "BrBG") ) )

colorsBr <- c("white",colorsBrBG2[12:22], "lightgrey")
# the cutoffs to used in the color scale at log-scale
brk.log <- c(seq(-3, 0,by=0.25),0.25)

####### calculate average occurcence
#Find the predicted maps for the different species
lista.ras <- Sys.glob(paste(rastermappath,"*linear.prob.*",sep="/"))
predvar<-stack(lista.ras)
# calculate the mean probability of presence for all species, as a rasterlayer.
mean.predvar <- mean(predvar, na.rm=T)
plot(mean.predvar)
newfile <- paste(rastermappath,"mean.prob.all.species.tif",sep="")
print(paste("newfile: ",newfile))
logpredvar <- log10(mean(predvar, na.rm=T) + 0.001)
if (file.exists(newfile)) {
    print(paste("Mean probability raster already exists:", newfile))
} else {
    print(paste("Writing mean probability raster to:", newfile))
    writeRaster(mean.predvar, filename= newfile, format = "GTiff", suffix='.tif', overwrite=TRUE)
    writelograster <- F
    if(writelograster){
      newfile <- paste(rastermappath,"/log 10 of mean.prob.all.species",sep="")
      writeRaster(logpredvar, filename= newfile, format = "GTiff", suffix='.tif', overwrite=TRUE)
    }
}

############################ ############################ ############################
############################ plot mean with same colors as before ############################
############################ ############################ ############################
#Plot mean of lin
colorsBrBG2 <- rev(divPalette(n=19, name = c( "BrBG") ) )
plot(seq(1:25), seq(1:25), col=colorsBrBG2, pch=16)
colorsBr <- c("white",colorsBrBG2[11:19])

map2 <- mean.predvar
map2[is.na(map2)]<- 1.005
colorsBr <- c(colorsBr, "lightgrey","lightgrey")
brk <- c(seq(0, 1,by=0.1),1.05)

logpredvar[is.na(mean(bar))]<- 0.25
xlim=c(-180,180)
ylim = c(-180,180)
mean.all.eur = paste(Plotpath,"mean.all.eur", ".png",sep="")
if (file.exists(mean.all.eur)) {
    print(paste("Mean probability plot already exists:", mean.all.eur))
} else {
    png(mean.all.eur ,  width = 180, height = 180, units = "mm", res=1200)
    plot(map2,col=colorsBr, breaks=brk)
    plot(shape2$geometry,add = TRUE, xlim=xlim, ylim=ylim, border = 1, lwd = 0.1)
    dev.off()
}


xlim=c(0,30)
ylim = c(50,70)
mean.all.Swe = paste(Plotpath,"mean.all.Swe", ".png",sep="")
if (file.exists(mean.all.Swe)) {
    print(paste("Mean probability plot already exists:", mean.all.Swe))
} else {
    png(mean.all.Swe,  width = 180, height = 180, units = "mm", res=1200)
    plot(map2,col=colorsBr, breaks=brk,xlim=xlim, ylim=ylim)
    plot(shape2$geometry,add = TRUE, xlim=xlim, ylim=ylim, border = 1, lwd = 0.1)
    dev.off()
}
# ############################ ############################ ############################
# ############################ Plot mean of log ############################
# ############################ ############################ ############################
logpredvar[is.na(mean(bar))]<- 0.25
mean.all.eur = paste(Plotpath,"logmean.all.eur", ".png",sep="")
if (file.exists(mean.all.eur)) {
    print(paste("Mean log10 probability plot already exists:", mean.all.eur))
} else {
    png(mean.all.eur ,  width = 180, height = 180, units = "mm", res=1200)
    plot(logpredvar,col=colorsBr, breaks=brk.log)
    plot(shape2$geometry,add = TRUE, xlim=xlim, ylim=ylim, border = 1, lwd = 0.1)
    dev.off()
}
logmean.all.Swe = paste(Plotpath,"logmean.all.Swe", ".png",sep="")
if (file.exists(logmean.all.Swe)) {
    print(paste("Mean log10 probability plot already exists:", logmean.all.Swe))
} else {
    png(logmean.all.Swe,  width = 180, height = 180, units = "mm", res=1200)
    plot(logpredvar,col=colorsBr, breaks=brk.log,xlim=xlim, ylim=ylim)
    plot(shape2$geometry,add = TRUE, xlim=xlim, ylim=ylim, border = 1, lwd = 0.1)
    dev.off()
}

### mean of log10 predvar
custom.pred <- log10(predvar[[1]]+0.001)
for(i in 2:length(names(predvar))){
  custom.pred <- custom.pred+ log10(predvar[[i]]+0.001)
}
plot(custom.pred)
newfile <- paste(rastermappath,"sum of log10 prob all.species",sep="")
if (file.exists(newfile)) {
    print(paste("Sum of log10 probability raster already exists:", newfile))
} else {
    writeRaster(custom.pred, filename= newfile, format = "GTiff", suffix='.tif', overwrite=TRUE)
}
custom.pred <- custom.pred/length(names(predvar))
plot(custom.pred)
newfile <- paste(rastermappath,"mean of log10 prob all.species",sep="")
if (file.exists(newfile)) {
    print(paste("Mean of log10 probability raster already exists:", newfile))
} else {
    writeRaster(custom.pred, filename= newfile, format = "GTiff", suffix='.tif', overwrite=TRUE)
}


#######################################
# customzed combined layers
##################################

########## get filenames for traffic data and load as rasterstack


newfile <- paste(rastermappath,"/mean.prob.all.species.tif",sep="")
mean.predvar <- raster(newfile)

lista.ras2 <- Sys.glob(paste(traffic_path,"*.tif",sep="/"))
print(paste("Found", length(lista.ras2), "traffic raster maps in", traffic_path))
if (length(lista.ras2) == 0) {
    print(paste("No traffic raster files found in:", traffic_path))
} else {
    print(paste("Using first traffic raster file:", lista.ras2[1]))
}

plotvar <- mean.predvar
shipvar <- raster(lista.ras2)
proj4string(shipvar)# chec projection
#shipvar <- projectRaster(shipvar,crs=CRS("+init=epsg:4326"))
#GDAL Message 1: +init=epsg:XXXX syntax is deprecated. It might return a CRS with a non-EPSG compliant axis order.

shipvar <- projectRaster(shipvar,crs=4326)

proj4string(shipvar)# chec projection

# create two layers of same size
plotvar.crop<- crop(plotvar,shipvar)
shipvar<- crop(shipvar,plotvar.crop)
plotvar2 <- raster::resample(plotvar.crop,shipvar, method = "bilinear")

extent(plotvar2)
extent(shipvar)
dim(plotvar2)
dim(shipvar)

shipvar2 <- shipvar
shipvar2[shipvar>1000] <- 1000

# # create the combined data layer sing some mathematical expression
combivar <- plotvar2 + shipvar2/1000
#combivar[combivar > 3 ]<- 3
par(mfrow = c(1,1))
plot(combivar)

combivar <- plotvar2 + shipvar2/5000

colorsBlue <- seqPalette(n=,24, name = c( "Blues") )
colorsorange <- seqPalette(n=,24, name = c( "Oranges") )
colorsdiv <- divPalette(n=,24, name = c( "BrBG") )

############################## make plot as in old project   ############
newfile <- paste(Plotpath,"combiplot.png",sep="")
if (file.exists(newfile)) {
    print(paste("Combined plot already exists:", newfile))
} else {
    png(newfile,  width = 360, height = 360, units = "mm", res=1200)

    par(mfrow=c(2,2))
    par(mar=c(0,2,6,8))

    plot(plotvar2, col= colorsBlue)
    title(main="mean risk")
    plot(shape2$geometry, add=T, lwd=0.1, col=NA)

    plot(shipvar2, col= colorsorange)
    title(main="traffic")
    plot(shape2$geometry, add=T, lwd=0.1, col=NA)

    new.plotvar2 <- plotvar2
    new.plotvar2[new.plotvar2>1.5]<-1.5
    plot(new.plotvar2, col= colorsBlue)
    #arg <- list(at=c(0,0), labels=c("NA","NA"))
    plot(shipvar2, add=T, col=colorsorange,alpha=0.5, legend=FALSE)
    title(main="superimposed")
    plot(shape2$geometry, add=T, lwd=0.1, col=NA)


    plot(combivar, col= colorsBlue)
    title(main="added by function")
    plot(shape2$geometry, add=T, lwd=0.1, col=NA)


    dev.off()
}


#######################################################################
#                              END                                         #####
#?#################################################################################
### pca on rasters
#https://www.rdocumentation.org/packages/RStoolbox/versions/0.2.6/topics/rasterPCA
library(ggplot2)
library(reshape2)
library(RStoolbox)
#Just a test to run PCA on the rastetrlayers...
## Run PCA
#rasterPCA(img, nSamples = NULL, nComp = nlayers(img), spca = FALSE, maskCheck = TRUE, ...)
set.seed(25)
rpc <- rasterPCA(predvar,nComp = 3,spca= T)

summary(rpc$model)
print(paste("summary of PCA model"),summary(rpc$model))
loadings(rpc$model)
print(paste("loadings of PCA model"),loadings(rpc$model))


ggRGB(rpc$map,1,2,3, stretch="lin", q=0)
  plots <- lapply(1:3, function(x) ggR(rpc$map, x, geom_raster = TRUE))
 plot(plots[[1]],xlim=xlim, ylim=ylim)
 plot(plots[[2]],xlim=xlim, ylim=ylim)
 plot(plots[[3]],xlim=xlim, ylim=ylim)
####### obtain oub prediction performance
 for(species in Data.table$species){
   file <-  paste(Modelpath,"RF.model.and.predictions.eur.wt.",species,".rda",sep="")
   #print(file)
   load(file)
   print("########################")
   print("")
   print(species)
   print(rf.output$"RF.selected"$"confusion")

 }

print("All done")
#######################################################################