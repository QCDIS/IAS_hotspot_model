#invasive.data.masterscript.r
setwd("~/Dokument/Projekt/HAV2025")
source("scripts/SEanalytics.functions2025.r")
require(raster)
#require(rgdal)
require(zip)
require(sf)
require(sp)
require(C50)
require(unix)
########################################################################################
### Ecological niche modelling for invasive species
### Gunnar Andersson, National Veterinary Institute, SVA, Sweden  gunnar.andersson@sva.se
### Prepared fro SE-analytics Mars 2020
########################################################################################

#suffixes <- c("",".no.chlora",".no.chlora.delim")
suffixes <- c(".marine.new.base")#c(".marine")
suffix <- suffixes[1]

stringsAsFactors = F
#Folder with the original rasterdata

path <- "~/Dokument/Projekt/HAV2025"

#Rasterpaths <- c("C:/Users/gunnar.andersson/OneDrive - SVA/SE analytics/Data/2019-11-01_09-49-26/Present_Benthic_Avg/Present_Benthic_Avg",
#                 "C:/Users/gunnar.andersson/OneDrive - SVA/SE analytics/Data/2019-11-01_09-49-26/Present_Surface/Present_Surface")

#Folder where the rasterstacks are stored
Stackpath <- paste(path, c("data/Biooracle.download/rasterstacks/baseline"), sep="/")

# Folder containg presence and absence data for each species and preudoabsence data used for all species.
SpeciespathRaw <- paste(path,"data/Indata2025/speciesIndata2025", sep ="/")
#Speciespath <- paste(path,"Data2022/RF.indata.no.chlora", sep ="/")
Speciespath <- paste(path,"/data/RF.indata",suffix, sep ="")
if (!dir.exists(Speciespath)) dir.create(Speciespath) 

#This folder contains files describing the cross validation scheme for each species
#Iterations.path <- paste(path,"Data2022/Iterations.no.chlora", sep ="/")
Iterations.path <- paste(path,"/data/Iterations",suffix, sep ="")
if (!dir.exists(Iterations.path)) dir.create(Iterations.path) 

#ROC.path <- paste(path,"Data2022/ROC.no.chlora", sep ="/")
ROC.path <- paste(path,"/data/ROC",suffix, sep ="")
if (!dir.exists(ROC.path)) dir.create(ROC.path) 

# where plots should be stored
#Plotpath <- paste(path, "Plots.no.chlora" ,sep="/")
Plotpath <- paste(path, "/Plots" ,suffix,sep="")
if (!dir.exists(Plotpath)) dir.create(Plotpath) 

#The "outpath" contains the .csv files with species data and extracted environmental variables
Outpath <- paste(path,"/data/RF.indata",suffix, sep ="")
if (!dir.exists(Outpath)) dir.create(Outpath) 

# the modelpath is where the models from random forest is saved. The files will also contain predictions from the cross validation,
# used as input when ROC curves are prpares
Modelpath <-  paste(path,"/data/Models",suffix, sep ="")
if (!dir.exists(Modelpath)) dir.create(Modelpath) 

#### traffic patn may not be used here
trafficpath <- paste(path,"/data/traffic layers/2016_All_shiptypes_AIS_Shipping_Density", sep ="")
#Modelpath <- "C:/Users/gunnar.andersson/OneDrive - SVA/SE analytics/Data/models/models apr 2020"

resultpath <-  paste(path, "/data/Results" ,suffix,sep="")
if (!dir.exists(resultpath)) dir.create(resultpath) 

selectionspath <-  paste(path, "/data/Results",suffix ,sep="")
if (!dir.exists(selectionspath)) dir.create(selectionspath) 


# This contains the predicted species distribuition maps as rda files. These are used input for the plots.
mappath <- paste(path,"/data/Maps",suffix, sep ="")
if (!dir.exists(mappath)) dir.create(mappath) 

#mappath <- "C:/Users/gunnar.andersson/OneDrive - SVA/SE analytics/maps/maps apr 2020"

# stores the prediction resutls as geotiff rasters. Both data for individual species and for the different measures
rastermappath <- paste(path,"/data/Rastermaps", suffix,sep ="")
if (!dir.exists(rastermappath)) dir.create(rastermappath) 


#Paths to folders where shapefiels are stored describing ICES areas. Used as overlays for plots only.
ICESpath <- paste(path,"data/ICES_areas" ,sep="/")
ICESecopath <- paste(path,"data/ICES_ecoregions" ,sep="/")
#ICESstatrectpath <- "C:/Users/gunnar.andersson/OneDrive - SVA/SE analytics/Data/ICES_StatRec_mapto_ICES_Areas"

# read in the region map used in the plots, and transfor to WGS84 - World Geodetic System 1984
# actually the map will not be needed until later. Readin git could be omitted here to save memory.
#shape2 <- readOGR(dsn="H:/R gruppen/Kartor/NUTS_shapefile",layer="nutsByHand")
#proj4string(shape2)#GRS80 EPSG:4019
#shape2 <- spTransform(shape2, CRS("+init=epsg:4326"))
#plot(shape2)# just to check that it appears right
#########################################################

# Read the
#Data.table <- read.csv2("data.table.csv",header=TRUE, sep=";",stringsAsFactors = F)
#Data.table <- read.csv2("data.table.apr2020.csv",header=TRUE, sep=";",stringsAsFactors = F)

#Data.table <- read.csv(paste(path,"/Data2022/data.table.first6.2022.csv", sep=""),header=TRUE, sep=",",stringsAsFactors = F)
#Data.table <- read.csv(paste(path,"/Data2022/data.table.2022.csv", sep=""),header=TRUE, sep=",",stringsAsFactors = F)
#Data.table <- read.csv(paste(path,"/Data2022/data.table.extra.download.2022.csv", sep=""),header=TRUE, sep=",",stringsAsFactors = F)
#Data.table <- read.csv(paste(path,"/data/data.table.feb22.csv", sep=""),header=TRUE, sep=",",stringsAsFactors = F)

#Data.table <- read.csv(paste(path,"/data/data.table.mars.some.removed.2024.csv", sep=""),header=TRUE, sep=",",stringsAsFactors = F)
#Data.table <- read.csv(paste(path,"/data/data.table.mars.filtered.2024.csv", sep=""),header=TRUE, sep=",",stringsAsFactors = F)
Data.table <- read.csv(paste(path,"/data/species.data/data.table.apr2025.new.baseline.csv", sep=""),header=TRUE, sep=",",stringsAsFactors = F)



which(duplicated(Data.table))

head(Data.table)

names(Data.table)[1]<- "species" # just to check

# select one species to try out the code. Not used in the loop.
Species = Data.table$species[2] # Ficopomatus enigmaticus  "Neogobius melanostomus"
#Define which stack to used when extracting environmental data-
# not using alternative rasterstacks
#Stack <- stack(paste(Stackpath,"/Biooracle.global2024", suffix,".TIFF", sep="")) #"globalStack.rda"#"globalStack.rda" or "europeStack.rda"
Stack <- stack(paste(Stackpath,"/Biooracle.filled.layers.global2025",".tif", sep="")) #"globalStack.rda"#"globalStack.rda" or "europeStack.rda"

#par(mfrow = c(2,2))
#    for(su in 1:3){
#suffix2 <- suffixes[su]
#Stack2 <- stack(paste(Stackpath,"/rasterstack.global.2022", suffix2,".TIFF", sep="")) #"globalStack.rda"#"globalStack.rda" or "europeStack.rda"
#plot(mean(Stack2))
#    }
#Stack2 <- stack(paste(Stackpath,"/rasterstack.Europe.2022.no.chlora.delim.TIFF", sep=""))
#plot(mean(Stack))

#rm(layernames)
#not using alternative stacks
#load(paste(Stackpath,"/layernames",suffix,".rda" ,sep =""))
load(paste(Stackpath,"/layernames",".rda" ,sep =""))
names(Stack) <- layernames

#Prevent that strings in csv files are red in as factors variables
stringsAsFactors= FALSE

stats <- c()
#for(Species in Data.table$species[1:2]){
  for(Species in Data.table$species[-1]){
    
  print(Species)
  # Read in present absent and pseudoabsent points, 
    #convert these points to spatial coordinates and extract environmental variables from rasterstack
 species.data.list <-  read.and.extract(data.table = Data.table,
                                   species = Species,
                                   stack = Stack, 
                                   speciespath = SpeciespathRaw,
                                   stackpath = Stackpath,
                                   plotpath = Plotpath,
                                   outpath = Outpath)
 species.data <- species.data.list[["complete.points"]]
 stats.temp <- c(Species,species.data.list[["stats"]])
 stats <- rbind(stats, stats.temp)
 head(species.data)
 ### Check that there are no erroneous entries in the files. 
 ## the most common synonyms for "present" and "absent" are identified and entries harmonized.
 #For entires where status cannot be determiend the line is removed
 # If data was checked in the earleir stage these lines should not find any mistakes
 my.data <- species.data
 print(paste(Species,paste(unique(my.data$occurrenceStatus) ))  )
 my.data$occurrenceStatus <- as.character(my.data$occurrenceStatus)
 present.synonyms <-  which(my.data$occurrenceStatus == "present"|my.data$occurrenceStatus == "Present"| my.data$occurrenceStatus == "established"| my.data$occurrenceStatus == "Established")
 my.data$occurrenceStatus[present.synonyms] <- "present"
 absent.synonyms <-  which(my.data$occurrenceStatus == "Absent"| my.data$occurrenceStatus == "")
 
 my.data$occurrenceStatus[absent.synonyms] <- "absent"
 remove <- which(!my.data$occurrenceStatus == "absent"  & !my.data$occurrenceStatus == "present")
 print(paste("remove",remove))# promt if lines are removed
 if(length(remove > 0)){
   my.data <- my.data[-remove,]                
 }
 species.data <- my.data
 # Save the species data with environmental variables in selected folder
 print(paste(Species, "npos=", length(which(species.data$occurrenceStatus == "present")),
                                      "nneg=", length(which(species.data$occurrenceStatus == "absent") )))
  write.csv(species.data, file = paste(Outpath,"/",Species,"_indata.csv", sep=""),row.names=F)
}
colnames(stats)<- c("Species","npos","nabs", "npos.unique", "nabs.unique", "npos.unique.complete", "nabs.unique.complete")
write.csv2(as.data.frame(stats),file =  paste(path,"/data/species.stats",suffix,".csv", sep=""))
#print(species)
###
stats <- as.data.frame(stats)
exclude <-  union(which(as.numeric(stats$npos.unique.complete) <6),
           which(as.numeric(stats$nabs.unique.complete) <6))
#krasch <- c(11,17,34)
#exclude <- union(exclude, krasch)
#stats[exclude,]
#stats <- read.csv2(file = "stats.mar 2025.csv", sep=",")
write.csv2(stats, file = paste("stats",suffix,"apr 2025.csv", sep=""))
save(exclude, file = paste("exclude",suffix,".apr2025.csv",sep=""))

Data.table$species[exclude]
#load(file = "excludemar2025.csv")
load(file= paste("exclude",suffix,".apr2025.csv",sep=""))
##################################################################
### prepare iterations
##################################################################


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
#for(Species in Data.table$species[-exclude]){
  for(Species in Data.table$species){
    
#for(Species in Data.table$species[-c(22,23,24,25,26)]){
    
  # Species <- Data.table$species[1] 
 # iterations.path
  print(Species)
  nsites <- split.data(species = Species, 
             indata.path = Outpath,
             iterations.path= Iterations.path)
#unique(all.site.iters.tot[[5]][[5]])
  if(nsites[1]<5 | nsites[2]<5){
    exclude2<- rbind(exclude2, c(Species, nsites[1],nsites[2]))
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
  #Species <-Data.table$species[5]){
  lista.rda<- Sys.glob(paste(Iterations.path,"*.rda",sep="/"))
 # for(Species in Data.table$species[-c(1:20, 22,23,24,25,26)]){
    for(Species in Data.table$species[-exclude]){
    #  for(Species in Data.table$species[-c(1:79)]){
        
    indata.path = Outpath
    lista.csv<- Sys.glob(paste(indata.path,"*.csv",sep="/"))
    my.data <- read.csv(lista.csv[grep(Species,lista.csv)],header=T)
    my.data$RANDOMVAR <- runif(length(my.data$ID),0,1)
    my.data$RANDOMVAR2 <- runif(length(my.data$ID),0,1)
    my.data$RANDOMVAR3 <- runif(length(my.data$ID),0,1)
    
    print(paste(Species, "n.abs=", length(which(my.data$occurrenceStatus == "absent")),
    ":n.pos=" ,length(which(my.data$occurrenceStatus == "present"))))
    unique(my.data$occurrenceStatus )
##    We are using occurrence iters in this project, not NUTS iters
     use.site.inters <- TRUE
    if(use.site.inters){
      load(lista.rda[
        intersect(grep("occurance.iters_",lista.rda),grep(Species,lista.rda))
      ])
    }
     # get info from iterations data structure
     nrep <- length(all.occurance.iters)
     CV.level <- length(all.occurance.iters[[1]])
     process <- seq(1,nrep*CV.level)
     rep <- sort(rep(seq(1,nrep),CV.level))
     iter <- rep(seq(1,CV.level),nrep)
     process.plan <- cbind(process,rep,iter)
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
     
     ## start parallell execution
     MCMCresult <- try(parLapply(cl,1:length(process.plan[,1]),
                             function(i) MCMC.process(.mydata =my.data, 
                                                      .iters =all.occurance.iters,
                                                      .process.plan = process.plan,
                                                      .process.ID = i)))
     stopCluster(cl)
     time.to.complete2 <- proc.time()-ptm2
     print(time.to.complete2)
   #  
  #  
     if(length(grep("Error",MCMCresult))>0){
       print("Error - fallback")
       myreturn <- list()
       for(prid in 1:25){
   myreturn[[prid]] <- MCMC.process(.mydata =my.data, 
                   .iters =all.occurance.iters,
                  .process.plan = process.plan,
                   .process.ID = prid)
     }
     MCMCresult <- myreturn
     }

     save(MCMCresult, file= paste(selectionspath,"/selected.vars.",Species,".rda",sep=""))
     rm(MCMCresult)
        gc()
  }
}

#Neogobius fluviatilis n.abs= 1469 :n.pos= 329"#
# Error in gzfile(file) : invalid 'description' argument
"Ponticola kessleri n.abs= 1469 :n.pos= 116"
#Error in gzfile(file) : invalid 'description' argument
#Corbicula fluminalis n.abs= 1469 :n.pos= 32"
#Error in gzfile(file) : invalid 'description' argument
#[1] "Dikerogammarus villosus n.abs= 1474 :n.pos= 43"
#Error in gzfile(file) : invalid 'description' argument
#[1] "Faxonius rusticus n.abs= 1469 :n.pos= 217"
#Error in gzfile(file) : invalid 'description' argument
##################################################################
## plot selected variables
##################################################################
#Plotpath
#selectionspath
#Species <-Data.table$species[1]){
lista.rda<- Sys.glob(paste(Iterations.path,"*.rda",sep="/"))
lista.selection <- Sys.glob(paste(selectionspath,"*.rda",sep="/"))
indata.path = Outpath
lista.csv<- Sys.glob(paste(indata.path,"*.csv",sep="/"))

#for(Species in Data.table$species){
  for(Species in Data.table$species[-exclude]){
    #  for(Species in Data.table$species[-exclude][-c(4,7,9,10,12,13)]){
    
    my.data <- read.csv(lista.csv[grep(Species,lista.csv)],header=T, sep=",")
  length(which(my.data$occurrenceStatus == "absent"))
  length(which(my.data$occurrenceStatus == "present"))
  
  load(lista.selection[intersect(grep(Species,lista.selection),
                                 grep("selected.vars",lista.selection))])
  all.attributes <- sort(MCMCresult[[1]]$RI$attribute.name)
  
  cutoffs <- mean(sapply(1:length(MCMCresult), function(x) MCMCresult[[x]]$cutoff))
  RI_cutoff <- mean( sapply(1:length(MCMCresult), function(x) MCMCresult[[x]]$RI$RI[  as.numeric(MCMCresult[[x]]$cutoff)] ) )
  
  selection <- as.data.frame(all.attributes)
  selection<- cbind(selection,
                    sapply(1:length(MCMCresult),function(i)
                      as.numeric(MCMCresult[[i]]$RI[ order(MCMCresult[[i]]$RI$attribute.name),c("RI")]
                      ))
)  
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
  plotname <- paste(Plotpath,"/","_",Species,"_seleced_variables",".png",sep="")
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
  title( Species)
  dev.off()# close plotfile
  
  # plot pdf var distribution
  pdf(paste(Plotpath,"/var_distribution",Species,".pdf", sep=""), paper = "a4r")# height= 11.7 , width= 16.6)#paper = "a4r"
 
  random <- grep("RANDOM",raw.labels )
  if(length(random) >0){ raw.labels <- raw.labels[-random]}
  par(mfrow= c(nrow,3))
  par(mar=c(4,2,2,0))
  for(l in 1:length(raw.labels)){
    # plot(c(),c(),xlim=c(0,1), ylim=c(0,1))
    #}
    
    my.var <- as.character(raw.labels[l])# names(my.data)[24]
    # my.var <- "bio_10"
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
    # i <- 1
    # temptab[which(chunks == i),]
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
    title(main= paste(Species,my.var), cex.main = 0.7)
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

##################################################################
  ## Run random forests. the function will prepare data and execute.
##################################################################
rm(Stack)§
rm(map)
rm(my.map)
rm(my.maps)
rm(mystack)
rm(my.data)
gc()
#rlimit_all()
#rlimit_as(1e12)  #increases to ~12GB
#rlimit_stack(100000000)

#for(Species in Data.table$species){
 for(Species in Data.table$species[-exclude]){
    #for(Species in Data.table$species[-c(7 ,9,10,16,18 ,20, 22,23,24,25,26)]){
  #Species <-  Data.table$species[6]
#for(Species in Data.table$species[-c(1:18, 20, 22,23,24,25,26)]){
  # the object "rf.output.list contains all resutls from random forests as a list object
  rf.output.list <- run.random.forests(species = Species ,  
                selvar = "all",
                indata.path = Outpath,
                iterations.path= Iterations.path)
  
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
  save(rf.output.cv, file= paste(Modelpath,"/RF.model.and.predictions.CV.eur.wt.",Species,".rda",sep=""))
  #save the random forest model for late prediction of maps
  save(rf.output, file= paste(Modelpath,"/RF.model.and.predictions.eur.wt.",Species,".rda",sep=""))
  #remove the large object generated to make space for the next species
  rm(rf.output.list)
  rm(rf.output.cv)
  rm(rf.output)
  gc()
  
}
####################
# try C50 rules
#for(Species in Data.table$species){
  for(Species in Data.table$species[-exclude]){
    
  #for(i in Data.table$species[-c(7 ,10,16,18, 20, 22,23,24,25,26)]){
#Species <-  i #Data.table$species[1]
lista.csv<- Sys.glob(paste(Speciespath,"*.csv",sep="/"))
my.data <- read.csv(lista.csv[grep(Species,lista.csv)],header=T)
head(my.data)
names(my.data)
#vars <- c("SST_.AMP", "SST_.MIN", "SST_.MAX", "SST_.MEAN", "CHLORA_.MAX",    "SALINITY")
#vars <- c("SST_.AMP", "SST_.MIN", "SST_.MAX", "SST_.MEAN",    "SALINITY")
vars <- names(my.data)[5:23]
  
rule_mod <- C5.0(x = my.data[, vars], y = as.factor(my.data$occurrenceStatus), rules = TRUE)

sink(paste(Modelpath,"/C50 rules.",Species,".txt",sep=""))
print(Species)
print("/n")
print(summary(rule_mod))
sink()
}
##################################################################
# calculate ROC curves and plot 
##################################################################

all.AUC <- c()
for(Species in Data.table$species[-exclude]){
  #for(Species in Data.table$species){
    
  #for(Species in Data.table$species[-c(7 ,9,10,16,18 ,20, 22,23,24,25,26)]){
  #Species <-  Data.table$species[2]
  load(paste(Modelpath,"/RF.model.and.predictions.CV.eur.wt.",Species,".rda",sep=""))
  load(paste(Iterations.path,"/occurance.iters_",Species,".rda",sep=""))
  
  lista.csv<- Sys.glob(paste(Outpath,"*.csv",sep="/"))
  my.data <- read.csv(lista.csv[grep(Species,lista.csv)],header=T)
  
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
 mean.AUC <-  plot.ROC(.ROC.path=ROC.path,
           .my.species = Species,
           .mean.ROC = mean.ROC,
           .all.ROC =all.ROC)
  #remove objects t make list for next run
    rm(all.occurance.iters)
    rm(rf.output.cv)
    rm(all.ROC)
    rm(mean.ROC)
    gc()
    all.AUC <- rbind(all.AUC, c(Species,mean.AUC))
}
write.csv2(as.data.frame(all.AUC),file =  paste(path,"/data/all.AUC",suffix,".csv", sep=""))
##################################################################
## predict maps
##################################################################
#plots
#load rasterlayers
#load(paste(Stackpath,"europeStack.rda", sep="/") )
#plot(var2[[1]])

#Stack <- stack(paste(Stackpath,"rasterstack.global.2022.TIFF", sep="/")) #"globalStack.rda"#"globalStack.rda" or "europeStack.rda"

#Stack <- stack(paste(Stackpath,"/Biooracle.global2024",".tif", sep="")) #"globalStack.rda"#"globalStack.rda" or "europeStack.rda"
#Stack <- stack(paste(Stackpath,"/Biooracle.Europe2024",".tif", sep="")) #"globalStack.rda"#"globalStack.rda" or "europeStack.rda"

dataset_scenarios <- c( "baseline" ,"ssp119" ,  "ssp126"  , "ssp245"   ,"ssp370" ,  "ssp460"  , "ssp585"  )
dec.vec <- c("", "dec50", "dec100")

dir <- paste(path, "data/Biooracle.download/", sep="/")


#for(sel.sen in c(1)){
  for(sel.sen in c(2,4,5,7)){
    
  scenario <- dataset_scenarios[[sel.sen]]
  for(dec in c[c(2,3)]){
 # for(dec in dec.vec[c(1)]){
    
    print(paste(scenario, dec))    
    # indir <- paste(dir,"datalayer.nc/",dataset_scenarios[[1]][1],"/",sep="")
    #  outdir <- paste(dir,"datalayer.nc/",scenario,dec,"/",sep="")
    Stackpath <- paste(dir,"rasterstacks/",scenario,dec,"/",sep="")
    #if (!dir.exists(Stackpath)) dir.create(Stackpath) 
   
     rastermappath <- paste(path,"/data/Rastermaps",scenario,dec,sep ="")
    if (!dir.exists(rastermappath)) dir.create(rastermappath) 
    
    
#Stack <- stack(paste(Stackpath,"Biooracle.filled.layers.Europe2025",".tif", sep="")) #Biooracle.filled.layers.global2025.tif
     Stack <- stack(paste(Stackpath,"Biooracle.filled.layers.global2025",".tif", sep=""))

#rm(layernames)
#load(paste(Stackpath,"/layernames",suffix,".rda" ,sep =""))
Stackpath.base <- paste(path, c("data/Biooracle.download/rasterstacks/baseline"), sep="/")

load(paste(Stackpath.base,"/layernames",".rda" ,sep =""))

#load(paste(Stackpath,"layernames.filled",".rda" ,sep =""))

print("load layenames")
names(Stack) <- layernames
bar <- Stack# stack(var2)
#plot(bar[[1]])

# Species  <- Data.table$species[1]

#call the function to predict the species distribution at raster level.
#The function return the prediction as a raster object but also make plots as .png 
# the lines between png() and dev.off() may be removed/inactivated if the png plots are not wanted.
#for(Species in Data.table$species){
  for(Species in Data.table$species[-exclude]){
    print(Species)
#  for(Species in Data.table$species[-c(7 ,9,10,16,18 ,20, 22,23,24,25,26)]){
    #predict.maps <- function(species,indata.path,modelpath){

    map <-  predict.maps(species = Species,
            modelpath = Modelpath
            )
    # save predicted probabilities as rasterfiles.
    #plot(map)
    #newfile <- paste(rastermappath,"/linear.prob.",Species,sep="")
    newfile <- paste(rastermappath,"/linear.prob.global.",Species,sep="")
    #print(paste(file, newfile , sep="###########"))
    writeRaster(map, filename= newfile, format = "GTiff", suffix='.tif', overwrite=TRUE)
    
}
}}#end scenario ,dec.vec
######################################################################################################################################
#####################################################################################################################################
## ###################################################################plot maps
##################################################################
######################################################################################################################################

#restore path 
#rastermappath <- paste(path,"/data/Rastermaps", suffix,sep ="")
rastermappath <- paste(path,"/data/Rastermaps",scenarios[1],dec.vec[1],sep ="")


#colorsBrBG <- rev(divPalette(n=12, name = c( "BrBG") ) )
colorsBrBG2 <- rev(divPalette(n=19, name = c( "BrBG") ) )

plot(seq(1:25), seq(1:25), col=colorsBrBG2, pch=16)
colorsBr <- c("white",colorsBrBG2[11:19])
colorsBlue <- seqPalette(n=,12, name = c( "Blues") ) 
#shape.ecoregions <- readOGR(dsn=ICESecopath,
#                            layer="ICES_ecoregions_20171207_erase_ESRI",encoding="UTF-8")
shape.ecoregions <- st_read(paste(ICESecopath,"/ICES_ecoregions_20171207_erase_ESRI.shp",sep=""), )
#plot(shape.ecoregions$geometry)

brk <- c(seq(0, 1,by=0.1),1.05)
# temporary solution
shape2 <-shape.ecoregions
xlim=c(0,30)# for swe
ylim = c(50,70)#for swe
#xlim=c(-180,180)
#ylim=c(-90,90)

#lista.ras <- Sys.glob(paste(rastermappath,"*linear.prob.*",sep="/"))
lista.ras <- Sys.glob(paste(rastermappath,"*linear.prob.global*",sep="/"))
 
predvar<-stack(lista.ras)
plot(mean(predvar))

# Plotpath <- "~/Dokument/Projekt/HAV2025/Plots.marine.new.base.global"
if(!dir.exists(Plotpath)){dir.create(Plotpath)}
  
for(Species in Data.table$species[-exclude]){
    #plot.maps <- function(.map,.plotpath, .colors,.brk,.shape){
  
  rastermap <- lista.ras[grep(Species,lista.ras)]
  
  plot.maps(.Species=Species,
            indata.path = Outpath,
            .rastermap=rastermap,
            .plotpath=Plotpath,
            .colors=colorsBr,
            .brk=brk,
            .shape=shape2,
            xlim=xlim,
            ylim= ylim,
            .cex1 = 0.1,
            .cex2 = 0.25,
            .lwd1 = 0.2,
            .lwd2 = 0.3
            )
 
  #save(map, file= paste(mappath,"/predicted.map.",Species,".rda",sep=""))
 #plot(map)
}
######################################################################################################################################
###################################################plot futire scenarios#########################################################
Plotpath2 <- gsub( ".marine","ssp",Plotpath)
if(!dir.exists(Plotpath2)){dir.create(Plotpath2)}
  
listatiff <- c()
scenarios <- c( "baseline" ,"ssp119" ,  "ssp126"  , "ssp245"   ,"ssp370" ,  "ssp460"  , "ssp585"  )
dec.vec <- c("", "dec50", "dec100")
for(Scenario in scenarios){
  for (dec in dec.vec){
    
rastermappath <- paste(path,"/data/Rastermaps",Scenario,dec,sep ="")
if(dir.exists(rastermappath)){
  print(paste("checking",rastermappath))
lista.tif.temp<- Sys.glob(paste(rastermappath,"/","*.tif",sep=""))
listatiff <- c(listatiff,lista.tif.temp)
}else{print(paste("no dir",rastermappath))}
  }}


for(sel.sen in c(2,3,5,6,7)){
  scenario <- scenarios[[sel.sen]]
  
for(Species in Data.table$species[-exclude]){
#Species <- Data.table$species[-exclude][10]



#get files for scenario
#sel.sen <- 4

  
print(paste(scenario, dec))    
    # indir <- paste(dir,"datalayer.nc/",dataset_scenarios[[1]][1],"/",sep="")
    #  outdir <- paste(dir,"datalayer.nc/",scenario,dec,"/",sep="")

cases <- union(grep("baseline",listatiff),grep(scenario,listatiff))
my.rasterfiles<- listatiff[ intersect(grep(Species,listatiff), cases)  ]

my.stack <- try(stack(my.rasterfiles))
r <- 1

my.layernames <- sapply(1:length(my.rasterfiles), function(r)
  strsplit(strsplit( my.rasterfiles[r],"Rastermaps")[[1]][2],"/")[[1]][1])
names(my.stack)<- my.layernames

my.change <- my.stack[[3]]- my.stack[[1]]

if(length(grep("Error",my.stack[1]))<1 ){
  print(paste(names(my.stack), "stack success"))
  par(mfrow=c(3,1))
#  plot(my.stack)
  #success <- TRUE
}else{
  print(paste(names(my.stack), "stack fail"))
  break
}    #########################
colorsBrBG2 <- rev(divPalette(n=19, name = c( "BrBG") ) )
colorsBrBG3 <- rev(divPalette(n=10, name = c( "BrBG") ) )

#plot(seq(1:25), seq(1:25), col=colorsBrBG2, pch=16)
colorsBr <- c("white",colorsBrBG2[11:19])
colorsBlue <- seqPalette(n=,12, name = c( "Blues") ) 
brk <- c(seq(0, 1,by=0.1),1.05)
brk2 <- c(seq(-1, 1,by=0.1),1.05)

png(paste(Plotpath2,"/",Species, ".",scenario,".png",sep=""),  width = 180, height = 180, units = "mm", res=1200)

par(mfrow=c(2,2))
par(mar=c(2,2,2,6))

for(i in 1:3){
  map<-my.stack[[i]]
  # plot(map)
  
  map2 <- map
  map2[is.na(map2)]<- 1.005
  my.colors <- c(colorsBr,"lightgrey")
  
  plot(map2, col=my.colors, breaks = brk)
  if(i == 1){
    title(main = Species)}else{
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
######################################################################################################################################
### ###################################################################Plot logmaps
#####################################################################################################################################
######################################################################################################################################
#Species  <- Data.table$species[1]

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
plot(mean(predvar))
for(Species in Data.table$species){
  #for(Species in Data.table$species[-exclude]){
    
    rastermap <- lista.ras[grep(Species,lista.ras)]
    map<-raster(rastermap)
 #map[is.na(mean(bar))]<- 0.000001
  lista.csv<- Sys.glob(paste(Outpath,"*.csv",sep="/"))
  my.data <- read.csv(lista.csv[grep(Species,lista.csv)],header=T)
  
  # prepare a rasterlayer with the 10-log of the predicted probability of presence.
  # An arbitrary small number is added, as log10(0) is not defined
  logmap <- log10(map+0.001) #log10(max(map , 0.001, na.rm=T))
  
  # an arbitrary large number is inserted at locations where presence is undefined, in other words land.
 logmap[is.na(mean(bar))]<- 0.24
 # save(logmap, file= paste(mappath,"/predicted.map.log.",Species,".rda",sep="") )
  
 # map[is.na(map)]<- -0.1
 
 # plot maps at Eurpean and Swedish scale this presence points added.
 # Plot without indicating xlim and ylim. This gives a plot area defined by the raster extent
  png(paste(Plotpath,"/",Species, ".logprob.Eur.png",sep=""),  width = 180, height = 180, units = "mm", res=1200)
  plot(logmap,col=colorsBr, breaks=brk.log)
  plot(shape2$geometry,add = TRUE, xlim=xlim, ylim=ylim, border = 1, lwd = 0.1)
  points(my.data$Lon[which(my.data$occurrenceStatus == "absent")],
         my.data$Lat[which(my.data$occurrenceStatus == "absent")], col=4, pch=1, cex=0.3, lwd=0.2)
  points(my.data$Lon[which(my.data$occurrenceStatus == "present")],
         my.data$Lat[which(my.data$occurrenceStatus == "present")], col=2, pch=1, cex=0.3, lwd=0.2)
  dev.off()
  ### same for Sweden
  # make the plot using the xlim and ylim defined earlier.
  png(paste(Plotpath,"/",Species, ".logprob.Swe.png",sep=""),  width = 180, height = 180, units = "mm", res=1200)
  plot(logmap,col=colorsBr, breaks=brk.log , xlim=xlim, ylim=ylim )
  plot(shape2$geometry,add = TRUE, xlim=xlim, ylim=ylim, border = 1, lwd = 0.1)
  points(my.data$Lon[which(my.data$occurrenceStatus == "absent")],
         my.data$Lat[which(my.data$occurrenceStatus == "absent")], col=4, pch=1, cex=1, lwd=0.25)
  points(my.data$Lon[which(my.data$occurrenceStatus == "present")],
         my.data$Lat[which(my.data$occurrenceStatus == "present")], col=2, pch=1, cex=1, lwd=0.25)
    dev.off()
  
  
  
  # save predicted probabilities as rasterfiles.
 # newfile <- paste(rastermappath,"/linear.prob.",Species,sep="")
  #print(paste(file, newfile , sep="###########"))
 # writeRaster(map, filename= newfile, format = "GTiff", suffix='.tif', overwrite=TRUE)
  
    writelograster <- F
    if(writelograster){
 newfile <- paste(rastermappath,"/log.prob.",Species,sep="")
  #print(paste(file, newfile , sep="###########"))
  writeRaster(logmap, filename= newfile, format = "GTiff", suffix='.tif', overwrite=TRUE)
    }
  rm(logmap)
  rm(map)
  gc()
}
######################################################################################################################################

###################################################################
### #################################################Plot combined maps
##################################################################
######################################################################################################################################
# define color palettes. White and lightgrey to incidate values outside range. 
colorsBrBG <- c(rev(divPalette(n=12, name = c( "BrBG") ) ),"lightgrey","lightgrey","lightgrey")
colorsBlue <- c(seqPalette(n=,12, name = c( "Blues") ),"lightgrey","lightgrey","lightgrey" )
colorsBrBG2 <- rev(divPalette(n=22, name = c( "BrBG") ) )
#plot(seq(1:25), seq(1:23), col=colorsBrBG2, pch=16)
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
newfile <- paste(rastermappath,"/mean.prob.all.species.tif",sep="")
writeRaster(mean.predvar, filename= newfile, format = "GTiff", suffix='.tif', overwrite=TRUE)

logpredvar <- log10(mean(predvar, na.rm=T) + 0.001)
writelograster <- F
if(writelograster){
  newfile <- paste(rastermappath,"/log 10 of mean.prob.all.species",sep="")
  writeRaster(logpredvar, filename= newfile, format = "GTiff", suffix='.tif', overwrite=TRUE)
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

#plot(seq(1:25), seq(1:23), col=colorsBrBG2, pch=16)
colorsBr <- c(colorsBr, "lightgrey","lightgrey")
brk <- c(seq(0, 1,by=0.1),1.05)

logpredvar[is.na(mean(bar))]<- 0.25
xlim=c(-180,180)
ylim = c(-180,180)
png(paste(Plotpath,"/mean.all.eur", ".png",sep=""),  width = 180, height = 180, units = "mm", res=1200)
plot(map2,col=colorsBr, breaks=brk)
plot(shape2$geometry,add = TRUE, xlim=xlim, ylim=ylim, border = 1, lwd = 0.1)

dev.off()
xlim=c(0,30)
ylim = c(50,70)
png(paste(Plotpath,"/mean.all.Swe", ".png",sep=""),  width = 180, height = 180, units = "mm", res=1200)
plot(map2,col=colorsBr, breaks=brk,xlim=xlim, ylim=ylim)
plot(shape2$geometry,add = TRUE, xlim=xlim, ylim=ylim, border = 1, lwd = 0.1)

dev.off()
############################ ############################ ############################ 
############################ Plot mean of log ############################ 
############################ ############################ ############################ 
logpredvar[is.na(mean(bar))]<- 0.25
png(paste(Plotpath,"/logmean.all.eur", ".png",sep=""),  width = 180, height = 180, units = "mm", res=1200)
plot(logpredvar,col=colorsBr, breaks=brk.log)
plot(shape2$geometry,add = TRUE, xlim=xlim, ylim=ylim, border = 1, lwd = 0.1)

dev.off()
png(paste(Plotpath,"/logmean.all.Swe", ".png",sep=""),  width = 180, height = 180, units = "mm", res=1200)
plot(logpredvar,col=colorsBr, breaks=brk.log,xlim=xlim, ylim=ylim)
plot(shape2$geometry,add = TRUE, xlim=xlim, ylim=ylim, border = 1, lwd = 0.1)

dev.off()

### mean of log10 predvar
custom.pred <- log10(predvar[[1]]+0.001)
for(i in 2:length(names(predvar))){
  custom.pred <- custom.pred+ log10(predvar[[i]]+0.001)
}
plot(custom.pred)
newfile <- paste(rastermappath,"/sum of log10 prob all.species",sep="")
writeRaster(custom.pred, filename= newfile, format = "GTiff", suffix='.tif', overwrite=TRUE)
custom.pred <- custom.pred/length(names(predvar))
plot(custom.pred)
newfile <- paste(rastermappath,"/mean of log10 prob all.species",sep="")
writeRaster(custom.pred, filename= newfile, format = "GTiff", suffix='.tif', overwrite=TRUE)


#######################################
# customzed combined layers
################################## 

########## get filenames for traffic data and load as rasterstack


newfile <- paste(rastermappath,"/mean.prob.all.species.tif",sep="")
mean.predvar <- raster(newfile)

lista.ras2 <- Sys.glob(paste(trafficpath,"*.tif",sep="/"))

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

# create the combined data layer sing some mathematical expression
combivar <- plotvar2 + shipvar2/1000
#combivar[combivar > 3 ]<- 3
par(mfrow = c(1,1))
plot(combivar)

combivar <- plotvar2 + shipvar2/5000

colorsBlue <- seqPalette(n=,24, name = c( "Blues") )
colorsorange <- seqPalette(n=,24, name = c( "Oranges") ) 
colorsdiv <- divPalette(n=,24, name = c( "BrBG") )

####
#par(mfrow = c(1,3))
#par(mar=rep(1,4))
#plot(plotvar.crop)# the risk map
#plot(shipvar2)# the traffic map
#plot(combivar)

###
#par(mfrow = c(1,1))
#plot(plotvar.crop, col= colorsBlue)
#arg <- list(at=c(0,0), labels=c("NA","NA"))
#plot(shipvar2, add=T, col=colorsorange,alpha=0.5, legend=FALSE)
#plot(shape2,add = TRUE, xlim=xlim, ylim=ylim, border = 1, lwd = 0.1)


############################## make plot as in old project   ############
newfile <- paste(Plotpath,"/combiplot.png",sep="")
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
loadings(rpc$model)


ggRGB(rpc$map,1,2,3, stretch="lin", q=0)
  plots <- lapply(1:3, function(x) ggR(rpc$map, x, geom_raster = TRUE))
 plot(plots[[1]],xlim=xlim, ylim=ylim)
 plot(plots[[2]],xlim=xlim, ylim=ylim) 
 plot(plots[[3]],xlim=xlim, ylim=ylim)
####### obtain oub prediction performance
 for(Species in Data.table$species){
   file <-  paste(Modelpath,"/RF.model.and.predictions.eur.wt.",Species,".rda",sep="")
   #print(file)
   load(file)
   print("########################")
   print("")
   print(Species)
   print(rf.output$"RF.selected"$"confusion")
   
 }
