##require(rgdal)
#require(CoordinateCleaner)
#require(speciesgeocodeR)# installed via zip-tar frpm CRAN
#require(sf)
#require(raster)
#path_home
path <- "~/Dokument/Projekt/HAV2025/"
#file <- "Data2022/Corbicula fluminea GBIF/occurrence.txt"

#file <- "Data2022/GBIF alla arter 0429558-210914110416597/GBIF alla arter 0429558-210914110416597.csv"
#download 1, freshwater file <- "data/speciesdata/0010903-240202131308920.csv"
#file <- "data/speciesdata/0013730-240202131308920.csv"

#file <- "data/speciesdata/0011294-240216155721649.csv" #2 arter 22a feb
#file <- "data/species.data/0018230-250310093411724.csv"# Arter från Matthias Mars 2024

#install.packages("unix") 
#library(unix)
#rlimit_all()
#rlimit_as(1e48) 

require(data.table)

#file <- "data/species.data/0015088-250402121839773.csv"# bakgrundstaxa från Matthias Apr 2025 bara från 2023
file <- "data/species.data/0015306-250402121839773.csv"# bakgrundstaxa från Matthias Apr 2025 

input <-  fread(paste(path,file,sep=""), data.table=FALSE)#load it in as a data.frame 

length(which(!is.na(input$depth)))



head(input)
dim(input)

input <- input[which(!is.na(input$depth)),]
unique.records <- which(!duplicated(input[,c("family","occurrenceStatus" ,"locality")]))
length(unique.records)

head(input)
types <- unique(input$basisOfRecord)
for(type in types){
print(paste(type, "antal ",length(which(input$basisOfRecord == type))))
}
for(taxon in unique(input$taxonKey)){
  print(paste(taxon, "antal ",length(which(input$taxonKey == taxon))))
}
length(unique(input$locality))

unique(input$order)
my.sample <- c()
for(order in unique(input$order)){
  x <- sample(which(input$order == order),min(length(which(input$order == order)),10000), replace = F)
  my.sample <- c(x, my.sample)
}
length(my.sample)
reduced.input <- input[my.sample,]
write.csv(reduced.input, file="data/species.data/background_2023_n317714.csv",row.names = F)