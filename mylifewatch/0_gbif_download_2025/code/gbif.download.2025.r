library(dplyr)
library(readr)  
library(rgbif) # for occ_download

#1

#gbif.org credentials
user=args$gbif_username
pwd=args$gbif_password
email=args$email


file_url <- "/mnt/inputs/NIS_list_combined_Mar2025_v2.csv"
download_path <- "/mnt/outputs/"
key <- "0010903-240202131308920"
download_file = paste(download_path,key,".zip",sep ="" )


########################################
gbif_taxon_keys <- 
  readr::read_delim(file_url, delim =",",na = c("", "NA"), comment = "",   col_names = TRUE,skip_empty_rows = TRUE)%>%
  pull("Taxon name") %>% # use fewer names if you want to just test 
  name_backbone_checklist()  %>% # match to backbone
  filter(!matchType == "NONE") %>% # get matched names
  pull(usageKey) # get the gbif taxonkeys

gbif_taxon_keys <- gbif_taxon_keys

occ_download(
  pred_in("taxonKey", gbif_taxon_keys),
  format = "SIMPLE_CSV",
  user=user,
  pwd=pwd,
  email=email
)

###########################
d <- occ_download_get(key,path=download_path) %>%
  occ_download_import()