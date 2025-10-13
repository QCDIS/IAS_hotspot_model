library(dplyr)
library(readr)  
library(rgbif) # for occ_download

#1

#gbif.org credentials
user=args$gbif_username
pwd=args$gbif_password
email=args$email
key <- args$key

inputs_path = "/mnt/inputs/"
outputs_path <- "/mnt/outputs/"
download_file = paste0(outputs_path,key,".zip")
nis_list_path <- paste0(inputs_path,"NIS_list_combined_Mar2025_v2.csv")
dest_path =  paste0(nis_list_path, ".zip")
print(paste("nis_list_path:", nis_list_path))
print(paste("dest_path:", dest_path))
download_zip_data_if_not_present_and_unzip(
    data_path = nis_list_path,
    data_url = args$nis_list_url,
    dest_path = inputs_path
    )



########################################
gbif_taxon_keys <- 
  readr::read_delim(nis_list_path, delim =",",na = c("", "NA"), comment = "",   col_names = TRUE,skip_empty_rows = TRUE)%>%
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
d <- occ_download_get(key,path=outputs_path) %>%
  occ_download_import()