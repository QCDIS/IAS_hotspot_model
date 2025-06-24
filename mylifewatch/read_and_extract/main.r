# Load the functions/SEanalytics.functions.r

source("functions/SEanalytics.functions.r")

# Call the read.and.extract function
read.and.extract("/mnt/inputs/data_table",
                 "/mnt/inputs/species_to_use",
                 "/mnt/inputs/rasterstacks",
                 "/mnt/outputs/plots",
                 "/mnt/outputs/outpath")

