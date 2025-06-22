#'  Read package in ./req/pacakges.txt and install depends on the repo
#'  if its a binary package (e.g binary:re-cran-remotes) its going to use
#'  apt-get install, if not is going to find the remotes function and use it
#'  (e.g. github:jeroen/jsonlite -> remotes::install_github(jeroen/jsonlite))

library(remotes)
library(doParallel)
FILENAME <- "./req/packages.txt"

# read package
packages <- readLines(FILENAME)
# take package types c('cran', 'binary', ...)
package_types <- unique(sapply(strsplit(packages, ":"), function(x) x[1]))
# create a list of package by type list(cran=c('dplyr', 'data.table', ...))
packages_by_types <- sapply(package_types, function(x) {
    prefix <- paste0("^", x, ":", collapse = "")
    gsub(prefix, "", packages[grepl(prefix, packages)])
}, simplify = F)
# check package with version are only sourced in cran
package_version <- strsplit(packages[grepl('-', packages)], '-')
check_is_version <- sapply(package_version, function(x){x[2]})
package_version <- package_version[grepl("^[[:digit:]|.]+$",check_is_version)]
for (i in seq_along(package_version)){
  package = package_version[[i]][1]
  if(!grepl('^cran', tolower(package))){
    packge_name = unlist(strsplit(package, ':'))[2]
    stop(paste0('Cannot install package ', packge_name,
           'with specific version if source is not CRAN'))
  }
}
# if some package versioned in cran add to packages_by_types
if(length(package_version)>0){
  packages_by_types[['version']] <- package_version
}

# install in parallel packages_by_types
foreach(type=tolower(names(packages_by_types))) %dopar% {
  packages <- packages_by_types[[type]]
  if (type == "binary") {
    cmd <- paste0(c("apt-get install -y -qq ", packages), collapse = " ")
    system(cmd)
  }else{
    library(remotes)
    if(type=='version'){
      for (i in seq_along(packages)){
        pck_name = unlist(strsplit(packages[[i]][1], ':'))[2]
        pck_version = packages[[i]][2]
        install_version(package = pck_name, version = pck_version, Ncpus = -1)
      }
    }else{
      fun <- match.fun(paste0("install_", type, collapse = ""))
      fun(packages, Ncpus = -1)
    }
  }
}
