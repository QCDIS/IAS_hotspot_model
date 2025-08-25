FROM ghcr.io/qcdis/habitat_suitability_model:v0.2.1

RUN apt-get update && apt-get install -y \
                               gdal-bin \
                               proj-bin \
                               libgdal-dev \
                               libproj-dev \
                               r-cran-fbasics


COPY functions/ /wrp/functions

# Install devtools::install_github("azizka/speciesgeocodeR")
RUN Rscript -e "install.packages('devtools')"
RUN Rscript -e "devtools::install_github('azizka/speciesgeocodeR')"
RUN Rscript -e "devtools::install_github('bio-oracle/biooracler')"
RUN Rscript -e "install.packages('rgbif')"
RUN Rscript -e "install.packages('CoordinateCleaner')"
RUN Rscript -e "install.packages('speciesgeocodeR')"
RUN Rscript -e "install.packages('sf')"
RUN Rscript -e "install.packages('raster')"
RUN Rscript -e "install.packages('terra')"
RUN Rscript -e "install.packages('dplyr')"
RUN Rscript -e "install.packages('readr')"
RUN Rscript -e "install.packages('ggplot2')"
RUN Rscript -e "install.packages('maps')"
RUN Rscript -e "install.packages('glue')"
RUN Rscript -e "install.packages('randomForest')"
RUN Rscript -e "install.packages('fBasics')"
RUN Rscript -e "install.packages('rnaturalearth')"
RUN Rscript -e "install.packages('rnaturalearthdata')"
RUN Rscript -e "install.packages('ncdf4')"
RUN Rscript -e "install.packages('RNetCDF')"
RUN Rscript -e "install.packages('readxl')"
RUN Rscript -e "install.packages('parallel')"