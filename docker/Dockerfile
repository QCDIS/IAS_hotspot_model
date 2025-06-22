FROM ghcr.io/qcdis/habitat_suitability_model:v0.2.1

RUN apt-get update && apt-get install -y \
                               gdal-bin \
                               proj-bin \
                               libgdal-dev \
                               libproj-dev \
                               r-cran-fbasics


COPY functions/ /wrp/functions