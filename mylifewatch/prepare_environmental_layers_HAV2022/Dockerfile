FROM ias_hotspot_model

# Copy dependencies and install them
COPY req req
RUN Rscript req/exe.R

WORKDIR wrp

# copy the folder code and the main.r code
ADD code code

RUN ls -lah


ADD main.r .

# ENTRYPOINT ["Rscript", "--vanilla", "main.r"]
