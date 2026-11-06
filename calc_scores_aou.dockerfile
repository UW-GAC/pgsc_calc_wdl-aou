FROM rocker/tidyverse:4.6.1

RUN cd /usr/local  && \
    git clone -b aou-updates https://github.com/UW-GAC/pgsc_calc_wdl-aou.git

RUN Rscript -e 'install.packages("optparse")'
