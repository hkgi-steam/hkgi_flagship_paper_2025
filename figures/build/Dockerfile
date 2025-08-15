# Use the official Jupyter R notebook image as base
FROM quay.io/jupyter/r-notebook:r-4.4.1

# Switch to root to install system dependencies (if needed)
USER root

# Install Arial font
# RUN apt-get update && \
#     apt-get install -y --no-install-recommends \
#     debconf-utils && \
#     echo "ttf-mscorefonts-installer msttcorefonts/accepted-mscorefonts-eula select true" | debconf-set-selections && \
#     apt-get install -y ttf-mscorefonts-installer && \
#     fc-cache -fv

RUN R -e "install.packages('remotes', repos='http://cran.rstudio.com/'); \
          remotes::install_version('circlize', version='0.4.16', repos='http://cran.rstudio.com/'); \
          remotes::install_version('cowplot', version='1.2.0', repos='http://cran.rstudio.com/'); \
          remotes::install_version('ggbreak', version='0.1.5', repos='http://cran.rstudio.com/')"

RUN conda install -c conda-forge \
    r-ggpattern==1.1.4 \
    r-ggrepel==0.9.6 \
    pandas==2.3.1 \
    numpy==2.3.2 \
    plotly==6.3.0

# Install additional Python packages if needed
# RUN pip install --no-cache-dir package1 package2