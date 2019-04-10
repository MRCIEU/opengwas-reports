FROM continuumio/miniconda3:4.5.12
ENV MRBASE_REPORT_ENV docker

# bcftools
RUN apt-get update && apt-get install -y make gcc zlib1g-dev libbz2-dev lzma-dev lzma liblzma-dev fonts-noto
RUN curl -SL https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2 \
  | tar -xvj \
  && bcftools-1.9/configure \
  && make install -C bcftools-1.9

# Configure conda
# RUN conda update --yes conda
COPY ./env/docker.yml docker.yml
RUN conda env create -f docker.yml
ENV PATH /opt/conda/envs/mrbase-report/bin:$PATH

# Update env
RUN /bin/bash -c "exec Rscript -e \"install.packages('argparse', repos = 'http://cran.rstudio.com')\""

# Copy files
RUN mkdir -p /home/mrbase-report
COPY . /home/mrbase-report/
ENV PATH /home/mrbase-report:$PATH

# Add watcher
# RUN wget -O /home/mrbase-report/watcher.py https://raw.githubusercontent.com/MRCIEU/bgc-upload-orchestrator/master/watcher.py?token=AB1fTMa-e8fsZrhTfgrH2VEnYtpvjtCBks5cgZIdwA%3D%3D && chmod 775 /home/mrbase-report/watcher.py
RUN mkdir -p /home/bin
ADD watcher.py /home/bin
RUN chmod 775 /home/bin/watcher.py

WORKDIR /home/mrbase-report
# Keep container persistent
CMD tail -f /dev/null
