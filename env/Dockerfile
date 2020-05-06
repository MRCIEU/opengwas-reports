FROM continuumio/miniconda3:4.7.10-alpine
USER root

RUN apk add --no-cache bash font-noto tzdata git
SHELL ["bash", "-c"]

# Configure conda
# RUN conda update --yes conda
COPY ./env/environment.yml ./
RUN /opt/conda/bin/conda env create -f environment.yml
ENV PATH /opt/conda/envs/ieu-gwas-report/bin:$PATH

# Codebase
RUN mkdir -p /home/ieu-gwas-report
WORKDIR /home/ieu-gwas-report

# Launch app
# CMD ["Rscript", "render_gwas_report.R"]
