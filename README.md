# IEU GWAS QC report

Pipeline report, metrics, and plots.

<p float="centre">
  <img src="assets/mrbase-gwas-report.png" width="400" />
  <img src="assets/mrbase-meta-report.png" width="400" />
</p>

**Table of Contents**

- [Report module for IEU GWAS](#report-module-for-ieu-gwas)
- [Local usage](#local-usage)
    - [Setting up environment](#setting-up-environment)
    - [GWAS report](#gwas-report)
        - [Example](#example)
    - [Meta report](#meta-report)
        - [Example](#example-1)
- [Deployment](#deployment)
    - [Setting up environment](#setting-up-environment-1)
    - [GWAS report](#gwas-report-1)


# Setting up

## Local conda

```
# Create conda environment
$ conda env create -f env/environment.yml

$ conda activate ieu-gwas-report
```

## Docker

```
# Build docker image
docker build -t ieu-gwas-report -f env/Dockerfile .

# Produce report
input_dir=/path/to/input_dir
ref_dir=/path/to/ref_dir
docker run --rm --name ieu-gwas-report \
  -v ${input_dir}:/input_dir \
  -v ${ref_dir}:/ref_data \
  -v $(pwd):/home/ieu-gwas-report \
  ieu-gwas-report Rscript render_gwas_report.R \
  --refdata /ref_data /input_dir/data.vcg.gz
```

# Usage

## GWAS report

```
$ Rscript render_gwas_report.R --help
usage: render_gwas_report.R [-h] [--refdata REFDATA] [-j N_CORES] [--id ID] [--output_dir OUTPUT_DIR] [--show]
                            [--reuse] [--no-report]
                            input

Generate report for a GWAS pipeline.

optional arguments:
  -h, --help            show this help message and exit
  --show                If True, show the report after it is generated [default: False]
  --reuse               If True, reuse processed files [default: False]
  --no-report           If True, only do processing and not rmarkdown report [default: False]

required arguments:
  input                 Input data file, path/to/file (e.g. gwas-files/IEU-a-2/IEU-a-2.vcf.gz)

Override config.yml:
  --refdata REFDATA     reference data e.g. 1000 genomes vcf/bcf annotation
                        file, path/to/file
  -j N_CORES, --n_cores N_CORES
                        Number of cores to use for multiprocessing.
  --id ID               ID of the GWAS, by default is the base name of the input file.
  --output_dir OUTPUT_DIR
                        Directory to store outputs, by default is the same to input.
```

## Meta report

```
```
