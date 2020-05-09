# IEU GWAS QC report

Pipeline report, metrics, and plots.

<p float="centre">
  <img src="assets/mrbase-gwas-report.png" width="400" />
  <img src="assets/mrbase-meta-report.png" width="400" />
</p>

## Setting up

Either use conda to create environment

```
# Create conda environment
conda env create -f env/environment.yml
conda activate ieu-gwas-report
```

or run inside a docker container

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

## Usage

### GWAS report

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

### Meta report

To do...


## Generate reports for IGD data

First clone into the `run_igd` branch of the repository

```
git clone git@github.com:MRCIEU/mrbase-report-module.git
cd mrbase-report-module
git checkout run_igd
```

Generate the reports for a large number of existing datasets. First create `igd_config.json` file specifying working directories e.g.:

```json
{
  "igd_dir": "/mnt/storage/private/mrcieu/research/scratch/IGD/data/public",
  "out_dir": "/mnt/storage/private/mrcieu/research/scratch/IGD/data/dev/igd-reports",
  "reference": "/mnt/storage/private/mrcieu/research/mr-eve/vcf-reference-datasets/1000g/1kg_v3_nomult.bcf"
}
```

Use Snakemake to orchestrate the jobs. This will submit the jobs to the cluster. Ideally you would run this command from inside a `screen` session. 

```
snakemake -prk \
-j 200 \
--cluster-config bc4-cluster.json \
--cluster "sbatch \
  --job-name={cluster.name} \
  --partition={cluster.partition} \
  --nodes={cluster.nodes} \
  --ntasks-per-node={cluster.ntask} \
  --cpus-per-task={cluster.ncpu} \
  --time={cluster.time} \
  --mem={cluster.mem} \
  --output={cluster.output} \
  --parsable" \
--cluster-status ./slurm_status.py
```

**Note** that sometimes the conda environment messes things up. In that case just run without the `--use-conda` flag, making sure that all the required R packages as listed in `env/environemnt.yml` are installed.

