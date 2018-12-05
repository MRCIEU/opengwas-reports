# Report module for MR-Base

Pipeline report, metrics, and plots.

Usage

```shell
$  Rscript render_gwas_report.R --help
usage: render_gwas_report.R [-h] [--gwas_id GWAS_ID] [--input INPUT] [--metadata METADATA] [-s]

optional arguments:
  -h, --help           show this help message and exit
  --gwas_id GWAS_ID    Directory with the associated gwas_id [default: example]
  --input INPUT        Input bcf file [default: harmonised.bcf]
  --metadata METADATA  metadata json file: [default harmonised.json]
  -s, --show           If True, show the report after it is generated
```
