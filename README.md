# Report module for MR-Base

Pipeline report, metrics, and plots.

Usage

```shell
$ rscript render_gwas_report.R --help
usage: render_gwas_report.R [-h] [-g GWAS_ID] [-i INPUT] [-m METADATA] [-q QUIETLY]

optional arguments:
  -h, --help            show this help message and exit
  -g GWAS_ID, --gwas_id GWAS_ID
                        Directory with the associated gwas_id [default: example]
  -i INPUT, --input INPUT
                        Input bcf file [default: harmonised.bcf]
  -m METADATA, --metadata METADATA
                        metadata json file: [default harmonised.json]
  -q QUIETLY, --quietly QUIETLY
                        If not quietly, will open the report file. [default: True]
```
