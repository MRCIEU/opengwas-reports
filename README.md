# Report module for MR-Base

Pipeline report, metrics, and plots.

![dummy_example](./assets/dummy_example.png)

# Usage

```
 $ Rscript render_gwas_report.R --help
usage: render_gwas_report.R [-h] [--gwas_id GWAS_ID] [--input INPUT]
                            [--metadata METADATA] [-s]

optional arguments:
  -h, --help           show this help message and exit
  --gwas_id GWAS_ID    Directory with the associated gwas_id [default:
                       example]
  --input INPUT        Input bcf file [default: harmonised.bcf]
  --metadata METADATA  metadata json file: [default harmonised.json]
  -s, --show           If True, show the report after it is generated
                       [default: False]
```

## A dummy example:

Suppose we have a GWAS QC run named `example` as the ID of this run,
the harmonised bcf files named `harmonised.bcf` and `harmonised.bcf.csi`,
and the metadata file named as `harmonised.json`.
Create a file structure to `./gwas_input` as:

>     .
>     ├── gwas_input
>     │   └── example
>     │       ├── harmonised.bcf
>     │       ├── harmonised.bcf.csi
>     │       ├── harmonised.json
>     │       └── report_example.html

Build the report for this QC run as

```
 $ Rscript render_gwas_report.R --gwas_id example --input harmonised.bcf --metadata harmonised.json
```

If the build is successful, the report will be generated as `./gwas_input/${example}/report_${example}.html` with `example` as the `gwas_id`:

>     .
>     ├── gwas_input
>     │   └── example
>     │       ├── harmonised.bcf
>     │       ├── harmonised.bcf.csi
>     │       ├── harmonised.json
>     │       ├── intermediate
>     │       │   ├── report_query.tsv
>     │       │   └── rmd_intermediate_files
>     │       └── report_example.html
