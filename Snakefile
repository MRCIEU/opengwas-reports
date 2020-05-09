import os
import re
import subprocess
import json

with open("igd-config.json", "r") as f:
  config = json.load(f)

igd = config['igd_dir']
out = config['out_dir']
reference = config['reference']

os.makedirs(out + '/job_reports', exist_ok=True)

ID = [x.strip() for x in [y for y in os.listdir(igd)] if 'eqtl-a' not in x]



rule all:
	input:
		expand("{out}/{id}/{id}_report.html", out=out, id=ID)

rule report:
	output:
		"{out}/{id}/{id}_report.html"
	conda:
		"env/environment.yml"
	shell:
		"""
mkdir -p {out}/{wildcards.id}
cp {igd}/{wildcards.id}/clump.txt {out}/{wildcards.id}
cp {igd}/{wildcards.id}/ldsc.txt.log {out}/{wildcards.id}
cp {igd}/{wildcards.id}/{wildcards.id}.json {out}/{wildcards.id}
/usr/bin/time -v Rscript render_gwas_report.R \
  --refdata {reference} \
  --output_dir {out}/{wildcards.id} \
  --n_cores 1 \
  --reuse \
  {igd}/{wildcards.id}/{wildcards.id}.vcf.gz
rm {out}/{wildcards.id}/clump.txt
rm {out}/{wildcards.id}/ldsc.txt.log
rm {out}/{wildcards.id}/{wildcards.id}.json
		"""
