#!/bin/bash

# Bins
declare -a bins=('conda' 'plink' 'bcftools')
# Patterns to grep in $(conda env list)
declare -a conda_env_patterns=('mrbase-report' 'ldsc ')
# Dependent paths
declare -a dependent_paths=(
  $(pwd)/ref_data
  $(pwd)/ref_data/1kg_v3_nomult.bcf
  $(pwd)/ref_data/1kg_v3_nomult.bcf.csi
  $(pwd)/ref_data/study-table-13-02-19.tsv
  $(pwd)/ref_data/eur_w_ld_chr
  $(pwd)/ref_data/ld_files
  $(pwd)/ref_data/snplist.gz
)

# bins
echo
echo "# Check existence of conda"
check_bin() {
  if [ -x "$(command -v $1)" ]; then
    echo "GOOD: $1"
  else
    echo "BAD: $1"
  fi
}

for bin in "${bins[@]}"; do
  check_bin $bin
done

# conda envs
echo
echo "# Check existence of conda envs"
conda_envs=$(conda env list)
check_conda_envs() {
  if grep -q -e $1 <<< $conda_envs; then
    echo "GOOD: $1"
  else
    echo "BAD: $1"
  fi
}

for env in "${conda_env_patterns[@]}"; do
  check_conda_envs $env
done

# dependency
check_data() {
  if [ -d $1 ] || [ -f $1 ]; then
    echo "GOOD: $1"
  else
    echo "BAD: $1"
  fi
}

echo
echo "# Check existence of dependency data"
for path in "${dependent_paths[@]}"; do
  check_data $path
done
