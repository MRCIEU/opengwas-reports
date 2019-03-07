#!/bin/bash
# Example usage
# bash symlink-gwas-files.sh --orig_dir /path/to/gwas_parent --target_dir /path/to/target_dir

# Call getopt to validate the provided input.
options=$(getopt -o n --long orig_dir:,target_dir:,dryrun -- "$@")
eval set -- "$options"

dryrun=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --orig_dir)
      orig_dir=$2
      shift 2
      ;;
    --target_dir)
      target_dir=$2
      shift 2
      ;;
    -n|--dryrun)
      dryrun=1
      shift
      ;;
    *)
      break
      ;;
  esac
done

echo "args: "
echo "- orig_dir: ${orig_dir}"
echo "- target_dir: ${target_dir}"
echo "- dryrun: ${dryrun}"

# Find all id directories
id_dirs=${orig_dir}/*/
echo "Number of id directories found" $(echo $id_dirs | tr ' ' '\n' | wc -l)

# mkdir $orig_dir
if [[ ${dryrun} -eq 0 ]] && [[ ! -d ${target_dir} ]]; then
  mkdir -p ${target_dir}
fi

# symlink $target_dir/$id/data.{bcf,bcf.csi} to $orig_dir/$id/data.{bcf,bcf.csi}
for id_dir in ${id_dirs}; do
  id=$(echo $id_dir | sed -e "s@/\$@@" -e "s@${orig_dir}/@@")
  echo $id
  if [[ -f $id_dir/data.bcf ]] && [[ -f $id_dir/data.bcf.csi ]]; then
    if [[ ${dryrun} -eq 0 ]]; then
      mkdir -p $target_dir/$id
      ln -vs $id_dir/data.bcf $target_dir/$id/data.bcf
      ln -vs $id_dir/data.bcf.csi $target_dir/$id/data.bcf.csi
    else
      echo "ln -vs $id_dir/data.bcf $target_dir/$id/data.bcf"
      echo "ln -vs $id_dir/data.bcf.csi $target_dir/$id/data.bcf.csi"
    fi
  fi
done
