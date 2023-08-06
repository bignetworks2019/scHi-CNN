#!/bin/bash

input_dir="sample-data"
cool_dir="sample-cool-out"
resol=100000
chr_len="hg19.chrom.sizes"
genome="hg19"

mkdir -p ${cool_dir}

for file in "${input_dir}"/*; do
  echo $file
  cooler cload pairs "${chr_len}":"$resol" "$file" "$file".cool --assembly "$genome" -c1 2 -p1 3 -c2 4 -p2 5
done

mv "${input_dir}"/*.cool "${cool_dir}"/

