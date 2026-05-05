#!/bin/bash

RSCAPE=~/bio_project/rscape_v2.6.4/bin/R-scape
BASE=/mnt/c/Users/timod/OneDrive/Documentos/School/UNI/MA1/SEM2/bio/bioinformatics_project_2526/RQ3_2

for noise in 1 3 5 10 25 50 75
do
  for seed in 1 2 3 4 5
  do
    cd "$BASE/$noise/$seed" || exit

    infile="RF00162_${noise}_${seed}.sto"
    outfile="RF00162_${noise}_${seed}_out.txt"

    echo "Running $infile..."

    $RSCAPE -s --cacofold "$infile" > "$outfile"

    echo "Done: $outfile"
  done
done

echo "All runs complete."