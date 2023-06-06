#!/bin/bash

input=$(echo SRR1039508 SRR1039509 SRR1039512	SRR1039513 SRR1039516 SRR1039517 SRR1039520 SRR1039521)

for files in ${input[@]}
do 

fastq-dump --split-3 --outdir '/home/berrit.kievith/daur2/les1/fastq_files/.' --gzip $files


done





