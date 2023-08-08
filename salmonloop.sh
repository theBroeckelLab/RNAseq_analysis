#!/bin/bash
#List of base name of files to loop
filelist=$1
#folder location where the fastq for fileslist names are stored
fastq=$2

while read -r line;
do
	salmon quant -i salmon_index_human -l A -1 ${fastq}/${line}_R1.fastq.gz -2 ${fastq}/${line}_R2.fastq.gz -p 32 -o quants/${line}_quant

done < $filelist
