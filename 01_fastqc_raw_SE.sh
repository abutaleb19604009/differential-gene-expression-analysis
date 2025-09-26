#!/bin/bash

sampleIDs=("dmso.br1" "dmso.br2" "lbet151.br1" "lbet151.br2" "sgc0946.br1" "sgc0946.br2")

myDir=$(pwd)

for sampleID in "${sampleIDs[@]}"
do

	mkdir -p ${myDir}/../output/${sampleID}/fastqc/raw/	

	## run fastqc on the raw fastq

	fastqc -o ${myDir}/../output/${sampleID}/fastqc/raw/ ${myDir}/../input/${sampleID}.fastq.gz

done


