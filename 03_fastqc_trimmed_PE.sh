#!/bin/bash

sampleIDs=("dmso.br1" "dmso.br2" "lbet151.br1" "lbet151.br2" "sgc0946.br1" "sgc0946.br2")

myDir=$(pwd)

for sampleID in "${sampleIDs[@]}"
do

	mkdir -p ${myDir}/../output/${sampleID}/fastqc/trimmed/

	## run fastqc on the trimmed fastq

	fastqc -o ${myDir}/../output/${sampleID}/fastqc/trimmed/ ${myDir}/../output/${sampleID}/trimmomatic/${sampleID}.trim_1.fq.gz ${myDir}/../output/${sampleID}/trimmomatic/${sampleID}.trim_2.fq.gz

done


