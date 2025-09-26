#!/bin/bash

sampleIDs=("dmso.br1" "dmso.br2" "lbet151.br1" "lbet151.br2" "sgc0946.br1" "sgc0946.br2")
contamFile="TruSeq3-SE.fa"

myDir=$(pwd)

## loop over each sample in turn

for sampleID in "${sampleIDs[@]}"
do

	## make a directory for each sample under the ../output/ directory

	mkdir -p ${myDir}/../output/${sampleID}/trimmomatic

	## run the trimmomatic command

	trimmomatic SE -phred33 ${myDir}/../input/${sampleID}.fastq.gz ${myDir}/../output/${sampleID}/trimmomatic/${sampleID}.trim.fq.gz ILLUMINACLIP:${myDir}/../resources/${contamFile}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15:MINLEN:26

done


