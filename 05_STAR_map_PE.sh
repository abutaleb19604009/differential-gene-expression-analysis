#!/bin/bash

sampleIDs=("dmso.br1" "dmso.br2" "lbet151.br1" "lbet151.br2" "sgc0946.br1" "sgc0946.br2")

myDir=$(pwd)

for sampleID in "${sampleIDs[@]}"
do

	mkdir -p ${myDir}/../output/${sampleID}/star/

	## run the star mapping command

	STAR --readFilesCommand gunzip -c --outSAMunmapped Within KeepPairs --outMultimapperOrder Random --outSAMmultNmax 1 --runThreadN 4 --runMode alignReads --quantMode GeneCounts --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ${myDir}/../output/${sampleID}/star/${sampleID}.onemap. --genomeDir ${myDir}/../resources/ --readFilesIn ${myDir}/../output/${sampleID}/trimmomatic/${sampleID}.trim_1.fq.gz ${myDir}/../output/${sampleID}/trimmomatic/${sampleID}.trim_2.fq.gz

done


