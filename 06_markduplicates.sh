#!/bin/bash

sampleIDs=("dmso.br1" "dmso.br2" "lbet151.br1" "lbet151.br2" "sgc0946.br1" "sgc0946.br2")

myDir=$(pwd)

for sampleID in "${sampleIDs[@]}"
do

        mkdir -p ${myDir}/../output/${sampleID}/markduplicates/

        ## run the markduplicate and samtools sort commands

	picard MarkDuplicates I=${myDir}/../output/${sampleID}/star/${sampleID}.onemap.Aligned.sortedByCoord.out.bam O=${myDir}/../output/${sampleID}/markduplicates/${sampleID}.markdup.bam M=${myDir}/../output/${sampleID}/markduplicates/${sampleID}.metrics.markdup.txt REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT

	samtools index ${myDir}/../output/${sampleID}/markduplicates/${sampleID}.markdup.bam

done

