#!/bin/bash

sampleIDs=("dmso.br1" "dmso.br2" "lbet151.br1" "lbet151.br2" "sgc0946.br1" "sgc0946.br2")

myDir=$(pwd)

for sampleID in "${sampleIDs[@]}"
do
        mkdir -p ${myDir}/../output/${sampleID}/featurecounts/

	## run featurecounts

	cd ${myDir}/../output/${sampleID}/featurecounts/ && featureCounts -p -F GTF -t exon -g gene_id -a ${myDir}/../resources/Homo_sapiens.GRCh38.109.22.gtf -o ${myDir}/../output/${sampleID}/featurecounts/${sampleID}.markdup.featurecount ${myDir}/../output/${sampleID}/markduplicates/${sampleID}.markdup.bam

done

