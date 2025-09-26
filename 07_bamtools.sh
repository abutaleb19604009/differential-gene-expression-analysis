#!/bin/bash

sampleIDs=("dmso.br1" "dmso.br2" "lbet151.br1" "lbet151.br2" "sgc0946.br1" "sgc0946.br2")

myDir=$(pwd)

for sampleID in "${sampleIDs[@]}"
do
        mkdir -p ${myDir}/../output/${sampleID}/bamtools/

        ## run bamtools 

        bamtools stats -in ${myDir}/../output/${sampleID}/markduplicates/${sampleID}.markdup.bam > ${myDir}/../output/${sampleID}/bamtools/${sampleID}.markdup.dupstats.txt

done

