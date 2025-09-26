#!/bin/bash

myDir=$(pwd)

STAR --runThreadN 4 --runMode genomeGenerate --genomeDir ${myDir}/../resources/ --genomeFastaFiles ${myDir}/../resources/Homo_sapiens.GRCh38.dna_sm.chromosome.22.fa --sjdbGTFfile ${myDir}/../resources/Homo_sapiens.GRCh38.109.22.gtf --sjdbOverhang 49


