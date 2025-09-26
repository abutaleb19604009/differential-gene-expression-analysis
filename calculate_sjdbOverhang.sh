#!/bin/bash

SEQLEN=$(gzip -cd ../input/dmso.br1_1.fastq.gz | head -2 | tail -1 | wc -c)
SEQLEN1=$(($SEQLEN - 2))

echo "$SEQLEN1"
