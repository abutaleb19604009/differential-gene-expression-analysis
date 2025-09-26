#!/bin/bash

myDir=$(pwd)

## run multiqc

multiqc ../output -n multiQC
