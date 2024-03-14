#!/usr/bin/env bash

MY_PWD=$(pwd)


# Run the nucMACC pipeline

/home/tools/nextflow run nucMACC \
--analysis 'nucMACC' \
--csvInput $MY_PWD'/data/samples_nucMACC.csv' \
--outDir Run_out \
--genomeIdx $MY_PWD'/data/Bowtie2Index/genome' \
--genomeSize 162367812 \
--genome $MY_PWD'/data/genome.fa' \
--bamEntry \
--TSS $MY_PWD'/data/genes.gtf' 
