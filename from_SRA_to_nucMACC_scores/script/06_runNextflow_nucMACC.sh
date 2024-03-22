#!/usr/bin/env bash
MY_PWD=$1
cd $MY_PWD


# Run the nucMACC pipeline

nextflow run nucMACC \
--analysis 'nucMACC' \
--csvInput $MY_PWD'/data/samples_nucMACC.csv' \
--outDir Run_out \
--genomeIdx $MY_PWD'/data/Bowtie2Index/genome' \
--genomeSize 162367812 \
--genome $MY_PWD'/data/genome.fa' \
--bamEntry \
--TSS $MY_PWD'/data/genes.gtf'
