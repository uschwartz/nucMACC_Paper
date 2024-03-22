#!/usr/bin/env bash
MY_PWD=$1
cd $MY_PWD


# Run QC for MNase-seq samples

nextflow run nucMACC \
--analysis 'MNaseQC' \
--csvInput $MY_PWD'/data/samples.csv' \
--outDir QC_out \
--genomeIdx $MY_PWD'/data/Bowtie2Index/genome' \
--genomeSize 162367812 \
--publishBamFlt \
--TSS $MY_PWD'/data/genes.gtf'
