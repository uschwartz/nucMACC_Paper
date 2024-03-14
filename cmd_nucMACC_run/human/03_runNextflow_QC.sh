#!/usr/bin/env bash

projectDir=/home/$USER/nucMACC_Paper
AnalysisDir=/home/$USER/nucMACC_Paper/data/nucMACC_H2AZ
AnnoDir=/home/$USER/nucMACC_Paper/data

cd $AnalysisDir

nextflow run nucMACC \
--analysis 'MNaseQC' \
--csvInput $projectDir'/script/addData/samples_WC.csv' \
--outDir $AnalysisDir'/WC_QC_out' \
--genomeIdx $AnnoDir'/bowtie2idx/bt2_' \
--genomeSize 2913022398 \
--publishBamFlt \
--TSS $AnnoDir'/GRCh38/protein_coding.gtf' \
--blacklist $AnnoDir'/GRCh38/blacklist_comp_ENCFF356LFX_MT_rDNA_chr_nonCanonical.bed'
