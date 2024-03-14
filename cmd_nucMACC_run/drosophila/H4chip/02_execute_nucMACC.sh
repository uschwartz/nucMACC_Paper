#!/usr/bin/env bash

AnalysisDir=/home/$USER/nucMACC_Paper/data
AnnoDir=/home/$USER/nucMACC_Paper/data/Bowtie2Index

cd $AnalysisDir

nextflow run nucMACC \
--csvInput $AnalysisDir'/script/H4chip/additionalFiles/H4ChIP_input.csv' \
--outDir $AnalysisDir'/H4_ChIP/' \
--genomeIdx $AnnoDir'/genome' \
--genomeSize 162367812 \
--genome $AnalysisDir'/genome.fa' \
--blacklist $AnalysisDir'/dm3-blacklistChromosomes.bed' 
