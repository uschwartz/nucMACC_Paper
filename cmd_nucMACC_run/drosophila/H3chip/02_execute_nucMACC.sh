#!/usr/bin/env bash

AnalysisDir=/home/$USER/nucMACC_Paper
AnnoDir=/home/$USER/nucMACC_Paper/data/Bowtie2Index

cd $AnalysisDir

nextflow run nucMACC \
--csvInput $AnalysisDir'/script/H3chip/additionalFiles/H3ChIP_input.csv' \
--outDir $AnalysisDir'/H3_ChIP' \
--genomeIdx $AnnoDir'/genome' \
--genomeSize 162367812 \
--genome $AnnoDir'/data/genome.fa' \
--blacklist $AnalysisDir'/data/dm3-blacklistChromosomes.bed' -resume
