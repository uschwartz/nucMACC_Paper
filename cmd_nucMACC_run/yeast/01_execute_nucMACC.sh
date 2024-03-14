#!/usr/bin/env bash

AnalysisDir=/home/$USER/nucMACC_Paper/data/MNase_Yeast
AnnoDir=$AnalysisDir'/annotation'
Input=$AnalysisDir'/data'
OutDir=$AnalysisDir'/results'

cd $AnalysisDir

nextflow run  nucMACC \
--csvInput $Input'/input_pooled.csv' \
--outDir $OutDir \
--genomeIdx $AnnoDir'/Bowtie2Index/genome' \
--genomeSize 12100000 \
--genome $AnnoDir'/Bowtie2Index/genome.fa' \
--TSS $AnnoDir'/Yeast_TSS.bed' \
--blacklist $AnnoDir'/blacklistChromosomes.bed' \
