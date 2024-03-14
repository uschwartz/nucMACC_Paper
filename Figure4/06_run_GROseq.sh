#!/usr/bin/env bash

projectDir=/home/$USER/nucMACC_Paper/data
AnalysisDir=/home/$USER/nucMACC_Paper/data/R001_nucMacc
cd $AnalysisDir

outDir=$AnalysisDir/manuscript_figures/Fig4/data/GROseq/raw

nextflow run RNAseq/main.nf  \
    --fastqPath $outDir \
	--outPath $AnalysisDir"/manuscript_figures/Fig4/GROseq/" \
	--STARidxPath $projectDir"/Drosophila_melanogaster_UCSC_dm3/STARidx" \
	--gtfPath projectDir"/Drosophila_melanogaster_UCSC_dm3/" \
	--gtfFile protein_coding.gtf --strandness unstranded
