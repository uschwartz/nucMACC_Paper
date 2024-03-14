#!/usr/bin/env bash

ProjectDir=/home/$USER/nucMACC_Paper/
DataDir=/home/$USER/nucMACC_Paper/data/S009_Rap1_rDNA
AnalysisDir=$ProjectDir/05_NETseq_mRNAseq

cd $ProjectDir

nextflow run uschwartz/RNAseq -r main  \
	--fastqPath $DataDir/data/NET_and_mRNA_seq/raw/fastqs/pooled_fastq \
	--outPath $AnalysisDir \
	--STARidxPath $ProjectDir/data/yeast/STARidx \
	--gtfPath $ProjectDir/data/yeast/ \
	--gtfFile all_genes.gtf  \
    -w $AnalysisDir/work
