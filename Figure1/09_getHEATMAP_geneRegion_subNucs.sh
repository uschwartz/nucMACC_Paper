#!/usr/bin/env bash

AnalysisDir=/home/$USER/nucMACC_Paper

OutDir=$AnalysisDir"/manuscript_figures/Fig1/"
DataDir=$OutDir"/data"


########
mkdir -p $OutDir"/heatmaps/subNucs_geneRegion/"
cd $OutDir"/heatmaps/subNucs_geneRegion/"

######## scale regions ########################

### BED protein coding genes #############

computeMatrix scale-regions -S $DataDir"/bigwigs_subNucs//"*".bigwig" \
 -R $DataDir"/expression/dm3.refGene_coding.bed" \
 -o "computeMatrix2plot_bed.txt.gz" \
 --outFileNameMatrix "computeMatrix2txt_bed.txt.gz" \
 --outFileSortedRegions "computeMatrix_geneList_bed.bed" \
 -b 500 -a 500 --smartLabels -p 10


###
plotHeatmap -m "computeMatrix2plot_bed.txt.gz" \
     -out 'DefaultHeatmap_bed.png' \
     --outFileSortedRegions 'sortedRegions_Heatmap_bed.txt'\
     --outFileNameMatrix 'values_Heatmap_bed.txt'
