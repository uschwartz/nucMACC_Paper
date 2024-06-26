#!/usr/bin/env bash

AnalysisDir=/home/$USER/nucMACC_Paper

OutDir=$AnalysisDir"/manuscript_figures/Fig1/"
DataDir=$OutDir"/data"

########
mkdir -p $OutDir"/heatmaps/monoNucs_TES/"
cd $OutDir"/heatmaps/monoNucs_TES/"

###compute matrix
computeMatrix reference-point -S $DataDir"/bigwigs_monoNucs//"*".bigwig" \
 -R $DataDir"/expression/dm3.refGene.gtf" \
 --referencePoint TES \
 -o "computeMatrix2plot_all.txt.gz" \
 --outFileNameMatrix "computeMatrix2txt_all.txt.gz" \
 --outFileSortedRegions "computeMatrix_geneList_all.bed" \
 -b 1500 -a 1500 --smartLabels -p 10


###
plotHeatmap -m "computeMatrix2plot_all.txt.gz" \
     -out 'DefaultHeatmap_all.png' \
     --outFileSortedRegions 'sortedRegions_Heatmap_all.txt'\
     --outFileNameMatrix 'values_Heatmap_all.txt'
