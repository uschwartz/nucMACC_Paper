#!/usr/bin/env bash

AnalysisDir=/home/$USER/nucMACC_Paper
DataDir=$OutDir"/data"

OutDir=$AnalysisDir"/manuscript_figures/Fig1/"
mkdir -p $OutDir


mkdir -p $OutDir"/heatmaps/monoNucs_occu/"
cd $OutDir"/heatmaps/monoNucs_occu/"


######## H4
H4pooled=$AnalysisDir/H4_ChIP_v2.1_featureCOUNTS/RUN/05_MONO-NUCS_PROFILE/pooled_monoNucs_profile.bw

###compute matrix
computeMatrix reference-point -S $H4pooled \
 -R $DataDir"/expression/dm3.refGene.gtf" \
 --referencePoint TSS \
 -o "computeMatrix2plot_H4.txt.gz" \
 --outFileNameMatrix "computeMatrix2txt_H4.txt.gz" \
 --outFileSortedRegions "computeMatrix_geneList_H4.bed" \
 -b 1500 -a 1500 --smartLabels -p 10


###
plotHeatmap -m "computeMatrix2plot_H4.txt.gz" \
     -out 'DefaultHeatmap_H4.png' \
     --outFileSortedRegions 'sortedRegions_Heatmap_H4.txt'\
     --outFileNameMatrix 'values_Heatmap_H4.txt'

######## H3
H3pooled=$AnalysisDir/H3_ChIP_v2_featureCOUNTS/RUN/05_MONO-NUCS_PROFILE/pooled_monoNucs_profile.bw

###compute matrix
computeMatrix reference-point -S   $H3pooled \
      -R $DataDir"/expression/dm3.refGene.gtf" \
      --referencePoint TSS \
      -o "computeMatrix2plot_H3.txt.gz" \
      --outFileNameMatrix "computeMatrix2txt_H3.txt.gz" \
      --outFileSortedRegions "computeMatrix_geneList_H3.bed" \
      -b 1500 -a 1500 --smartLabels -p 10

###
plotHeatmap -m "computeMatrix2plot_H3.txt.gz" \
          -out 'DefaultHeatmap_H3.png' \
          --outFileSortedRegions 'sortedRegions_Heatmap_H3.txt'\
          --outFileNameMatrix 'values_Heatmap_H3.txt'
