#!/usr/bin/env bash

AnalysisDir=/home/$USER/nucMACC_Paper
OutDir=$AnalysisDir"/manuscript_figures/Fig1/"
DataDir=$OutDir"/data"

mkdir -p $OutDir"/profiles/monoNucs_TSS_profile/"
cd $OutDir"/profiles/monoNucs_TSS_profile/"


######## H4
H4bw=$AnalysisDir/H4_ChIP_v2.1_featureCOUNTS/RUN/05_MONO-NUCS_PROFILE/

###compute matrix
computeMatrix reference-point -S $H4bw"/H4"*".bw" \
 -R $DataDir"/expression/dm3.refGene.gtf" \
 --referencePoint TSS \
 -o "computeMatrix2plot_H4.txt.gz" \
 --outFileNameMatrix "computeMatrix2txt_H4.txt.gz" \
 --outFileSortedRegions "computeMatrix_geneList_H4.bed" \
 -b 1500 -a 1500 --smartLabels -p 10


###
plotProfile -m "computeMatrix2plot_H4.txt.gz" \
     -out 'DefaultHeatmap_H4.png' \
     --outFileSortedRegions 'sortedRegions_Profile_H4.txt'\
     --outFileNameData 'values_Profile_H4.txt'




######## H3
H3bw=$AnalysisDir/H3_ChIP_v2_featureCOUNTS/RUN/05_MONO-NUCS_PROFILE/

###compute matrix
computeMatrix reference-point -S   $H3bw"/H3"*".bw" \
      -R $DataDir"/expression/dm3.refGene.gtf" \
      --referencePoint TSS \
      -o "computeMatrix2plot_H3.txt.gz" \
      --outFileNameMatrix "computeMatrix2txt_H3.txt.gz" \
      --outFileSortedRegions "computeMatrix_geneList_H3.bed" \
      -b 1500 -a 1500 --smartLabels -p 10

###
plotProfile -m "computeMatrix2plot_H3.txt.gz" \
          -out 'DefaultHeatmap_H3.png' \
          --outFileSortedRegions 'sortedRegions_Profile_H3.txt'\
          --outFileNameData 'values_Profile_H3.txt'
