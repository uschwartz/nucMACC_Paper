#!/usr/bin/env bash

AnalysisDir=/home/$USER/nucMACC_Paper/data/R001_nucMacc
OutDir=$AnalysisDir"/manuscript_figures/Fig3/"
DataDir=$OutDir"/data"

mkdir -p $OutDir"/DNaseHS/"
cd $OutDir"/DNaseHS/"


######## DNaseHS
DNase_bw=$AnalysisDir"/manuscript_figures/Fig2/data/"
beds=$AnalysisDir"/H4_ChIP_v2.1_featureCOUNTS/RUN/10_sub-nucMACC"

###compute matrix 
computeMatrix reference-point -S $DNase_bw"/DNase_S2DHS-seq_r1.bw" \
 -R  $beds"/unStable_subNucs.bed" $beds"/nonCanonical_subNucs.bed" \
 --referencePoint center \
 -o "computeMatrix2plot_DNaseHS.txt.gz" \
 --outFileNameMatrix "computeMatrix2txt_DNaseHS.txt.gz" \
 -b 3000 -a 3000 --smartLabels -p 10


### plot DNase profile
plotProfile -m "computeMatrix2plot_DNaseHS.txt.gz" \
     -out 'DefaultProfile.pdf' \
     --outFileNameData 'values_Profile_DNaseHS.txt'
