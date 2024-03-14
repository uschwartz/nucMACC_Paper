#!/usr/bin/env bash

ProjectDir=/home/$USER/nucMACC_Paper/
AnalysisDir=$ProjectDir/NETseq
DataDir=$AnalysisDir"/data"
OutDir=$AnalysisDir"/Pausing_Profile_sub"
AnnoDir=$ProjectDir/data/MNase_Yeast/annotation/TSS_TTS/sub

mkdir -p $OutDir

#### Unscale downstream

computeMatrix scale-regions -S $DataDir"/NETseq2.bw" \
      -R  $AnnoDir"/NDR_TSS_TTS.bed" \
      $AnnoDir"/unstable_TSS_TTS.bed" \
      -o $OutDir"/computeMatrix2plot_5scale.txt.gz" \
      --outFileSortedRegions $OutDir"/sortedRegions.bed"\
      --unscaled5prime 150 \
      -m 2000 --skipZeros \
      -b 500 -a 500 --smartLabels -p 5

plotHeatmap -m $OutDir"/computeMatrix2plot_5scale.txt.gz" \
      -out $OutDir"/DefaultHeat_5scale.pdf" \
      --outFileSortedRegions $OutDir"/sortedRegions_Heatmap.bed" \
      --outFileNameMatrix $OutDir'/values_heat_5scale.txt'
