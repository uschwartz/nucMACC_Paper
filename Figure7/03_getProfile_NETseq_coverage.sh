#!/usr/bin/env bash

AnalysisDir=~/Analysis/NETseq/
DataDir=$AnalysisDir"/data"
OutDir=$AnalysisDir"NETseq_Profile_hyper_hypo_neu"

mkdir -p $OutDir

#### Unscale downstream

computeMatrix reference-point -S $DataDir"/NETseq_coverage.bw" \
      -R  $DataDir"/hyper_TSS.bed" \
       $DataDir"/hypo_TSS.bed" \
      -o $OutDir"/computeMatrix2plot_5scale.txt.gz" \
      --outFileSortedRegions $OutDir"/sortedRegions.bed"\
      -b 300 -a 300 --smartLabels -p 5

plotHeatmap -m $OutDir"/computeMatrix2plot_5scale.txt.gz" \
     -out $OutDir"/DefaultHeat_5scale.pdf" \
     --outFileSortedRegions $OutDir"/sortedRegions_Heatmap.bed" \
     --outFileNameMatrix $OutDir'/values_heat_5scale.txt'
