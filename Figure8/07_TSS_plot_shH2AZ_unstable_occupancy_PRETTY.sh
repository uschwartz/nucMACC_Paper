#!/usr/bin/env bash

AnalysisDir="/home/$USER/nucMACC_Paper/data/nucMACC_H2AZ"
bw=$AnalysisDir"/TSSplots/"

mkdir -p $AnalysisDir/TSSplots/shH2AZ_unstable_PRETTY/data
cd $AnalysisDir/TSSplots/shH2AZ_unstable_PRETTY/data


########### nuccMACC scores #################
bwOccshH2AZ=/Volumes/PromisePegasus/_Research_/nucMACC_H2AZ/shH2AZ_nucMACC_out/RUN/01_MONO-NUCS_PROFILE
bwOccWC=/Volumes/PromisePegasus/_Research_/nucMACC_H2AZ/WC_nucMACC_out/RUN/01_MONO-NUCS_PROFILE
bwOccH2AZ=/Volumes/PromisePegasus/_Research_/nucMACC_H2AZ/H2AZ_nucMACC_out/RUN/01_MONO-NUCS_PROFILE

computeMatrix reference-point \
  -S $bwOccshH2AZ"/pooled"*".bw" \
    $bwOccWC"/pooled"*".bw" \
    $bwOccH2AZ"/pooled"*".bw" \
    -R 'unstable_TSS_MANsort.bed' \
    --referencePoint TSS  \
    --skipZeros \
    --missingDataAsZero \
    -o "computeMatrix2plot_occ.txt.gz" \
    -b 900 -a 900 --smartLabels -p 12


 plotHeatmap -m "computeMatrix2plot_occ.txt.gz" \
               -out "../heatmap_occ.png" \
               --sortRegions "keep" \
               --outFileNameMatrix "values_Heatmap_occ.txt"
