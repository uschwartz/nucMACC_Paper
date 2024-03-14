#!/usr/bin/env bash

ProjectDir=/home/$USER/nucMACC_Paper/data/GRCh38
AnalysisDir="/home/$USER/nucMACC_Paper/data/nucMACC_H2AZ"
bw=$AnalysisDir"/TSSplots/"

mkdir -p $AnalysisDir/TSSplots/shH2AZ_unstable_PRETTY/data
cd $AnalysisDir/TSSplots/shH2AZ_unstable_PRETTY/data


########### (sub-)nuccMACC scores #################

computeMatrix reference-point \
  -S $bw"/shH2AZ/data/shH2AZ_sub-nucMACC_scores.bw" \
    $bw"/WC/data/WC_sub-nucMACC_scores.bw" \
    -R $bw'/shH2AZ/data/unstable_TSS_shH2AZ.bed' \
    --referencePoint TSS  \
    --skipZeros \
    --missingDataAsZero \
    -o "computeMatrix2plot_sub-nucMACC_scores.txt.gz" \
    -b 900 -a 900 --smartLabels -p 12


 plotHeatmap -m "computeMatrix2plot_sub-nucMACC_scores.txt.gz" \
               -out "../heatmap_sub-nucMACC_scores.png" \
               --outFileSortedRegions "sortedRegions_Heatmap_sub_bed.txt" \
               --outFileNameMatrix "values_Heatmap_sub-nucMACC.txt"
















########### get scores in wt ####################
bwWC=$AnalysisDir"/TSSplots/WC/data/"

computeMatrix reference-point \
  -S $bwWC"/"*"scores.bw"  \
    -R $bw'/unstable_TSS_shH2AZ.bed' \
    --referencePoint TSS  \
    --skipZeros \
     --missingDataAsZero \
    -o "computeMatrix2plot_scores_WC.txt.gz" \
    -b 900 -a 900 --smartLabels -p 12


 plotHeatmap -m "computeMatrix2plot_scores_WC.txt.gz" \
               -out "../heatmap_scores_WC.png"


########### get scores in H2AZ ####################
bwH2AZ=$AnalysisDir"/TSSplots/H2AZ/data/"

computeMatrix reference-point \
 -S $bwH2AZ"/"*"scores.bw"  \
   -R $bw'/unstable_TSS_shH2AZ.bed' \
   --referencePoint TSS  \
   --skipZeros \
    --missingDataAsZero \
   -o "computeMatrix2plot_scores_H2AZ.txt.gz" \
   -b 900 -a 900 --smartLabels -p 12


plotHeatmap -m "computeMatrix2plot_scores_H2AZ.txt.gz" \
              -out "../heatmap_scores_H2AZ.png"

################### get occupancy H2AZ ##############
bwOccH2AZ=/Volumes/PromisePegasus/_Research_/nucMACC_H2AZ/H2AZ_nucMACC_out/RUN/01_MONO-NUCS_PROFILE

computeMatrix reference-point \
 -S $bwOccH2AZ"/"*".bw"  \
   -R $bw'/unstable_TSS_shH2AZ.bed' \
   --referencePoint TSS  \
   --skipZeros \
    --missingDataAsZero \
   -o "computeMatrix2plot_scores_H2AZocc.txt.gz" \
   -b 900 -a 900 --smartLabels -p 12

plotHeatmap -m "computeMatrix2plot_scores_H2AZocc.txt.gz" \
                 -out "../heatmap_scores_H2AZocc.png"

 ################### get occupancy WC ##############
 bwOccWC=/Volumes/PromisePegasus/_Research_/nucMACC_H2AZ/WC_nucMACC_out/RUN/01_MONO-NUCS_PROFILE

 computeMatrix reference-point \
  -S $bwOccWC"/"*".bw"  \
    -R $bw'/unstable_TSS_shH2AZ.bed' \
    --referencePoint TSS  \
    --skipZeros \
     --missingDataAsZero \
    -o "computeMatrix2plot_scores_WCocc.txt.gz" \
    -b 900 -a 900 --smartLabels -p 12

 plotHeatmap -m "computeMatrix2plot_scores_WCocc.txt.gz" \
                  -out "../heatmap_scores_WCocc.png"

################### get occupancy H2AZ ##############
bwOccshH2AZ=/Volumes/PromisePegasus/_Research_/nucMACC_H2AZ/shH2AZ_nucMACC_out/RUN/01_MONO-NUCS_PROFILE

computeMatrix reference-point \
-S $bwOccshH2AZ"/"*".bw"  \
-R $bw'/unstable_TSS_shH2AZ.bed' \
--referencePoint TSS  \
--skipZeros \
--missingDataAsZero \
-o "computeMatrix2plot_scores_shH2AZocc.txt.gz" \
-b 900 -a 900 --smartLabels -p 12

plotHeatmap -m "computeMatrix2plot_scores_shH2AZocc.txt.gz" \
           -out "../heatmap_scores_shH2AZocc.png"
