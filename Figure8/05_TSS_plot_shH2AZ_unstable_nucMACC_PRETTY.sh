AnalysisDir="/Volumes/PromisePegasus/_Research_/nucMACC_H2AZ"
bw=$AnalysisDir"/TSSplots/"

mkdir -p $AnalysisDir/TSSplots/shH2AZ_unstable_PRETTY/data
cd $AnalysisDir/TSSplots/shH2AZ_unstable_PRETTY/data


########### nuccMACC scores #################

computeMatrix reference-point \
  -S $bw"/shH2AZ/data/shH2AZ_nucMACC_scores.bw" \
    $bw"/WC/data/WC_nucMACC_scores.bw" \
    -R 'unstable_TSS_MANsort.bed' \
    --referencePoint TSS  \
    --skipZeros \
    --missingDataAsZero \
    -o "computeMatrix2plot_nucMACC_scores.txt.gz" \
    -b 900 -a 900 --smartLabels -p 12


 plotHeatmap -m "computeMatrix2plot_nucMACC_scores.txt.gz" \
               -out "../heatmap_nucMACC_scores.png" \
               --sortRegions "keep" \
               --outFileNameMatrix "values_Heatmap_nucMACC.txt"
