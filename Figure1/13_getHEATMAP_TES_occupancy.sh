AnalysisDir=~/Analysis/R001_nucMacc
DataDir=$OutDir"/data"

OutDir=$AnalysisDir"/manuscript_figures/Fig1/"
mkdir -p $OutDir"/heatmaps/monoNucs_occu_TES/"
cd $OutDir"/heatmaps/monoNucs_occu_TES/"
######## H4
H4pooled=$AnalysisDir/H4_ChIP_v2.1_featureCOUNTS/RUN/05_MONO-NUCS_PROFILE/pooled_monoNucs_profile.bw

###compute matrix
computeMatrix reference-point -S $H4pooled \
 -R $DataDir"/expression/dm3.refGene.gtf" \
 --referencePoint TES \
 -o "computeMatrix2plot_H4.txt.gz" \
 --outFileNameMatrix "computeMatrix2txt_H4.txt.gz" \
 --outFileSortedRegions "computeMatrix_geneList_H4.bed" \
 -b 1500 -a 1500 --smartLabels -p 10


###
plotHeatmap -m "computeMatrix2plot_H4.txt.gz" \
     -out 'DefaultHeatmap_H4.png' \
     --outFileSortedRegions 'sortedRegions_Heatmap_H4.txt'\
     --outFileNameMatrix 'values_Heatmap_H4.txt'
