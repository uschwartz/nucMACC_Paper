AnalysisDir=~/Analysis/R001_nucMacc

OutDir=$AnalysisDir"/manuscript_figures/Fig1/"
DataDir=$OutDir"/data"


########
mkdir -p $OutDir"/heatmaps/monoNucs_geneRegion/"
cd $OutDir"/heatmaps/monoNucs_geneRegion/"

######## scale regions ########################

### BED protein coding genes #############

computeMatrix scale-regions -S $DataDir"/bigwigs_monoNucs//"*".bigwig" \
 -R $DataDir"/expression/dm3.refGene_coding.bed" \
 -o "computeMatrix2plot_bed.txt.gz" \
 --outFileNameMatrix "computeMatrix2txt_bed.txt.gz" \
 --outFileSortedRegions "computeMatrix_geneList_bed.bed" \
 -b 500 -a 500 --smartLabels -p 10


###
plotHeatmap -m "computeMatrix2plot_bed.txt.gz" \
     -out 'DefaultHeatmap_bed.png' \
     --outFileSortedRegions 'sortedRegions_Heatmap_bed.txt'\
     --outFileNameMatrix 'values_Heatmap_bed.txt'
