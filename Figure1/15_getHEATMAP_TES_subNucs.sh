AnalysisDir=~/Analysis/R001_nucMacc

OutDir=$AnalysisDir"/manuscript_figures/Fig1/"
DataDir=$OutDir"/data"

########
mkdir -p $OutDir"/heatmaps/subnucs_TES/"
cd $OutDir"/heatmaps/subnucs_TES/"

###compute matrix
computeMatrix reference-point -S $DataDir"/bigwigs_subNucs//"*".bigwig" \
 -R $DataDir"/expression/dm3.refGene.gtf" \
 --referencePoint TES \
 -o "computeMatrix2plot_all.txt.gz" \
 --outFileNameMatrix "computeMatrix2txt_all.txt.gz" \
 --outFileSortedRegions "computeMatrix_geneList_all.bed" \
 -b 1500 -a 1500 --smartLabels -p 10


###
plotHeatmap -m "computeMatrix2plot_all.txt.gz" \
     -out 'DefaultHeatmap_all.png' \
     --outFileSortedRegions 'sortedRegions_Heatmap_all.txt'\
     --outFileNameMatrix 'values_Heatmap_all.txt'
