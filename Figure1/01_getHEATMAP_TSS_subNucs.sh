AnalysisDir=~/Analysis/R001_nucMacc

OutDir=$AnalysisDir"/manuscript_figures/Fig1/"
mkdir -p $OutDir

DataDir=$OutDir"/data"
mkdir -p $DataDir"/bigwigs_subNucs/"

cd $DataDir"/bigwigs_subNucs/"

## convert bedgraph to bigwig
#H3
#remove overlaps
Rscript $AnalysisDir/script/manuscript_figures/Fig1/convert_ovrlp_bedGraph.R $AnalysisDir/H3_ChIP_v2_featureCOUNTS/RUN/10_sub-nucMACC/sub-nucMACC_scores.bedGraph ./H3_sub-nucMACC_scores_disjoint.bedGraph

bedGraphToBigWig H3_sub-nucMACC_scores_disjoint.bedGraph \
~/Annotation/Drosophila_melanogaster_UCSC_dm3/chr_sizes.txt \
 H3_sub-nucMACC_scores.bigwig

#H4
#remove overlaps
Rscript $AnalysisDir/script/manuscript_figures/Fig1/convert_ovrlp_bedGraph.R $AnalysisDir/H4_ChIP_v2.1_featureCOUNTS/RUN/10_sub-nucMACC/sub-nucMACC_scores.bedGraph ./H4_sub-nucMACC_scores_disjoint.bedGraph

bedGraphToBigWig H4_sub-nucMACC_scores_disjoint.bedGraph \
~/Annotation/Drosophila_melanogaster_UCSC_dm3/chr_sizes.txt \
 H4_sub-nucMACC_scores.bigwig


rm *.bedGraph
########
mkdir -p $OutDir"/heatmaps/subnucs/"
cd $OutDir"/heatmaps/subnucs/"

###compute matrix
computeMatrix reference-point -S $DataDir"/bigwigs_subNucs//"*".bigwig" \
 -R $DataDir"/expression/dm3.refGene.gtf" \
 --referencePoint TSS \
 -o "computeMatrix2plot_all.txt.gz" \
 --outFileNameMatrix "computeMatrix2txt_all.txt.gz" \
 --outFileSortedRegions "computeMatrix_geneList_all.bed" \
 -b 1500 -a 1500 --smartLabels -p 10


###
plotHeatmap -m "computeMatrix2plot_all.txt.gz" \
     -out 'DefaultHeatmap_all.png' \
     --outFileSortedRegions 'sortedRegions_Heatmap_all.txt'\
     --outFileNameMatrix 'values_Heatmap_all.txt'
