#!/usr/bin/env bash

AnalysisDir=/home/$USER/nucMACC_Paper

OutDir=$AnalysisDir"/manuscript_figures/Fig1/"
mkdir -p $OutDir

DataDir=$OutDir"/data"
mkdir -p $DataDir"/bigwigs_subNucs/"

cd $DataDir"/bigwigs_subNucs/"

OUT=$DataDir"/bigwigs_subNucs"

## convert bedgraph to bigwig
#H3
#remove overlaps
Rscript $AnalysisDir/script/manuscript_figures/Fig1/convert_ovrlp_bedGraph.R $AnalysisDir/H3_ChIP/RUN/10_sub-nucMACC/sub-nucMACC_scores.bedGraph $OUT/H3_sub-nucMACC_scores_disjoint.bedGraph

bedGraphToBigWig H3_sub-nucMACC_scores_disjoint.bedGraph \
$AnalysisDir/data/Drosophila_melanogaster_UCSC_dm3/chr_sizes.txt \
$AnalysisDir/data/H3_sub-nucMACC_scores.bigwig

#H4
#remove overlaps
Rscript $AnalysisDir/script/manuscript_figures/Fig1/convert_ovrlp_bedGraph.R $AnalysisDir/data/H4_ChIP/RUN/10_sub-nucMACC/sub-nucMACC_scores.bedGraph $AnalysisDir/data/H4_ChIP/H4_sub-nucMACC_scores_disjoint.bedGraph

bedGraphToBigWig H4_sub-nucMACC_scores_disjoint.bedGraph \
$AnalysisDir/data/Drosophila_melanogaster_UCSC_dm3/chr_sizes.txt \
$AnalysisDir/data/H4_ChIP/H4_sub-nucMACC_scores.bigwig

rm $OUT/*.bedGraph

########
mkdir -p $OutDir"/heatmaps/subnucs/"
cd $OutDir"/heatmaps/subnucs/"

###compute matrix for the heatmap
computeMatrix reference-point -S $DataDir"/bigwigs_subNucs//"*".bigwig" \
 -R $DataDir"/expression/dm3.refGene.gtf" \
 --referencePoint TSS \
 -o "computeMatrix2plot_all.txt.gz" \
 --outFileNameMatrix "computeMatrix2txt_all.txt.gz" \
 --outFileSortedRegions "computeMatrix_geneList_all.bed" \
 -b 1500 -a 1500 --smartLabels -p 10

### plot heatmap
plotHeatmap -m $OUT"/computeMatrix2plot_all.txt.gz" \
     -out $OUT'/DefaultHeatmap_all.png' \
     --outFileSortedRegions $OUT'/sortedRegions_Heatmap_all.txt'\
     --outFileNameMatrix $OUT'/values_Heatmap_all.txt'
