AnnoPath="/Users/admin/Annotation/GRCh38/"
AnalysisDir="/Volumes/PromisePegasus/_Research_/nucMACC_H2AZ"
genomINFO=$AnalysisDir"/script/addData/chrom_Sizes.txt"

mkdir -p $AnalysisDir/TSSplots/shH2AZ/data
cd $AnalysisDir/TSSplots/shH2AZ/data/

### HAZ
bw=/Volumes/PromisePegasus/_Research_/nucMACC_H2AZ/shH2AZ_nucMACC_out/RUN/01_MONO-NUCS_PROFILE

computeMatrix reference-point \
 -S $bw"/"*".bw"  \
   -R $AnnoPath'/nextflow/protein_coding.gtf' \
   --referencePoint TSS  \
   -o "computeMatrix2plot_mono.txt.gz" \
   -b 900 -a 900 --smartLabels -p 12

plotHeatmap -m "computeMatrix2plot_mono.txt.gz" \
        -out "../heatmap_cl7.png" \
        --kmeans 7

### hypo and hypo
mono=/Volumes/PromisePegasus/_Research_/nucMACC_H2AZ/shH2AZ_nucMACC_out/RUN/05_nucMACC/

bedtools merge -i $mono"/hyperAcc_monoNucs.bed" | bedtools genomecov -bga \
 -i - \
  -g $genomINFO \
  | sort -k1,1 -k2,2n  >"hyperAcc_monoNucs.bedgraph"

bedGraphToBigWig "hyperAcc_monoNucs.bedgraph" \
  $genomINFO \
  hyperAcc_monoNucs.bw


  bedtools merge -i $mono"/hypoAcc_monoNucs.bed" | bedtools genomecov -bga \
   -i - \
    -g $genomINFO \
    | sort -k1,1 -k2,2n  >"hypoAcc_monoNucs.bedgraph"

  bedGraphToBigWig "hypoAcc_monoNucs.bedgraph" \
    $genomINFO  \
    hypoAcc_monoNucs.bw

computeMatrix reference-point \
 -S *"monoNucs.bw"  \
   -R $AnnoPath'/nextflow/protein_coding.gtf' \
   --referencePoint TSS  \
   -o "computeMatrix2plot_mono_bed.txt.gz" \
   -b 900 -a 900 --smartLabels -p 5

plotHeatmap -m "computeMatrix2plot_mono_bed.txt.gz" \
           -out "../heatmap_mono_bed.png"

############## sub ########
sub=/Volumes/PromisePegasus/_Research_/nucMACC_H2AZ/shH2AZ_nucMACC_out/RUN/06_sub-nucMACC/


bedtools merge -i $sub"/unStable_subNucs.bed" | bedtools genomecov -bga \
 -i - \
  -g $genomINFO \
  | sort -k1,1 -k2,2n  >"unStable_subNucs.bedgraph"

bedGraphToBigWig "unStable_subNucs.bedgraph" \
  $genomINFO \
  unStable_subNucs.bw


  bedtools merge -i $sub"/nonCanonical_subNucs.bed" | bedtools genomecov -bga \
   -i - \
    -g $genomINFO \
    | sort -k1,1 -k2,2n  >"nonCanonical_subNucs.bedgraph"

  bedGraphToBigWig "nonCanonical_subNucs.bedgraph" \
    $genomINFO \
    nonCanonical_subNucs.bw

computeMatrix reference-point \
 -S *"subNucs.bw"  \
   -R $AnnoPath'/nextflow/protein_coding.gtf' \
   --referencePoint TSS  \
   -o "computeMatrix2plot_sub_bed.txt.gz" \
   -b 900 -a 900 --smartLabels -p 12

plotHeatmap -m "computeMatrix2plot_sub_bed.txt.gz" \
           -out "../heatmap_sub_bed.png" \
           --outFileSortedRegions "sortedRegions_Heatmap_sub_bed.txt" \
           --outFileNameMatrix "values_Heatmap_sub_bed.txt"




########### (sub-)nuccMACC scores #################

#remove overlaps
Rscript ../../../script/X_figure_plot/convert_ovrlp_bedGraph.R $mono"/nucMACC_scores.bedGraph" ./shH2AZ_nucMACC_scores_disjoint.bedGraph

bedGraphToBigWig shH2AZ_nucMACC_scores_disjoint.bedGraph \
../../../script/addData/chrom_Sizes.txt \
 shH2AZ_nucMACC_scores.bw

## sub
Rscript ../../../script/X_figure_plot/convert_ovrlp_bedGraph.R $sub"/sub-nucMACC_scores.bedGraph" ./shH2AZ_sub-nucMACC_scores_disjoint.bedGraph

bedGraphToBigWig shH2AZ_sub-nucMACC_scores_disjoint.bedGraph \
../../../script/addData/chrom_Sizes.txt \
 shH2AZ_sub-nucMACC_scores.bw


computeMatrix reference-point \
  -S *"scores.bw"  \
    -R $AnnoPath'/nextflow/protein_coding.gtf' \
    --referencePoint TSS  \
    --skipZeros \
     --missingDataAsZero \
    -o "computeMatrix2plot_scores.txt.gz" \
    -b 900 -a 900 --smartLabels -p 12


 plotHeatmap -m "computeMatrix2plot_scores.txt.gz" \
               -out "../heatmap_scores.png"
