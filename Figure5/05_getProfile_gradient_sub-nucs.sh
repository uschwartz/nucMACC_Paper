AnalysisDir=~/Analysis/R001_nucMacc
OutDir=$AnalysisDir"/manuscript_figures/Fig3/ChIP_atlas"
DataDir=$AnalysisDir"/manuscript_figures/Fig2/ChIP_atlas/data"

cd $OutDir/data
mkdir sub-nucMACC_coverage

## unstable
unstable=$AnalysisDir/H4_ChIP_v2.1_featureCOUNTS/RUN/10_sub-nucMACC/unStable_subNucs.bed

bedtools merge -i $unstable | bedtools genomecov -bga \
 -i - \
  -g ~/Annotation/Drosophila_melanogaster_UCSC_dm3/chr_sizes.txt \
  | sort -k1,1 -k2,2n  >sub-nucMACC_coverage/unstable_pos.bedgraph

bedGraphToBigWig sub-nucMACC_coverage/unstable_pos.bedgraph \
 ~/Annotation/Drosophila_melanogaster_UCSC_dm3/chr_sizes.txt \
  sub-nucMACC_coverage/unstable_pos.bigwig


## nonCanonical
noncan=$AnalysisDir/H4_ChIP_v2.1_featureCOUNTS/RUN/10_sub-nucMACC/nonCanonical_subNucs.bed

bedtools merge -i $noncan | bedtools genomecov -bga \
 -i - \
  -g ~/Annotation/Drosophila_melanogaster_UCSC_dm3/chr_sizes.txt \
  | sort -k1,1 -k2,2n  >sub-nucMACC_coverage/noncan_pos.bedgraph

bedGraphToBigWig sub-nucMACC_coverage/noncan_pos.bedgraph \
 ~/Annotation/Drosophila_melanogaster_UCSC_dm3/chr_sizes.txt \
  sub-nucMACC_coverage/noncan_pos.bigwig

##bigwigs
H4bw=$AnalysisDir/H4_ChIP_v2.1_featureCOUNTS/RUN/06_SUB-NUCS_PROFILE/
unstableBw=$OutDir/data/sub-nucMACC_coverage/unstable_pos.bigwig
noncanBw=$OutDir/data/sub-nucMACC_coverage/noncan_pos.bigwig

subnucMACCbw=$AnalysisDir"/manuscript_figures/Fig1/data/bigwigs_subNucs/H4_sub-nucMACC_scores.bigwig"


######## AbdA
mkdir -p $OutDir"/profiles/AbdA"
cd $OutDir"/profiles/AbdA"

###compute matrix
computeMatrix reference-point -S $H4bw"/H4"*".bw" \
 -R $DataDir"/AbdA_HA/GSE101554_S2-AbdA_HA_ChIP_peaks.bed" \
 --referencePoint center \
 -o "computeMatrix2plot_H4.txt.gz" \
 --outFileNameMatrix "computeMatrix2txt_H4.txt.gz" \
 --outFileSortedRegions "computeMatrix_peakList_H4.bed" \
 -b 1500 -a 1500 --smartLabels -p 10



plotProfile -m "computeMatrix2plot_H4.txt.gz" \
     -out 'DefaultHeatmap.png' \
     --outFileSortedRegions 'sortedRegions_Profile.txt'\
     --outFileNameData 'values_Profile.txt' \
     --perGroup

########## M1BP
mkdir -p $OutDir"/profiles/M1BP"
cd $OutDir"/profiles/M1BP"


computeMatrix reference-point -S $H4bw"/H4"*".bw" \
      -R $DataDir"/M1BP/GSE101554_S2_M1BP_ChIP_peaks.bed" \
      --referencePoint center \
      -o "computeMatrix2plot_H4.txt.gz" \
      --outFileNameMatrix "computeMatrix2txt_H4.txt.gz" \
      --outFileSortedRegions "computeMatrix_peakList_H4.bed" \
      -b 1500 -a 1500 --smartLabels -p 10



plotProfile -m "computeMatrix2plot_H4.txt.gz" \
  -out 'DefaultHeatmap.png' \
  --outFileSortedRegions 'sortedRegions_Profile.txt'\
  --outFileNameData 'values_Profile.txt' \
  --perGroup



  ##########   Heatmap  nucMACC scores ##########################

   computeMatrix reference-point -S $subnucMACCbw \
        -R $DataDir"/M1BP/GSE101554_S2_M1BP_ChIP_peaks.bed" \
        --referencePoint center \
        -o "computeMatrix2plot_subnucMACCscore.txt.gz" \
        --outFileSortedRegions "computeMatrix_peakList_subnucMACCscore.bed" \
        -b 1500 -a 1500 --smartLabels -p 10

  plotHeatmap -m "computeMatrix2plot_subnucMACCscore.txt.gz" \
        -out 'subnucMACCscore_Heatmap.png' \
        --outFileSortedRegions 'sortedRegions_Heatmap_subnucMACCscore.txt' \
        --outFileNameMatrix 'values_Heatmap_subnucMACCscores.txt'


  ###########################
  computeMatrix reference-point -S $unstableBw $noncanBw \
      -R $DataDir"/M1BP/GSE101554_S2_M1BP_ChIP_peaks.bed" \
      --referencePoint center \
      -o "computeMatrix2plot_subnucMACC.txt.gz" \
      --outFileNameMatrix "computeMatrix2txt_subnucMACC.txt.gz" \
      --outFileSortedRegions "computeMatrix_peakList_subnucMACC.bed" \
      -b 1500 -a 1500 --smartLabels -p 10


  plotHeatmap -m "computeMatrix2plot_subnucMACC.txt.gz" \
  -out 'subnucMACC_Heatmap.png' \
  --outFileSortedRegions 'sortedRegions_Heatmap.txt' \
  --outFileNameMatrix 'values_Heatmap.txt'


  ############ M1BP study

  mkdir -p $OutDir"/profiles/study"
  cd $OutDir"/profiles/study"


  computeMatrix reference-point -S $H4bw"/H4"*".bw" \
        -R $DataDir"/study/"*".bed" \
        --referencePoint center \
        -o "computeMatrix2plot_H4.txt.gz" \
        --outFileNameMatrix "computeMatrix2txt_H4.txt.gz" \
        --outFileSortedRegions "computeMatrix_peakList_H4.bed" \
        -b 1500 -a 1500 --smartLabels -p 10



  plotProfile -m "computeMatrix2plot_H4.txt.gz" \
    -out 'DefaultProfile.png' \
    --outFileSortedRegions 'sortedRegions_Profile.txt'\
    --outFileNameData 'values_Profile.txt' \
    --perGroup



  ######################





########## ZKSCAN3
mkdir -p $OutDir"/profiles/ZKSCAN3"
cd $OutDir"/profiles/ZKSCAN3"

/Users/admin/Tools/ucsctools/liftOver $DataDir"/ZKSCAN3/GSE149116_ZKSCAN3_S2_peaks.bed" \
 $DataDir"/dm6ToDm3.over.chain.gz" \
 $DataDir"/ZKSCAN3/GSE149116_ZKSCAN3_S2_peaks_dm3.bed" \
 $DataDir"/ZKSCAN3/unMapped"


 computeMatrix reference-point -S $H4bw"/H4"*".bw" \
       -R $DataDir"/ZKSCAN3/GSE149116_ZKSCAN3_S2_peaks_dm3.bed" \
       --referencePoint center \
       -o "computeMatrix2plot_H4.txt.gz" \
       --outFileNameMatrix "computeMatrix2txt_H4.txt.gz" \
       --outFileSortedRegions "computeMatrix_peakList_H4.bed" \
       -b 1500 -a 1500 --smartLabels -p 10



 plotProfile -m "computeMatrix2plot_H4.txt.gz" \
   -out 'DefaultHeatmap.png' \
   --outFileSortedRegions 'sortedRegions_Profile.txt'\
   --outFileNameData 'values_Profile.txt' \
   --perGroup

######  SuVar3-9

cd $DataDir"/SUvar3-9/"



mkdir -p $OutDir"/profiles/SUvar3-9"
cd $OutDir"/profiles/SUvar3-9"


computeMatrix reference-point -S $H4bw"/H4"*".bw" \
 -R $DataDir"/SUvar3-9/SUvar3-9.bed" \
 --referencePoint center \
 -o "computeMatrix2plot_H4.txt.gz" \
 --outFileNameMatrix "computeMatrix2txt_H4.txt.gz" \
 --outFileSortedRegions "computeMatrix_peakList_H4.bed" \
 -b 1500 -a 1500 --smartLabels -p 10


 plotProfile -m "computeMatrix2plot_H4.txt.gz" \
   -out 'DefaultHeatmap.png' \
   --outFileSortedRegions 'sortedRegions_Profile.txt'\
   --outFileNameData 'values_Profile.txt' \
   --perGroup


######  MSL-1

mkdir -p $OutDir"/profiles/MSL-1"
cd $OutDir"/profiles/MSL-1"


computeMatrix reference-point -S $H4bw"/H4"*".bw" \
      -R $DataDir"/MSL-1/MSL-1.bed" \
      --referencePoint center \
      -o "computeMatrix2plot_H4.txt.gz" \
      --outFileNameMatrix "computeMatrix2txt_H4.txt.gz" \
      --outFileSortedRegions "computeMatrix_peakList_H4.bed" \
      -b 1500 -a 1500 --smartLabels -p 10


plotProfile -m "computeMatrix2plot_H4.txt.gz" \
  -out 'DefaultHeatmap.png' \
  --outFileSortedRegions 'sortedRegions_Profile.txt'\
  --outFileNameData 'values_Profile.txt' \
  --perGroup

  ######  HP2


  mkdir -p $OutDir"/profiles/HP2"
  cd $OutDir"/profiles/HP2"


  computeMatrix reference-point -S $H4bw"/H4"*".bw" \
        -R $DataDir"/HP2/HP2.bed" \
        --referencePoint center \
        -o "computeMatrix2plot_H4.txt.gz" \
        --outFileNameMatrix "computeMatrix2txt_H4.txt.gz" \
        --outFileSortedRegions "computeMatrix_peakList_H4.bed" \
        -b 1500 -a 1500 --smartLabels -p 10


  plotProfile -m "computeMatrix2plot_H4.txt.gz" \
    -out 'DefaultHeatmap.png' \
    --outFileSortedRegions 'sortedRegions_Profile.txt'\
    --outFileNameData 'values_Profile.txt' \
    --perGroup
