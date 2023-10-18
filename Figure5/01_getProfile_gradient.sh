AnalysisDir=~/Analysis/R001_nucMacc
OutDir=$AnalysisDir"/manuscript_figures/Fig2/ChIP_atlas"
DataDir=$OutDir"/data"


cd $DataDir
mkdir nucMACC_coverage

## hyper acc
hyperacc=$AnalysisDir/H4_ChIP_v2.1_featureCOUNTS/RUN/09_nucMACC/hyperAcc_monoNucs.bed

bedtools merge -i $hyperacc | bedtools genomecov -bga \
 -i - \
  -g ~/Annotation/Drosophila_melanogaster_UCSC_dm3/chr_sizes.txt \
  | sort -k1,1 -k2,2n  >nucMACC_coverage/hyperAcc_pos.bedgraph

bedGraphToBigWig nucMACC_coverage/hyperAcc_pos.bedgraph \
 ~/Annotation/Drosophila_melanogaster_UCSC_dm3/chr_sizes.txt \
  nucMACC_coverage/hyperAcc_pos.bigwig


## hypoacc
hypoacc=$AnalysisDir/H4_ChIP_v2.1_featureCOUNTS/RUN/09_nucMACC/hypoAcc_monoNucs.bed

bedtools merge -i $hypoacc | bedtools genomecov -bga \
 -i - \
  -g ~/Annotation/Drosophila_melanogaster_UCSC_dm3/chr_sizes.txt \
  | sort -k1,1 -k2,2n  >nucMACC_coverage/hypoAcc_pos.bedgraph

bedGraphToBigWig nucMACC_coverage/hypoAcc_pos.bedgraph \
 ~/Annotation/Drosophila_melanogaster_UCSC_dm3/chr_sizes.txt \
  nucMACC_coverage/hypoAcc_pos.bigwig

##bigwigs
H4bw=$AnalysisDir/H4_ChIP_v2.1_featureCOUNTS/RUN/05_MONO-NUCS_PROFILE/
hyperAccBw=$DataDir/nucMACC_coverage/hyperAcc_pos.bigwig
hypoAccBw=$DataDir/nucMACC_coverage/hypoAcc_pos.bigwig

nucMACCbw=$AnalysisDir"/manuscript_figures/Fig1/data/bigwigs_monoNucs/H4_nucMACC_scores.bigwig"

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

 computeMatrix reference-point -S $nucMACCbw \
      -R $DataDir"/M1BP/GSE101554_S2_M1BP_ChIP_peaks.bed" \
      --referencePoint center \
      -o "computeMatrix2plot_nucMACCscore.txt.gz" \
      --outFileSortedRegions "computeMatrix_peakList_nucMACCscore.bed" \
      -b 1500 -a 1500 --smartLabels -p 10

plotHeatmap -m "computeMatrix2plot_nucMACCscore.txt.gz" \
      -out 'nucMACCscore_Heatmap.png' \
      --outFileSortedRegions 'sortedRegions_Heatmap_nucMACCscore.txt' \
      --outFileNameMatrix 'values_Heatmap_nucMACCscores.txt'


###########################
computeMatrix reference-point -S $hyperAccBw $hypoAccBw \
    -R $DataDir"/M1BP/GSE101554_S2_M1BP_ChIP_peaks.bed" \
    --referencePoint center \
    -o "computeMatrix2plot_nucMACC.txt.gz" \
    --outFileNameMatrix "computeMatrix2txt_nucMACC.txt.gz" \
    --outFileSortedRegions "computeMatrix_peakList_nucMACC.bed" \
    -b 1500 -a 1500 --smartLabels -p 10


plotHeatmap -m "computeMatrix2plot_nucMACC.txt.gz" \
-out 'nucMACC_Heatmap.png' \
--outFileSortedRegions 'sortedRegions_Heatmap.txt' \
--outFileNameMatrix 'values_Heatmap.txt'



############################


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

for entry in *.gff3; do

echo $entry
name=$(echo $entry| cut -d'.' -f 1)
awk -v OFS='\t' 'BEGIN{n=5} NR>5 {print $1,$4,$5,$9,$6,$7 }' $entry |  sed 's/^/chr/' >$name".bed"

done


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

cd $DataDir"/MSL-1/"

for entry in *.gff3; do

echo $entry
name=$(echo $entry| cut -d'.' -f 1)
awk -v OFS='\t' 'BEGIN{n=5} NR>5 {print $1,$4,$5,$9,$6,$7 }' $entry |  sed 's/^/chr/' >$name".bed"


done


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

  ##########   Heatmap  nucMACC scores ##########################

   computeMatrix reference-point -S $nucMACCbw \
        -R $DataDir"/MSL-1/MSL-1.bed" \
        --referencePoint center \
        -o "computeMatrix2plot_nucMACCscore.txt.gz" \
        --outFileSortedRegions "computeMatrix_peakList_nucMACCscore.bed" \
        -b 1500 -a 1500 --smartLabels -p 10

  plotHeatmap -m "computeMatrix2plot_nucMACCscore.txt.gz" \
        -out 'nucMACCscore_Heatmap.png' \
        --outFileSortedRegions 'sortedRegions_Heatmap_nucMACCscore.txt' \
        --outFileNameMatrix 'values_Heatmap_nucMACCscores.txt'


  ###########################
  computeMatrix reference-point -S $hyperAccBw $hypoAccBw \
      -R $DataDir"/MSL-1/MSL-1.bed" \
      --referencePoint center \
      -o "computeMatrix2plot_nucMACC.txt.gz" \
      --outFileNameMatrix "computeMatrix2txt_nucMACC.txt.gz" \
      --outFileSortedRegions "computeMatrix_peakList_nucMACC.bed" \
      -b 1500 -a 1500 --smartLabels -p 10


  plotHeatmap -m "computeMatrix2plot_nucMACC.txt.gz" \
  -out 'nucMACC_Heatmap.png' \
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


  ######  HP2

  cd $DataDir"/HP2/"

  for entry in *.gff3; do

  echo $entry
  name=$(echo $entry| cut -d'.' -f 1)
  awk -v OFS='\t' 'BEGIN{n=5} NR>5 {print $1,$4,$5,$9,$6,$7 }' $entry |  sed 's/^/chr/' >$name".bed"


  done


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
