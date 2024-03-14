#!/usr/bin/env bash

projectDir=/home/$USER/nucMACC_Paper/data
AnalysisDir=/home/$USER/nucMACC_Paper/data/R001_nucMacc
OutDir=$AnalysisDir"manuscript_figures/Fig5/ChIP_atlas"

DataDir=$OutDir"/data"

# create new directory for the analysis
cd $DataDir
mkdir nucMACC_coverage

# positions of hyper-accesssible nucleosomes
hyperacc=$AnalysisDir/H4_ChIP/RUN/09_nucMACC/hyperAcc_monoNucs.bed

bedtools merge -i $hyperacc | bedtools genomecov -bga \
      -i - -g $DataDir/data/Drosophila_melanogaster_UCSC_dm3/chr_sizes.txt \
      | sort -k1,1 -k2,2n  > $DataDir/nucMACC_coverage/hyperAcc_pos.bedgraph

bedGraphToBigWig $DataDir/nucMACC_coverage/hyperAcc_pos.bedgraph \
      $projectDir/data/Drosophila_melanogaster_UCSC_dm3/chr_sizes.txt \
      $DataDir/nucMACC_coverage/hyperAcc_pos.bigwig

# positions of hypo-accesssible nucleosomes
hypoacc=$AnalysisDir/H4_ChIP/RUN/09_nucMACC/hypoAcc_monoNucs.bed

bedtools merge -i $hypoacc | bedtools genomecov -bga \
      -i - \
      -g $DataDir/data/Drosophila_melanogaster_UCSC_dm3/chr_sizes.txt \
      | sort -k1,1 -k2,2n  > $DataDir/nucMACC_coverage/hypoAcc_pos.bedgraph

bedGraphToBigWig $DataDir/nucMACC_coverage/hypoAcc_pos.bedgraph \
      $DataDir/data/Drosophila_melanogaster_UCSC_dm3/chr_sizes.txt \
      $DataDir/nucMACC_coverage/hypoAcc_pos.bigwig

# bigwigs of H4 MNase samples
H4bw=$AnalysisDir/H4_ChIP/RUN/05_MONO-NUCS_PROFILE/
hyperAccBw=$DataDir/nucMACC_coverage/hyperAcc_pos.bigwig
hypoAccBw=$DataDir/nucMACC_coverage/hypoAcc_pos.bigwig

nucMACCbw=$AnalysisDir"/manuscript_figures/Fig1/data/bigwigs_monoNucs/H4_nucMACC_scores.bigwig"

#1. AbdA
mkdir -p $OutDir"/profiles/AbdA"
cd $OutDir"/profiles/AbdA"

# Heatmap of the H4 signal on AbdA binding sites
computeMatrix reference-point -S $H4bw"/H4"*".bw" \
      -R $DataDir"/AbdA_HA/GSE101554_S2-AbdA_HA_ChIP_peaks.bed" \
      --referencePoint center \
      -o $DataDir"/computeMatrix2plot_H4.txt.gz" \
      --outFileNameMatrix $DataDir"/computeMatrix2txt_H4.txt.gz" \
      --outFileSortedRegions $DataDir"/computeMatrix_peakList_H4.bed" \
      -b 1500 -a 1500 --smartLabels -p 10

plotProfile -m $DataDir"/computeMatrix2plot_H4.txt.gz" \
      -out $DataDir'/DefaultHeatmap.png' \
      --outFileSortedRegions $DataDir'/sortedRegions_Profile.txt'\
      --outFileNameData $DataDir'/values_Profile.txt' \
      --perGroup

#2. M1BP
mkdir -p $OutDir"/profiles/M1BP"
cd $OutDir"/profiles/M1BP"

# Heatmap of the H4 signal on M1BP TF binding sites
computeMatrix reference-point -S $H4bw"/H4"*".bw" \
      -R $DataDir"/M1BP/GSE101554_S2_M1BP_ChIP_peaks.bed" \
      --referencePoint center \
      -o $DataDir"/computeMatrix2plot_H4.txt.gz" \
      --outFileNameMatrix $DataDir"/computeMatrix2txt_H4.txt.gz" \
      --outFileSortedRegions $DataDir"/computeMatrix_peakList_H4.bed" \
      -b 1500 -a 1500 --smartLabels -p 10

plotProfile -m $DataDir"/computeMatrix2plot_H4.txt.gz" \
      -out $DataDir'/DefaultHeatmap.png' \
      --outFileSortedRegions $DataDir'/sortedRegions_Profile.txt'\
      --outFileNameData $DataDir'/values_Profile.txt' \
      --perGroup

# Heatmap nucMACC scores on M1BP TF binding sites  
computeMatrix reference-point -S $nucMACCbw \
      -R $DataDir"/M1BP/GSE101554_S2_M1BP_ChIP_peaks.bed" \
      --referencePoint center \
      -o $DataDir"/computeMatrix2plot_nucMACCscore.txt.gz" \
      --outFileSortedRegions $DataDir"/computeMatrix_peakList_nucMACCscore.bed" \
      -b 1500 -a 1500 --smartLabels -p 10

plotHeatmap -m $DataDir"/computeMatrix2plot_nucMACCscore.txt.gz" \
      -out $DataDir'/nucMACCscore_Heatmap.png' \
      --outFileSortedRegions $DataDir'/sortedRegions_Heatmap_nucMACCscore.txt' \
      --outFileNameMatrix $DataDir'/values_Heatmap_nucMACCscores.txt'

# Heatmap nucMACC on M1BP TF binding sites
# for hypo and hyper accessbile nucleosome positions 
computeMatrix reference-point -S $hyperAccBw $hypoAccBw \
      -R $DataDir"/M1BP/GSE101554_S2_M1BP_ChIP_peaks.bed" \
      --referencePoint center \
      -o $DataDir"/computeMatrix2plot_nucMACC.txt.gz" \
      --outFileNameMatrix $DataDir"/computeMatrix2txt_nucMACC.txt.gz" \
      --outFileSortedRegions $DataDir"/computeMatrix_peakList_nucMACC.bed" \
      -b 1500 -a 1500 --smartLabels -p 10

plotHeatmap -m $DataDir"/computeMatrix2plot_nucMACC.txt.gz" \
-out $DataDir'/nucMACC_Heatmap.png' \
--outFileSortedRegions $DataDir'/sortedRegions_Heatmap.txt' \
--outFileNameMatrix $DataDir'/values_Heatmap.txt'

#3. ZKSCAN3
mkdir -p $OutDir"/profiles/ZKSCAN3"
cd $OutDir"/profiles/ZKSCAN3"

/home/$USER/tools/liftOver $DataDir"/ZKSCAN3/GSE149116_ZKSCAN3_S2_peaks.bed" \
      $DataDir"/dm6ToDm3.over.chain.gz" \
      $DataDir"/ZKSCAN3/GSE149116_ZKSCAN3_S2_peaks_dm3.bed" \
      $DataDir"/ZKSCAN3/unMapped"

# Heatmap nucMACC on ZKSCAN3 binding sites
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

#4. SuVar3-9

cd $DataDir"/SUvar3-9/"

for entry in *.gff3; do

echo $entry
name=$(echo $entry| cut -d'.' -f 1)
awk -v OFS='\t' 'BEGIN{n=5} NR>5 {print $1,$4,$5,$9,$6,$7 }' $entry |  sed 's/^/chr/' >$name".bed"
done

mkdir -p $OutDir"/profiles/SUvar3-9"
cd $OutDir"/profiles/SUvar3-9"

# Heatmap nucMACC on SuVar3-9 binding sites
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

#5. MSL-1

cd $DataDir"/MSL-1/"

for entry in *.gff3; do
echo $entry
name=$(echo $entry| cut -d'.' -f 1)
awk -v OFS='\t' 'BEGIN{n=5} NR>5 {print $1,$4,$5,$9,$6,$7 }' $entry |  sed 's/^/chr/' >$name".bed"
done

mkdir -p $OutDir"/profiles/MSL-1"
cd $OutDir"/profiles/MSL-1"

# Heatmap nucMACC on MSL-1 binding sites
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

# Heatmap nucMACC on MSL-1 binding sites
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

# Heatmap nucMACC on MSL-1 binding sites
# for hypo and hyper accessbile nucleosome positions 
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

#6.  M1BP study

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

#7. HP2

cd $DataDir"/HP2/"

for entry in *.gff3; do
echo $entry
name=$(echo $entry| cut -d'.' -f 1)
awk -v OFS='\t' 'BEGIN{n=5} NR>5 {print $1,$4,$5,$9,$6,$7 }' $entry |  sed 's/^/chr/' >$name".bed"
done

mkdir -p $OutDir"/profiles/HP2"
cd $OutDir"/profiles/HP2"

# Heatmap nucMACC on HP2 binding sites
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
