AnalysisDir=~/Analysis/R001_nucMacc
OutDir=$AnalysisDir"/manuscript_figures/Fig2/"


mkdir -p $OutDir"/DNaseHS/"
cd $OutDir"/DNaseHS/"


######## DNaseHS
DNase_bw=$OutDir"/data/"
beds=$AnalysisDir"/H4_ChIP_v2.1_featureCOUNTS/RUN/09_nucMACC"

###compute matrix
computeMatrix reference-point -S $DNase_bw"/DNase_S2DHS-seq_r1.bw" \
 -R  $beds"/hyperAcc_monoNucs.bed" $beds"/hypoAcc_monoNucs.bed" \
 --referencePoint center \
 -o "computeMatrix2plot_DNaseHS.txt.gz" \
 --outFileNameMatrix "computeMatrix2txt_H4.txt.gz" \
 -b 3000 -a 3000 --smartLabels -p 10


###
plotProfile -m "computeMatrix2plot_DNaseHS.txt.gz" \
     -out 'DefaultProfile.pdf' \
     --outFileNameData 'values_Profile_DNaseHS.txt'
