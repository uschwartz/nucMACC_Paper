AnalysisDir='/Users/admin/Analysis/R001_nucMacc/manuscript_figures/Fig7'
OutDir=$AnalysisDir"/in_vitro_Kubik/"

mkdir -p $OutDir

cd $AnalysisDir
############################ SGD ############################

mkdir $OutDir"/SGD/"

computeMatrix reference-point -S data/in_vitro_bigwigs/subNucs/SGD-Replicate3_subNucs_profile.bw \
data/in_vitro_bigwigs/monoNucs/SGD-Replicate3_monoNucs_profile.bw \
 -R data/Kubik/TSS_NDR_Kubik.bed data/Kubik/TSS_unstable_Kubik.bed \
 --referencePoint TSS \
 -o $OutDir"/SGD/cmpt_mx.gz" \
 -b 1000 -a 1000 --smartLabels -p 10

plotProfile -m $OutDir"/SGD/cmpt_mx.gz" \
      -out $OutDir'/SGD/DefaultHeatmap_median.png' \
      --outFileNameData $OutDir'/SGD/values_Profile_median.txt' \
      --averageType "median"

plotProfile -m $OutDir"/SGD/cmpt_mx.gz" \
            -out $OutDir'/SGD/DefaultHeatmap.png' \
            --outFileNameData $OutDir'/SGD/values_Profile.txt'

plotHeatmap -m $OutDir"/SGD/cmpt_mx.gz" \
      -out $OutDir'/SGD/DefaultHeatmap_HEAT.png' \
      --outFileNameMatrix $OutDir'/SGD/values_Heatmap.txt'



############################ INO80 ############################

mkdir $OutDir"/INO80/"

computeMatrix reference-point -S data/in_vitro_bigwigs/subNucs/INO80_subNucs_profile.bw \
data/in_vitro_bigwigs/monoNucs/INO80_monoNucs_profile.bw \
-R data/Kubik/TSS_NDR_Kubik.bed data/Kubik/TSS_unstable_Kubik.bed  \
--referencePoint TSS \
-o $OutDir"/INO80/cmpt_mx.gz" \
-b 1000 -a 1000 --smartLabels -p 10

plotProfile -m $OutDir"/INO80/cmpt_mx.gz" \
-out $OutDir'/INO80/DefaultHeatmap_median.png' \
--outFileNameData $OutDir'/INO80/values_Profile_median.txt' \
--averageType "median"

plotProfile -m $OutDir"/INO80/cmpt_mx.gz" \
      -out $OutDir'/INO80/DefaultHeatmap.png' \
      --outFileNameData $OutDir'/INO80/values_Profile.txt'

plotHeatmap -m $OutDir"/INO80/cmpt_mx.gz" \
-out $OutDir'/INO80/DefaultHeatmap_HEAT.png' \
--outFileNameMatrix $OutDir'/INO80/values_Heatmap.txt'


############################ INO80_ISW2_RSC_ISW1a ############################

mkdir $OutDir"/INO80_ISW2_RSC_ISW1a/"

computeMatrix reference-point -S data/in_vitro_bigwigs/subNucs/INO80_ISW2_RSC_ISW1a_subNucs_profile.bw \
data/in_vitro_bigwigs/monoNucs/INO80_ISW2_RSC_ISW1a_monoNucs_profile.bw \
-R data/Kubik/TSS_NDR_Kubik.bed data/Kubik/TSS_unstable_Kubik.bed  \
--referencePoint TSS \
-o $OutDir"/INO80_ISW2_RSC_ISW1a/cmpt_mx.gz" \
-b 1000 -a 1000 --smartLabels -p 10

plotProfile -m $OutDir"/INO80_ISW2_RSC_ISW1a/cmpt_mx.gz" \
-out $OutDir'/INO80_ISW2_RSC_ISW1a/DefaultHeatmap_median.png' \
--outFileNameData $OutDir'/INO80_ISW2_RSC_ISW1a/values_Profile_median.txt' \
--averageType "median"

plotProfile -m $OutDir"/INO80_ISW2_RSC_ISW1a/cmpt_mx.gz" \
      -out $OutDir'/INO80_ISW2_RSC_ISW1a/DefaultHeatmap.png' \
      --outFileNameData $OutDir'/INO80_ISW2_RSC_ISW1a/values_Profile.txt'

plotHeatmap -m $OutDir"/INO80_ISW2_RSC_ISW1a/cmpt_mx.gz" \
-out $OutDir'/INO80_ISW2_RSC_ISW1a/DefaultHeatmap_HEAT.png' \
--outFileNameMatrix $OutDir'/INO80_ISW2_RSC_ISW1a/values_Heatmap.txt'



############################ INO80_ISW2_RSC_ISW1a_ISW1b_Chd1 ############################

mkdir $OutDir"/INO80_ISW2_RSC_ISW1a_ISW1b_Chd1/"

computeMatrix reference-point -S data/in_vitro_bigwigs/subNucs/INO80_ISW2_RSC_ISW1a_ISW1b_Chd1_subNucs_profile.bw \
data/in_vitro_bigwigs/monoNucs/INO80_ISW2_RSC_ISW1a_ISW1b_Chd1_monoNucs_profile.bw \
-R data/Kubik/TSS_NDR_Kubik.bed data/Kubik/TSS_unstable_Kubik.bed  \
--referencePoint TSS \
-o $OutDir"/INO80_ISW2_RSC_ISW1a_ISW1b_Chd1/cmpt_mx.gz" \
-b 1000 -a 1000 --smartLabels -p 10

plotProfile -m $OutDir"/INO80_ISW2_RSC_ISW1a_ISW1b_Chd1/cmpt_mx.gz" \
-out $OutDir'/INO80_ISW2_RSC_ISW1a_ISW1b_Chd1/DefaultHeatmap_median.png' \
--outFileNameData $OutDir'/INO80_ISW2_RSC_ISW1a_ISW1b_Chd1/values_Profile_median.txt' \
--averageType "median"

plotProfile -m $OutDir"/INO80_ISW2_RSC_ISW1a_ISW1b_Chd1/cmpt_mx.gz" \
      -out $OutDir'/INO80_ISW2_RSC_ISW1a_ISW1b_Chd1/DefaultHeatmap.png' \
      --outFileNameData $OutDir'/INO80_ISW2_RSC_ISW1a_ISW1b_Chd1/values_Profile.txt'

plotHeatmap -m $OutDir"/INO80_ISW2_RSC_ISW1a_ISW1b_Chd1/cmpt_mx.gz" \
-out $OutDir'/INO80_ISW2_RSC_ISW1a_ISW1b_Chd1/DefaultHeatmap_HEAT.png' \
--outFileNameMatrix $OutDir'/INO80_ISW2_RSC_ISW1a_ISW1b_Chd1/values_Heatmap.txt'
