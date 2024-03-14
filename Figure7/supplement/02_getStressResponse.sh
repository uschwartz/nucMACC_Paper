#!/usr/bin/env bash

AnalysisDir='/Users/admin/Analysis/R001_nucMacc/manuscript_figures/Fig7'
OutDir=$AnalysisDir"/RSCprofiles"
mkdir -p $OutDir
cd $AnalysisDir

#1. RSC3

mkdir $OutDir"/RSC3/"

computeMatrix reference-point -S data/Rossi_bigwigs/Rsc3_rep2_ChIP-exo.bw \
data/Rossi_bigwigs/Rsc3_rep4_ChIP-exo.bw \
 -R data/subNuc_TSS/NDR.bed data/subNuc_TSS/unstable.bed \
 --referencePoint TSS \
 -o $OutDir"/RSC3/cmpt_mx.gz" \
 -b 500 -a 500 --smartLabels -p 10

plotProfile -m $OutDir"/RSC3/cmpt_mx.gz" \
      -out $OutDir'/RSC3/DefaultHeatmap_median.png' \
      --outFileNameData $OutDir'/RSC3/values_Profile_median.txt' \
      --averageType "median"

plotProfile -m $OutDir"/RSC3/cmpt_mx.gz" \
            -out $OutDir'/RSC3/DefaultHeatmap.png' \
            --outFileNameData $OutDir'/RSC3/values_Profile.txt'

plotHeatmap -m $OutDir"/RSC3/cmpt_mx.gz" \
      -out $OutDir'/RSC3/DefaultHeatmap_HEAT.png' \
      --outFileNameMatrix $OutDir'/RSC3/values_Heatmap.txt'

#2. RSC9

mkdir $OutDir"/RSC9/"

computeMatrix reference-point -S data/Rossi_bigwigs/Rsc9_rep3_ChIP-exo.bw \
data/Rossi_bigwigs/Rsc9_rep2_ChIP-exo.bw \
-R data/subNuc_TSS/NDR.bed data/subNuc_TSS/unstable.bed \
--referencePoint TSS \
-o $OutDir"/RSC9/cmpt_mx.gz" \
-b 500 -a 500 --smartLabels -p 10

plotProfile -m $OutDir"/RSC9/cmpt_mx.gz" \
  -out $OutDir'/RSC9/DefaultHeatmap_median.png' \
  --outFileNameData $OutDir'/RSC9/values_Profile_median.txt' \
  --averageType "median"

plotProfile -m $OutDir"/RSC9/cmpt_mx.gz" \
        -out $OutDir'/RSC9/DefaultHeatmap.png' \
        --outFileNameData $OutDir'/RSC9/values_Profile.txt'


#3. RSC1

mkdir $OutDir"/RSC1/"

computeMatrix reference-point -S data/Rossi_bigwigs/Rsc1_rep3_ChIP-exo.bw \
data/Rossi_bigwigs/Rsc1_rep2_ChIP-exo.bw \
-R data/subNuc_TSS/NDR.bed data/subNuc_TSS/unstable.bed \
--referencePoint TSS \
-o $OutDir"/RSC1/cmpt_mx.gz" \
-b 500 -a 500 --smartLabels -p 10

plotProfile -m $OutDir"/RSC1/cmpt_mx.gz" \
  -out $OutDir'/RSC1/DefaultHeatmap_median.png' \
  --outFileNameData $OutDir'/RSC1/values_Profile_median.txt' \
  --averageType "median"

plotProfile -m $OutDir"/RSC1/cmpt_mx.gz" \
        -out $OutDir'/RSC1/DefaultHeatmap.png' \
        --outFileNameData $OutDir'/RSC1/values_Profile.txt'

#4. RSC58

mkdir $OutDir"/RSC58/"

computeMatrix reference-point -S data/Rossi_bigwigs/Rsc58_rep3_ChIP-exo.bw \
data/Rossi_bigwigs/Rsc58_rep1_ChIP-exo.bw \
-R data/subNuc_TSS/NDR.bed data/subNuc_TSS/unstable.bed \
--referencePoint TSS \
-o $OutDir"/RSC58/cmpt_mx.gz" \
-b 500 -a 500 --smartLabels -p 10

plotProfile -m $OutDir"/RSC58/cmpt_mx.gz" \
  -out $OutDir'/RSC58/DefaultHeatmap_median.png' \
  --outFileNameData $OutDir'/RSC58/values_Profile_median.txt' \
  --averageType "median"

plotProfile -m $OutDir"/RSC58/cmpt_mx.gz" \
        -out $OutDir'/RSC58/DefaultHeatmap.png' \
        --outFileNameData $OutDir'/RSC58/values_Profile.txt'
