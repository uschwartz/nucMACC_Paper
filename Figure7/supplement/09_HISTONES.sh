#!/usr/bin/env bash

AnalysisDir=/home/$USER/nucMACC_Paper/data/R001_nucMacc/manuscript_figures/Fig7
OutDir=$AnalysisDir"/histones"

mkdir -p $OutDir
cd $AnalysisDir

#1. H3

mkdir $OutDir"/H3"

computeMatrix reference-point -S data/histones/H3_rep2_MNase-ChIP-seq.bw \
    -R data/subNuc_TSS/NDR.bed data/subNuc_TSS/unstable.bed \
    --referencePoint TSS \
    -o $OutDir"/H3/cmpt_mx.gz" \
    -b 350 -a 350 --smartLabels -p 10

plotHeatmap -m $OutDir"/H3/cmpt_mx.gz" \
    -out $OutDir'/H3/DefaultHeatmap_HEAT_col_blue4_short.pdf' \
    --outFileNameMatrix $OutDir'/H3/values_Heatmap.txt' \
    --colorList '#00008B,yellow' --zMax 4 --heatmapHeight 12 \
    --heatmapWidth 5

plotHeatmap -m $OutDir"/H3/cmpt_mx.gz" \
        -out $OutDir'/H3/DefaultHeatmap_HEAT_col_blue4_short.pdf' \
        --outFileNameMatrix $OutDir'/H3/values_Heatmap.txt' \
        --colorList '#00008B,yellow' --zMax 4 --heatmapHeight 12 \
        --heatmapWidth 5

#2. H3K4me3

mkdir $OutDir"/H3K4me3"

computeMatrix reference-point -S data/histones/H3K4me3_rep2_MNase-ChIP-seq.bw \
    -R data/subNuc_TSS/NDR.bed data/subNuc_TSS/unstable.bed \
    --referencePoint TSS \
    -o $OutDir"/H3K4me3/cmpt_mx.gz" \
    -b 350 -a 350 --smartLabels -p 10

plotHeatmap -m $OutDir"/H3K4me3/cmpt_mx.gz" \
    -out $OutDir'/H3K4me3/DefaultHeatmap_HEAT_col_blue4_short.pdf' \
    --outFileNameMatrix $OutDir'/H3K4me3/values_Heatmap.txt' \
    --colorList '#00008B,yellow' --zMax 4 --heatmapHeight 12 \
    --heatmapWidth 5

plotHeatmap -m $OutDir"/H3K4me3/cmpt_mx.gz" \
        -out $OutDir'/H3K4me3/DefaultHeatmap_HEAT_col_blue4_short.pdf' \
        --outFileNameMatrix $OutDir'/H3K4me3/values_Heatmap.txt' \
        --colorList '#00008B,yellow' --zMax 8 --heatmapHeight 12 \
        --heatmapWidth 5

#3. H2B

mkdir $OutDir"/H2B"

computeMatrix reference-point -S data/histones/H2B_rep1_MNase-ChIP-seq.bw \
    -R data/subNuc_TSS/NDR.bed data/subNuc_TSS/unstable.bed \
    --referencePoint TSS \
    -o $OutDir"/H2B/cmpt_mx.gz" \
    -b 350 -a 350 --smartLabels -p 10

plotHeatmap -m $OutDir"/H2B/cmpt_mx.gz" \
    -out $OutDir'/H2B/DefaultHeatmap_HEAT_col_blue4_short.pdf' \
    --outFileNameMatrix $OutDir'/H2B/values_Heatmap.txt' \
    --colorList '#00008B,yellow' --zMax 4 --heatmapHeight 12 \
    --heatmapWidth 5

plotHeatmap -m $OutDir"/H2B/cmpt_mx.gz" \
        -out $OutDir'/H2B/DefaultHeatmap_HEAT_col_blue4_short.png' \
        --outFileNameMatrix $OutDir'/H2B/values_Heatmap.txt' \
        --colorList '#00008B,yellow' --zMax 4 --heatmapHeight 12 \
        --heatmapWidth 5

#4. H2AZ 

mkdir $OutDir"/H2AZ"

computeMatrix reference-point -S data/histones/Htz1_rep1_ChIP-exo.bw \
    -R data/subNuc_TSS/NDR.bed data/subNuc_TSS/unstable.bed \
    --referencePoint TSS \
    -o $OutDir"/H2AZ/cmpt_mx.gz" \
    -b 350 -a 350 --smartLabels -p 10

plotHeatmap -m $OutDir"/H2AZ/cmpt_mx.gz" \
    -out $OutDir'/H2AZ/DefaultHeatmap_HEAT_col_blue4_short.pdf' \
    --outFileNameMatrix $OutDir'/H2AZ/values_Heatmap.txt' \
    --colorList '#00008B,yellow'  --heatmapHeight 12 \
    --heatmapWidth 5 --zMax 23

plotHeatmap -m $OutDir"/H2AZ/cmpt_mx.gz" \
    -out $OutDir'/H2AZ/DefaultHeatmap_HEAT_col_blue4_short.png' \
    --outFileNameMatrix $OutDir'/H2AZ/values_Heatmap.txt' \
    --colorList '#00008B,yellow'  --heatmapHeight 12 \
    --heatmapWidth 5 --zMax 23
