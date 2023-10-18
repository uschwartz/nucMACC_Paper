AnalysisDir=~/Analysis/R001_nucMacc/manuscript_figures/Fig4/GROseq/
OutDir=$AnalysisDir"TSSprofile"


mkdir -p $OutDir
cd $AnalysisDir"/alignment"

samtools merge groseq_merge.bam SRR1503602_dupmark.bam SRR1503603_dupmark.bam
samtools index groseq_merge.bam
bamCoverage -b groseq_merge.bam -o profiles/groseq_merge.bw \
--binSize 10 -p 3 --normalizeUsing CPM



#### Unscale downstream

computeMatrix scale-regions -S "profiles/groseq_merge.bw" \
      -R  "../../data/TSSgroups/TSS_unstable.bed" \
       "../../data/TSSgroups/TSS_NDR.bed" \
      -o $OutDir"/computeMatrix2plot_5scale.txt.gz" \
      --unscaled5prime 200 \
      -m 2000 --skipZeros \
      -b 500 -a 500 --smartLabels -p 5



plotHeatmap -m $OutDir"/computeMatrix2plot_5scale.txt.gz" \
     -out $OutDir"/DefaultHeat_5scale.pdf" \
     --outFileNameMatrix $OutDir'/values_heat_5scale.txt'
