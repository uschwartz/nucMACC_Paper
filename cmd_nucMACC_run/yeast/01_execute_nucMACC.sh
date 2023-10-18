AnalysisDir=/Users/mac-pro3/Analysis/MNase_Yeast/H4_IP
AnnoDir=/Users/mac-pro3/Analysis/MNase_Yeast/annotation
Input=$AnalysisDir'/data'
OutDir=$AnalysisDir'/results'


cd $AnalysisDir

nextflow run  ~/Analysis/nucMACC \
--csvInput $Input'/input_pooled.csv' \
--outDir $OutDir \
--genomeIdx $AnnoDir'/Bowtie2Index/genome' \
--genomeSize 12100000 \
--genome $AnnoDir'/Bowtie2Index/genome.fa' \
--TSS $AnnoDir'/Yeast_TSS.bed' \
--blacklist $AnnoDir'/blacklistChromosomes.bed' \
-w ./work_H4_IP -resume
