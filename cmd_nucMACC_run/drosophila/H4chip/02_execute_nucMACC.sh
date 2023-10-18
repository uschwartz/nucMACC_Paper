AnalysisDir=/Users/admin/Analysis/R001_nucMacc
AnnoDir=/Users/admin/Annotation/Drosophila_melanogaster_UCSC_dm3/Bowtie2Index


##run script CHANGED FEATURECOUNTS

cd $AnalysisDir

nextflow run  ~/00_scripts/nextflow/nucMACC \
--csvInput './script/H4chip/additionalFiles/H4ChIP_input.csv' \
--outDir 'H4_ChIP_v2.1_featureCOUNTS' \
--genomeIdx $AnnoDir'/genome' \
--genomeSize 162367812 \
--genome $AnnoDir'/genome.fa' \
--blacklist $AnalysisDir'/data/Annotation/dm3-blacklistChromosomes.bed' \
-w ./work_H4_v2.1 -resume
