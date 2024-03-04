cd /Volumes/PromisePegasus/_Research_/nucMACC_H2AZ
AnnoPath="/Users/admin/Annotation/GRCh38/"


nextflow run ~/00_scripts/nextflow/nucMACC \
--analysis 'MNaseQC' \
--csvInput 'script/addData/samples_WC.csv' \
--outDir WC_QC_out \
--genomeIdx $AnnoPath'/bowtie2idx/bt2_' \
--genomeSize 2913022398 \
--publishBamFlt \
--TSS $AnnoPath'/nextflow/protein_coding.gtf' \
--blacklist $AnnoPath'/blacklist_comp_ENCFF356LFX_MT_rDNA_chr_nonCanonical.bed'
