cd /Volumes/PromisePegasus/_Research_/nucMACC_H2AZ
AnnoPath="/Users/admin/Annotation/GRCh38/"



nextflow run ~/00_scripts/nextflow/nucMACC \
--analysis 'nucMACC' \
--csvInput 'script/addData/samples_WC_nucMACC.csv' \
--outDir WC_nucMACC_out \
--genomeSize  2913022398 \
--genome $AnnoPath'/GRCh38.primary_assembly_plus_humRibosomal.fa' \
--bamEntry \
--TSS $AnnoPath'/nextflow/protein_coding.gtf' \
--blacklist $AnnoPath'/blacklist_comp_ENCFF356LFX_MT_rDNA_chr_nonCanonical.bed'
-resume
