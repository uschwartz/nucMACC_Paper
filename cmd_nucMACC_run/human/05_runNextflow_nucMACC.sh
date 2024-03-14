#!/usr/bin/env bash

projectDir=/home/$USER/nucMACC_Paper
AnalysisDir=/home/$USER/nucMACC_Paper/data/nucMACC_H2AZ
AnnoDir=/home/$USER/nucMACC_Paper/data

nextflow run nucMACC \
--analysis 'nucMACC' \
--csvInput $AnalysisDir'/samples_WC_nucMACC.csv' \
--outDir $AnalysisDir'/WC_nucMACC_out' \
--genomeSize  2913022398 \
--genome $AnnoDir'/GRCh38/GRCh38.primary_assembly_plus_humRibosomal.fa' \
--bamEntry \
--TSS $AnnoDir'/GRCh38/protein_coding.gtf' \
--blacklist $AnnoDir'/GRCh38/blacklist_comp_ENCFF356LFX_MT_rDNA_chr_nonCanonical.bed'
