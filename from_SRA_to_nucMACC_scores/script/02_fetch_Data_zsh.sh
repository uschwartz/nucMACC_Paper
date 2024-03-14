#!/usr/bin/env bash

MY_PWD=$1
cd $MY_PWD"/data"

# Declare an associative array
typeset -A sra_map

# Populate the array with SRA IDs and corresponding sample names
sra_map[SRR3383177]="S2_H4_ChIP_100U_MNase_rep1"
sra_map[SRR3383178]="S2_H4_ChIP_25U_MNase_rep1"
sra_map[SRR3383179]="S2_H4_ChIP_6.25U_MNase_rep1"
sra_map[SRR3383180]="S2_H4_ChIP_1.5U_MNase_rep1"
sra_map[SRR3383181]="S2_H4_ChIP_100U_MNase_rep2"
sra_map[SRR3383182]="S2_H4_ChIP_25U_MNase_rep2"
sra_map[SRR3383183]="S2_H4_ChIP_6.25U_MNase_rep2"
sra_map[SRR3383184]="S2_H4_ChIP_1.5U_MNase_rep2"

# Iterate over the SRA IDs
mkdir sra
for sra_id in "${!sra_map[@]}"; do
    prefetch -O sra "$sra_id"
done

# Create new folder and dump the FASTQ files inside
mkdir raw_data

# Rename files from SRA ID to the name used in the publication
for sra_id in "${!sra_map[@]}"; do
    # Run fastq-dump command
    fastq-dump "sra/${sra_id}" --split-files --gzip
    
    # Move the generated files to the raw_data directory
    mv "${sra_id}_1.fastq.gz" "raw_data/${sra_map[$sra_id]}_1.fastq.gz"
    mv "${sra_id}_2.fastq.gz" "raw_data/${sra_map[$sra_id]}_2.fastq.gz"
done

