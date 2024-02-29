#!/bin/zsh
MY_PWD=$1
cd $MY_PWD"/data"

# Declare an associative array
typeset -A sra_map

# Populate the array
sra_map[SRR3383177]="S2_H4_ChIP_100U_MNase_rep1"
sra_map[SRR3383178]="S2_H4_ChIP_25U_MNase_rep1"
sra_map[SRR3383179]="S2_H4_ChIP_6.25U_MNase_rep1"
sra_map[SRR3383180]="S2_H4_ChIP_1.5U_MNase_rep1"
sra_map[SRR3383181]="S2_H4_ChIP_100U_MNase_rep2"
sra_map[SRR3383182]="S2_H4_ChIP_25U_MNase_rep2"
sra_map[SRR3383183]="S2_H4_ChIP_6.25U_MNase_rep2"
sra_map[SRR3383184]="S2_H4_ChIP_1.5U_MNase_rep2"

# Iterate over the array keys and print them
mkdir sra
for sra_id in ${(k)sra_map}; do;
    prefetch  -O sra $sra_id
done


mkdir raw_data
for sra_id in ${(k)sra_map}; do
fastq-dump sra/${sra_id} --split-files --gzip
mv "${sra_id}_1.fastq.gz" "raw_data/${sra_map[$sra_id]}_1.fastq.gz"
mv "${sra_id}_2.fastq.gz" "raw_data/${sra_map[$sra_id]}_2.fastq.gz"
done

fastq-dump sra/SRR3383184 --split-files --gzip
mv "SRR3383184_1.fastq.gz" "raw_data/S2_H4_ChIP_1.5U_MNase_rep2_1.fastq.gz"
mv "SRR3383184_2.fastq.gz" "raw_data/S2_H4_ChIP_1.5U_MNase_rep2_2.fastq.gz"
